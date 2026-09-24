import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication

import numpy as np
from collections import defaultdict
import h5py

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
    output_name = settings["Parameters"]["output_name"].GetString()
    postprocess_parameters = settings["Parameters"]["postprocess_parameters"]

    if model_part.IsDistributed():
        raise Exception("Distributed model parts are not supported yet.")
    else:
        return GiDHDF5OutputProcess(model_part, output_name, postprocess_parameters)

    

class GiDHDF5OutputProcess(KratosMultiphysics.OutputProcess):

    #names for geometry types
    type_names = {
        0: 'Kratos_generic_type',
        1: 'Kratos_Hexahedra3D20', 
        2: 'Kratos_Hexahedra3D27', 
        3: 'Kratos_Hexahedra3D8', 
        4: 'Kratos_Prism3D15', 
        5: 'Kratos_Prism3D6', 
        6: 'Kratos_Pyramid3D13', 
        7: 'Kratos_Pyramid3D5', 
        8: 'Kratos_Quadrilateral2D4', 
        9: 'Kratos_Quadrilateral2D8', 
        10: 'Kratos_Quadrilateral2D9', 
        11: 'Kratos_Quadrilateral3D4', 
        12: 'Kratos_Quadrilateral3D8', 
        13: 'Kratos_Quadrilateral3D9', 
        14: 'Kratos_Tetrahedra3D10', 
        15: 'Tetrahedra', #'Kratos_Tetrahedra3D4', 
        16: 'Triangle', 
        17: 'Kratos_Triangle2D6', 
        18: 'Kratos_Triangle2D10', 
        19: 'Kratos_Triangle2D15', 
        20: 'Triangle', #'Kratos_Triangle3D3', 
        21: 'Kratos_Triangle3D6', 
        22: 'Kratos_Line2D2', 
        23: 'Kratos_Line2D3', 
        24: 'Kratos_Line2D4', 
        25: 'Kratos_Line2D5', 
        26: 'Kratos_Line3D2', 
        27: 'Kratos_Line3D3', 
        28: 'Kratos_Point2D', 
        29: 'Kratos_Point3D', 
        30: 'Kratos_Sphere3D1', 
        31: 'Kratos_Nurbs_Curve', 
        32: 'Kratos_Nurbs_Surface', 
        33: 'Kratos_Nurbs_Volume', 
        34: 'Kratos_Nurbs_Curve_On_Surface', 
        35: '???', 
        36: 'Kratos_Surface_In_Nurbs_Volume', 
        37: 'Kratos_Brep_Curve', 
        38: 'Kratos_Brep_Surface', 
        39: '???', 
        40: 'Kratos_Brep_Curve_On_Surface', 
        41: '???', 
        42: 'Kratos_Quadrature_Point_Geometry', 
        43: 'Kratos_Coupling_Geometry'
        }
    
    
    defaults_results_file_configuration = KratosMultiphysics.Parameters('''
                {
                    "output_control_type": "step",
                    "output_interval": 1.0,
                    "nodal_results": [],
                    "nodal_nonhistorical_results": [],
                    "nodal_flags_results": [],
                    "elemental_conditional_flags_results": [],
                    "gauss_point_results": [],
                    "dtype":"float32"
                }
            ''')

        
    def __init__(self, model_part, output_name, parameters):
        super().__init__()
        self.model_part = model_part
        self.parameters = parameters

        self.parameters["result_file_configuration"].ValidateAndAssignDefaults(self.defaults_results_file_configuration)

        self.filename = output_name

        self.list_of_variables_with_history = self._GetListOfVariables(self.parameters["result_file_configuration"]["nodal_results"])
        self.list_of_variables_with_no_history = self._GetListOfVariables(self.parameters["result_file_configuration"]["nodal_nonhistorical_results"])
        self.list_of_flags = self._GetListOfVariables(self.parameters["result_file_configuration"]["nodal_flags_results"])
        self.list_of_gauss_point_variables = self._GetListOfVariables(self.parameters["result_file_configuration"]["gauss_point_results"])
        if self.parameters["result_file_configuration"]["dtype"].GetString() == "float32":
            self.dtype=np.float32
        else:
            self.dtype=np.float64
        self.int_type = np.int32 #by default. we will check if this is sufficient and eventually upgrade

        controller_settings = KratosMultiphysics.Parameters("""{}""")
        controller_settings.AddString("model_part_name", model_part.FullName())
        if self.parameters["result_file_configuration"].Has("output_control_type"): controller_settings.AddValue("output_control_type", self.parameters["result_file_configuration"]["output_control_type"])
        if self.parameters["result_file_configuration"].Has("output_interval"): controller_settings.AddValue("output_interval", self.parameters["result_file_configuration"]["output_interval"])
        self.controller = KratosMultiphysics.OutputController(model_part.GetModel(), controller_settings)

        if(controller_settings["output_control_type"].GetString() == "step"):
            self.time_label_var = KratosMultiphysics.STEP
        else:
            self.time_label_var = KratosMultiphysics.TIME

        self.printing_index = 1 #unique identity of a given output step
        
        self.f = h5py.File(self.filename, 'w', libver="latest")
        self.f.swmr_mode = True #makes it possible to read the file while it is being written to by another process

        dt = np.dtype('S5')
        self.f.attrs.create('GiD Post Results File', data="1.1", dtype=dt)

    def IsOutputStep(self):
        return self.controller.Evaluate()

    def PrintOutput(self):
        
        if(self.IsOutputStep()):
            time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
            self._WriteNodalSolutionStepResults(step, time, self.list_of_variables_with_history)
            # self._WriteNodalSolutionStepResults(self.printing_index,step, time, self.list_of_variables_with_history)
            # self._WriteNodalSolutionStepResults(self.printing_index,step, time, self.list_of_variables_with_history)
            # self._WriteNodalSolutionStepResults(self.printing_index,step, time, self.list_of_variables_with_history)
            self.printing_index += 1
            self.f.flush() #flush to disk to make sure the results are written before the next output step

        # Schedule next output
        self.controller.Update()


    def _GetListOfVariables(self,var_names):
        list_of_variables_with_history = []
        for var_name in var_names.values():
            variable = KratosMultiphysics.KratosGlobals.GetVariable(var_name.GetString())
            list_of_variables_with_history.append(variable)
        return list_of_variables_with_history


    def __categorize_geometries_by_geometry_type(self, geometries_container):
        """
        Categorizes geometries from a geometry container into their respective
        geometry types based on the KratosGeometryType enum.

        Parameters
        ----------
        geometries_container : iterable
            A container of geometry objects (e.g., ModelPart.Geometries,
            geometry list from submodelpart).

        Returns
        -------
        list
            A list of KratosMultiphysics.GeometryContainerType, where each index corresponds to a geometry type as defined in KratosMultiphysics.GeometryData_KratosGeometryType enum. Each
        """
        #TODO: make a utility to do this operation in c++
        #initialize "categorized" with empty lists for each geometry type
        categorized = []
        for i in range(len(KratosMultiphysics.GeometryData_KratosGeometryType.__members__)):
            categorized.append(KratosMultiphysics.GeometryContainerType())

        for geometry in geometries_container:
            geometry_type = geometry.GetGeometryType()
            categorized[geometry_type][geometry.Id] = geometry
        
        return categorized

    def ExecuteInitialize(self):
        self.categorized_geometries = self.__categorize_geometries_by_geometry_type(self.model_part.Geometries)

        # #TODO: make it efficient by using a tensor adaptor for this
        self.node_ids = []
        for node in self.model_part.Nodes:
            self.node_ids.append(node.Id)
        
        max_id = max(self.node_ids)
        if max_id > np.iinfo(np.int32).max:
            self.int_type = np.int64 #note that this is int32 by defaut

        self.node_ids = np.array(self.node_ids, dtype=self.int_type)

        self.__WriteMesh()
        
        self.base_results_group = self.f.create_group("Results")
        #self.base_mesh_group.attrs['ModelPartName'] = self.model_part.Name



    def __WriteMesh(self):
        
        base_mesh_group = self.f.create_group("Meshes")
        

        mesh_id = 1
        for i in range(len(self.categorized_geometries)):
            item = self.categorized_geometries[i]
            type_name = self.type_names[int(KratosMultiphysics.GeometryData_KratosGeometryType(i))]

            if len(item) != 0:
                mesh_group = base_mesh_group.create_group(str(mesh_id))

                nodes_group = mesh_group.create_group("Coordinates") #leave it empty, as coordinates are stored in the "Coordinates" group at the mesh level
                if mesh_id  ==1:
                    ##print node ids and coordinates
                    #get coordinates
                    ta = KratosMultiphysics.TensorAdaptors.NodePositionTensorAdaptor(self.model_part.Nodes, KratosMultiphysics.Configuration.Current)
                    ta.Check()
                    ta.CollectData()
                    xyz = ta.data

                    nodes_group.create_dataset('1', data=self.node_ids,compression="szip",compression_opts=("nn", 8))
                    nodes_group.create_dataset('2', data=xyz[:,0],compression="szip",compression_opts=("nn", 8))
                    nodes_group.create_dataset('3', data=xyz[:,1],compression="szip",compression_opts=("nn", 8))
                    nodes_group.create_dataset('4', data=xyz[:,2],compression="szip",compression_opts=("nn", 8))
                
                
                conn_group = mesh_group.create_group("Elements")

                # Get connectivity using TensorAdaptor
                connectivity_adaptor = KratosMultiphysics.TensorAdaptors.ConnectivityIdsTensorAdaptor(item)
                connectivity_adaptor.CollectData()
                connectivity = connectivity_adaptor.data  # numpy array: shape (num_geometries, num_nodes_per_elem)

                dtype_dim = np.dtype('S1')
                dtype_type = np.dtype('S50')
                dtype_nnode = np.dtype('S5')
                mesh_group.attrs.create('Dimension', data="3", dtype=dtype_dim)
                mesh_group.attrs.create('ElemType', data=type_name, dtype=dtype_type)
                mesh_group.attrs.create('Name', data=type_name, dtype=dtype_type)
                mesh_group.attrs.create('Nnode', data=str(connectivity.shape[1]), dtype=dtype_nnode)

                #TODO: make this efficient
                ids = []
                for geom in item:
                    ids.append(geom.Id)
                ids = np.array(ids, dtype=self.int_type)

                conn_group.create_dataset('1', data=ids,compression="szip",compression_opts=("nn", 8))
                for i in range(connectivity.shape[1]):
                    conn_group.create_dataset(str(i+2), data=connectivity[:, i],compression="szip",compression_opts=("nn", 8))
                conn_group.create_dataset(str(connectivity.shape[1]+2), data=np.zeros(connectivity.shape[0], dtype=self.int_type),compression="szip",compression_opts=("nn", 8)) #dummy dataset for element partition, as we are not supporting distributed meshes yet
                mesh_id += 1
        self.f.flush() #flush to disk to make sure the mesh is written before results are written

    def WriteResultsAttribute(self, group, variable):
        dtype_analysis = np.dtype('S6')
        dtype_name = np.dtype('S50')
        dtype_step = np.dtype('S10')
        dtype_time = np.dtype('S20')
        dtype_single = np.dtype('S1')
        dtype_location = np.dtype('S7')
        dtype_result_type = np.dtype('S6')
        
        group.attrs.create('Analysis', data="Kratos", dtype=dtype_analysis)
        group.attrs.create('Name', data=variable.Name(), dtype=dtype_name)
        group.attrs.create('Step', data=str(self.model_part.ProcessInfo[self.time_label_var]), dtype=dtype_step)
        #group.attrs.create('Time', data=str(self.model_part.ProcessInfo[KratosMultiphysics.TIME]), dtype=dtype_time)

        if(type(variable) == KratosMultiphysics.DoubleVariable):
            group.attrs.create('NumComponents', data="1", dtype=dtype_single)
            group.attrs.create('ResultLocation', data="OnNodes", dtype=dtype_location)
            group.attrs.create('ResultType', data="Scalar", dtype=dtype_result_type)
            group.attrs.create('Component 1', data=variable.Name(), dtype=dtype_name)
        elif(type(variable) == KratosMultiphysics.Array1DVariable3):
            group.attrs.create('NumComponents', data="3", dtype=dtype_single)
            group.attrs.create('ResultLocation', data="OnNodes", dtype=dtype_location)
            group.attrs.create('ResultType', data="Vector", dtype=dtype_result_type)
            group.attrs.create('Component 1', data=variable.Name() + "_X", dtype=dtype_name)
            group.attrs.create('Component 2', data=variable.Name() + "_Y", dtype=dtype_name)
            group.attrs.create('Component 3', data=variable.Name() + "_Z", dtype=dtype_name)

    def _WriteNodalSolutionStepResults(self, step, time, list_of_variables_with_history):
        """
        Writes nodal solution step results for specified variables to the HDF5 file. The results are stored in a group named 'NodalResults', with each variable having its own subgroup containing the variable values for all nodes.

        Parameters
        ----------
        list_of_variables_with_history : list of KratosMultiphysics.Variable
            A list of KratosMultiphysics.Variable objects representing the nodal variables to be written to the HDF5 file.
        """

        for variable in list_of_variables_with_history:
            var_name = variable.Name()

            ta = KratosMultiphysics.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes,variable,0)
            ta.CollectData()
            var_values = ta.data

            nodal_results_group = self.base_results_group.create_group(str(self.printing_index)) #dataset is identified by an unique integer "printing_index"

            self.WriteResultsAttribute(nodal_results_group, variable)
            nodal_results_group.create_dataset(str(1), data=self.node_ids,compression="szip",compression_opts=("nn", 8))
            if(len(var_values.shape) == 1):
                nodal_results_group.create_dataset(str(2), data=var_values,compression="szip",compression_opts=("nn", 8))
            else:
                for i in range(var_values.shape[1]):
                    nodal_results_group.create_dataset(str(i+2), data=var_values[:,i],compression="szip",compression_opts=("nn", 8))

            self.printing_index += 1


 
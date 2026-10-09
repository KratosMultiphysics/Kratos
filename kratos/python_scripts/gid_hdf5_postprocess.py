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
        1: 'Hexahedra', #Kratos_Hexahedra3D20', 
        2: 'Hexahedra', #Kratos_Hexahedra3D27', 
        3: 'Hexahedra', #Kratos_Hexahedra3D8', 
        4: 'Prism', #Kratos_Prism3D15', 
        5: 'Prism', #Kratos_Prism3D6', 
        6: 'Prism', #Kratos_Pyramid3D13', 
        7: 'Pyramid', #Kratos_Pyramid3D5', 
        8: 'Quadrilateral', #Kratos_Quadrilateral2D4', 
        9: 'Quadrilateral', #Kratos_Quadrilateral2D8', 
        10: 'Quadrilateral', #Kratos_Quadrilateral2D9', 
        11: 'Quadrilateral', #Kratos_Quadrilateral3D4', 
        12: 'Quadrilateral', #Kratos_Quadrilateral3D8', 
        13: 'Quadrilateral', #Kratos_Quadrilateral3D9', 
        14: 'Tetrahedra', #Kratos_Tetrahedra3D10', 
        15: 'Tetrahedra', #'Kratos_Tetrahedra3D4', 
        16: 'Triangle', 
        17: 'Triangle', #Kratos_Triangle2D6', 
        18: 'Triangle', #Kratos_Triangle2D10', 
        19: 'Triangle', #Kratos_Triangle2D15', 
        20: 'Triangle', #'Kratos_Triangle3D3', 
        21: 'Triangle', #Kratos_Triangle3D6', 
        22: 'Line', #Kratos_Line2D2', 
        23: 'Line', #Kratos_Line2D3', 
        24: 'Line', #Kratos_Line2D4', 
        25: 'Line', #Kratos_Line2D5', 
        26: 'Line', #Kratos_Line3D2', 
        27: 'Line', #Kratos_Line3D3', 
        28: 'Point', #Kratos_Point2D', 
        29: 'Point', #Kratos_Point3D', 
        30: 'Point', 
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
         43: 'Kratos_Coupling_Geometry',
         44: 'Kratos_Quadrature_Point_Curve_Geometry',
         45: 'Kratos_Quadrature_Point_Curve_On_Surface_Geometry',
         46: 'Kratos_Quadrature_Point_Surface_In_Volume_Geometry'
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
        self.list_of_flags = self._GetListOfFlags(self.parameters["result_file_configuration"]["nodal_flags_results"])
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
            self._WriteNodalNonHistoricalResults(step, time, self.list_of_variables_with_no_history)
            self._WriteNodalFlagsResults(step, time, self.list_of_flags)
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

    def _GetListOfFlags(self,flag_names):
        list_of_flags = []
        for flag_name in flag_names.values():
            name = flag_name.GetString()
            flag = getattr(KratosMultiphysics, name, None)
            if not isinstance(flag, KratosMultiphysics.Flags):
                available_flags = sorted(n for n in dir(KratosMultiphysics) if isinstance(getattr(KratosMultiphysics, n, None), KratosMultiphysics.Flags))
                raise Exception(f"Flag {name} is unknown. Available flags: {available_flags}")
            list_of_flags.append((name, flag))
        return list_of_flags


    def __categorize_entities_by_geometry_type(self):
        """
        Categorizes the model part entities by geometry type using
        ContainerClassificationUtilities.

        If the model part contains elements, the elements are classified and
        the mesh is printed keyed by element id. Otherwise the standalone
        geometries of the model part are classified and printed keyed by
        geometry id. A geometry without an explicitly assigned id raises an
        error.

        Returns
        -------
        list
            A list of containers where each index corresponds to a geometry
            type as defined in the KratosMultiphysics.GeometryData_KratosGeometryType
            enum. Each container holds the entities of that type.
        """
        classifier = KratosMultiphysics.ContainerClassificationUtilities()
        if len(self.model_part.Elements) > 0:
            return classifier.ClassifyByGeometryType(self.model_part.Elements)
        else:
            for node in self.model_part.Geometries:
                if node.IsIdSelfAssigned() or node.IsIdGeneratedFromString():
                    raise Exception(f"Geometry with id {node.Id} does not have an explicitly assigned id. "
                                    "Assign an id to all geometries before printing the mesh with the GiD HDF5 postprocess.")
            return classifier.ClassifyByGeometryType(self.model_part.Geometries)


    def ExecuteInitialize(self):
        self.categorized_entities = self.__categorize_entities_by_geometry_type()
        
        self.node_ids = self.model_part.Nodes.IdsList()
        
        max_id = max(self.node_ids)
        if max_id > np.iinfo(np.int32).max:
            self.int_type = np.int64 #note that this is int32 by defaut

        self.node_ids = np.array(self.node_ids, dtype=self.int_type)

        self.__WriteMesh()
        
        self.base_results_group = self.f.create_group("Results")

    def _CreateDataset(self, group, name, data):
        """Create a 1-D HDF5 dataset, applying szip compression only when the
        fastest-varying axis is at least as long as the szip block size (8).
        szip cannot compress shorter blocks (e.g. a few elements in a small
        mesh), so those datasets are stored uncompressed.
        """
        if data.shape[-1] >= 8:
            group.create_dataset(name, data=data, compression="szip", compression_opts=("nn", 8))
        else:
            group.create_dataset(name, data=data)

    def __WriteMesh(self):
        
        base_mesh_group = self.f.create_group("Meshes")
        

        mesh_id = 1
        for i in range(len(self.categorized_entities)):
            item = self.categorized_entities[i]
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

                    self._CreateDataset(nodes_group, '1', self.node_ids)
                    self._CreateDataset(nodes_group, '2', xyz[:,0])
                    self._CreateDataset(nodes_group, '3', xyz[:,1])
                    self._CreateDataset(nodes_group, '4', xyz[:,2])
                
                
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

                ids = item.IdsList()

                self._CreateDataset(conn_group, '1', ids)
                for i in range(connectivity.shape[1]):
                    self._CreateDataset(conn_group, str(i+2), connectivity[:, i])
                self._CreateDataset(conn_group, str(connectivity.shape[1]+2), np.zeros(connectivity.shape[0], dtype=self.int_type)) #dummy dataset for element partition, as we are not supporting distributed meshes yet
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

    def WriteFlagResultsAttribute(self, group, flag_name):
        dtype_analysis = np.dtype('S6')
        dtype_name = np.dtype('S50')
        dtype_step = np.dtype('S10')
        dtype_single = np.dtype('S1')
        dtype_location = np.dtype('S7')
        dtype_result_type = np.dtype('S6')

        group.attrs.create('Analysis', data="Kratos", dtype=dtype_analysis)
        group.attrs.create('Name', data=flag_name, dtype=dtype_name)
        group.attrs.create('Step', data=str(self.model_part.ProcessInfo[self.time_label_var]), dtype=dtype_step)
        group.attrs.create('NumComponents', data="1", dtype=dtype_single)
        group.attrs.create('ResultLocation', data="OnNodes", dtype=dtype_location)
        group.attrs.create('ResultType', data="Scalar", dtype=dtype_result_type)
        group.attrs.create('Component 1', data=flag_name, dtype=dtype_name)

    def _WriteNodalResultsDatasets(self, nodal_results_group, var_values):
        self._CreateDataset(nodal_results_group, str(1), self.node_ids)
        if(len(var_values.shape) == 1):
            self._CreateDataset(nodal_results_group, str(2), var_values)
        else:
            for i in range(var_values.shape[1]):
                self._CreateDataset(nodal_results_group, str(i+2), var_values[:,i])

    def _WriteNodalSolutionStepResults(self, step, time, list_of_variables_with_history):
        """
        Writes nodal solution step results for specified variables to the HDF5 file. The results are stored in a group named 'NodalResults', with each variable having its own subgroup containing the variable values for all nodes.

        Parameters
        ----------
        list_of_variables_with_history : list of KratosMultiphysics.Variable
            A list of KratosMultiphysics.Variable objects representing the nodal variables to be written to the HDF5 file.
        """

        for variable in list_of_variables_with_history:
            ta = KratosMultiphysics.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes,variable,0)
            ta.CollectData()

            nodal_results_group = self.base_results_group.create_group(str(self.printing_index)) #dataset is identified by an unique integer "printing_index"

            self.WriteResultsAttribute(nodal_results_group, variable)
            self._WriteNodalResultsDatasets(nodal_results_group, ta.data)

            self.printing_index += 1

    def _WriteNodalNonHistoricalResults(self, step, time, list_of_variables_with_no_history):
        """
        Writes nodal non-historical results for specified variables to the HDF5 file, using the current (non-historical) storage of each node.
        """

        for variable in list_of_variables_with_no_history:
            ta = KratosMultiphysics.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes,variable)
            ta.CollectData()

            nodal_results_group = self.base_results_group.create_group(str(self.printing_index))

            self.WriteResultsAttribute(nodal_results_group, variable)
            self._WriteNodalResultsDatasets(nodal_results_group, ta.data)

            self.printing_index += 1

    def _WriteNodalFlagsResults(self, step, time, list_of_flags):
        """
        Writes nodal flags results to the HDF5 file. Each flag is written as an integer (0/1) scalar field over the nodes.
        """

        for flag_name, flag in list_of_flags:
            ta = KratosMultiphysics.TensorAdaptors.FlagsTensorAdaptor(self.model_part.Nodes,flag)
            ta.CollectData()

            nodal_results_group = self.base_results_group.create_group(str(self.printing_index))

            self.WriteFlagResultsAttribute(nodal_results_group, flag_name)
            self._WriteNodalResultsDatasets(nodal_results_group, np.array(ta.data, dtype=np.float32)) #note that we are converting to float32 as required by the output format

            self.printing_index += 1
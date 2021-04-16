# Import system python
import os

import numpy as np # 
import h5py
import KratosMultiphysics as Kratos
import KratosMultiphysics.ConvectionDiffusionApplication as CD
class WriteConvergenceNodalErrorToHdf5:
    def __init__(self, destination_model_part, velocity_error, concentration_error):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        test_number -- It is the number of the mesh that is computing
        """
        
        file_path = os.getcwd()
        dir_path = os.path.dirname(file_path)
        self.file_path = dir_path+'/Norouzi_mesh_2_MeshConvergenceL2ErrorNorm.hdf5'
        self.model_part_fluid = destination_model_part  
        self.group_name = str(self.model_part_fluid.ProcessInfo[Kratos.STEP])
        self.velocity_error = velocity_error
        self.concentration_error = concentration_error

        num_nodes = self.model_part_fluid.NumberOfNodes()
        node_id_array = [0]*num_nodes
        coordinate_array = [0]*num_nodes
        velocity_array = [0]*num_nodes
        temperature_array = [0]*num_nodes
        vectorial_error_array = [0]*num_nodes
        scalar_error_array = [0]*num_nodes
        # np.zeros(num_nodes)
        for node_i, node in enumerate(self.model_part_fluid.Nodes):
            node_id_array[node_i] = node.Id
            velocity_array[node_i] = node.GetSolutionStepValue(Kratos.VELOCITY)
            temperature_array[node_i] = node.GetSolutionStepValue(Kratos.TEMPERATURE)
            coordinate_array[node_i] = node
            vectorial_error_array[node_i] = node.GetSolutionStepValue(CD.VECTORIAL_MESH_ERROR)
            scalar_error_array[node_i] = node.GetSolutionStepValue(CD.SCALAR_MESH_ERROR)
        
        self.node_id = node_id_array
        self.velocities = velocity_array
        self.temperatures = temperature_array
        self.coordinates = coordinate_array
        self.vectorial_error = vectorial_error_array
        self.scalar_error = scalar_error_array
        
        for Element in self.model_part_fluid.Elements:
            self.element_size = Element.GetGeometry().Length()
            break
     
        with h5py.File(self.file_path, 'a') as f:
                self.WriteDataToFile(file_or_group = f,
                            names = ['COORDINATES','VELOCITIES','TEMPERATURES','ID','VELOCITY_NODAL_ERROR', 'CONCENTRATION_NODAL_ERROR', 'VELOCITY_ERROR', 'CONCENTRATION_ERROR'],
                            data = [self.coordinates, self.velocities, self.temperatures, self.node_id, self.vectorial_error, self.scalar_error, self.velocity_error, self.concentration_error])

    def WriteDataToFile(self, file_or_group, names, data):
        if self.group_name in file_or_group:
            file_or_group['/'].__delitem__(self.group_name)
            self.sub_group = file_or_group.create_group(self.group_name)
        else:
            self.sub_group = file_or_group.create_group(self.group_name)
        self.sub_group.attrs['element_size'] = str(self.element_size)
        self.sub_group.attrs['n_elements'] = str(len(self.model_part_fluid.Elements))
        self.sub_group.attrs['time'] = str(self.model_part_fluid.ProcessInfo[Kratos.TIME])
        for name, datum in zip(names, data):
            if name in file_or_group:
                file_or_group.__delitem__(name)
        for name, datum in zip(names, data):
            self.sub_group.create_dataset(name = name, data = datum)
        file_or_group.close()


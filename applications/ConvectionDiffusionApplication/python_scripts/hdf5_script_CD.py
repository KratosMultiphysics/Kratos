# Import system python
import os

import numpy as np # 
import h5py
import KratosMultiphysics as Kratos
#Might be a process
# class ErrorProjectionPostProcessTool(Kratos.Process):
class ErrorProjectionPostProcessTool(object):
    def __init__(self, file_name, fluid_sub_modelpart):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        test_number -- It is the number of the mesh that is computing
        """
        # self.parameters = Kratos.Parameters( """
        # {
        #     "file_name": "sp_data.hdf5",
        #     "target_porosity" : 0.3,
        #     "probe_height" : 0.2032
        # }  """ )

        self.file_path = file_name
        self.model_part_fluid = fluid_sub_modelpart  
        self.group_name = str(self.model_part_fluid.ProcessInfo[Kratos.STEP])

        num_nodes = self.model_part_fluid.NumberOfNodes()
        node_id_array = [0]*num_nodes
        coordinate_array = [0]*num_nodes
        velocity_array = [0]*num_nodes
        temperature_array = [0]*num_nodes
        # np.zeros(num_nodes)
        for node_i, node in enumerate(self.model_part_fluid.Nodes):
            node_id_array[node_i] = node.Id
            velocity_array[node_i] = node.GetSolutionStepValue(Kratos.VELOCITY)
            temperature_array[node_i] = node.GetSolutionStepValue(Kratos.TEMPERATURE)
            coordinate_array[node_i] = node
        
        self.node_id = node_id_array
        self.velocities = velocity_array
        self.temperatures = temperature_array
        self.coordinates = coordinate_array

        
        for Element in self.model_part_fluid.Elements:
            self.element_size = Element.GetGeometry().Length()
            break
     
        with h5py.File(self.file_path, 'a') as f:
                self.WriteDataToFile(file_or_group = f,
                            names = ['COORDINATES', 'VELOCITIES', 'TEMPERATURES', 'ID'],
                            data = [self.coordinates, self.velocities, self.temperatures, self.node_id])

    def WriteDataToFile(self, file_or_group, names, data):
        if self.group_name in file_or_group:
            file_or_group['/'].__delitem__(self.group_name)
            self.sub_group = file_or_group.create_group(self.group_name)
        else:
            self.sub_group = file_or_group.create_group(self.group_name)
        # self.sub_group.attrs['element_size'] = str(self.element_size)
        self.sub_group.attrs['n_elements'] = str(len(self.model_part_fluid.Elements))
        # self.sub_group.attrs['delta_time'] = str(self.model_part_fluid.ProcessInfo[Kratos.DELTA_TIME])
        self.sub_group.attrs['time'] = str(self.model_part_fluid.ProcessInfo[Kratos.TIME])
        for name, datum in zip(names, data):
            if name in file_or_group:
                file_or_group.__delitem__(name)
        for name, datum in zip(names, data):
            self.sub_group.create_dataset(name = name, data = datum)
        file_or_group.close()
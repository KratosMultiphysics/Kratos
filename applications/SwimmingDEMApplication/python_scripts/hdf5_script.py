# Import system python
import os

import numpy as np
import h5py
import KratosMultiphysics as Kratos
import KratosMultiphysics.SwimmingDEMApplication as Dem

class ErrorProjectionPostProcessTool(object):
    def __init__(self, test_number):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        test_number -- It is the number of the mesh that is computing
        """
        self.parameters = Kratos.Parameters( """
        {
            "file_name": "sp_data.hdf5",
            "target_porosity" : 0.3,
            "probe_height" : 0.2032
        }  """ )

        self.time = []
        self.v_error = []
        self.p_error = []
        self.av_mod_error = []
        self.problem_path = os.getcwd()
        self.file_path = os.path.join(str(self.problem_path),self.parameters["file_name"].GetString())
        self.dtype = np.float64
        self.group_name = str(test_number)

    def WriteData(self, error_model_part, velocity_error_projected, pressure_error_projected):
        self.error_model_part = error_model_part
        mod_error = 0.0
        iterator = 0
        self.element_size = []
        for node in self.error_model_part.Nodes:
            mod_error += np.sqrt(node.GetSolutionStepValue(Dem.ERROR_X)**2 + node.GetSolutionStepValue(Dem.ERROR_Y)**2 + node.GetSolutionStepValue(Dem.ERROR_Z)**2)
        self.av_mod_error.append(mod_error/ len(self.error_model_part.Elements))

        for Element in self.error_model_part.Elements:
            self.element_size.append(Element.GetGeometry().Length())
            iterator += 1
        #self.max_element = np.sum(self.element_size)/iterator
        self.max_element = np.amax(self.element_size)
        self.time.append(self.error_model_part.ProcessInfo[Kratos.TIME])
        self.v_error.append(velocity_error_projected)
        self.p_error.append(pressure_error_projected)
        with h5py.File(self.file_path, 'a') as f:
                self.WriteDataToFile(file_or_group = f,
                            names = ['TIME', 'V_ERROR', 'P_ERROR', 'AVERAGE_ERROR'],
                            data = [self.time, self.v_error, self.p_error, self.av_mod_error])

    def WriteDataToFile(self, file_or_group, names, data):
        overwrite_previous = True
        if self.group_name in file_or_group:
            if overwrite_previous:
                file_or_group['/'].__delitem__(self.group_name)
                self.sub_group = file_or_group.create_group(self.group_name)
            else:
                self.sub_group = file_or_group['/' + self.group_name]
        self.sub_group.attrs['element_size'] = str(self.max_element)
        self.sub_group.attrs['n_elements'] = str(len(self.error_model_part.Elements))
        for name, datum in zip(names, data):
            if name in file_or_group:
                file_or_group.__delitem__(name)
        for name, datum in zip(names, data):
            self.sub_group.create_dataset(name = name, data = datum)
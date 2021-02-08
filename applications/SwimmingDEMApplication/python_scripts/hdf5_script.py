# Import system python
import os

import numpy as np
import h5py
import KratosMultiphysics as Kratos

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

        for Element in self.error_model_part.Elements:
            self.element_size = Element.GetGeometry().Length()
            break

        self.time.append(self.error_model_part.ProcessInfo[Kratos.TIME])
        self.v_error.append(velocity_error_projected)
        self.p_error.append(pressure_error_projected)
        with h5py.File(self.file_path, 'a') as f:
                self.WriteDataToFile(file_or_group = f,
                            names = ['TIME', 'V_ERROR', 'P_ERROR'],
                            data = [self.time, self.v_error, self.p_error])

    def WriteDataToFile(self, file_or_group, names, data):
        if self.group_name in file_or_group:
            file_or_group['/'].__delitem__(self.group_name)
            self.sub_group = file_or_group.create_group(self.group_name)
        else:
            self.sub_group = file_or_group.create_group(self.group_name)
        self.sub_group.attrs['element_size'] = str(self.element_size)
        self.sub_group.attrs['n_elements'] = str(len(self.error_model_part.Elements))
        self.sub_group.attrs['delta_time'] = str(self.error_model_part.ProcessInfo[Kratos.DELTA_TIME])
        for name, datum in zip(names, data):
            if name in file_or_group:
                file_or_group.__delitem__(name)
        for name, datum in zip(names, data):
            self.sub_group.create_dataset(name = name, data = datum)
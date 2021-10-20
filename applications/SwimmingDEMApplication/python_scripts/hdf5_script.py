# Import system python
import os

import numpy as np
import h5py
import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.FluidDynamicsApplication as Fluid

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
        self.mean_iteration = []
        self.problem_path = os.getcwd()
        self.file_path = os.path.join(str(self.problem_path),self.parameters["file_name"].GetString())
        self.dtype = np.float64
        self.group_name = str(test_number)

    def WriteData(self, error_model_part, velocity_error_projected, pressure_error_projected, projection_type, model_type, subscale_type):
        self.error_model_part = error_model_part

        self.projection_type = projection_type
        self.model_type = model_type
        self.subscale_type = subscale_type

        for Element in self.error_model_part.Elements:
            self.element_size = Element.GetGeometry().Length()
            break

        iterations = 0.0
        for Element in self.error_model_part.Elements:
            iterations += Element.GetValue(Fluid.ADJOINT_FLUID_SCALAR_1)

        self.mean_iteration.append(iterations/len(self.error_model_part.Elements))
        self.time.append(self.error_model_part.ProcessInfo[Kratos.TIME])
        self.v_error.append(velocity_error_projected)
        self.p_error.append(pressure_error_projected)

        with h5py.File(self.file_path, 'a') as f:
                self.WriteDataToFile(file_or_group = f,
                            names = ['TIME', 'V_ERROR', 'P_ERROR', 'MEAN_ITERATION'],
                            data = [self.time, self.v_error, self.p_error, self.mean_iteration])

    def WriteDataToFile(self, file_or_group, names, data):
        if self.group_name in file_or_group:
            file_or_group['/'].__delitem__(self.group_name)
            self.sub_group = file_or_group.create_group(self.group_name)
        else:
            self.sub_group = file_or_group.create_group(self.group_name)
        self.sub_group.attrs['element_size'] = str(self.element_size)
        self.sub_group.attrs['n_elements'] = str(len(self.error_model_part.Elements))
        self.sub_group.attrs['delta_time'] = str(self.error_model_part.ProcessInfo[Kratos.DELTA_TIME])
        self.sub_group.attrs['projection_type'] = str(self.projection_type)
        self.sub_group.attrs['model_type'] = str(self.model_type)
        self.sub_group.attrs['subscale_type'] = str(self.subscale_type)

        for name, datum in zip(names, data):
            if name in file_or_group:
                file_or_group.__delitem__(name)
        for name, datum in zip(names, data):
            self.sub_group.create_dataset(name = name, data = datum)

class ParticleDragForcePostProcessTool(object):
    def __init__(self):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.parameters = Kratos.Parameters( """
        {
            "file_name": "sp_data.hdf5",
            "target_porosity" : 0.3,
            "probe_height" : 0.2032
        }  """ )


        self.problem_path = os.getcwd()
        self.file_path = os.path.join(str(self.problem_path),self.parameters["file_name"].GetString())
        self.dtype = np.float64

    def WriteData(self, spheres_model_part, fluid_model_part):
        self.spheres_model_part = spheres_model_part
        self.group_name = str(fluid_model_part.ProcessInfo[Kratos.TIME])

        self.id_particle = []
        self.position_x = []
        self.position_y = []
        self.position_z = []
        self.radius = []
        self.drag_force_x = []
        self.drag_force_y = []
        self.drag_force_z = []
        self.y = []
        self.porosity = []
        self.slip_velocity_x = []
        self.slip_velocity_y = []
        self.slip_velocity_z = []

        for node in self.spheres_model_part.Nodes:
            self.id_particle.append(node.Id)
            self.position_x.append(node.X)
            self.position_y.append(node.Y)
            self.position_z.append(node.Z)
            self.radius.append(node.GetSolutionStepValue(Kratos.RADIUS,0))
            self.drag_force_x.append(node.GetSolutionStepValue(Kratos.DRAG_FORCE_X))
            self.drag_force_y.append(node.GetSolutionStepValue(Kratos.DRAG_FORCE_Y))
            self.drag_force_z.append(node.GetSolutionStepValue(Kratos.DRAG_FORCE_Z))
            self.porosity.append(node.GetSolutionStepValue(Kratos.FLUID_FRACTION_PROJECTED))
            self.slip_velocity_x.append(node.GetSolutionStepValue(Kratos.SLIP_VELOCITY_X))
            self.slip_velocity_y.append(node.GetSolutionStepValue(Kratos.SLIP_VELOCITY_Y))
            self.slip_velocity_z.append(node.GetSolutionStepValue(Kratos.SLIP_VELOCITY_Z))
            self.y.append(node.GetSolutionStepValue(DEM.IMPACT_WEAR))

        #self.time.append(self.spheres_model_part.ProcessInfo[Kratos.TIME])
        with h5py.File(self.file_path, 'a') as f:
                self.WriteDataToFile(file_or_group = f,
                            names = ['id_particle', 'position_x', 'position_y', 'position_z', 'radius', 'drag_force_x', 'drag_force_y',
                                    'drag_force_z', 'y', 'porosity', 'slip_velocity_x', 'slip_velocity_y', 'slip_velocity_z'],
                            data = [self.id_particle, self.position_x, self.position_y, self.position_z, self.radius,
                                    self.drag_force_x, self.drag_force_y, self.drag_force_z, self.y, self.porosity,
                                    self.slip_velocity_x, self.slip_velocity_y, self.slip_velocity_z])

    def WriteDataToFile(self, file_or_group, names, data):
        if self.group_name in file_or_group:
            file_or_group['/'].__delitem__(self.group_name)
        self.sub_group = file_or_group.create_group(self.group_name)
        column_name = np.dtype({'names': names, 'formats':[(float)]*len(names)})
        data_array = np.rec.fromarrays(data, dtype = column_name)

        self.sub_group.create_dataset(name = self.group_name, data = data_array)
# Import system python
import os

import numpy as np
import h5py
import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.SwimmingDEMApplication as SDEM
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
        self.v_L2_error = []
        self.p_L2_error = []
        self.p_H1_error = []
        self.v_H1_error = []
        self.av_mod_error = []
        self.n_iterations = []
        self.problem_path = os.getcwd()
        self.file_path = os.path.join(str(self.problem_path),self.parameters["file_name"].GetString())
        self.dtype = np.float64
        self.group_name = str(test_number)

    def WriteData(self, error_model_part, velocity_L2_error_norm, pressure_L2_error_norm, velocity_H1_error_seminorm, pressure_H1_error_seminorm, projection_type, model_type, subscale_type, porosity_mean, n_iterations, max_iteration, relax_alpha, lowest_alpha, damkohler, omega,reynolds_number):
        self.error_model_part = error_model_part
        self.projection_type = projection_type
        self.model_type = model_type
        self.subscale_type = subscale_type
        self.porosity_mean = porosity_mean
        self.max_iteration = max_iteration
        self.relax_alpha = relax_alpha
        self.lowest_alpha = lowest_alpha
        self.damkohler_number = damkohler
        self.omega = omega

        element_size = [Element.GetGeometry().Length() for Element in self.error_model_part.Elements]
        self.element_size = np.max(element_size)

        self.reynolds_number = reynolds_number

        self.n_iterations.append(n_iterations)
        self.time.append(self.error_model_part.ProcessInfo[Kratos.TIME])
        self.v_L2_error.append(velocity_L2_error_norm)
        self.p_L2_error.append(pressure_L2_error_norm)
        self.v_H1_error.append(velocity_H1_error_seminorm)
        self.p_H1_error.append(pressure_H1_error_seminorm)

        with h5py.File(self.file_path, 'a') as f:
                self.WriteDataToFile(file_or_group = f,
                            names = ['TIME', 'V_L2_ERROR', 'P_L2_ERROR', 'V_H1_ERROR', 'P_H1_ERROR', 'N_ITERATIONS'],
                            data = [self.time, self.v_L2_error, self.p_L2_error, self.v_H1_error, self.p_H1_error, self.n_iterations])

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
        self.sub_group.attrs['porosity_mean'] = str(self.porosity_mean)
        self.sub_group.attrs['max_iteration'] = str(self.max_iteration)
        self.sub_group.attrs['relaxation_alpha'] = str(self.relax_alpha)
        self.sub_group.attrs['lowest_alpha'] = str(self.lowest_alpha)
        self.sub_group.attrs['damkohler_number'] = str(self.damkohler_number)
        self.sub_group.attrs['omega'] = str(self.omega)
        self.sub_group.attrs['reynolds_number'] = str(self.reynolds_number)

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

class ForceBalance(object):
    def __init__(self):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.parameters = Kratos.Parameters( """
        {
            "file_name": "sp_data",
            "target_porosity" : 0.3,
            "probe_height" : 0.2032
        }  """ )


        self.problem_path = os.getcwd()
        self.dtype = np.float64
        self.hydrodynamic_force_x = []
        self.hydrodynamic_force_y = []
        self.hydrodynamic_force_z = []
        self.hydrodynamic_reaction_x = []
        self.hydrodynamic_reaction_y = []
        self.hydrodynamic_reaction_z = []
        self.hydrodynamic_reaction_x = []
        self.hydrodynamic_reaction_y = []
        self.hydrodynamic_reaction_z = []
        self.projected_hydrodynamic_reaction_x = []
        self.projected_hydrodynamic_reaction_y = []
        self.projected_hydrodynamic_reaction_z = []
        self.error = []
        self.time = []
        self.velocity = []

    def WriteData(self, spheres_model_part, fluid_model_part, time):

        self.delta_time = round(fluid_model_part.ProcessInfo.GetValue(Kratos.DELTA_TIME),3)

        element_size = [Element.GetGeometry().Length() for Element in fluid_model_part.Elements]
        self.element_size = np.max(element_size)
        self.group_name = str(self.element_size)

        self.file_path = os.path.join(str(self.problem_path),self.parameters["file_name"].GetString() + "_" + str(self.delta_time) + ".hdf5")

        hydrodynamic_force_x = hydrodynamic_force_y = hydrodynamic_force_z = hydrodynamic_reaction_x = hydrodynamic_reaction_y = hydrodynamic_reaction_z = 0.0
        projected_hydrodynamic_reaction_x = projected_hydrodynamic_reaction_y = projected_hydrodynamic_reaction_z = velocity = 0.0

        for node in spheres_model_part.Nodes:
            hydrodynamic_force_x += node.GetSolutionStepValue(Kratos.HYDRODYNAMIC_FORCE_X)
            hydrodynamic_force_y += node.GetSolutionStepValue(Kratos.HYDRODYNAMIC_FORCE_Y)
            hydrodynamic_force_z += node.GetSolutionStepValue(Kratos.HYDRODYNAMIC_FORCE_Z)
            projected_hydrodynamic_reaction_x += node.GetSolutionStepValue(SDEM.HYDRODYNAMIC_REACTION_PROJECTED_X)
            projected_hydrodynamic_reaction_y += node.GetSolutionStepValue(SDEM.HYDRODYNAMIC_REACTION_PROJECTED_Y)
            projected_hydrodynamic_reaction_z += node.GetSolutionStepValue(SDEM.HYDRODYNAMIC_REACTION_PROJECTED_Z)
            velocity += node.GetSolutionStepValue(Kratos.VELOCITY_Z)

        for node in fluid_model_part.Nodes:
            hydrodynamic_reaction_x += node.GetSolutionStepValue(Kratos.HYDRODYNAMIC_REACTION_X)
            hydrodynamic_reaction_y += node.GetSolutionStepValue(Kratos.HYDRODYNAMIC_REACTION_Y)
            hydrodynamic_reaction_z += node.GetSolutionStepValue(Kratos.HYDRODYNAMIC_REACTION_Z)

        error = np.abs(np.sqrt(hydrodynamic_force_x**2 + hydrodynamic_force_y**2 + hydrodynamic_force_z**2)-np.sqrt(hydrodynamic_reaction_x**2 + hydrodynamic_reaction_y**2 + hydrodynamic_reaction_z**2))/np.sqrt(hydrodynamic_force_x**2 + hydrodynamic_force_y**2 + hydrodynamic_force_z**2)

        self.time.append(time)
        self.hydrodynamic_force_x.append(hydrodynamic_force_x)
        self.hydrodynamic_force_y.append(hydrodynamic_force_y)
        self.hydrodynamic_force_z.append(hydrodynamic_force_z)
        self.hydrodynamic_reaction_x.append(hydrodynamic_reaction_x)
        self.hydrodynamic_reaction_y.append(hydrodynamic_reaction_y)
        self.hydrodynamic_reaction_z.append(hydrodynamic_reaction_z)
        self.projected_hydrodynamic_reaction_x.append(projected_hydrodynamic_reaction_x)
        self.projected_hydrodynamic_reaction_y.append(projected_hydrodynamic_reaction_y)
        self.projected_hydrodynamic_reaction_z.append(projected_hydrodynamic_reaction_z)
        self.error.append(error)
        self.velocity.append(velocity)

        with h5py.File(self.file_path, 'a') as f:
                self.WriteDataToFile(file_or_group = f,
                            names = ['HYDRODYNAMIC_FORCE_X','HYDRODYNAMIC_FORCE_Y','HYDRODYNAMIC_FORCE_Z','HYDRODYNAMIC_REACTION_X','HYDRODYNAMIC_REACTION_Y', 'HYDRODYNAMIC_REACTION_Z', 'HYDRODYNAMIC_REACTION_PROJECTED_X', 'HYDRODYNAMIC_REACTION_PROJECTED_Y', 'HYDRODYNAMIC_REACTION_PROJECTED_Z', 'ERROR', 'VELOCITY', 'TIME'],
                            data = [self.hydrodynamic_force_x ,self.hydrodynamic_force_y, self.hydrodynamic_force_z, self.hydrodynamic_reaction_x, self.hydrodynamic_reaction_y, self.hydrodynamic_reaction_z, self.projected_hydrodynamic_reaction_x, self.projected_hydrodynamic_reaction_y, self.projected_hydrodynamic_reaction_z, self.error, self.velocity, self.time])

    def WriteDataToFile(self, file_or_group, names, data):
        if self.group_name in file_or_group:
            file_or_group['/'].__delitem__(self.group_name)
            self.sub_group = file_or_group.create_group(self.group_name)
        else:
            self.sub_group = file_or_group.create_group(self.group_name)

        self.sub_group.attrs['element_size'] = str(self.element_size)
        self.sub_group.attrs['delta_time'] = str(self.delta_time)

        for name, datum in zip(names, data):
            if name in file_or_group:
                file_or_group.__delitem__(name)
        for name, datum in zip(names, data):
            self.sub_group.create_dataset(name = name, data = datum)

class Pressure(object):
    def __init__(self,test_number):
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
        self.dtype = np.float64
        self.pressure = []
        self.velocity_x = []
        self.velocity_y = []
        self.time = []
        self.file_path = os.path.join(str(self.problem_path),self.parameters["file_name"].GetString())
        self.group_name = str(test_number)

    def WriteData(self, fluid_model_part, model_type, alpha_min, reynolds_number):

        self.model_type = model_type
        self.alpha_min = alpha_min
        self.reynolds_number = reynolds_number
        time = fluid_model_part.ProcessInfo[Kratos.TIME]

        for node in fluid_model_part.Nodes:
            if node.X == 6.15 and node.Y == 4.0:
                pressure = node.GetSolutionStepValue(Kratos.PRESSURE)
                vel_x = node.GetSolutionStepValue(Kratos.VELOCITY_X)
                vel_y = node.GetSolutionStepValue(Kratos.VELOCITY_Y)
                break

        self.pressure.append(pressure)
        self.velocity_x.append(vel_x)
        self.velocity_y.append(vel_y)
        self.time.append(time)

        with h5py.File(self.file_path, 'a') as f:
                self.WriteDataToFile(file_or_group = f,
                            names = ['PRESSURE', 'VELOCITY_X', 'VELOCITY_Y', 'TIME'],
                            data = [self.pressure, self.velocity_x, self.velocity_y, self.time])
    def WriteDataToFile(self, file_or_group, names, data):
        if self.group_name in file_or_group:
            file_or_group['/'].__delitem__(self.group_name)
            self.sub_group = file_or_group.create_group(self.group_name)
        else:
            self.sub_group = file_or_group.create_group(self.group_name)

        self.sub_group.attrs['model_type'] = str(self.model_type)
        self.sub_group.attrs['alpha_min'] = str(self.alpha_min)
        self.sub_group.attrs['Reynolds_number'] = str(self.reynolds_number)

        for name, datum in zip(names, data):
            if name in file_or_group:
                file_or_group.__delitem__(name)
        for name, datum in zip(names, data):
            self.sub_group.create_dataset(name = name, data = datum)

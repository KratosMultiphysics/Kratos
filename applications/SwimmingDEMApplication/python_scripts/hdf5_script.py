# Import system python
import os

import numpy as np
import h5py
import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.FluidDynamicsApplication as Fluid
from KratosMultiphysics.SwimmingDEMApplication import *

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

    def WriteData(self, error_model_part, velocity_L2_error_norm, pressure_L2_error_norm, velocity_H1_error_seminorm, pressure_H1_error_seminorm, projection_type, model_type, subscale_type, porosity_mean, n_iterations, max_iteration, relax_alpha, lowest_alpha, damkohler, omega,reynolds_number, viscous_projector):
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
        if viscous_projector == 0:
            self.viscous_projector = 'Identity'
        elif viscous_projector == 1:
            self.viscous_projector = 'Symmetric'
        elif viscous_projector == 2:
            self.viscous_projector = 'SymmetricAndDeviatoric'

        elem_size = [Element.GetGeometry().Length() for Element in self.error_model_part.Elements]

        self.element_size = np.max(elem_size)

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
        self.sub_group.attrs['viscous_projector'] = str(self.viscous_projector)

        for name, datum in zip(names, data):
            if name in file_or_group:
                file_or_group.__delitem__(name)
        for name, datum in zip(names, data):
            self.sub_group.create_dataset(name = name, data = datum)

class NodalErrorPostProcessTool(object):
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
        self.n_iterations = []
        self.problem_path = os.getcwd()
        self.file_path = os.path.join(str(self.problem_path),self.parameters["file_name"].GetString())
        self.dtype = np.float64
        self.group_name = str(test_number)

    def WriteData(self, fluid_model_part, reynolds_number, model_type, n_iterations, max_iteration):
        self.fluid_model_part = fluid_model_part
        self.model_type = model_type
        nodal_error = total_area = 0.0
        for node in fluid_model_part.Nodes:
            nodal_error += node.GetSolutionStepValue(Kratos.NODAL_AREA) * (node.GetSolutionStepValue(EXACT_VELOCITY_X)-node.GetSolutionStepValue(Kratos.VELOCITY_X))**2+(node.GetSolutionStepValue(EXACT_VELOCITY_Y)-node.GetSolutionStepValue(Kratos.VELOCITY_Y))**2+(node.GetSolutionStepValue(EXACT_VELOCITY_Z)-node.GetSolutionStepValue(Kratos.VELOCITY_Z))**2
            total_area += node.GetSolutionStepValue(Kratos.NODAL_AREA)

        nodal_error = nodal_error/total_area

        elem_size = [Element.GetGeometry().Length() for Element in self.fluid_model_part.Elements]

        self.element_size = np.max(elem_size)
        self.max_iteration = max_iteration
        self.reynolds_number = reynolds_number

        self.n_iterations.append(n_iterations)
        self.time.append(self.fluid_model_part.ProcessInfo[Kratos.TIME])
        self.v_L2_error.append(nodal_error)

        with h5py.File(self.file_path, 'a') as f:
                self.WriteDataToFile(file_or_group = f,
                            names = ['TIME', 'V_L2_ERROR', 'N_ITERATIONS'],
                            data = [self.time, self.v_L2_error, self.n_iterations])

    def WriteDataToFile(self, file_or_group, names, data):
        if self.group_name in file_or_group:
            file_or_group['/'].__delitem__(self.group_name)
            self.sub_group = file_or_group.create_group(self.group_name)
        else:
            self.sub_group = file_or_group.create_group(self.group_name)
        self.sub_group.attrs['element_size'] = str(self.element_size)
        self.sub_group.attrs['n_elements'] = str(len(self.fluid_model_part.Elements))
        self.sub_group.attrs['delta_time'] = str(self.fluid_model_part.ProcessInfo[Kratos.DELTA_TIME])
        self.sub_group.attrs['model_type'] = str(self.model_type)
        self.sub_group.attrs['max_iteration'] = str(self.max_iteration)
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

class EnergyAnalytics(object):
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
        self.fluid_kinetic_energy_x = []
        self.fluid_kinetic_energy_y = []
        self.fluid_kinetic_energy_z = []
        self.particle_kinetic_energy_x = []
        self.particle_kinetic_energy_z = []
        self.particle_total_kinetic_energy = []
        self.front_position = []
        self.time = []
        self.group_name = str(1)

    def WriteData(self, spheres_model_part, fluid_model_part, time):
        self.spheres_model_part = spheres_model_part
        self.fluid_model_part = fluid_model_part

        fke_x = fke_y = fke_z = pke_x = pke_z = t_pke = 0.0
        max_node_x = nodal_area_upper = nodal_area_bottom = 0.0

        solid_fraction_cont = solid_fraction_disc = 0.0


        for node in self.fluid_model_part.Nodes:
            fke_x += 0.5*node.GetSolutionStepValue(Kratos.NODAL_AREA)*node.GetSolutionStepValue(Kratos.FLUID_FRACTION)*node.GetSolutionStepValue(Kratos.DENSITY)*node.GetSolutionStepValue(Kratos.VELOCITY_X)**2
            fke_y += 0.5*node.GetSolutionStepValue(Kratos.NODAL_AREA)*node.GetSolutionStepValue(Kratos.FLUID_FRACTION)*node.GetSolutionStepValue(Kratos.DENSITY)*node.GetSolutionStepValue(Kratos.VELOCITY_Y)**2
            fke_z += 0.5*node.GetSolutionStepValue(Kratos.NODAL_AREA)*node.GetSolutionStepValue(Kratos.FLUID_FRACTION)*node.GetSolutionStepValue(Kratos.DENSITY)*node.GetSolutionStepValue(Kratos.VELOCITY_Z)**2
            solid_fraction_cont += node.GetSolutionStepValue(Kratos.NODAL_AREA) - node.GetSolutionStepValue(Kratos.NODAL_AREA) * node.GetSolutionStepValue(Kratos.FLUID_FRACTION)

        for node in self.spheres_model_part.Nodes:
            pke_x += 0.5*self.spheres_model_part.GetProperties()[1][DEM.PARTICLE_DENSITY]*4/3*np.pi*node.GetSolutionStepValue(Kratos.RADIUS)**3*node.GetSolutionStepValue(Kratos.VELOCITY_X)**2
            pke_z += 0.5*self.spheres_model_part.GetProperties()[1][DEM.PARTICLE_DENSITY]*4/3*np.pi*node.GetSolutionStepValue(Kratos.RADIUS)**3*node.GetSolutionStepValue(Kratos.VELOCITY_Z)**2
            t_pke += 0.5*self.spheres_model_part.GetProperties()[1][DEM.PARTICLE_DENSITY]*4/3*np.pi*node.GetSolutionStepValue(Kratos.RADIUS)**3*np.sqrt(node.GetSolutionStepValue(Kratos.VELOCITY_X)**2 + node.GetSolutionStepValue(Kratos.VELOCITY_Z)**2 + node.GetSolutionStepValue(Kratos.VELOCITY_Y)**2)**2
            if node.X >= max_node_x and node.Z != 0.0025:
                max_node_x = node.X
            solid_fraction_disc += 4/3*np.pi*node.GetSolutionStepValue(Kratos.RADIUS)**3

        print('projected_volume', solid_fraction_cont)
        print('real_volume', solid_fraction_disc)
        print('difference', solid_fraction_cont - solid_fraction_disc)
        self.time.append(time)
        self.particle_kinetic_energy_x.append(pke_x)
        self.particle_kinetic_energy_z.append(pke_z)
        self.fluid_kinetic_energy_x.append(fke_x)
        self.fluid_kinetic_energy_y.append(fke_y)
        self.fluid_kinetic_energy_z.append(fke_z)
        self.particle_total_kinetic_energy.append(t_pke)
        self.front_position.append(max_node_x)

        with h5py.File(self.file_path, 'a') as f:
                self.WriteDataToFile(file_or_group = f,
                            names = ['FLUID_KINETIC_ENERGY_X', 'FLUID_KINETIC_ENERGY_Y','FLUID_KINETIC_ENERGY_Z','PARTICLES_KINETIC_ENERGY_X','PARTICLES_KINETIC_ENERGY_Z', 'PARTICLES_TOTAL_KINETIC_ENERGY', 'FRONT_POSITION', 'TIME'],
                            data = [self.fluid_kinetic_energy_x, self.fluid_kinetic_energy_y, self.fluid_kinetic_energy_z, self.particle_kinetic_energy_x, self.particle_kinetic_energy_z, self.particle_total_kinetic_energy, self.front_position, self.time])

    def WriteDataToFile(self, file_or_group, names, data):
        if self.group_name in file_or_group:
            file_or_group['/'].__delitem__(self.group_name)
            self.sub_group = file_or_group.create_group(self.group_name)
        else:
            self.sub_group = file_or_group.create_group(self.group_name)

        for name, datum in zip(names, data):
            if name in file_or_group:
                file_or_group.__delitem__(name)
        for name, datum in zip(names, data):
            self.sub_group.create_dataset(name = name, data = datum)

class KineticEnergy(object):
    def __init__(self, test_number):
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
        self.fluid_kinetic_energy = []
        self.time = []
        self.group_name = str(test_number)

    def WriteData(self, fluid_model_part, model_type, amean, time):
        self.fluid_model_part = fluid_model_part
        self.model_type = model_type
        self.alpha_mean = amean

        fke = total_volume = 0.0

        for node in self.fluid_model_part.Nodes:
            fke += 0.5*node.GetSolutionStepValue(Kratos.NODAL_AREA)*node.GetSolutionStepValue(Kratos.FLUID_FRACTION)*node.GetSolutionStepValue(Kratos.DENSITY)*(node.GetSolutionStepValue(Kratos.VELOCITY_X)**2+node.GetSolutionStepValue(Kratos.VELOCITY_Y)**2 + node.GetSolutionStepValue(Kratos.VELOCITY_Z)**2)
            total_volume += node.GetSolutionStepValue(Kratos.NODAL_AREA)*node.GetSolutionStepValue(Kratos.FLUID_FRACTION)

        self.time.append(time)
        self.fluid_kinetic_energy.append(fke/total_volume)

        with h5py.File(self.file_path, 'a') as f:
                self.WriteDataToFile(file_or_group = f,
                            names = ['FLUID_KINETIC_ENERGY', 'TIME'],
                            data = [self.fluid_kinetic_energy, self.time])

    def WriteDataToFile(self, file_or_group, names, data):
        if self.group_name in file_or_group:
            file_or_group['/'].__delitem__(self.group_name)
            self.sub_group = file_or_group.create_group(self.group_name)
        else:
            self.sub_group = file_or_group.create_group(self.group_name)

        self.sub_group.attrs['model_type'] = str(self.model_type)
        self.sub_group.attrs['alpha_mean'] = str(self.alpha_mean)

        for name, datum in zip(names, data):
            if name in file_or_group:
                file_or_group.__delitem__(name)
        for name, datum in zip(names, data):
            self.sub_group.create_dataset(name = name, data = datum)

class AveragingVariablesPostProcessTool(object):

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
        self.dtype = np.float64
        self.file_path = os.path.join(str(self.problem_path),self.parameters["file_name"].GetString())
        n_layers = 100
        self.velocity_layers = dict(("average_velocity_layer_" + str(i), []) for i in range(1, n_layers + 1))
        self.velocity_gradient_z_layers = dict(("average_velocity_gradient_z_layer_" + str(i), []) for i in range(1, n_layers + 1))
        self.granular_temperature_layers = dict(("granular_temperature_layer_" + str(i), []) for i in range(1, n_layers + 1))
        self.averaged_stress_layers = dict(("averaged_stress_layer_" + str(i), []) for i in range(1, n_layers + 1))
        self.averaged_packing_density_layers = dict(("averaged_packing_density_" + str(i), []) for i in range(1, n_layers + 1))
        self.time = []

        self.average_utility = AveragingVariablesUtility()

    def WriteData(self, spheres_model_part):
        time = spheres_model_part.ProcessInfo[Kratos.TIME]
        n_layers = 100
        n_sublayers = 11
        layer_width = 0.0004
        bottom = 10*layer_width
        plane_area = 40*layer_width*32*layer_width
        self.time.append(time)

        for i in range(1, n_layers + 1):
            centroid_layer_i = bottom + (i-0.5)*layer_width
            layer_averaged_velocity = Kratos.Vector(3)
            layer_averaged_dv_dz = Kratos.Vector(3)
            layer_averaged_particle_stress = Kratos.Matrix(3,3)
            layer_granular_temperature = Kratos.Matrix(3,3)
            layer_packing_density = self.average_utility.AverageVariables(spheres_model_part, centroid_layer_i, layer_width, n_sublayers, plane_area, layer_averaged_velocity, layer_averaged_dv_dz, layer_averaged_particle_stress, layer_granular_temperature)
            kinetic_stress = -spheres_model_part.GetProperties()[1][DEM.PARTICLE_DENSITY]*layer_packing_density*np.array(layer_granular_temperature)
            layer_averaged_stress = np.array(layer_averaged_particle_stress) + kinetic_stress
            self.velocity_layers["average_velocity_layer_" + str(i)].append(layer_averaged_velocity)
            self.velocity_gradient_z_layers["average_velocity_gradient_z_layer_" + str(i)].append(layer_averaged_dv_dz)
            self.granular_temperature_layers["granular_temperature_layer_" + str(i)].append(np.array(layer_granular_temperature).flatten())
            self.averaged_stress_layers["averaged_stress_layer_" + str(i)].append(layer_averaged_stress.flatten())
            self.averaged_packing_density_layers["averaged_packing_density_" + str(i)].append(layer_packing_density)
            self.group_name = str(i)
            with h5py.File(self.file_path, 'a') as f:
                if self.group_name in f:
                    f['/'].__delitem__(self.group_name)
                    self.sub_group = f.create_group(self.group_name)
                else:
                    self.sub_group = f.create_group(self.group_name)

                self.AverageVelocityToFile(['AVERAGED_VELOCITY'], self.velocity_layers["average_velocity_layer_" + str(i)])
                self.AverageGradientZVelocityToFile(['AVERAGED_GRADIENT_Z_VELOCITY'], self.velocity_gradient_z_layers["average_velocity_gradient_z_layer_" + str(i)])
                self.GranularTemperatureToFile(['GRANULAR_TEMPERATURE'], self.granular_temperature_layers["granular_temperature_layer_" + str(i)])
                self.StressTensorToFile(['STRESS'], self.averaged_stress_layers["averaged_stress_layer_" + str(i)])
                self.PackingDensityToFile(['PACKING_DENSITY'],self.averaged_packing_density_layers["averaged_packing_density_" + str(i)])
                self.TimeToFile(['TIME'], self.time)

    def AverageVelocityToFile(self, group_name, data):
        column_names = ['U_x', 'U_y', 'U_z']
        self.WriteTensorDataToFile(group_name, column_names, data)

    def AverageGradientZVelocityToFile(self, group_name, data):
        column_names = ['DU_x', 'DU_y', 'DU_z']
        self.WriteTensorDataToFile(group_name, column_names, data)

    def GranularTemperatureToFile(self, group_name, data):
        column_names = ['T_xx', 'T_xy', 'T_xz', 'T_yx', 'T_yy', 'T_yz', 'T_zx', 'T_zy', 'T_zz']
        self.WriteTensorDataToFile(group_name, column_names, data)

    def StressTensorToFile(self, group_name, data):
        column_names = ['S_xx', 'S_xy', 'S_xz', 'S_yx', 'S_yy', 'S_yz', 'S_zx', 'S_zy', 'S_zz']
        self.WriteTensorDataToFile(group_name, column_names, data)

    def PackingDensityToFile(self, group_name, data):
        self.WriteScalarDataToFile(group_name, data)

    def TimeToFile(self, group_name, data):
        self.WriteScalarDataToFile(group_name, data)

    def WriteTensorDataToFile(self, group_name, column_names, data):
        column_name = np.dtype({'names': column_names, 'formats':[(float)]*len(column_names)})
        data_array = []
        for i in range(0,len(data)):
            data_array.append(np.rec.fromarrays(data[i], dtype = column_name))
        self.sub_group.create_dataset(name = group_name[0], data = data_array)

    def WriteScalarDataToFile(self, group_name, data):
        data_array = []
        data_array.append(data)
        self.sub_group.create_dataset(name = group_name[0], data = data_array)

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
        self.sub_group.attrs['reynolds_number'] = str(self.reynolds_number)

        for name, datum in zip(names, data):
            if name in file_or_group:
                file_or_group.__delitem__(name)
        for name, datum in zip(names, data):
            self.sub_group.create_dataset(name = name, data = datum)
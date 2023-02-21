# Import system python
import os

import numpy as np
import h5py
import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as Dem
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

def Norm(vector):
    return np.sqrt(sum(v**2 for v in vector))

class SandProductionPostProcessTool():
    def __init__(self, fem_model_part, dem_model_part, test_number):

        if not test_number:
            return

        self.parameters = Kratos.Parameters( """
        {
            "file_name": "sp_data.hdf5",
            "target_porosity" : 0.3,
            "probe_height" : 0.2032
        }  """ )

        if test_number == 1: # CTW16
            test_id = "CTW16"
        elif test_number == 2: # CTW10
            test_id = "CTW10"
        else: # Blind test
            test_id = "BlindTest"

        self.parameters.AddEmptyValue("test_id")
        self.parameters["test_id"].SetString(test_id)

        self.dem_model_part = dem_model_part
        self.fem_model_part = fem_model_part

        self.post_process_utility = DemFem.PostProcessUtilities(self.dem_model_part)

        self.internal_radius, self.external_radius, self.interface_radius, self.thickness, self.volume, self.average_radius = self.GetProbeDimensions()
        for elem in self.dem_model_part.Elements:
            self.density = elem.Properties[Dem.PARTICLE_DENSITY]
            break
        self.porosity = self.CalculatePorosity()
        self.problem_path = os.getcwd()
        self.file_path = os.path.join(str(self.problem_path),self.parameters["file_name"].GetString())
        self.target_porosity = self.parameters["target_porosity"].GetDouble()
        self.real_probe_height = self.parameters["probe_height"].GetDouble()

        self.dtype = np.float64
        self.compression_type = 'gzip'
        self.number_of_readings = 0
        self.last_time = - float('inf')

        # This assumes annulus centered at the origin
        with h5py.File(self.file_path, 'w') as f:
            f.attrs['test_id'] = self.parameters["test_id"].GetString()
            f.attrs['internal_radius'] = self.internal_radius
            f.attrs['external_radius'] = self.external_radius
            f.attrs['interface_radius'] = self.interface_radius
            f.attrs['thickness'] = self.thickness
            f.attrs['volume'] = self.volume
            f.attrs['real_probe_height'] = self.real_probe_height
            f.attrs['target_porosity'] = self.target_porosity
            f.attrs['porosity'] = self.porosity
            f.attrs['density'] = self.density

    def CalculateAverageRadius(self):
        radii_array = np.array([node.GetSolutionStepValue(Kratos.RADIUS) for node in self.dem_model_part.Nodes])
        return sum(radii_array)/len(radii_array)

    def GetProbeDimensions(self):
        average_radius = self.CalculateAverageRadius()
        internal_radius = min(Norm([node.X, node.Y]) for node in self.dem_model_part.Nodes) - average_radius
        external_radius = max(Norm([node.X, node.Y]) for node in self.fem_model_part.Nodes)
        interface_radius = max(Norm([node.X, node.Y]) for node in self.dem_model_part.Nodes) + average_radius
        min_z = min(node.Z for node in self.fem_model_part.Nodes)
        max_z = max(node.Z for node in self.fem_model_part.Nodes)
        thickness = max_z - min_z
        area = np.pi * (external_radius**2 - internal_radius**2)
        volume = thickness * area
        return internal_radius, external_radius, interface_radius, thickness, volume, average_radius

    def CalculatePorosity(self):
        volume = 0.
        for node in self.dem_model_part.Nodes:
            volume += 4/3 * np.pi * node.GetSolutionStepValue(Kratos.RADIUS) ** 3
        volume_dem_part = np.pi * self.thickness * (self.interface_radius**2 - self.internal_radius**2)
        porosity = volume / volume_dem_part
        return porosity

    def WriteData(self):
        time = self.dem_model_part.ProcessInfo[Kratos.TIME]
        name = str(self.number_of_readings)
        with h5py.File(self.file_path, 'r+') as f:
            f.create_group(name = name)
            f[name].attrs['time'] = time
        column_shape = (len(self.dem_model_part.Nodes), )
        self.ids_array = np.array([node.Id for node in self.dem_model_part.Nodes])
        self.radii_array = np.array([node.GetSolutionStepValue(Kratos.RADIUS) for node in self.dem_model_part.Nodes])
        self.x_array = np.array([node.X for node in self.dem_model_part.Nodes])
        self.y_array = np.array([node.Y for node in self.dem_model_part.Nodes])
        self.z_array = np.array([node.Z for node in self.dem_model_part.Nodes])
        is_sticky_list = []
        self.post_process_utility.GetStickyStatus(is_sticky_list)
        self.is_sticky_array = np.array(is_sticky_list)
        initial_continuum_bonds_list = []
        self.post_process_utility.GetInitialContinuumBonds(initial_continuum_bonds_list)
        self.initial_continuum_bonds_array = np.array(initial_continuum_bonds_list)
        current_continuum_bonds_list = []
        self.post_process_utility.GetCurrentContinuumBonds(current_continuum_bonds_list)
        self.current_continuum_bonds_array = np.array(current_continuum_bonds_list)

        if not self.last_time == time:
            with h5py.File(self.file_path, 'r+') as f:
                f.create_dataset(name + '/id', compression = self.compression_type, shape = column_shape, dtype = self.dtype) # TODO: change type as int
                f[name + '/id'][:] = self.ids_array[:]
                f.create_dataset(name + '/radius', compression = self.compression_type, shape = column_shape, dtype = self.dtype)
                f[name + '/radius'][:] = self.radii_array[:]
                f.create_dataset(name + '/x', compression = self.compression_type, shape = column_shape, dtype = self.dtype)
                f[name + '/x'][:] = self.x_array[:]
                f.create_dataset(name + '/y', compression = self.compression_type, shape = column_shape, dtype = self.dtype)
                f[name + '/y'][:] = self.y_array[:]
                f.create_dataset(name + '/z', compression = self.compression_type, shape = column_shape, dtype = self.dtype)
                f[name + '/z'][:] = self.z_array[:]
                f.create_dataset(name + '/is_sticky', compression = self.compression_type, shape = column_shape, dtype = self.dtype)
                f[name + '/is_sticky'][:] = self.is_sticky_array[:]
                f.create_dataset(name + '/initial_continuum_bonds', compression = self.compression_type, shape = column_shape, dtype = self.dtype)
                f[name + '/initial_continuum_bonds'][:] = self.initial_continuum_bonds_array[:]
                f.create_dataset(name + '/current_continuum_bonds', compression = self.compression_type, shape = column_shape, dtype = self.dtype)
                f[name + '/current_continuum_bonds'][:] = self.current_continuum_bonds_array[:]

        self.last_time = time
        self.number_of_readings += 1

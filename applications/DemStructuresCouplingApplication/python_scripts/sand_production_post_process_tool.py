import numpy as np
import h5py
import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as Dem
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem

def Norm(vector):
    return np.sqrt(sum(v**2 for v in vector))

class SandProductionPostProcessTool(object):
    def __init__(self, parameters, fem_model_part, dem_model_part):

        self.parameters = parameters
        self.dem_model_part = dem_model_part
        self.fem_model_part = fem_model_part
        self.internal_radius, self.external_radius, self.interface_radius, self.thickness, self.volume = self.GetProbeDimensions()
        for elem in self.dem_model_part.Elements:
            self.density = elem.Properties[Dem.PARTICLE_DENSITY]
            break
        # self.density = self.dem_model_part.Elements(0).GetProperties()[PARTICLE_DENSITY]
        self.porosity = self.CalculatePorosity()
        self.file_path = self.parameters["file_path"].GetString() + 'sp_data.hdf5'
        self.target_porosity = self.parameters["target_porosity"].GetDouble()
        self.real_probe_height = self.parameters["probe_height"].GetDouble()

        self.dtype = np.float64
        self.compression_type = 'gzip'
        self.number_of_readings = 0
        self.last_time = - float('inf')

        # This assumes annulus centered at the origin
        print(self.file_path)
        with h5py.File(self.file_path, 'w') as f:
            f.attrs['test_id'] = self.parameters["test_id"].GetString()
            f.attrs['internal_radius'] = self.internal_radius
            f.attrs['external_radius'] = self.external_radius
            f.attrs['thickness'] = self.thickness
            f.attrs['volume'] = self.volume
            f.attrs['real_probe_height'] = self.real_probe_height
            f.attrs['target_porosity'] = self.target_porosity
            f.attrs['porosity'] = self.porosity

    def GetProbeDimensions(self):
        internal_radius = min(Norm([node.X, node.Y, node.Z]) for node in self.dem_model_part.Nodes)
        external_radius = max(Norm([node.X, node.Y, node.Z]) for node in self.fem_model_part.Nodes)
        interface_radius = max(Norm([node.X, node.Y, node.Z]) for node in self.dem_model_part.Nodes)
        min_z = min(node.Z for node in self.fem_model_part.Nodes)
        max_z = max(node.Z for node in self.fem_model_part.Nodes)
        thickness = max_z - min_z
        area = np.pi * (external_radius**2 - internal_radius**2)
        volume = thickness * area
        return internal_radius, external_radius, interface_radius, thickness, volume

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
        shape = (len(self.dem_model_part.Nodes), 7)
        column_shape = (len(self.dem_model_part.Nodes), )
        self.ids_array = np.array([node.Id for node in self.dem_model_part.Nodes])
        self.radii_array = np.array([node.GetSolutionStepValue(Kratos.RADIUS) for node in self.dem_model_part.Nodes])

        if not self.last_time == time:
            with h5py.File(self.file_path, 'r+') as f:
                f.create_dataset(name + '/id', compression = self.compression_type, shape = column_shape, dtype = self.dtype)
                f[name + '/id'][:] = self.ids_array[:]
                f.create_dataset(name + '/radius', compression = self.compression_type, shape = column_shape, dtype = self.dtype)
                f[name + '/radius'][:] = self.radii_array[:]

        self.last_time = time
        self.number_of_readings += 1

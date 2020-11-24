'''An inlet boundary condition process for KratosMultiphysics

license: license.txt
'''

__all__ = ['Factory', 'ImposeWindInletProcess']

from collections import namedtuple
from collections import Mapping
#from math import isclose
def isclose(a, b, rel_tol=1e-9, abs_tol=0.):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

from math import floor
from math import ceil
from math import log

import numpy as np
from scipy.special import comb
import matplotlib.pyplot as plt

import h5py
from pyevtk.hl import imageToVTK

import KratosMultiphysics
from KratosMultiphysics import TIME
from KratosMultiphysics import DELTA_TIME
from KratosMultiphysics import VELOCITY_X
from KratosMultiphysics import VELOCITY_Y
from KratosMultiphysics import VELOCITY_Z
from KratosMultiphysics import Logger

from time import time

from GaussianRandomField import *
from CovarianceKernels import VonKarmanCovariance, MannCovariance
from GenerateWind import GenerateWind

##DONE: added support for power law
##DONE: read all mean profile parameters from json file
##DONE: can choose between power law and logarithmic profile

class Parameters(Mapping):

    def __init__(self, kratos_parameters):
        self._kratos_parameters = kratos_parameters

    def __getitem__(self, key):
        param = self._kratos_parameters[key]
        if param.IsDouble():
            value = param.GetDouble()
        elif param.IsInt():
            value = param.GetInt()
        elif param.IsString():
            value = param.GetString()
        else:
            value = param
        return value

    def __iter__(self):
        yield from self._kratos_parameters.keys()

    def __len__(self):
        return self._kratos_parameters.size()

        
def Factory(settings, Model):
    return ImposeWindInletProcess(Model, Parameters(settings['Parameters']))


Extent = namedtuple('Extent', ['lower', 'upper'])


class RegularGrid1D:

    @property
    def lower_bound(self):
        return self.extent.lower

    @property
    def upper_bound(self):
        return self.extent.upper

    def __init__(self, start_pos, length, size):
        self.extent = Extent(start_pos, start_pos+length)
        self.size = size
        self.step_size = (self.upper_bound-self.lower_bound) / (self.size-1)

    def __getitem__(self, index):
        return self.lower_bound + self.step_size*index

    def __len__(self):
        return self.size

    def floor_index(self, coord):
        if isclose(coord, self.lower_bound, abs_tol=1e-4):
            # Guard against negative index due to floating point representation.
            return 0
        local_coord = coord - self.lower_bound
        return int(floor(local_coord / self.step_size))

    def ceil_index(self, coord):
        if isclose(coord, self.upper_bound, abs_tol=1e-4):
            return len(self) - 1
        local_coord = coord - self.lower_bound
        return int(ceil(local_coord / self.step_size))

    def has(self, coord):
        return (self.floor_index(coord) >= 0
                and self.ceil_index(coord) < len(self))


IndexSpan = namedtuple('IndexSpan', ['begin', 'end'])


class SubGrid1D:

    @property
    def lower_bound(self):
        return self.grid[self.span.begin]

    @property
    def upper_bound(self):
        return self.grid[self.span.end]

    @property
    def step_size(self):
        return self.grid.step_size

    @property
    def gslice(self):
        return slice(self.span.begin, self.span.end+1)

    def __init__(self, grid, start_pos, end_pos):
        self.grid = grid
        start_index = grid.floor_index(start_pos)
        end_index = grid.ceil_index(end_pos)
        end_index = max(end_index, start_index+1)
        self.span = IndexSpan(start_index, end_index)
        
    def __getitem__(self, index):
        return self.grid[self.span.begin + index]

    def __len__(self):
        return self.span.end - self.span.begin + 1

    def has(self, coord):
        return (self.floor_index(coord) >= 0
                and self.ceil_index(coord) < len(self))

    def floor_index(self, coord):
        if isclose(coord, self.lower_bound, abs_tol=1e-4):
            return 0
        return self.grid.floor_index(coord) - self.span.begin

    def ceil_index(self, coord):
        if isclose(coord, self.upper_bound, abs_tol=1e-4):
            return len(self) - 1
        return self.grid.ceil_index(coord) - self.span.begin


Grid3D = namedtuple('Grid3D', ['x', 'y', 'z'])

class InletPanel3D:

    def __init__(self, grid, data):
        self.grid = grid
        self.data = data
        self.dy = self.grid.y.step_size
        self.dz = self.grid.z.step_size
        self.y0 = self.grid.y.lower_bound
        self.z0 = self.grid.z.lower_bound

    def update(self, pos):
        i = self.grid.x.floor_index(pos)
        tx = (pos-self.grid.x[i]) / self.grid.x.step_size
        data_0 = self.data[i, self.grid.y.gslice, self.grid.z.gslice]
        data_1 = self.data[i+1, self.grid.y.gslice, self.grid.z.gslice]
        self.cut_data = (1.0 - tx) * data_0 + tx * data_1

    def interpolate(self, node):
        # xi and eta are local mesh cell coordinates in the interval [-1, 1].
        eps = 1e-8 # to counteract roundoff error which appears when performing modular arithmetic with doubles
        xi = ((node.Y - self.y0 + eps) % self.dy - eps) / self.dy
        eta = ((node.Z - self.z0 + eps) % self.dz - eps) / self.dz
        xi = 2.0 * (xi-0.5)
        eta = 2.0 * (eta-0.5)
        # Interpolate using bilinear shape functions.
        weights = np.array([
            0.25 * (1.0-xi) * (1.0-eta),
            0.25 * (1.0+xi) * (1.0-eta),
            0.25 * (1.0+xi) * (1.0+eta),
            0.25 * (1.0-xi) * (1.0+eta)
        ])
        j = self.grid.y.floor_index(node.Y)
        k = self.grid.z.floor_index(node.Z)
        if max(abs(weights[1:])) < 1e-11:
            return (weights[0] * self.cut_data[j, k])
        elif max(abs(weights[1:3])) < 1e-11:
            return (weights[0] * self.cut_data[j, k] + weights[3] * self.cut_data[j, k+1])
        elif max(abs(weights[2:])) < 1e-11:
            return (weights[0] * self.cut_data[j, k] + weights[1] * self.cut_data[j+1, k])
        else:
            return (
                weights[0] * self.cut_data[j, k]
                + weights[1] * self.cut_data[j+1, k]
                + weights[2] * self.cut_data[j+1, k+1]
                + weights[3] * self.cut_data[j, k+1]
            )


def weak_min(nodes, key): return key(min(nodes, key=key))


def weak_max(nodes, key): return key(max(nodes, key=key))


def get_extent(nodes, key):
    lower_bound = weak_min(nodes, key=key)
    upper_bound = weak_max(nodes, key=key)
    return Extent(lower_bound, upper_bound)


class LogMeanProfile:

    def __init__(self, friction_velocity, roughness_height, bulk_wind_speed, dim):

        self.friction_velocity = friction_velocity
        self.roughness_height = roughness_height
        self.bulk_wind_speed = bulk_wind_speed
        self.dim = dim
    
    def get_height(self,node):
        if self.dim == 2: 
            height = node.Y
        elif self.dim == 3: 
            height = node.Z
        return height
            
    def wind_speed(self, node):
        return (self.friction_velocity / 0.41
                * log((self.get_height(node) + self.roughness_height) / self.roughness_height))

class ImposeWindInletProcess:

    @property
    def inlet_nodes(self):
        return self.model_part.Nodes

    def __init__(self, Model, settings):
        for name, value in settings.items():
            setattr(self, name, value)
        self.model_part = Model[self.inlet_model_part_name]
        self.x0 = 0.0

        if len(self.inlet_nodes) > 0:
            # define dimensions of the inlet domain
            y_extent = get_extent(self.inlet_nodes, key=lambda node: node.Y)
            z_extent = get_extent(self.inlet_nodes, key=lambda node: node.Z)
            setattr(self, 'ly', y_extent.upper - y_extent.lower)
            setattr(self, 'lz', z_extent.upper - z_extent.lower)
            self.mean_profile = self.CreateMeanProfile()
            # x_extent = Extent(0, self.time_interval_length * self.mean_profile.bulk_wind_speed)
            # setattr(self, 'lx', x_extent.upper - x_extent.lower)
            self.inlet_position = 0.0
            
            # define wind field and map
            grid_dimensions = [self.lx, self.ly, self.lz]
            self.wind = GenerateWind(self.friction_velocity, self.reference_height, grid_dimensions, self.wind_grid_levels, self.seed)
            self.block_num=0
            self.blend_region, self.mappers = self.Create3DMappers(self.wind)
            
            # total_wind = self.wind.total_wind
            # plt.imshow(total_wind[:,0,:,0])
            # plt.show()

    def ExecuteInitialize(self):
        for node in self.inlet_nodes:
            for var, _ in self.mappers:
                node.Fix(var)
        
    def ExecuteInitializeSolutionStep(self):
        self.UpdateInletPosition()
        Logger.PrintInfo('ImposeWindInletProcess',
                         'inlet position = %e' % self.inlet_position)
        self.AssignVelocity()
        self.ApplyRamp()
        
    def CreateMeanProfile(self, dim=3):
        reference_height = self.reference_height
        # lz = self.lz
        umean = self.umean
        roughness_height = self.roughness_height
        self.bulk_wind_speed = umean
        # self.bulk_wind_speed = umean * (log(lz/roughness_height) - 1.0) / log(reference_height/roughness_height)
        self.friction_velocity = umean * 0.41 / log((reference_height + roughness_height) / roughness_height)
        return LogMeanProfile(self.friction_velocity, roughness_height, self.bulk_wind_speed, dim)

    def smoothstep(self, x, x_min=0.0, x_max=1.0, N=3):
        x = np.clip((x - x_min) / (x_max - x_min), 0.0, 1.0)

        result = 0
        for n in range(0, N + 1):
            result += comb(N + n, n) * comb(2 * N + 1, N - n) * (-x) ** n

        result *= x ** (N + 1)

        return result

    def Create3DMappers(self, wind, old_blend_region=None):
        output_wind_field = False

        wind_field = wind()
        new_blend_region = wind.blend_region
        blend_num = wind.blend_num

        if old_blend_region is not None:
            
            # diff = wind_field[:blend_num,...] - old_blend_region.copy()
            # fig, ax = plt.subplots(nrows=1)
            # im = ax.imshow(diff[0,...,0])
            # fig.colorbar(im, fraction=.1, orientation='horizontal')
            # plt.show()

            for i in range(0,blend_num):
                factor = (i+1)/(blend_num+1)
                factor = self.smoothstep(factor)
                wind_field[i,...] = (1-factor)* old_blend_region[i,...] + factor*wind_field[i,...]
        
        if output_wind_field is True:
            self.wind_field = wind_field
            self.ExportToVTK()
            self.block_num += 1
        
        nx, ny, nz = wind_field[...,0].shape
        x_grid = RegularGrid1D(self.x0, self.lx, nx)
        y_grid = RegularGrid1D(self.y0, self.ly, ny)
        z_grid = RegularGrid1D(self.z0, self.lz, nz)
        y_extent = get_extent(self.inlet_nodes, key=lambda node: node.Y)
        z_extent = get_extent(self.inlet_nodes, key=lambda node: node.Z)
        y_subgrid = SubGrid1D(y_grid, y_extent.lower, y_extent.upper)
        z_subgrid = SubGrid1D(z_grid, z_extent.lower, z_extent.upper)
        grid = Grid3D(x_grid, y_subgrid, z_subgrid)
        mappers = ((VELOCITY_X, InletPanel3D(grid, wind_field[...,0])),
                   (VELOCITY_Y, InletPanel3D(grid, wind_field[...,1])),
                   (VELOCITY_Z, InletPanel3D(grid, wind_field[...,2])))
        return new_blend_region, mappers

    def UpdateInletPosition(self):
        dt = self.model_part.ProcessInfo[DELTA_TIME]
        self.inlet_position += dt * self.mean_profile.bulk_wind_speed

        #If the end of the time block is reached, generate a new wind block and reset map
        if self.inlet_position >= self.x0 + self.lx:
            self.x0 += self.lx
            self.blend_region, self.mappers = self.Create3DMappers(self.wind, self.blend_region)
            # total_wind = self.wind.total_wind
            # plt.imshow(total_wind[:,0,:,0])
            # plt.show()


    def AssignVelocity(self):
        for var, mapper in self.mappers:
            mapper.update(self.inlet_position)
            if var == VELOCITY_X:
                for node in self.inlet_nodes:
                    vel = self.mean_profile.wind_speed(node)
                    vel += mapper.interpolate(node)
                    node.SetSolutionStepValue(var, vel)
            else:
                for node in self.inlet_nodes:
                    vel = mapper.interpolate(node)
                    node.SetSolutionStepValue(var, vel)

    def ApplyRamp(self):
        time = self.model_part.ProcessInfo[TIME]
        if time < self.ramp_time:
            scal = time / self.ramp_time
            for node in self.inlet_nodes:
                for var, _ in self.mappers:
                    vel = node.GetSolutionStepValue(var)
                    node.SetSolutionStepValue(var, scal * vel)
    
    def Check(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass
    
    def ExecuteFinalizeSolutionStep(self):
        pass
    
    def ExecuteBeforeOutputStep(self):
        pass
    
    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass

    def ExportToVTK(self):
        FileName = 'WindField_' + str(self.block_num)
        h = 1/self.wind_field.shape[0]
        spacing = (h,h,h)

        wind_field_ = tuple([np.copy(self.wind_field[...,i], order='C') for i in range(3)])

        cellData = {'grid': np.zeros(self.wind_field[...,0].shape), 'velocity': wind_field_}
        imageToVTK(FileName, cellData = cellData, spacing=spacing)

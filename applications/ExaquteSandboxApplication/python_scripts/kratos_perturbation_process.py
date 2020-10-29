"""
An initial condition process for KratosMultiphysics
license: license.txt
"""

__all__ = ['Factory', 'ImposePerturbedInitialConditionProcess']

from collections import namedtuple
from collections import Mapping

from math import floor, ceil

import numpy as np

import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication
from  KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
import KratosMultiphysics.MappingApplication
import KratosMultiphysics.ExaquteSandboxApplication

from KratosMultiphysics.ExaquteSandboxApplication.GenerateCN import GenerateCN

def isclose(a, b, rel_tol=1e-9, abs_tol=0.):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def _DotProduct(A,B):
    return sum(i[0]*i[1] for i in zip(A, B))

Extent = namedtuple('Extent', ['lower', 'upper'])

IndexSpan = namedtuple('IndexSpan', ['begin', 'end'])

Grid3D = namedtuple('Grid3D', ['x', 'y', 'z'])
Grid2D = namedtuple('Grid2D', ['x', 'z'])


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
    return ImposePerturbedInitialConditionProcess(Model,Parameters(settings['Parameters']))


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


class DomainPanel3D:

    def __init__(self, grid, data):
        self.grid = grid
        self.data = data
        self.dx = self.grid.x.step_size
        self.x0 = self.grid.x.lower_bound
        self.dy = self.grid.y.step_size
        self.y0 = self.grid.y.lower_bound
        self.dz = self.grid.z.step_size
        self.z0 = self.grid.z.lower_bound

    def interpolate(self, node):
        # xi and gamma and eta are local mesh cell coordinates in the interval [-1, 1].
        xi = ((node.X - self.x0) % self.dx) / self.dx
        gamma = ((node.Y - self.y0) % self.dy) / self.dy
        eta = ((node.Z - self.z0) % self.dz) / self.dz
        xi = 2.0 * (xi-0.5)
        gamma = 2.0 * (gamma-0.5)
        eta = 2.0 * (eta-0.5)
        # Interpolate using trilinear shape functions.
        weights = (
            0.125 * (1.0-xi) * (1.0-gamma) * (1.0-eta), # - - -
            0.125 * (1.0+xi) * (1.0-gamma) * (1.0-eta), # + - -
            0.125 * (1.0+xi) * (1.0+gamma) * (1.0-eta), # + + -
            0.125 * (1.0+xi) * (1.0+gamma) * (1.0+eta), # + + +
            0.125 * (1.0+xi) * (1.0-gamma) * (1.0+eta), # + - +
            0.125 * (1.0-xi) * (1.0-gamma) * (1.0+eta), # - - +
            0.125 * (1.0-xi) * (1.0+gamma) * (1.0-eta), # - + -
            0.125 * (1.0-xi) * (1.0+gamma) * (1.0+eta)  # - + +
        )
        i = self.grid.x.floor_index(node.X)
        j = self.grid.y.floor_index(node.Y)
        k = self.grid.z.floor_index(node.Z)
        return (
            weights[0] * self.data[i,j,k]
            + weights[1] * self.data[i+1,j,k]
            + weights[2] * self.data[i+1,j+1,k]
            + weights[3] * self.data[i+1,j+1,k+1]
            + weights[4] * self.data[i+1,j,k+1]
            + weights[5] * self.data[i,j,k+1]
            + weights[6] * self.data[i,j+1,k]
            + weights[7] * self.data[i,j+1,k+1]
        )


class DomainPanel2D:

    def __init__(self, grid, data):
        self.grid = grid
        self.data = data
        self.dx = self.grid.x.step_size
        self.x0 = self.grid.x.lower_bound
        self.dz = self.grid.z.step_size
        self.z0 = self.grid.z.lower_bound

    def interpolate(self, node):
        # xi and eta are local mesh cell coordinates in the interval [-1, 1].
        xi = ((node.X - self.x0) % self.dx) / self.dx
        eta = ((node.Y - self.z0) % self.dz) / self.dz
        xi = 2.0 * (xi-0.5)
        eta = 2.0 * (eta-0.5)
        # Interpolate using bilinear shape functions.
        weights = (
            0.25 * (1.0-xi) * (1.0-eta), # - -
            0.25 * (1.0+xi) * (1.0-eta), # + -
            0.25 * (1.0+xi) * (1.0+eta), # + +
            0.25 * (1.0-xi) * (1.0+eta)  # - +
        )
        j = self.grid.x.floor_index(node.X)
        k = self.grid.z.floor_index(node.Y)
        return (
            weights[0] * self.data[j, k]
            + weights[1] * self.data[j+1, k]
            + weights[2] * self.data[j+1, k+1]
            + weights[3] * self.data[j, k+1]
        )


class ImposePerturbedInitialConditionAnalysisStage(ConvectionDiffusionAnalysis):
    """
    This analysis stage class solves the Poisson problem,
    required by the correlated initial condition to apply no penetrability into the bluff body.
    """
    def __init__(self,input_model,input_parameters):
        super(ImposePerturbedInitialConditionAnalysisStage,self).__init__(input_model,input_parameters)
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_LAPLACIAN) # for storing \nabla TEMPERATURE
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_LAPLACIAN_RATE) # for storing P(u_{cn})
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_X_GRADIENT) # for storing VELOCITY_X gradient
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_Y_GRADIENT) # for storing VELOCITY_Y gradient
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_Z_GRADIENT) # for storing VELOCITY_Z gradient

    def ApplyBoundaryConditions(self):
        super(ImposePerturbedInitialConditionAnalysisStage,self).ApplyBoundaryConditions()
        structure_model_part = self._GetSolver().main_model_part.GetSubModelPart(self.project_parameters["problem_data"]["structure_model_part"].GetString())
        KratosMultiphysics.MortarUtilities.ComputeNodesMeanNormalModelPart(structure_model_part, True)
        main_model_part_name = self.project_parameters["problem_data"]["model_part_name"].GetString()

        if self.project_parameters["problem_data"].Has("load_velocity_field"):
            self.LoadVelocityField()
        else: # required if a velocity field from file is not loaded
            for node in self.model.GetModelPart(main_model_part_name).Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_LAPLACIAN_RATE, node.GetSolutionStepValue(KratosMultiphysics.VELOCITY))

        # set forcing term
        self.ComputeVelocityGradients() # of c*P(u_{cn}) + u_T (which is stored in VELOCITY)
        for node in self.model.GetModelPart(main_model_part_name).Nodes:
            divergence = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X_GRADIENT)[0] + \
                node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y_GRADIENT)[1] + \
                node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z_GRADIENT)[2]
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,-1*divergence) # forcing term = - \nabla \cdot (c*P(u_{cn}) + u_T)

        # set boundary condition
        for node in structure_model_part.Nodes:
            perturbed_velocity = self.project_parameters["problem_data"]["penalty_coefficient"].GetDouble() * node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_LAPLACIAN_RATE) # c*P(u_{cn})
            normal = node.GetSolutionStepValue(KratosMultiphysics.NORMAL)
            computed_flux = - _DotProduct(normal,perturbed_velocity)
            node.SetSolutionStepValue(KratosMultiphysics.FACE_HEAT_FLUX,computed_flux) # boundary term = - c*P(u_{cn}) \cdot NORMAL

    def ComputeVelocityGradients(self):
        KratosMultiphysics.ComputeNodalGradientProcess(self._GetSolver().main_model_part,KratosMultiphysics.VELOCITY_X,KratosMultiphysics.VELOCITY_X_GRADIENT).Execute()
        KratosMultiphysics.ComputeNodalGradientProcess(self._GetSolver().main_model_part,KratosMultiphysics.VELOCITY_Y,KratosMultiphysics.VELOCITY_Y_GRADIENT).Execute()
        KratosMultiphysics.ComputeNodalGradientProcess(self._GetSolver().main_model_part,KratosMultiphysics.VELOCITY_Z,KratosMultiphysics.VELOCITY_Z_GRADIENT).Execute()

    def LoadVelocityField(self):
        main_model_part_name = self.project_parameters["problem_data"]["model_part_name"].GetString()
        main_model_part = self.model.GetModelPart(main_model_part_name)
        penalty_coeff = self.project_parameters["problem_data"]["penalty_coefficient"].GetDouble()
        with open(self.project_parameters["problem_data"]["load_velocity_field"].GetString()) as dat_file:
            lines=dat_file.readlines()
            for line, node in zip(lines,main_model_part.Nodes):
                velocity = KratosMultiphysics.Vector(3,0.0)
                velocity[0] = float(line.split(' ')[0])
                velocity[1] = float(line.split(' ')[1])
                velocity[2] = float(line.split(' ')[2])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_LAPLACIAN_RATE, node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)) # P(u_{cn}) stored in VELOCITY_LAPLACIAN_RATE
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, penalty_coeff*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY) + velocity) # c*P(u_{cn}) + u_T stored in VELOCITY

    def RestoreBoundaryValues(self):
        for submdpa in self.project_parameters["solver_settings"]["processes_sub_model_part_list"].GetStringArray():
            for node in self._GetSolver().main_model_part.GetSubModelPart(submdpa).Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_LAPLACIAN, -1*node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_LAPLACIAN_RATE)) # recall VELOCITY_LAPLACIAN_RATE is P(u_{cn}) (without penalty coefficient)

    def FinalizeSolutionStep(self):
        super(ImposePerturbedInitialConditionAnalysisStage,self).FinalizeSolutionStep()
        KratosMultiphysics.ComputeNodalGradientProcess(self._GetSolver().main_model_part,KratosMultiphysics.TEMPERATURE,KratosMultiphysics.VELOCITY_LAPLACIAN).Execute()
        self.RestoreBoundaryValues()
        # save new velocity
        penalty_coeff = self.project_parameters["problem_data"]["penalty_coefficient"].GetDouble()
        for node in self._GetSolver().main_model_part.Nodes:
            velocity_field_loaded_and_correlated_noise = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY) # u_T + c*P(u_{cn})
            gradient_solution = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_LAPLACIAN) # no-boundary: \nabla T ; boundary: - \nabla T
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,velocity_field_loaded_and_correlated_noise + penalty_coeff*gradient_solution)


class ImposePerturbedInitialConditionProcess(KratosMultiphysics.Process):
    """
    Process generating divergence-free correlated noise, imposing no-penetrability condition on bluff bodies,
    preserving boundary conditions on boundaries and (if required) adding the contribution of a manually-loaded velocity field.
    The generated initial field reads
        u_0 = u_T + c*P(u_{cn}) + c*\nabla potential ,
    where u_T is the velocity field we load, c a penalty coefficient, and P(u_{cn}) and \nabla potential are discussed next.
    @details The process first projects the correlated noise u_{cn} into the model part. The projected correlated noise is denoted as P(u_{cn}).
    Then, the associated poisson problem for the unknown TEMPERATURE is solved:
        - \nabla \cdot \nabla TEMPERATURE = \nabla \cdot (c*P(u_{cn}) + u_{T})  on \Omega
        c*P(u_{cn}) \cdot NORMAL = - \nabla TEMPERATURE \cdot NORMAL  on \partial \Omega_{bluff body}
        TEMPERATURE = 0  on \Omega_{outlet}
    Then, the tangential velocity is set to zero on the boundaries \partial \Omega_{bluff body} and \partial \Omega_{external domain} by doing:
        P(u_{cn}) + \nabla TEMPERATURE = 0 on \Gamma
    The correlated noise satisfying the no-penetrability and solenoidal conditions and respecting the boundary conditions is: VELOCITY = u_T + c*P(u_{cn}) + c*\nabla TEMPERATURE.

    Parameters:
    model    : Kratos Model
    settings : Kratos Parameters
    """

    @property
    def model_part_nodes(self):
        return self.model_part.Nodes

    def __init__(self,model,settings):
        super().__init__()
        for name, value in settings.items():
            setattr(self, name, value)
        self.model = model
        self.model_part = model.GetModelPart(self.model_part_name)
        if len(self.model_part_nodes) > 0:
            # Mappers are only valid during the lifetime of the file
            dim = len(self.correlated_noise_generator["grid_dimensions"].GetVector())
            if self.correlated_noise_generator.Has("grid_shape"):
                self.vector_field_generator = GenerateCN(
                    corrlen = self.correlated_noise_generator["correlated_length"].GetDouble(),
                    grid_dimensions = self.correlated_noise_generator["grid_dimensions"].GetVector(),
                    grid_shape = self.correlated_noise_generator["grid_shape"].GetInt(),
                    ndim=dim)
            else:
                self.vector_field_generator = GenerateCN(
                    corrlen = self.correlated_noise_generator["correlated_length"].GetDouble(),
                    grid_dimensions = self.correlated_noise_generator["grid_dimensions"].GetVector(),
                    ndim=dim)
            if dim == 3:
                self.mappers = self.Create3DMappers(self.vector_field_generator)
            else:
                self.mappers = self.Create2DMappers(self.vector_field_generator)

    def Create3DMappers(self, vector_field_generator):
        if hasattr(self, 'seed'):
            vf = vector_field_generator(seed=self.seed)
        else:
            vf = vector_field_generator()
        # normalize generated correlated noise
        vf = vf/np.linalg.norm(vf)
        nx, ny, nz = vf.shape[:-1]
        x_grid = RegularGrid1D(0., self.lx, nx)
        y_grid = RegularGrid1D(self.y0, self.ly, ny)
        z_grid = RegularGrid1D(self.z0, self.lz, nz)
        grid = Grid3D(x_grid, y_grid, z_grid)
        # mappers become invalid after file_ is destructed
        mappers = ((KratosMultiphysics.VELOCITY_X, DomainPanel3D(grid, vf[...,0])),
                   (KratosMultiphysics.VELOCITY_Y, DomainPanel3D(grid, vf[...,1])),
                   (KratosMultiphysics.VELOCITY_Z, DomainPanel3D(grid, vf[...,2])))
        return mappers

    def Create2DMappers(self, vector_field_generator):
        vf = vector_field_generator(seed=self.seed)
        # normalize generated correlated noise
        vf = vf/np.linalg.norm(vf)
        nx, nz = vf.shape[:-1]
        x_grid = RegularGrid1D(self.x0, self.lx, nx)
        z_grid = RegularGrid1D(self.z0, self.lz, nz)
        grid = Grid2D(x_grid, z_grid)
        # mappers become invalid after file_ is destructed.
        mappers = ((KratosMultiphysics.VELOCITY_X, DomainPanel2D(grid, vf[...,0])),
                   (KratosMultiphysics.VELOCITY_Y, DomainPanel2D(grid, vf[...,1])))
        return mappers

    def ExecuteInitialize(self):
        self.MapCorrelatedNoiseToModelPart()
        self.SolveAssociatedPoissonProblem()

    def MapCorrelatedNoiseToModelPart(self):
        for var, mapper in self.mappers:
            for node in self.model_part_nodes:
                vel = mapper.interpolate(node)
                node.SetSolutionStepValue(var, vel)

    def SolveAssociatedPoissonProblem(self):
        default_parameters = self.GetDefaultParametersAnalysisStage()
        self.poisson_parameters["problem_data"].ValidateAndAssignDefaults(default_parameters["problem_data"])
        self.poisson_parameters["solver_settings"].ValidateAndAssignDefaults(default_parameters["solver_settings"])
        self.poisson_parameters["processes"]["constraints_process_list"][0]["Parameters"].ValidateAndAssignDefaults(default_parameters["processes"]["constraints_process_list"][0]["Parameters"])
        model = KratosMultiphysics.Model()
        simulation = ImposePerturbedInitialConditionAnalysisStage(model,self.poisson_parameters)
        simulation.Initialize() # required before mapping
        mapper = KratosMultiphysics.MappingApplication.MapperFactory.CreateMapper(self.model.GetModelPart(self.poisson_parameters["problem_data"]["model_part_name"].GetString()),simulation._GetSolver().main_model_part,self.mapper_parameters)
        mapper.Map(KratosMultiphysics.VELOCITY,KratosMultiphysics.VELOCITY)
        simulation.RunSolutionLoop()
        simulation.Finalize()
        mapper.InverseMap(KratosMultiphysics.VELOCITY,KratosMultiphysics.VELOCITY)

    def GetDefaultParametersAnalysisStage(self):
        default_parameters = KratosMultiphysics.Parameters( """ {
            "problem_data"             : {
                "model_part_name" : "specify_model_part_name",
                "domain_size"     : -1,
                "parallel_type"   : "OpenMP",
                "time_step"       : 1.1,
                "start_time"      : 0.0,
                "end_time"        : 1.0,
                "echo_level"      : 0,
                "penalty_coefficient" : 1.0,
                "load_velocity_field" : "specify_velocity_field",
                "structure_model_part" : "specify_structure_model_part"
            },
            "solver_settings": {
                "model_part_name" : "specify_model_part_name",
                "domain_size" : -1,
                "solver_type" : "stationary",
                "echo_level": 0,
                "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "specify_mdpa_path"
                },
                "material_import_settings"           : {
                    "materials_filename" : "specify_poisson_materials_path"
                },
                "time_stepping" : {
                    "time_step": 1.1
                },
                "element_replace_settings" : {
                    "element_name" : "LaplacianElement",
                    "condition_name" : "ThermalFace"
                },
                "problem_domain_sub_model_part_list": [],
                "processes_sub_model_part_list": [],
                "auxiliary_variables_list" : []
            },
            "processes" : {
                "constraints_process_list" : [{
                    "python_module" : "assign_scalar_variable_process",
                    "kratos_module" : "KratosMultiphysics",
                    "Parameters"    : {
                        "model_part_name" : "specify_any_submodelpart",
                        "variable_name"   : "TEMPERATURE",
                        "constrained"     : true,
                        "value"           : 0.0,
                        "interval"        : [0.0,"End"]
                    }
                }]
            }
        } """ )
        return default_parameters

    def Check(self):
        pass

    def ExecuteInitializeSolutionStep(self):
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
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# -*- coding: utf-8 -*-
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
CheckForPreviousImport()


def AddVariables(model_part, config):

    
    thermal_settings = ConvectionDiffusionSettings()
    thermal_settings.SetUnknownVariable(globals()[config.unknown_variable])
    thermal_settings.SetDensityVariable(globals()[config.density_variable])
    thermal_settings.SetDiffusionVariable(globals()[config.diffusion_variable])
    #thermal_settings.SetVolumeSourceVariable(globals()[config.volume_source_variable])
    #thermal_settings.SetSurfaceSourceVariable(globals()[config.surface_source_variable])
    thermal_settings.SetMeshVelocityVariable(globals()[config.mesh_velocity_variable])
    #thermal_settings.SetProjectionVariable(globals()[config.projection_variable])
    thermal_settings.SetSpecificHeatVariable(globals()[config.specific_heat_variable])
    thermal_settings.SetVelocityVariable(globals()[config.velocity_variable])

    model_part.AddNodalSolutionStepVariable(thermal_settings.GetUnknownVariable())
    model_part.AddNodalSolutionStepVariable(thermal_settings.GetDensityVariable())
    model_part.AddNodalSolutionStepVariable(thermal_settings.GetDiffusionVariable())
    #model_part.AddNodalSolutionStepVariable(thermal_settings.GetVolumeSourceVariable())
    #model_part.AddNodalSolutionStepVariable(thermal_settings.GetSurfaceSourceVariable())
    model_part.AddNodalSolutionStepVariable(thermal_settings.GetMeshVelocityVariable());
    #model_part.AddNodalSolutionStepVariable(thermal_settings.GetProjectionVariable());
    model_part.AddNodalSolutionStepVariable(thermal_settings.GetSpecificHeatVariable());
    model_part.AddNodalSolutionStepVariable(thermal_settings.GetVelocityVariable());
    print("variables for the convection diffusion solver added correctly")


def AddDofs(model_part, config):
    thermal_settings = ConvectionDiffusionSettings()
    thermal_settings.SetUnknownVariable(globals()[config.unknown_variable])
    model_part.AddNodalSolutionStepVariable(thermal_settings.GetUnknownVariable());

    for node in model_part.Nodes:
        node.AddDof(thermal_settings.GetUnknownVariable());
    print("DOFs for the convection diffusion solver added correctly")


class ConvectionDiffusionSolver:

    def __init__(self, model_part, domain_size, my_settings):

        self.settings = my_settings
        self.thermal_settings = ConvectionDiffusionSettings()

        self.thermal_settings = ConvectionDiffusionSettings()
        self.thermal_settings.SetUnknownVariable(globals()[self.settings.unknown_variable])
        self.thermal_settings.SetDensityVariable(globals()[self.settings.density_variable])
        self.thermal_settings.SetDiffusionVariable(globals()[self.settings.diffusion_variable])
        self.thermal_settings.SetVolumeSourceVariable(globals()[self.settings.volume_source_variable])
        self.thermal_settings.SetSurfaceSourceVariable(globals()[self.settings.surface_source_variable])
        self.thermal_settings.SetMeshVelocityVariable(globals()[self.settings.mesh_velocity_variable])
        self.thermal_settings.SetProjectionVariable(globals()[self.settings.projection_variable])
        self.thermal_settings.SetSpecificHeatVariable(globals()[self.settings.specific_heat_variable])
        self.thermal_settings.SetVelocityVariable(globals()[self.settings.velocity_variable])

        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part, number_of_avg_elems, number_of_avg_nodes)

        self.model_part = model_part
        self.domain_size = domain_size

        # assignation of parameters to be used
        self.time_order = 2;
        self.prediction_order = 1;
        self.ReformDofAtEachIteration = False;

        self.echo_level = 0

        # definition of the solvers
        pILUPrecond = ILU0Preconditioner()
        self.linear_solver = BICGSTABSolver(1e-6, 5000, pILUPrecond)

    def Initialize(self):
        (self.neighbour_search).Execute()
        self.model_part.ProcessInfo

#	self.settings = my_settings

        self.thermal_settings.SetUnknownVariable(globals()[self.settings.unknown_variable])
        self.thermal_settings.SetDensityVariable(globals()[self.settings.density_variable])
        self.thermal_settings.SetDiffusionVariable(globals()[self.settings.diffusion_variable])
        self.thermal_settings.SetVolumeSourceVariable(globals()[self.settings.volume_source_variable])
        self.thermal_settings.SetSurfaceSourceVariable(globals()[self.settings.surface_source_variable])
        self.thermal_settings.SetMeshVelocityVariable(globals()[self.settings.mesh_velocity_variable])
        self.thermal_settings.SetProjectionVariable(globals()[self.settings.projection_variable])
        self.thermal_settings.SetSpecificHeatVariable(globals()[self.settings.specific_heat_variable])
        self.thermal_settings.SetVelocityVariable(globals()[self.settings.velocity_variable])

        (self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, self.thermal_settings)

        self.solver = ResidualBasedConvectionDiffusionStrategy(self.model_part, self.linear_solver, self.ReformDofAtEachIteration, self.time_order, self.prediction_order)
        (self.solver).SetEchoLevel(self.echo_level)
        print("finished initialization of the fluid strategy")

    def Solve(self):
        if(self.ReformDofAtEachIteration):
            (self.neighbour_search).Execute()

        self.thermal_settings = ConvectionDiffusionSettings()
        self.thermal_settings.SetUnknownVariable(globals()[self.settings.unknown_variable])
        self.thermal_settings.SetDensityVariable(globals()[self.settings.density_variable])
        self.thermal_settings.SetDiffusionVariable(globals()[self.settings.diffusion_variable])
        self.thermal_settings.SetVolumeSourceVariable(globals()[self.settings.volume_source_variable])
        self.thermal_settings.SetSurfaceSourceVariable(globals()[self.settings.surface_source_variable])
        self.thermal_settings.SetMeshVelocityVariable(globals()[self.settings.mesh_velocity_variable])
        self.thermal_settings.SetProjectionVariable(globals()[self.settings.projection_variable])
        self.thermal_settings.SetSpecificHeatVariable(globals()[self.settings.specific_heat_variable])
        self.thermal_settings.SetVelocityVariable(globals()[self.settings.velocity_variable])

        (self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, self.thermal_settings)

        (self.solver).Solve()


#
#
def CreateSolver(model_part, config):
    conv_diff_solver = ConvectionDiffusionSolver(model_part, config.domain_size, config)

    # default settings
    if(hasattr(config, "domain_size")):
        conv_diff_solver.domain_size = config.domain_size
    if(hasattr(config, "ReformDofAtEachIteration")):
        conv_diff_solver.ReformDofAtEachIteration = config.ReformDofAtEachIteration
    if(hasattr(config, "echo_level")):
        conv_diff_solver.echo_level = config.echo_level

    # linear solver settings
    import linear_solver_factory
    if(hasattr(config, "linear_solver_config")):
        conv_diff_solver.linear_solver = linear_solver_factory.ConstructSolver(config.linear_solver_config)

    return conv_diff_solver

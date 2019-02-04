from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.RadiationApplication import *
from KratosMultiphysics.ExternalSolversApplication import *


def AddVariables(model_part, config):
    #
    radiation_settings = RadiationSettings()
    radiation_settings.SetUnknownVariable(globals()[config.unknown_variable])
    #radiation_settings.SetDensityVariable(globals()[config.density_variable])
    #radiation_settings.SetMeshVelocityVariable(globals()[config.mesh_velocity_variable])

    #
    model_part.AddNodalSolutionStepVariable(radiation_settings.GetUnknownVariable())
    #model_part.AddNodalSolutionStepVariable(radiation_settings.GetDensityVariable())
    #model_part.AddNodalSolutionStepVariable(radiation_settings.GetMeshVelocityVariable())

    #model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    #model_part.AddNodalSolutionStepVariable(VELOCITY);
    #model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    #model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    #model_part.AddNodalSolutionStepVariable(ENTHALPY);


def AddDofs(model_part, config):
    radiation_settings = RadiationSettings()
    radiation_settings.SetUnknownVariable(globals()[config.unknown_variable])
    model_part.AddNodalSolutionStepVariable(radiation_settings.GetUnknownVariable());

    for node in model_part.Nodes:
        node.AddDof(radiation_settings.GetUnknownVariable());

class ConvectionDiffusionrSolver:

    def __init__(self,model_part,domain_size,my_settings):

        self.settings = my_settings
        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)

        self.model_part = model_part
        self.domain_size = domain_size

        #assignation of parameters to be used
        self.time_order = 1;
        self.prediction_order = 1;
        self.ReformDofAtEachIteration = False;
        self.max_iter = 15;
        self.toll = 1e-9;

        self.echo_level = 0

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        self.linear_solver =  SuperLUSolver()


    def Initialize(self):

        (self.neighbour_search).Execute()
        self.model_part.ProcessInfo
        self.radiation_settings = RadiationSettings()
        self.radiation_settings.SetUnknownVariable(globals()[self.settings.unknown_variable])
        #self.radiation_settings.SetDensityVariable(globals()[self.settings.density_variable])
        #self.radiation_settings.SetMeshVelocityVariable(globals()[self.settings.mesh_velocity_variable])

        (self.model_part.ProcessInfo).SetValue(RADIATION_SETTINGS,self.radiation_settings)

        self.solver = ResidualBasedConvectionDiffusionrStrategyNonLinear(self.model_part,self.linear_solver,self.ReformDofAtEachIteration,self.time_order,self.max_iter,self.toll)
        (self.solver).SetEchoLevel(self.echo_level)


    def Solve(self):
        if(self.ReformDofAtEachIteration == True):
            (self.neighbour_search).Execute()
            (self.solver).Clear()


        self.radiation_settings = RadiationSettings()
        self.radiation_settings.SetUnknownVariable(globals()[self.settings.unknown_variable])
        #self.radiation_settings.SetDensityVariable(globals()[self.settings.density_variable])
        #self.radiation_settings.SetMeshVelocityVariable(globals()[self.settings.mesh_velocity_variable])

        (self.model_part.ProcessInfo).SetValue(RADIATION_SETTINGS,self.radiation_settings)

        (self.solver).Solve()

        if(self.ReformDofAtEachIteration == True):
            (self.solver).Clear()

def CreateSolver(model_part, config):
    conv_diffr_solver = ConvectionDiffusionrSolver(model_part, config.domain_size, config)
    print (conv_diffr_solver)
    #sssssssss
    # default settings
    if(hasattr(config, "time_order")):
        conv_diffr_solver.time_order = config.time_order
    if(hasattr(config, "domain_size")):
        conv_diffr_solver.domain_size = config.domain_size
    if(hasattr(config, "prediction_order")):
        conv_diffr_solver.prediction_order = config.prediction_order
    if(hasattr(config, "ReformDofAtEachIteration")):
        conv_diffr_solver.ReformDofAtEachIteration = config.ReformDofAtEachIteration
    if(hasattr(config, "echo_level")):
        conv_diffr_solver.echo_level = config.echo_level
    if(hasattr(config, "max_ite")):
        conv_diffr_solver.max_iter = config.max_iter
    if(hasattr(config, "toll")):
        conv_diffr_solver.toll = config.toll

    # linear solver settings
    import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
    if(hasattr(config, "linear_solver_config")):
        conv_diffr_solver.linear_solver = linear_solver_factory.ConstructSolver(config.linear_solver_config)

    return conv_diffr_solver


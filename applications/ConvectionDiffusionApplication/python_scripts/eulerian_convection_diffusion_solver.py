from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ConvectionDiffusionApplication import *

#implementation for string (py) settings
def AddVariables(model_part, py_settings=None):
    #if python(string) settings are given, we copy all the settings to the model part.
    #otherwise, we assume the user has already defined the settings in the model part.
    if py_settings is not None:
        #we must copy the settings to a c++ object, from now on we will only use the model part, input (py) settings will be ignored in the adddofs, initialize, etc
        thermal_settings = ConvectionDiffusionSettings() #c++ settings
        if hasattr(py_settings, "velocity_variable"):
           print(py_settings.velocity_variable)
           thermal_settings.SetVelocityVariable(eval(py_settings.velocity_variable))
        if hasattr(py_settings, "mesh_velocity_variable"):
            thermal_settings.SetMeshVelocityVariable(eval(py_settings.mesh_velocity_variable))
        if hasattr(py_settings, "diffusion_variable"):
            thermal_settings.SetDiffusionVariable(eval(py_settings.diffusion_variable))
        if hasattr(py_settings, "unknown_variable"):
            thermal_settings.SetUnknownVariable(eval(py_settings.unknown_variable))
        if hasattr(py_settings, "specific_heat_variable"):
            thermal_settings.SetSpecificHeatVariable(eval(py_settings.specific_heat_variable))
        if hasattr(py_settings, "density_variable"):
            thermal_settings.SetDensityVariable(eval(py_settings.density_variable))
        #and now we save it in the model part.
        (model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, thermal_settings)
    #now we can add the variable, as the name suggest
    if (model_part.ProcessInfo).Has(CONVECTION_DIFFUSION_SETTINGS):
        thermal_settings  = (model_part.ProcessInfo).GetValue(CONVECTION_DIFFUSION_SETTINGS)
        if thermal_settings.IsDefinedVelocityVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetVelocityVariable())
        if thermal_settings.IsDefinedMeshVelocityVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetMeshVelocityVariable())
        if thermal_settings.IsDefinedDiffusionVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetDiffusionVariable())
        if thermal_settings.IsDefinedUnknownVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetUnknownVariable())
        if thermal_settings.IsDefinedSpecificHeatVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetSpecificHeatVariable())
        if thermal_settings.IsDefinedDensityVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetDensityVariable())
        if thermal_settings.IsDefinedVolumeSourceVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetVolumeSourceVariable())
    else:
        raise ValueError('CONVECTION_DIFFUSION_SETTINGS not defined in the model part!')

def AddDofs(model_part, settings=None):
    thermal_settings  = (model_part.ProcessInfo).GetValue(CONVECTION_DIFFUSION_SETTINGS)
    for node in model_part.Nodes:
        node.AddDof(thermal_settings.GetUnknownVariable());
    print("user should include the variables as needed by nodes")


def CreateSolver(model_part, config):
    convection_solver = EulerianConvectionDiffusionSolver(model_part, config.domain_size)
    # linear solver settings
    import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
    if(hasattr(config, "convection_linear_solver_config")):
        self.linear_solver = linear_solver_factory.ConstructSolver(
            config.convection_linear_solver_config)

    return convection_solver

class EulerianConvectionDiffusionSolver:
    def __init__(self, model_part, domain_size):
        self.model_part = model_part
        self.domain_size = domain_size

        # assignation of parameters to be used
        self.ReformDofAtEachIteration = False;

        self.scalar_var_convected = 1

        # definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        self.linear_solver = BICGSTABSolver(1e-9, 5000, pDiagPrecond)

        # linear solver settings
        #import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        #if(hasattr(config, "convection_linear_solver_config")):
        #    self.linear_solver = linear_solver_factory.ConstructSolver(
        #        config.convection_linear_solver_config)

        model_part.ProcessInfo[TIME_INTEGRATION_THETA] = 0.5 #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)

    def Initialize(self):
        # convection diffusion tool
        self.convection_solver = ResidualBasedEulerianConvectionDiffusionStrategy(self.model_part,self.linear_solver,self.ReformDofAtEachIteration,self.domain_size)
        print("Finished Initialize")

    def SetEchoLevel(self, level):
        self.convection_solver.SetEchoLevel(level)

    def Solve(self):
        self.convection_solver.Solve()

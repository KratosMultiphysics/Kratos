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
        if (thermal_settings).IsDefinedVelocityVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetVelocityVariable())
        if (thermal_settings).IsDefinedMeshVelocityVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetMeshVelocityVariable())
        if (thermal_settings).IsDefinedDiffusionVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetDiffusionVariable())
        if (thermal_settings).IsDefinedUnknownVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetUnknownVariable())
        if (thermal_settings).IsDefinedSpecificHeatVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetSpecificHeatVariable())
        if (thermal_settings).IsDefinedDensityVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetDensityVariable())
        if (thermal_settings).IsDefinedProjectionVariable():
            model_part.AddNodalSolutionStepVariable(thermal_settings.GetProjectionVariable())

    else:
        raise ValueError('CONVECTION_DIFFUSION_SETTINGS not defined in the model part!')
    #now we add the specific variables needed for the PFEM-2 convection
    model_part.AddNodalSolutionStepVariable(PROJECTED_SCALAR1);

def AddDofs(model_part, settings=None):
    thermal_settings  = (model_part.ProcessInfo).GetValue(CONVECTION_DIFFUSION_SETTINGS)
    for node in model_part.Nodes:
        node.AddDof(thermal_settings.GetUnknownVariable());
    print("user should include the variables as needed by nodes")


def CreateSolver(model_part, config):
    convection_solver = PFEM2ConvectionDiffusionSolver(model_part, config.domain_size)
    # linear solver settings
    import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
    if(hasattr(config, "convection_linear_solver_config")):
        self.linear_solver = linear_solver_factory.ConstructSolver(
            config.convection_linear_solver_config)

    return convection_solver

class BFECCConvectionDiffusionSolver:
    def __init__(self, model_part, domain_size):
        self.model_part = model_part
        self.domain_size = domain_size

        # assignation of parameters to be used
        self.ReformDofAtEachIteration = False;

        self.scalar_var_convected = 1

        # definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        self.linear_solver = BICGSTABSolver(1e-9, 5000, pDiagPrecond)

        #for the bfecc
        #mount the search structure
        self.locator = BinBasedFastPointLocator2D(self.model_part)
        self.locator.UpdateSearchDatabase()
        #self.locator.UpdateSearchDatabaseAssignedSize(0.01)


    def Initialize(self):
        # diffusion only tool
        self.diffusion_solver = ResidualBasedSemiEulerianConvectionDiffusionStrategy(self.model_part,self.linear_solver,self.ReformDofAtEachIteration,self.domain_size)


        self.thermal_settings  = (self.model_part.ProcessInfo).GetValue(CONVECTION_DIFFUSION_SETTINGS)
        self.unknown_var = (self.thermal_settings).GetUnknownVariable()
        self.projection_var = (self.thermal_settings).GetProjectionVariable()
        self.velocity_var = (self.thermal_settings).GetVelocityVariable()

        #now tools for pfem2:
        self.VariableUtils = VariableUtils()

        #construct the utility to move the points
        if self.domain_size ==2:
            self.bfecc_utility = BFECCConvection2D(self.locator)
        else:
            self.bfecc_utility = BFECCConvection3D(self.locator)


        print("Finished Initialize")

    def Solve(self):
        substepping  = 10.0
        (self.VariableUtils).CopyScalarVar(self.unknown_var,self.projection_var,self.model_part.Nodes)
        self.bfecc_utility.CopyScalarVarToPreviousTimeStep(self.model_part,self.projection_var)
        self.bfecc_utility.BFECCconvect(self.model_part,self.projection_var,self.velocity_var,substepping)
        #self.bfecc_utility.ResetBoundaryConditions(self.model_part,self.unknown_var)
        #self.bfecc_utility.CopyScalarVarToPreviousTimeStep(self.model_part,self.unknown_var)
        #we only solve the mesh problem if there is diffusion or heat sources. otherwise->pure convection problem
        if (self.thermal_settings).IsDefinedDiffusionVariable() or (self.thermal_settings).IsDefinedSurfaceSourceVariable() or (self.thermal_settings).IsDefinedVolumeSourceVariable():
              (self.diffusion_solver).Solve()



def CreateSolver(model_part, config):
    convection_solver = BFECCConvectionDiffusionSolver(model_part, config.domain_size)
    # linear solver settings
    import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
    if(hasattr(config, "convection_linear_solver_config")):
        self.linear_solver = linear_solver_factory.ConstructSolver(
            config.convection_linear_solver_config)

    return convection_solver


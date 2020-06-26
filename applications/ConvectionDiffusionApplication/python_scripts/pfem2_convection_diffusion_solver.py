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

    else:
        raise ValueError('CONVECTION_DIFFUSION_SETTINGS not defined in the model part!')
    #now we add the specific variables needed for the PFEM-2 convection
    model_part.AddNodalSolutionStepVariable(PROJECTED_SCALAR1);
    model_part.AddNodalSolutionStepVariable(DELTA_SCALAR1);
    model_part.AddNodalSolutionStepVariable(YP);
    model_part.AddNodalSolutionStepVariable(MEAN_SIZE);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    #model_part.AddNodalSolutionStepVariable(NORMAL);

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

class PFEM2ConvectionDiffusionSolver:
    def __init__(self, model_part, domain_size):
        self.model_part = model_part
        self.domain_size = domain_size

        self.thermal_settings  = (model_part.ProcessInfo).GetValue(CONVECTION_DIFFUSION_SETTINGS)
        self.unknown_var = (self.thermal_settings).GetUnknownVariable()

        # assignation of parameters to be used
        self.ReformDofAtEachIteration = False;

        self.scalar_var_convected = 1

        # definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        self.linear_solver = BICGSTABSolver(1e-9, 5000, pDiagPrecond)

        #for the pfem2
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part)
        (self.neighbour_search).Execute()
        self.neighbour_elements_search= FindElementalNeighboursProcess(model_part,domain_size,number_of_avg_elems)
        (self.neighbour_elements_search).Execute()
        ##calculate normals
        #self.normal_tools = BodyNormalCalculationUtils()
        #self.normal_tools.CalculateBodyNormals(self.model_part,self.domain_size);

    def Initialize(self):
        # diffusion only tool
        self.diffusion_solver = ResidualBasedSemiEulerianConvectionDiffusionStrategy(self.model_part,self.linear_solver,self.ReformDofAtEachIteration,self.domain_size)

        #now tools for pfem2:
        self.VariableUtils = VariableUtils()

        maximum_number_of_particles= 8*self.domain_size
        if self.domain_size==2:
            self.moveparticles = MoveParticleUtilityScalarTransport2D(self.model_part,maximum_number_of_particles)
        else:
            self.moveparticles = MoveParticleUtilityScalarTransport3D(self.model_part,maximum_number_of_particles)
        self.moveparticles.MountBin()


        print("Finished Initialize")

    def Solve(self):
        (self.moveparticles).CalculateVelOverElemSize();
        (self.moveparticles).MoveParticles();
        pre_minimum_number_of_particles=self.domain_size;
        (self.moveparticles).PreReseed(pre_minimum_number_of_particles);
        (self.moveparticles).TransferLagrangianToEulerian();
        #(self.VariableUtils).CopyScalarVar(PROJECTED_SCALAR1,self.unknown_var,self.model_part.Nodes)
        #(self.moveparticles).ResetBoundaryConditions()
        #(self.moveparticles).CopyScalarVarToPreviousTimeStep(self.unknown_var,self.model_part.Nodes)

        #we only solve the mesh problem if there is diffusion or heat sources. otherwise->pure convection problem
        if (self.thermal_settings).IsDefinedDiffusionVariable() or (self.thermal_settings).IsDefinedSurfaceSourceVariable() or (self.thermal_settings).IsDefinedVolumeSourceVariable():
             (self.diffusion_solver).Solve()

        (self.moveparticles).CalculateDeltaVariables();
        (self.moveparticles).CorrectParticlesWithoutMovingUsingDeltaVariables();
        post_minimum_number_of_particles=self.domain_size*2;
        (self.moveparticles).PostReseed(post_minimum_number_of_particles);


def CreateSolver(model_part, config):
    convection_solver = PFEM2ConvectionDiffusionSolver(model_part, config.domain_size)
    # linear solver settings
    import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
    if(hasattr(config, "convection_linear_solver_config")):
        self.linear_solver = linear_solver_factory.ConstructSolver(
            config.convection_linear_solver_config)

    return convection_solver


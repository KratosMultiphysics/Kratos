from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

## Import base class file
# from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_compressible_solver import NavierStokesCompressibleSolver

from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid

from Kratos import *

def CreateSolver(model, custom_settings):
    return NavierStokesBiphaseDGCompressibleExplicitSolver(model, custom_settings)

class NavierStokesBiphaseDGCompressibleExplicitSolver(NavierStokesCompressibleSolver):

    @classmethod
    def GetDefaultSettings(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "compressible_solver_from_defaults",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "two_element_test",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "maximum_iterations": 1,
            "echo_level": 1,
            "time_order": 2,  
            "compute_reactions": false,
            "reform_dofs_at_each_step" : true,
            "relative_tolerance" : 1e-3,
            "absolute_tolerance" : 1e-5,
            "linear_solver_settings"       : {
                "solver_type"         : "amgcl",
                "max_iteration"       : 200,
                "tolerance"           : 1e-7,
                "provide_coordinates" : false,
                "smoother_type"       : "ilu0",
                "krylov_type"         : "gmres",
                "coarsening_type"     : "aggregation",
                "scaling"             : true,
                "verbosity"           : 0
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "assign_neighbour_elements_to_conditions": false,
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01
            },
            "periodic": "periodic",
            "move_mesh_flag": false
        }""")

        default_settings.AddMissingParameters(super(NavierStokesBiphaseDGCompressibleExplicitSolver, cls).GetDefaultSettings())
        return default_settings

    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(NavierStokesBiphaseDGCompressibleExplicitSolver,self).__init__(model,custom_settings)

##        self.element_name = "CompressibleNavierStokes"
        self.element_name = "CompressibleBiphaseDGNavierStokesExplicit"
        self.condition_name = "Condition"
        self.min_buffer_size = 2

        ## Set the nodal properties flag
        self.element_has_nodal_properties = True

        print("Construction of NavierStokesBiphaseDGCompressibleExplicitSolver finished.")

    def AddVariables(self):     ## Che cosa fa questa funzione? Devono esistere da qualche parte le variabili che aggiungo?
        
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.MOMENTUM_RK4)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.MOMENTUM_RHS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.DENSITY_RHS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.DENSITY_RK4)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.DENSITY_GAS_RHS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.DENSITY_GAS_RK4)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.DENSITY_SOLID_RHS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.DENSITY_SOLID_RK4)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.TOTAL_ENERGY_RHS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.TOTAL_ENERGY_RK4)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
                
        
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.DENSITY_SOLID)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.DENSITY_GAS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENTUM)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TOTAL_ENERGY)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONDUCTIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SPECIFIC_HEAT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.HEAT_CAPACITY_RATIO)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.DENSITY_MATERIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.SPECIFIC_HEAT_MATERIAL)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H) ## ?
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA) ## ?
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)  #for momentum
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.REACTION_DENSITY)  #for momentum
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.REACTION_DENSITY_GAS)  #for momentum
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.REACTION_DENSITY_SOLID)  #for momentum
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.REACTION_ENERGY)  #for momentum

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE) ## ?
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.Y_WALL) ## ?
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.KINEMATIC_VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)

        # Post-process
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.MACH)  #for momentum
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.SOLID_CONCENTRATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.DYNAMIC_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.EXTRA_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.GAS_PRESSURE)
        
        print("Monolithic compressible fluid solver variables added correctly")

    def AddDofs(self):
        
        print("\n\nHere I am - AddDofs \n\n")
        
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_X, KratosMultiphysics.REACTION_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_Y, KratosMultiphysics.REACTION_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_Z, KratosMultiphysics.REACTION_Z, self.main_model_part)
        # KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DENSITY, KratosFluid.REACTION_DENSITY, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosFluid.DENSITY_GAS, KratosFluid.REACTION_DENSITY_GAS, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosFluid.DENSITY_SOLID, KratosFluid.REACTION_DENSITY_SOLID, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.TOTAL_ENERGY, KratosFluid.REACTION_ENERGY, self.main_model_part)

        print("\n\nHere I am - AddedDofs \n\n")

    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            print("ERROR: _GetAutomaticTimeSteppingUtility out of date")
            #self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        # Set the time discretization utility to compute the BDF coefficients
        time_order = self.settings["time_order"].GetInt()
        if time_order == 2:
            self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order) ## To develop for different RK schemes
        else:
            raise Exception("Only \"time_order\" equal to 2 is supported. Provided \"time_order\": " + str(time_order))

        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]


        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.computing_model_part, self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        neighbour_finder = FindNodalNeighboursProcess(self.computing_model_part,6, 6)
        neighbour_finder.Execute() ##at wish ... when it is needed

#        self.solver = KratosFluid.RungeKuttaStrategy(
#            self.GetComputingModelPart(),
#            self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
#            self.settings["compute_reactions"].GetBool(),
#            self.settings["reform_dofs_at_each_step"].GetBool(),
#            self.settings["move_mesh_flag"].GetBool())

        self.solver = KratosFluid.ExplicitEulerDGStrategy(
            self.GetComputingModelPart(),
            self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
            self.settings["compute_reactions"].GetBool(),
            self.settings["reform_dofs_at_each_step"].GetBool(),
            self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
       
        (self.solver).Initialize()
        

        # self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble()) # REMEMBER TO CHECK MY STAB CONSTANTS

        print ("Monolithic compressible explicit solver initialization finished.")


    def InitializeSolutionStep(self):
        (self.solver).InitializeSolutionStep()


    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            is_converged = self.solver.SolveSolutionStep()
            return is_converged
        else:
            return True
        

    def Solve(self):
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)
        (self.solver).Solve()

    # def PrepareModelPart(self):
    #     super(NavierStokesCompressibleSolver,self).PrepareModelPart()
    #     if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
    #         self._ExecuteAfterReading()

    def _ExecuteAfterReading(self):
        ## Replace element and conditions
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

        ## Check that the input read has the shape we like
        prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
        prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
        prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])

        check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()

    def _SetNodalProperties(self):
        # Get density from the properties of the first element
        # Set it as initial condition for the DENSITY unknown in all the nodes
        for el in self.main_model_part.Elements:
            rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            if rho <= 0.0:
                raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            break
        else:
            raise Exception("No fluid elements found in the main model part.")
        # Transfer the obtained properties to the nodes
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)



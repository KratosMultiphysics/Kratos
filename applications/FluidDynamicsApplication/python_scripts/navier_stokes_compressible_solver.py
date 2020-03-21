from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

## Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid

def CreateSolver(model, custom_settings):
    return NavierStokesCompressibleSolver(model, custom_settings)

class NavierStokesCompressibleSolver(FluidSolver):

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
            "maximum_iterations": 10,
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

        default_settings.AddMissingParameters(super(NavierStokesCompressibleSolver, cls).GetDefaultSettings())
        return default_settings

    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(NavierStokesCompressibleSolver,self).__init__(model,custom_settings)

        self.element_name = "CompressibleNavierStokes"
        self.condition_name = "Condition"
        self.min_buffer_size = 3

        ## Construct the linear solver
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        ## Set the element replace settings
        #self._SetCompressibleElementReplaceSettings()

        print("Construction of NavierStokesCompressibleSolver finished.")


    def AddVariables(self):

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENTUM)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TOTAL_ENERGY)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONDUCTIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SPECIFIC_HEAT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.HEAT_CAPACITY_RATIO)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H) ## ?
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA) ## ?
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)  #for momentum
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluid.REACTION_DENSITY)  #for momentum
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

        print("Monolithic compressible fluid solver variables added correctly")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_X, KratosMultiphysics.REACTION_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_Y, KratosMultiphysics.REACTION_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MOMENTUM_Z, KratosMultiphysics.REACTION_Z, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DENSITY, KratosFluid.REACTION_DENSITY, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.TOTAL_ENERGY, KratosFluid.REACTION_ENERGY, self.main_model_part)

    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            print("ERROR: _GetAutomaticTimeSteppingUtility out of date")
            #self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        # Set the time discretization utility to compute the BDF coefficients
        time_order = self.settings["time_order"].GetInt()
        if time_order == 2:
            self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
        else:
            raise Exception("Only \"time_order\" equal to 2 is supported. Provided \"time_order\": " + str(time_order))

        # Creating the solution strategy
        self.conv_criteria = KratosMultiphysics.ResidualCriteria(self.settings["relative_tolerance"].GetDouble(),
                                                                 self.settings["absolute_tolerance"].GetDouble())


        #(self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt()
        (self.conv_criteria).SetEchoLevel(3)

        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        rotation_utility = KratosFluid.CompressibleElementRotationUtility(domain_size,KratosMultiphysics.SLIP)
        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(rotation_utility)
        #time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme() # DOFs (4,5)


        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)


        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part,
                                                                            time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())


        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
        #(self.solver).SetEchoLevel(1)


        (self.solver).Initialize()


        # self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble()) # REMEMBER TO CHECK MY STAB CONSTANTS

        print ("Monolithic compressible solver initialization finished.")


    def InitializeSolutionStep(self):
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)
        (self.solver).InitializeSolutionStep()


    def Solve(self):
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)
        (self.solver).Solve()

    def PrepareModelPart(self):
        super(NavierStokesCompressibleSolver,self).PrepareModelPart()
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self._ExecuteAfterReading()

    def _ExecuteAfterReading(self):
        ## Replace element and conditions
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

        ## Check that the input read has the shape we like
        prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
        prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
        prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])

        check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()


    #def _SetCompressibleElementReplaceSettings(self):
        #domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        #self.settings.AddEmptyValue("element_replace_settings")

        #if(domain_size == 3):
            #self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
            #{
                #"element_name":"CompressibleNavierStokes3D4N",
                #"condition_name": "SurfaceCondition3D3N"
            #}
            #""")
        #elif(domain_size == 2):
            #self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
            #{
                #"element_name":"CompressibleNavierStokes2D3N",
                #"condition_name": "LineCondition2D2N"
            #}
            #""")
        #else:
            #raise Exception("Domain size is not 2 or 3!!")

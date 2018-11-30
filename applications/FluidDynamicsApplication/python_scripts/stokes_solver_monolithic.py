from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics as kratoscore
import KratosMultiphysics.FluidDynamicsApplication as cfd
import KratosMultiphysics.IncompressibleFluidApplication

# Check that KratosMultiphysics was imported in the main script
kratoscore.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return StokesSolver(main_model_part, custom_settings)

class StokesSolver:

    ##constructor. the constructor shall only take care of storing the settings
    ##and the pointer to the main_model part. This is needed since at the point of constructing the
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = kratoscore.Parameters("""
        {
            "solver_type": "stokes_solver_monolithic",
            "force_steady_state": false,
            "velocity_tolerance": 1e-3,
            "pressure_tolerance": 1e-2,
            "absolute_velocity_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-2,
            "maximum_iterations": 3,
            "echo_level": 1,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "dynamic_tau": 0.001,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "volume_model_part_name" : "volume_model_part",
            "skin_parts":[""],
            "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
            },
            "linear_solver_settings": {
                    "solver_type": "Super LU",
                    "max_iteration": 500,
                    "tolerance": 1e-9,
                    "scaling": false,
                    "verbosity": 1
            }
        }""")

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #construct the linear solvers
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])


        print("Construction of NavierStokesSolver_FractionalStep finished")

    def GetMinimumBufferSize(self):
        return 3;

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(kratoscore.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(kratoscore.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(kratoscore.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(kratoscore.BODY_FORCE) #TODO: decide if it is needed. if constant it could be passed in properties
        self.main_model_part.AddNodalSolutionStepVariable(kratoscore.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(kratoscore.NORMAL) #TODO: this variable is not strictly needed by the solver - may be needed by other utilities

        if(self.settings["compute_reactions"].GetBool()):
            self.main_model_part.AddNodalSolutionStepVariable(kratoscore.REACTION)
            self.main_model_part.AddNodalSolutionStepVariable(kratoscore.REACTION_WATER_PRESSURE)

        print("variables for the  monolithic stokes solver added correctly")


    def AddDofs(self):
        kratoscore.VariableUtils().AddDof(kratoscore.VELOCITY_X, kratoscore.REACTION_X,self.main_model_part)
        kratoscore.VariableUtils().AddDof(kratoscore.VELOCITY_Y, kratoscore.REACTION_Y,self.main_model_part)
        kratoscore.VariableUtils().AddDof(kratoscore.VELOCITY_Z, kratoscore.REACTION_Z,self.main_model_part)
        kratoscore.VariableUtils().AddDof(kratoscore.PRESSURE, kratoscore.REACTION_WATER_PRESSURE,self.main_model_part)


    def ImportModelPart(self):

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            #here it would be the place to import restart data if required
            kratoscore.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)

            self.settings.AddEmptyValue("element_replace_settings")
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "element_name":"StokesTwoFluid3D4N",
                    "condition_name": "StokesWallCondition3D"
                    }
                    """)
            elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                raise Exception("sorry 2D case not implemented")
            else:
                raise Exception("domain size is not 2 or 3")

            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()


            ##here we shall check that the input read has the shape we like
            aux_params = kratoscore.Parameters("{}")
            aux_params.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
            aux_params.AddValue("skin_parts",self.settings["skin_parts"])
            import check_and_prepare_model_process_fluid
            check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()


            ##here we must construct correctly the constitutive law
            print("don't forget constructing the constitutive law!!!!!!!!!!!!!!")


            #if needed here we shall generate the constitutive laws
            #import constitutive_law_python_utility as constitutive_law_utils
            #constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(main_model_part, self.settings["DomainSize"]);
            #constitutive_law.Initialize();
        else:
            raise Exception("other input options are not yet implemented")

        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )

        print ("model reading finished")

    def Initialize(self):
        compute_model_part = self.GetComputingModelPart()

        self.bdf_process = kratoscore.ComputeBDFCoefficientsProcess(compute_model_part,2)

        time_scheme = kratoscore.ResidualBasedIncrementalUpdateStaticScheme()

        convergence_criteria = KratosMultiphysics.FluidDynamicsApplication.VelPrCriteria(self.settings["velocity_tolerance"].GetDouble(),
                                                                    self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                                    self.settings["pressure_tolerance"].GetDouble(),
                                                                    self.settings["absolute_pressure_tolerance"].GetDouble()
                                                                    )

        (convergence_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        builder_and_solver = kratoscore.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        move_mesh_flag = False #user should NOT configure this
        self.fluid_solver = kratoscore.ResidualBasedNewtonRaphsonStrategy(
                                                                          compute_model_part,
                                                                          time_scheme,
                                                                          self.linear_solver,
                                                                          convergence_criteria,
                                                                          builder_and_solver,
                                                                          self.settings["maximum_iterations"].GetInt(),
                                                                          self.settings["compute_reactions"].GetBool(),
                                                                          self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                          move_mesh_flag
                                                                          )

        (self.fluid_solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        print("Construction stokes solver finished")

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass

    def SaveRestart(self):
        pass #one should write the restart file here

    def Solve(self):
        self.bdf_process.Execute()

        if(self.settings["force_steady_state"].GetBool()):
            bdf_vec = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.BDF_COEFFICIENTS]
            for i in range(len(bdf_vec)):
                bdf_vec[i] = 0.0
            self.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.BDF_COEFFICIENTS, bdf_vec)

        self.fluid_solver.Solve()

    def InitializeSolutionStep(self):
        self.bdf_process.Execute()

        if(self.settings["force_steady_state"].GetBool()):
            bdf_vec = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.BDF_COEFFICIENTS]
            for i in range(len(bdf_vec)):
                bdf_vec[i] = 0.0
            self.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.BDF_COEFFICIENTS, bdf_vec)

        self.fluid_solver.InitializeSolutionStep()

    def Predict(self):
        self.fluid_solver.Predict()

    def SolveSolutionStep(self):
        self.fluid_solver.SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()

    def SetEchoLevel(self, level):
        self.fluid_solver.SetEchoLevel(level)

    def Clear(self):
        self.fluid_solver.Clear()

    def Check(self):
        self.fluid_solver.Check()



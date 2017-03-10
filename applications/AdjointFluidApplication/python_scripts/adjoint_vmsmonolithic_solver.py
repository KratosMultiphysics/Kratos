from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.AdjointFluidApplication as AdjointFluidApplication
#import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return AdjointVMSMonolithicSolver(main_model_part, custom_settings)

class AdjointVMSMonolithicSolver:

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
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "adjoint_vmsmonolithic_solver",
            "scheme_settings" : {
                "scheme_type" : "bossak"
            },
            "objective_settings" : {
                "objective_type" : "drag"
            },
            "model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name"
            },
            "linear_solver_settings" : {
                "solver_type" : "AMGCL"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts"  : [""],
            "dynamic_tau" : 0.0,
            "oss_switch"  : 0,
            "echo_level"  : 0
        }""")

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        print("Construction of AdjointVMSMonolithicSolver finished")

    def GetMinimumBufferSize(self):
        return 2

    def AddVariables(self):

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_SENSITIVITY)
        
        print("variables for the adjoint fluid solver added correctly")

    def ImportModelPart(self):

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            #here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)

            ##here we shall check that the input read has the shape we like
            aux_params = KratosMultiphysics.Parameters("{}")
            aux_params.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
            aux_params.AddValue("skin_parts",self.settings["skin_parts"])

            #here we replace the dummy elements we read with proper elements
            self.settings.AddEmptyValue("element_replace_settings")
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "element_name": "VMSAdjointElement3D",
                    "condition_name": "Condition3D4N"
                    }
                    """)
            elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "element_name": "VMSAdjointElement2D",
                    "condition_name": "Condition2D3N"
                    }
                    """)
            else:
                raise Exception("domain size is not 2 or 3")

            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

            import check_and_prepare_model_process_fluid
            check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()

            #here we read the KINEMATIC VISCOSITY and DENSITY and we apply it to the nodes
            for el in self.main_model_part.Elements:
                rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
                kin_viscosity = el.Properties.GetValue(KratosMultiphysics.VISCOSITY)
                break

            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)

        else:
            raise Exception("Other input options are not yet implemented.")

        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )

        print ("Model reading finished.")


    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.ADJOINT_PRESSURE)
            node.AddDof(KratosMultiphysics.ADJOINT_VELOCITY_X)
            node.AddDof(KratosMultiphysics.ADJOINT_VELOCITY_Y)
            node.AddDof(KratosMultiphysics.ADJOINT_VELOCITY_Z)

        print("DOFs for the VMS adjoint solver added correctly.")


    def Initialize(self):

        self.computing_model_part = self.GetComputingModelPart()

        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if self.settings["objective_settings"]["objective_type"].GetString() == "drag":
            if (domain_size == 2):
                self.objective_function = AdjointFluidApplication.DragObjectiveFunction2D(self.settings["objective_settings"])
            elif (domain_size == 3):
                self.objective_function = AdjointFluidApplication.DragObjectiveFunction3D(self.settings["objective_settings"])
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))
        else:
            raise Exception("invalid objective_type: " + self.settings["objective_settings"]["objective_type"].GetString())

        if self.settings["scheme_settings"]["scheme_type"].GetString() == "bossak":
            self.time_scheme = AdjointFluidApplication.AdjointBossakScheme(self.settings["scheme_settings"], self.objective_function)
        elif self.settings["scheme_settings"]["scheme_type"].GetString() == "steady":
            self.time_scheme = AdjointFluidApplication.AdjointSteadyScheme(self.settings["scheme_settings"], self.objective_function)
        else:
            raise Exception("invalid scheme_type: " + self.settings["scheme_settings"]["scheme_type"].GetString())

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(self.main_model_part,
                                                                     self.time_scheme,
                                                                     self.linear_solver,
                                                                     builder_and_solver,
                                                                     False,
                                                                     False,
                                                                     False,
                                                                     False)

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.solver.Check()
        
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())

        print ("Adjoint solver initialization finished.")

    def GetComputingModelPart(self):
        # Get as computational model part the "volume_model_part_name" in the ProjectParameters Json string
        #~ return self.main_model_part.GetSubModelPart(self.settings["volume_model_part_name"].GetString())

        # Get as computational model part the submodelpart generated in CheckAndPrepareModelProcess
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass

    def SaveRestart(self):
        pass #one should write the restart file here

    def DivergenceClearance(self):
        pass
        
    def SolverInitialize(self):
        self.solver.Initialize()

    def SolverInitializeSolutionStep(self):
        self.solver.InitializeSolutionStep()

    def SolverPredict(self):
        self.solver.Predict()

    def SolverSolveSolutionStep(self):
        self.solver.SolveSolutionStep()

    def SolverFinalizeSolutionStep(self):
        self.solver.FinalizeSolutionStep()

    def Solve(self):
        self.solver.Solve()

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()

    def Check(self):
        self.solver.Check()

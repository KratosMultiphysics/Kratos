from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
## Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.MetisApplication as MetisApplication
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.AdjointFluidApplication as AdjointFluidApplication

## Checks that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return AdjointVMSMonolithicMPISolver(main_model_part, custom_settings)

import adjoint_vmsmonolithic_solver

class AdjointVMSMonolithicMPISolver(adjoint_vmsmonolithic_solver.AdjointVMSMonolithicSolver):

    def __init__(self, main_model_part, custom_settings):

        self.main_model_part = main_model_part

        # default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "trilinos_adjoint_vmsmonolithic_solver",
            "scheme_settings" : {
                "scheme_type": "bossak"
            },
            "objective_settings" : {
                "objective_type" : "drag"
            },
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "linear_solver_settings" : {
                "solver_type" : "MultiLevelSolver"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "dynamic_tau": 0.0,
            "oss_switch": 0,
            "echo_level": 0
        }""")

        # overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        # construct the linear solver
        import trilinos_linear_solver_factory
        self.trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        if KratosMPI.mpi.rank == 0:
            print("Construction of AdjointVMSMonolithicMPISolver finished.")


    def GetMinimumBufferSize(self):
        return 2


    def AddVariables(self):
        ## Add variables from the base class
        super(AdjointVMSMonolithicMPISolver, self).AddVariables()

        ## Add specific MPI variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        KratosMPI.mpi.world.barrier()

        if KratosMPI.mpi.rank == 0:
            print("Variables for the AdjointVMSMonolithicMPISolver added correctly in each processor.")


    def ImportModelPart(self):

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            # here we read the already existing partitions from the primal solution.
            input_filename = self.settings["model_import_settings"]["input_filename"].GetString()
            mpi_input_filename = input_filename + "_" + str(KratosMPI.mpi.rank)
            self.settings["model_import_settings"]["input_filename"].SetString(mpi_input_filename)
            KratosMultiphysics.ModelPartIO(mpi_input_filename).ReadModelPart(self.main_model_part)

            # here we shall check that the input read has the shape we like
            aux_params = KratosMultiphysics.Parameters("{}")
            aux_params.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
            aux_params.AddValue("skin_parts",self.settings["skin_parts"])

            # here we replace the dummy elements we read with proper elements
            self.settings.AddEmptyValue("element_replace_settings")
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                        "element_name": "VMSAdjointElement3D",
                        "condition_name": "Condition3D3N"
                    }
                    """)
            elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                        "element_name": "VMSAdjointElement2D",
                        "condition_name": "Condition2D2N"
                    }
                    """)
            else:
                raise Exception("domain size is not 2 nor 3")

            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

            import check_and_prepare_model_process_fluid
            check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()

            # here we read the KINEMATIC VISCOSITY and DENSITY and we apply it to the nodes
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

        MetisApplication.SetMPICommunicatorProcess(self.main_model_part).Execute()

        ParallelFillCommunicator = TrilinosApplication.ParallelFillCommunicator(self.main_model_part.GetRootModelPart())
        ParallelFillCommunicator.Execute()

        if KratosMPI.mpi.rank == 0 :
            print("MPI communicators constructed.")
            print ("MPI model reading finished.")


    def AddDofs(self):
        super(AdjointVMSMonolithicMPISolver, self).AddDofs()
        KratosMPI.mpi.world.barrier()

        if KratosMPI.mpi.rank == 0:
            print("DOFs for the AdjointVMSMonolithicMPISolver added correctly in all processors.")


    def Initialize(self):

        self.epetra_comm = TrilinosApplication.CreateCommunicator()

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
            self.time_scheme = TrilinosApplication.TrilinosAdjointBossakScheme(self.settings["scheme_settings"], self.objective_function)
        elif self.settings["scheme_settings"]["scheme_type"].GetString() == "steady":
            self.time_scheme = TrilinosApplication.TrilinosAdjointSteadyScheme(self.settings["scheme_settings"], self.objective_function)
        else:
            raise Exception("invalid scheme_type: " + self.settings["scheme_settings"]["scheme_type"].GetString())

        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
            guess_row_size = 20*4
        elif self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            guess_row_size = 10*3

        self.builder_and_solver = TrilinosApplication.TrilinosBlockBuilderAndSolver(self.epetra_comm,
                                                                               guess_row_size,
                                                                               self.trilinos_linear_solver)

        self.solver = TrilinosApplication.TrilinosLinearStrategy(self.main_model_part,
                                                                 self.time_scheme,
                                                                 self.trilinos_linear_solver,
                                                                 self.builder_and_solver,
                                                                 False,
                                                                 False,
                                                                 False,
                                                                 False)

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.solver.Check()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())

        print ("Monolithic MPI solver initialization finished.")


    def GetComputingModelPart(self):
        # get the submodelpart generated in CheckAndPrepareModelProcess
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")


    def Solve(self):
        (self.solver).Solve()

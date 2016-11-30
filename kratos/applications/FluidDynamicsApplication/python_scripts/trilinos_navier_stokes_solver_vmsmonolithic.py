from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI                          # MPI-python interface
import KratosMultiphysics.MetisApplication as KratosMetis           # Partitioning
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD     # Fluid dynamics application

## Checks that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesMPISolver_VMSMonolithic(main_model_part, custom_settings)

class NavierStokesMPISolver_VMSMonolithic:


    ##constructor. the constructor shall only take care of storing the settings
    ##and the pointer to the main_model part. This is needed since at the point of constructing the
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part
        
        ## Default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "trilinos_navier_stokes_solver_vmsmonolithic",
            "scheme_type": "Bossak",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "maximum_iterations": 10,
            "dynamic_tau": 0.0,
            "oss_switch": 0,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "time_order": 2,
            "compute_reactions": false,
            "divergence_clearance_steps": 0,
            "reform_dofs_at_each_step": false,
            "relative_velocity_tolerance": 1e-5,
            "absolute_velocity_tolerance": 1e-7,
            "relative_pressure_tolerance": 1e-5,
            "absolute_pressure_tolerance": 1e-7,
            "linear_solver_settings"        : {
                "solver_type" : "AztecSolver"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "alpha":-0.3,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "regularization_coef": 1000,
            "MoveMeshFlag": false,
            "use_slip_conditions": false,
            "turbulence_model": "None",
            "use_spalart_allmaras": false
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Construct the linear solver
        import new_trilinos_linear_solver_factory
        self.trilinos_linear_solver = new_trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        print("Construction of NavierStokesMPISolver_VMSMonolithic finished.")
        

    def GetMinimumBufferSize(self):
        return 3;

    def AddVariables(self):

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DIVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.Y_WALL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_TAU) 
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.OSS_SWITCH)  
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.M)           
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)          
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)    

        KratosMPI.mpi.world.barrier()
        
        if KratosMPI.mpi.rank == 0:
            print("Variables for the VMS fluid Trilinos solver added correctly.")

    def ImportModelPart(self):

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            input_filename = self.settings["model_import_settings"]["input_filename"].GetString()

            # Serial partition of the original .mdpa file
            if KratosMPI.mpi.rank == 0 :
                
                # Original .mdpa file reading
                model_part_io = KratosMultiphysics.ModelPartIO(input_filename)
                
                # Partition of the original .mdpa file
                number_of_partitions = KratosMPI.mpi.size # Number of partitions equals the number of processors
                domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
                verbosity = 1
                sync_conditions = True # Make sure that the condition goes to the same partition as the element is a face of
                partitioner = KratosMetis.MetisDivideHeterogeneousInputProcess(model_part_io, number_of_partitions , domain_size, verbosity, sync_conditions)
                partitioner.Execute() 
                
                print("Metis divide finished.")

            KratosMPI.mpi.world.barrier()
            
            ## Read the partitioned .mdpa files
            mpi_input_filename = input_filename + "_" + str(KratosMPI.mpi.rank)
            model_part_io = KratosMultiphysics.ModelPartIO(mpi_input_filename)
            model_part_io.ReadModelPart(self.main_model_part)
            
            ##here we shall check that the input read has the shape we like
            aux_params = KratosMultiphysics.Parameters("{}")
            aux_params.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
            aux_params.AddValue("skin_parts",self.settings["skin_parts"])

            #here we replace the dummy elements we read with proper elements
            self.settings.AddEmptyValue("element_replace_settings")
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""{
                                                                                              "element_name":"VMS3D4N",
                                                                                              "condition_name": "MonolithicWallCondition3D"
                                                                                             }""")
            elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""{
                                                                                              "element_name":"VMS2D3N",
                                                                                              "condition_name": "MonolithicWallCondition2D"
                                                                                             }""")
            else:
                raise Exception("Domain size is neither 2 nor 3!!")

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

            ## Construct and execute the MPICommunicator
            KratosMetis.SetMPICommunicatorProcess(self.main_model_part).Execute()
            
            ## Construct and execute the Parallel fill communicator
            ParallelFillCommunicator = KratosTrilinos.ParallelFillCommunicator(self.main_model_part.GetRootModelPart())
            ParallelFillCommunicator.Execute()
            ParallelFillCommunicator.PrintDebugInfo()
            
            print(self.main_model_part)
        else:
            raise Exception("Other input options are not yet implemented.")
            
        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )

        print ("Model reading finished.")
                

    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE)
            node.AddDof(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.REACTION_Z)

        KratosMPI.mpi.world.barrier()
        
        if KratosMPI.mpi.rank == 0:
            print("DOFs for the VMS Trilinos fluid solver added correctly.")


    def Initialize(self):

        ## Construct the communicator
        self.EpetraCommunicator = KratosTrilinos.CreateCommunicator()

        ## Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        ## Creating the Trilinos convergence criteria
        self.conv_criteria = KratosTrilinos.TrilinosUPCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                               self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                               self.settings["relative_pressure_tolerance"].GetDouble(),
                                                               self.settings["absolute_pressure_tolerance"].GetDouble(),
                                                               self.EpetraCommunicator)
        
        ## Creating the Trilinos time scheme
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            if self.settings["scheme_type"].GetString() == "Bossak":
                self.time_scheme = KratosTrilinos.TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(self.settings["alpha"].GetDouble(),
                                                                                                          self.settings["move_mesh_strategy"].GetInt(),
                                                                                                          self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                                                          KratosCFD.PATCH_INDEX)
            elif self.settings["scheme_type"].GetString() == "BDF":
                self.time_scheme = KratosTrilinos.TrilinosResidualBasedPredictorCorrectorBDFScheme(self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                                                   KratosCFD.PATCH_INDEX)
        else:
            self.time_scheme = KratosTrilinos.TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(self.settings["alpha"].GetDouble(),
                                                                                                      self.settings["move_mesh_strategy"].GetInt(),
                                                                                                      self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
            

        ## Set the guess_row_size (guess about the number of zero entries) for the Trilinos builder and solver 
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3:
            guess_row_size = 20*4
        elif self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            guess_row_size = 10*3

        ## Construct the Trilinos builder and solver
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.builder_and_solver = KratosTrilinos.TrilinosBlockBuilderAndSolverPeriodic(self.EpetraCommunicator,
                                                                                           guess_row_size,
                                                                                           self.trilinos_linear_solver,
                                                                                           KratosCFD.PATCH_INDEX)
        else:
            self.builder_and_solver = KratosTrilinos.TrilinosBlockBuilderAndSolver(self.EpetraCommunicator,
                                                                                   guess_row_size,
                                                                                   self.trilinos_linear_solver)

        ## Construct the Trilinos Newton-Raphson strategy
        self.solver = KratosTrilinos.TrilinosNewtonRaphsonStrategy(self.main_model_part,
                                                                            self.time_scheme,
                                                                            self.trilinos_linear_solver,
                                                                            self.conv_criteria,
                                                                            self.builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["MoveMeshFlag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
        
        self.solver.Initialize()
        self.solver.Check()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.OSS_SWITCH, self.settings["oss_switch"].GetInt())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.M, self.settings["regularization_coef"].GetDouble())

        print ("Monolithic solver initialization finished.")

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

        #~ if self.settings["divergence_clearance_steps"].GetInt() > 0:
            #~ print("Calculating divergence-free initial condition")
            #~ ## Initialize with a Stokes solution step
            #~ try:
                #~ import KratosMultiphysics.ExternalSolversApplication as KratosExternalSolvers
                #~ smoother_type = KratosExternalSolvers.AMGCLSmoother.DAMPED_JACOBI
                #~ solver_type = KratosExternalSolvers.AMGCLIterativeSolverType.CG
                #~ gmres_size = 50
                #~ max_iter = 200
                #~ tol = 1e-7
                #~ verbosity = 0
                #~ stokes_linear_solver = KratosExternalSolvers.AMGCLSolver(smoother_type,
                                                                         #~ solver_type,
                                                                         #~ tol,
                                                                         #~ max_iter,
                                                                         #~ verbosity,
                                                                         #~ gmres_size)
            #~ except:
                #~ pPrecond = DiagonalPreconditioner()
                #~ stokes_linear_solver = BICGSTABSolver(1e-9, 5000, pPrecond)

            #~ stokes_process = KratosCFD.StokesInitializationProcess(self.main_model_part,
                                                                   #~ stokes_linear_solver,
                                                                   #~ self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                                                   #~ KratosCFD.PATCH_INDEX)
            #~ ## Copy periodic conditions to Stokes problem
            #~ stokes_process.SetConditions(self.main_model_part.Conditions)
            #~ ## Execute Stokes process
            #~ stokes_process.Execute()
            #~ stokes_process = None

            #~ for node in self.main_model_part.Nodes:
                #~ node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0, 0.0)
                #~ node.SetSolutionStepValue(KratosMultiphysics.ACCELERATION_X, 0, 0.0)
                #~ node.SetSolutionStepValue(KratosMultiphysics.ACCELERATION_Y, 0, 0.0)
                #~ node.SetSolutionStepValue(KratosMultiphysics.ACCELERATION_Z, 0, 0.0)
#~ ##                vel = node.GetSolutionStepValue(VELOCITY)
#~ ##                for i in range(0,2):
#~ ##                    node.SetSolutionStepValue(VELOCITY,i,vel)

            #~ self.settings["divergence_clearance_steps"].SetInt(0)
            #~ print("Finished divergence clearance.")


    def SolverInitialize(self):
        self.DivergenceClearance()
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
        self.DivergenceClearance()
        self.solver.Solve()

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()

    def Check(self):
        self.solver.Check()

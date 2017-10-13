##################################################################
################### distributed_include.py   #####################
##################################################################
##### KRATOS Multiphysics                                    #####
##### include file for distributed-memory simulation         #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Institute for Structural Mechanics, RUB   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
##################################################################
##################################################################
import sys
import os
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import **
*if(strcmp(GenData(Structural_Application),"1")==0)
from KratosMultiphysics.StructuralApplication import **
*endif
*if(strcmp(GenData(External_Solvers_Application),"1")==0)
from KratosMultiphysics.ExternalSolversApplication import **
*endif
*if(strcmp(GenData(MKL_Solvers_Application),"1")==0)
from KratosMultiphysics.MKLSolversApplication import **
*endif
*if(strcmp(GenData(Mortar_Application),"1")==0)
from KratosMultiphysics.MortarApplication import **
*endif
*if(strcmp(GenData(AGMG_Solver_Application),"1")==0)
from KratosMultiphysics.AGMGSolverApplication import **
*endif
*if(strcmp(GenData(Multithreaded_Solvers_Application),"1")==0)
from KratosMultiphysics.MultithreadedSolversApplication import **
*endif
from KratosMultiphysics.mpi import **
from KratosMultiphysics.MetisApplication import **
from KratosMultiphysics.DistributedBuildersApplication import **
from KratosMultiphysics.TrilinosSolversApplication import **
from KratosMultiphysics.PetscSolversApplication import **
kernel = Kernel()   #defining kernel

##################################################################
##################################################################
class Model:
    def __init__( self, problem_name, path ):
        ##################################################################
        ## DEFINE MODELPART ##############################################
        ##################################################################
        self.model_part = ModelPart("*GenData(Simulation_Name)")
        self.path = path
        self.problem_name = problem_name
        ##################################################################
        ## ADD VARIABLES #################################################
        ##################################################################
*if(strcmp(GenData(analysis_type),"static")==0)
        import structural_solver_static
        structural_solver_static.AddVariables(self.model_part)
*endif
*if(strcmp(GenData(analysis_type),"quasi-static")==0)
        import structural_solver_dynamics
        structural_solver_dynamics.AddVariables(self.model_part)
*endif
*if(strcmp(GenData(analysis_type),"dynamic")==0)
        import structural_solver_dynamics
        structural_solver_dynamics.AddVariables(self.model_part)
*endif
        set_communicator_process = SetMPICommunicatorProcess(self.model_part)
        set_communicator_process.Execute() #this is important for the model_part to read the partition correctly
        ##################################################################
        ## READ MODELPART ################################################
        ##################################################################
        #reading a model
        write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_elements = WriteConditionsFlag.WriteConditions
        #write_elements = WriteConditionsFlag.WriteElementsOnly
*if(strcmp(GenData(Output_Format),"ASCII")==0)
        post_mode = GiDPostMode.GiD_PostAscii
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
        multi_file_flag = MultiFileFlag.MultipleFiles
*else
        multi_file_flag = MultiFileFlag.SingleFile
*endif
*else
        post_mode = GiDPostMode.GiD_PostBinary
        multi_file_flag = MultiFileFlag.MultipleFiles
*endif
        self.gid_io = StructuralGidIO( self.path+self.problem_name, post_mode, multi_file_flag, write_deformed_flag, write_elements )
        self.model_part_io = ModelPartIO(self.path+self.problem_name)
        self.read_partition_file = True
        number_of_partitions = mpi.size
        if number_of_partitions > 1:
            if self.read_partition_file == True:
                dimension = 3
                verbose = 0
                synchronize_conds = True
                model_part_io = ModelPartIO(self.path+self.problem_name)
                partitioning_process = MetisDivideHeterogeneousInputProcess(model_part_io, number_of_partitions, dimension, verbose, synchronize_conds)
                partitioning_process.Execute()
            self.model_part_io = ModelPartIO(self.path+self.problem_name+'_'+str(mpi.rank))
            self.model_part_io.ReadModelPart(self.model_part)
        else:
            self.model_part_io = ModelPartIO(self.path+self.problem_name)
            self.model_part_io.ReadModelPart(self.model_part)
        mpi.world.barrier()
        self.meshWritten = False
        if( mpi.rank == 0 ):
            self.mergefile = open( self.path+self.problem_name+"_merge_results.bch", 'w' )
            self.mergefile.write("Postprocess\n")
            self.mergefile.write("mescape\n \n")
        mpi.world.barrier()
*if(strcmp(GenData(Transfer_To_Node),"1")==0)
        variable_transfer_solver = MKLPardisoSolver()
        #variable_transfer_solver.SetEchoLevel(0)
        self.variable_transfer_util = VariableTransferUtility(variable_transfer_solver)
*endif
        ## READ DEACTIVATION FILE ########################################
        self.cond_file = open(self.path+self.problem_name+".mdpa",'r' )
        self.cond_activation_flags = []
        for line in self.cond_file:
            if "//ElementAssignment" in line:
                val_set = line.split(' ')
                self.model_part.Conditions[int(val_set[1])].SetValue( ACTIVATION_LEVEL, self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL) )
                #print( "assigning ACTIVATION_LEVEL of element: " +str(int(val_set[2])) + " to Condition: " + str(int(val_set[1])) + " as " + str(self.model_part.Elements[int(val_set[2])].GetValue(ACTIVATION_LEVEL)) )
        self.cond_file.close()

        #the buffer size should be set up here after the mesh is read for the first time
        self.model_part.SetBufferSize(2)

        ##################################################################
        ## ADD DOFS ######################################################
        ##################################################################
*if(strcmp(GenData(analysis_type),"static")==0)
        structural_solver_static.AddDofs(self.model_part)
*endif
*if(strcmp(GenData(analysis_type),"quasi-static")==0)
        structural_solver_dynamics.AddDofs(self.model_part)
*endif
*if(strcmp(GenData(analysis_type),"dynamic")==0)
        structural_solver_dynamics.AddDofs(self.model_part)
*endif
        ##################################################################
        ## INITIALISE SOLVER FOR PARTICULAR SOLUTION #####################
        ##################################################################
        # define space and parallel communicator
*if(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Epetra")==0)
        self.parallel_space = TrilinosEpetraSpace()
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Tpetra")==0)
        self.parallel_space = TrilinosTpetraSpace()
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"PETSc")==0)
        self.parallel_space = PETScSpace()
*endif
        self.comm = self.parallel_space.CreateCommunicator(mpi.world)

        # defining linear solver
*if(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Epetra")==0)
        self.structure_linear_solver = TrilinosEpetraSolver()
        # Please add your preferred solver for Trilinos Epetra backend
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Tpetra")==0)
        self.structure_linear_solver = TrilinosTpetraSolver()
        # Please add your preferred solver for Trilinos Tpetra backend
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"PETSc")==0)
        self.structure_linear_solver = PETScSolver()
*endif

        # defining builder_and_solver
*if(strcmp(GenData(Builder_And_Solver),"Residual_Based_Elimination_Builder_And_Solver_Deactivation")==0)
*if(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Epetra")==0)
        self.builder_and_solver = TrilinosEpetraResidualBasedEliminationBuilderAndSolverDeactivation(self.comm, self.structure_linear_solver)
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Tpetra")==0)
        self.builder_and_solver = TrilinosTpetraResidualBasedEliminationBuilderAndSolverDeactivation(self.comm, self.structure_linear_solver)
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"PETSc")==0)
        self.builder_and_solver = PETScResidualBasedEliminationBuilderAndSolverDeactivation(self.comm, self.structure_linear_solver)
*endif
*elseif(strcmp(GenData(Builder_And_Solver),"Residual_Based_Block_Builder_And_Solver_Deactivation")==0)
*if(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Epetra")==0)
        self.builder_and_solver = TrilinosEpetraResidualBasedBlockBuilderAndSolverDeactivation(self.comm, self.structure_linear_solver)
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Tpetra")==0)
        self.builder_and_solver = TrilinosTpetraResidualBasedBlockBuilderAndSolverDeactivation(self.comm, self.structure_linear_solver)
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"PETSc")==0)
        self.builder_and_solver = PETScResidualBasedBlockBuilderAndSolverDeactivation(self.comm, self.structure_linear_solver)
*endif
*endif

        # defining time scheme
*if(strcmp(GenData(analysis_type),"static")==0)
*if(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Epetra")==0)
        self.time_scheme = TrilinosEpetraResidualBasedIncrementalUpdateStaticScheme()
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Tpetra")==0)
        self.time_scheme = TrilinosTpetraResidualBasedIncrementalUpdateStaticScheme()
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"PETSc")==0)
        self.time_scheme = PETScResidualBasedIncrementalUpdateStaticScheme()
*endif
*endif
*if(strcmp(GenData(analysis_type),"quasi-static")==0)
*set var dissipationradius(real)=GenData(Dissipation_Radius,real)
*if(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Epetra")==0)
        self.time_scheme = TrilinosEpetraResidualBasedNewmarkScheme(*dissipationradius)
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Tpetra")==0)
        self.time_scheme = TrilinosTpetraResidualBasedNewmarkScheme(*dissipationradius)
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"PETSc")==0)
        self.time_scheme = PETScResidualBasedNewmarkScheme(*dissipationradius)
*endif
*endif
*if(strcmp(GenData(analysis_type),"dynamic")==0)
*set var dissipationradius(real)=GenData(Dissipation_Radius,real)
*if(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Epetra")==0)
        self.time_scheme = TrilinosEpetraResidualBasedNewmarkScheme(*dissipationradius)
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Tpetra")==0)
        self.time_scheme = TrilinosTpetraResidualBasedNewmarkScheme(*dissipationradius)
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"PETSc")==0)
        self.time_scheme = PETScResidualBasedNewmarkScheme(*dissipationradius)
*endif
*endif

        # defining convergence criteria
*set var reltol(real)=GenData(Relative_Tolerance,real)
*set var abstol(real)=GenData(Absolute_Tolerance,real)
*if(strcmp(GenData(Convergence_Criteria),"Displacement_Criteria")==0)
        # TODO: implement distributed DisplacementCriteria
*endif
*if(strcmp(GenData(Convergence_Criteria),"Multiphase_Flow_Criteria")==0)
*if(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Epetra")==0)
        self.conv_criteria = TrilinosEpetraMultiphaseFlowCriteria(*reltol, *abstol)
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"Trilinos_Tpetra")==0)
        self.conv_criteria = TrilinosTpetraMultiphaseFlowCriteria(*reltol, *abstol)
*elseif(strcmp(GenData(Linear_Algebra_BackEnd),"PETSc")==0)
        self.conv_criteria = PETScMultiphaseFlowCriteria(*reltol, *abstol)
*endif
*endif

        # defining solving strategy flags
        self.MaxNewtonRapshonIterations = *GenData(Max_Newton_Raphson_Iterations,int)
*if(strcmp(GenData(Reactions),"1")==0)
        self.CalculateReactionFlag = True
*else
        self.CalculateReactionFlag = False
*endif
        self.ReformDofSetAtEachStep = *GenData(Reform_DofSet_At_Each_Step)
        self.MoveMeshFlag = *GenData(Move_Mesh_Flag)
        self.EchoLevel = 0

        # defining solving strategy
        import distributed_strategies
        self.solver = distributed_strategies.ResidualBasedNewtonRaphsonStrategy(self.parallel_space, \
                        self.comm, \
                        self.model_part, \
                        self.time_scheme, \
                        self.structure_linear_solver, \
                        self.conv_criteria, \
                        self.builder_and_solver, \
                        self.MaxNewtonRapshonIterations, \
                        self.CalculateReactionFlag, \
                        self.ReformDofSetAtEachStep, \
                        self.MoveMeshFlag)
        self.solver.SetEchoLevel(self.EchoLevel)

*if(strcmp(GenData(Solving_Strategy),"Residual_Based_BFGS_Strategy")==0)
        # TODO: implement distributed ResidualBasedBFGSStrategy
*endif

    def InitializeModel( self ):
        ##################################################################
        ## INITIALISE CONSTITUTIVE LAWS ##################################
        ##################################################################
        #set material parameters
        append_manual_data = False
*loop materials
*if(strcmp(MatProp(ConstitutiveLaw),"Isotropic3D")==0)
        self.model_part.Properties[*MatNum].SetValue(DENSITY, *MatProp(Density,real) )
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
        print "Linear elastic model selected"
*elseif(strcmp(MatProp(ConstitutiveLaw),"PlaneStrain")==0)
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(THICKNESS, *MatProp(Thickness,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, PlaneStrain() )
*elseif(strcmp(MatProp(ConstitutiveLaw),"PlaneStress")==0)
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(THICKNESS, *MatProp(Thickness,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, PlaneStress() )
*elseif(strcmp(MatProp(ConstitutiveLaw),"Isotropic3Dshell")==0)
        self.model_part.Properties[*MatNum].SetValue(DENSITY, *MatProp(Density,real) )
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, Isotropic3D() )
        self.model_part.Properties[*MatNum].SetValue(THICKNESS, *MatProp(Thickness,real) )
        print "Linear elastic model selected"
*elseif(strcmp(MatProp(ConstitutiveLaw),"GroutingMortar")==0)
        self.model_part.Properties[*MatNum].SetValue(DENSITY, *MatProp(Density,real) )
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(PRIMARY_HYDRATION_TIME, *MatProp(prim_hyd_time,real) )
        self.model_part.Properties[*MatNum].SetValue(PRIMARY_HYDRATION_TIME_GRADIENT, *MatProp(gradient_prim_hyd_time,real) )
        self.model_part.Properties[*MatNum].SetValue(STIFFNESS_RATIO, *MatProp(E_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, GroutingMortar() )
        print "Grouting Mortar material selected"
*elseif(strcmp(MatProp(ConstitutiveLaw),"DruckerPrager_3D")==0)
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(COHESION, *MatProp(Cohesion,real) )
        self.model_part.Properties[*MatNum].SetValue(INTERNAL_FRICTION_ANGLE, *MatProp(Friction_angle,real) )
        self.model_part.Properties[*MatNum].SetValue(ISOTROPIC_HARDENING_MODULUS, *MatProp(Isotropic_hardening_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, DruckerPrager() )
*elseif(strcmp(MatProp(ConstitutiveLaw),"DruckerPrager_PlaneStrain")==0)
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(COHESION, *MatProp(Cohesion,real) )
        self.model_part.Properties[*MatNum].SetValue(INTERNAL_FRICTION_ANGLE, *MatProp(Friction_angle,real) )
        self.model_part.Properties[*MatNum].SetValue(ISOTROPIC_HARDENING_MODULUS, *MatProp(Isotropic_hardening_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, DruckerPrager() )
*elseif(strcmp(MatProp(ConstitutiveLaw),"DruckerPrager_PlaneStress")==0)
        print("Warning: DruckerPrager_PlaneStress is not implemented")
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(COHESION, *MatProp(Cohesion,real) )
        self.model_part.Properties[*MatNum].SetValue(INTERNAL_FRICTION_ANGLE, *MatProp(Friction_angle,real) )
        self.model_part.Properties[*MatNum].SetValue(ISOTROPIC_HARDENING_MODULUS, *MatProp(Isotropic_hardening_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, DruckerPrager() )
*elseif(strcmp(MatProp(ConstitutiveLaw),"IsotropicDamage3D")==0)
        self.model_part.Properties[*MatNum].SetValue(YOUNG_MODULUS, *MatProp(Compressive_Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(CONCRETE_YOUNG_MODULUS_C, *MatProp(Compressive_Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(CONCRETE_YOUNG_MODULUS_T, *MatProp(Tensile_Young_modulus,real) )
        self.model_part.Properties[*MatNum].SetValue(FC, *MatProp(Compressive_strength,real) )
        self.model_part.Properties[*MatNum].SetValue(FT, *MatProp(Tensile_strength,real) )
        self.model_part.Properties[*MatNum].SetValue(YIELD_STRESS, *MatProp(Yield_stress,real) )
        self.model_part.Properties[*MatNum].SetValue(FRACTURE_ENERGY, *MatProp(Fracture_energy,real) )
        self.model_part.Properties[*MatNum].SetValue(POISSON_RATIO, *MatProp(Poisson_ratio,real) )
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, IsotropicDamage3DImplex() )
*elseif(strcmp(MatProp(ConstitutiveLaw),"UserDefined")==0)
        append_manual_data = True
        self.model_part.Properties[*MatNum].SetValue(CONSTITUTIVE_LAW, DummyConstitutiveLaw() )
*else
        self.model_part.Properties[*MatNum].SetValue(0)
        # this constitutive law is not yet supported by the problemtype
*if(strcmp(GenData(analysis_type),"quasi-static")==0)
*set var rayleighalpha(real)=GenData(Rayleigh_Damping_Alpha,real)
*set var rayleighbeta(real)=GenData(Rayleigh_Damping_Beta,real)
        self.model_part.Properties[*MatNum].SetValue(RAYLEIGH_DAMPING_ALPHA, *rayleighalpha )
        self.model_part.Properties[*MatNum].SetValue(RAYLEIGH_DAMPING_BETA, *rayleighbeta )
*endif
*if(strcmp(GenData(analysis_type),"dynamic")==0)
*set var rayleighalpha(real)=GenData(Rayleigh_Damping_Alpha,real)
*set var rayleighbeta(real)=GenData(Rayleigh_Damping_Beta,real)
        self.model_part.Properties[*MatNum].SetValue(RAYLEIGH_DAMPING_ALPHA, *rayleighalpha )
        self.model_part.Properties[*MatNum].SetValue(RAYLEIGH_DAMPING_BETA, *rayleighbeta )
*endif
*end materials

        ##################################################################
        ## STORE LAYER SETS ##############################################
        ##################################################################
        model_layers = __import__(self.problem_name+"_layers")
        ## ELEMENTS on layers ############################################
        self.layer_sets = model_layers.ReadLayerSets()
        ## NODES on layers ###############################################
        self.layer_nodes_sets = model_layers.ReadLayerNodesSets()
        ##################################################################
        print "layer sets stored"
        ##################################################################
        ## STORE NODES CORRECTLY FOR CONDITIONS ##########################
        ##################################################################
        self.node_groups = model_layers.ReadNodeGroups()
        print "node groups stored"
        ##################################################################
        ## ACTIVATION ####################################################
        ##################################################################
        self.deac = DeactivationUtility()
        self.deac.Initialize( self.model_part )
        if(mpi.rank == 0):
            print "activation utility initialized"
        self.model_part.Check( self.model_part.ProcessInfo )
        if(mpi.rank == 0):
            print "model successfully initialized"

    def WriteOutput( self, time ):
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
        mpi_fn_step = 0.0001
        meshname = time+mpi_fn_step**(mpi.rank+1)
        self.gid_io.InitializeMesh( meshname )
        mesh = self.model_part.GetMesh()
        #self.gid_io.WriteNodeMesh( mesh )
        self.gid_io.WriteMesh( mesh )
        if(mpi.rank == 0):
            print("mesh written...")
        self.gid_io.FinalizeMesh()
*else
        if( self.meshWritten == False ):
            meshname = time+0.0001**mpi.rank
            self.gid_io.InitializeMesh( meshname )
            mesh = self.model_part.GetMesh()
            self.gid_io.WriteMesh( mesh )
            self.meshWritten = True
            self.gid_io.FinalizeMesh()
*endif
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
        self.gid_io.InitializeResults( meshname, mesh )
*else
        self.gid_io.InitializeResults( 0.0 + mpi_fn_step**(mpi.rank+1), mesh )
*endif
*if(strcmp(GenData(Displacements),"1")==0)
        if(mpi.rank == 0):
            print("write nodal displacements")
        self.gid_io.WriteNodalResults(DISPLACEMENT, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
*endif
*if(strcmp(GenData(Reactions),"1")==0)
        self.gid_io.WriteNodalResults(REACTION, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
*endif
*if(strcmp(GenData(PK2_Stresses),"1")==0)
        self.gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR, mesh, time)
*endif
*if(strcmp(GenData(Green_Lagrange_Strains),"1")==0)
        self.gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_STRAIN_TENSOR, self.model_part, time)
        self.gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR, self.model_part, time)
*endif
*if(strcmp(GenData(Stresses),"1")==0)
*if(strcmp(GenData(Transfer_To_Node),"1")==0)
        self.variable_transfer_util.TransferVariablesToNodes(self.model_part, STRESSES)
        self.gid_io.WriteNodalResults(STRESSES, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
*else
        self.gid_io.PrintOnGaussPoints(STRESSES, self.model_part, time)
*endif
*endif
*if(strcmp(GenData(Strains),"1")==0)
*if(strcmp(GenData(Transfer_To_Node),"1")==0)
        self.variable_transfer_util.TransferVariablesToNodes(self.model_part, STRAIN)
        self.gid_io.WriteNodalResults(STRAIN, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
*else
        self.gid_io.PrintOnGaussPoints(STRAIN, self.model_part, time)
*endif
*endif
*if(strcmp(GenData(Jack_Forces),"1")==0)
        self.gid_io.PrintOnGaussPoints(JACK_FORCE, self.model_part, time)
*endif
*if(strcmp(GenData(Insitu_Stress),"1")==0)
*if(strcmp(GenData(Transfer_To_Node),"1")==0)
        self.variable_transfer_util.TransferVariablesToNodes(self.model_part, PRESTRESS)
        self.gid_io.WriteNodalResults(PRESTRESS, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
*else
        self.gid_io.PrintOnGaussPoints(PRESTRESS, self.model_part, time)
*endif
*endif
*if(strcmp(GenData(Plastic_Strains),"1")==0)
*if(strcmp(GenData(Transfer_To_Node),"1")==0)
        self.variable_transfer_util.TransferVariablesToNodes(self.model_part, PLASTIC_STRAIN_VECTOR)
        self.gid_io.WriteNodalResults(PLASTIC_STRAIN_VECTOR, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
*else
        self.gid_io.PrintOnGaussPoints(PLASTIC_STRAIN_VECTOR, self.model_part, time)
*endif
*if(strcmp(GenData(Transfer_To_Node),"1")==0)
        self.variable_transfer_util.TransferVariablesToNodes(self.model_part, PLASTICITY_INDICATOR)
        self.gid_io.WriteNodalResults(PLASTICITY_INDICATOR, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
*else
        self.gid_io.PrintOnGaussPoints(PLASTICITY_INDICATOR, self.model_part, time)
*endif
*endif
*if(strcmp(GenData(Internal_Variables),"1")==0)
        for i in range(0, 13):
            self.gid_io.PrintOnGaussPoints(INTERNAL_VARIABLES, self.model_part, time, i)
*endif
*if(strcmp(GenData(Air_Pressure),"1")==0)
*if(strcmp(GenData(Transfer_To_Node),"1")==0)
*else
        self.gid_io.PrintOnGaussPoints(AIR_PRESSURE, self.model_part, time)
*endif
*endif
*if(strcmp(GenData(Water_Pressure),"1")==0)
*if(strcmp(GenData(Transfer_To_Node),"1")==0)
        self.variable_transfer_util.TransferVariablesToNodes(self.model_part, WATER_PRESSURE)
        self.gid_io.WriteNodalResults(WATER_PRESSURE, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
        self.variable_transfer_util.TransferVariablesToNodes(self.model_part, EXCESS_PORE_WATER_PRESSURE)
        self.gid_io.WriteNodalResults(EXCESS_PORE_WATER_PRESSURE, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
*else
        self.gid_io.PrintOnGaussPoints(WATER_PRESSURE, self.model_part, time)
        self.gid_io.PrintOnGaussPoints(EXCESS_PORE_WATER_PRESSURE, self.model_part, time)
*endif
*endif
*if(strcmp(GenData(Saturation),"1")==0)
        self.gid_io.PrintOnGaussPoints(SATURATION, self.model_part, time)
*endif
        self.gid_io.WriteNodalResults(PARTITION_INDEX, self.model_part.GetCommunicator().InterfaceMesh().Nodes, time, 0)
*if(strcmp(GenData(New_mesh_for_each_step),"1")==0)
        self.gid_io.FinalizeResults()
*endif
        if(mpi.rank == 0):
            meshname = time+mpi_fn_step
            self.mergefile.write("Files Read "+self.path+self.problem_name+"_"+str(meshname)+".post.bin\n")
            self.mergefile.write("mescape\n")
            for rank in range(2, mpi.size+1 ):
                meshname = time+mpi_fn_step**(rank)
                self.mergefile.write("Files Add "+self.path+self.problem_name+"_"+str(meshname)+".post.bin\n")
                self.mergefile.write("mescape\n")
            self.mergefile.write("Files SaveAll BinMeshesSets\n")
            self.mergefile.write(self.path+"merged_"+self.problem_name+"_"+ str(time)+".bin\n")
            self.mergefile.write("mescape\n")
            self.mergefile.write("\n")

    def FinalizeModel( self ):
        self.gid_io.CloseResultFile()
        if(mpi.rank == 0):
            self.mergefile.close()

    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
        self.deac.Reactivate( self.model_part, from_reac, to_reac )
        self.deac.Deactivate( self.model_part, from_deac, to_deac )
        self.model_part.CloneTimeStep(time)
        self.solver.Solve()
##################################################################

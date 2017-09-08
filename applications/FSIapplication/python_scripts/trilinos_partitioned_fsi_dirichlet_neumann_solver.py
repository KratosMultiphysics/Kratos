from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import utilities
import math
import python_solvers_wrapper_fluid            # Import the fluid Python solvers wrapper

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.MetisApplication as KratosMetis
import KratosMultiphysics.TrilinosApplication as KratosTrilinos
import KratosMultiphysics.ALEApplication as KratosALE
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Python files import 
import convergence_accelerator_factory          # Convergence accelerator factory
import trilinos_partitioned_fsi_base_solver     # Base class file


def CreateSolver(structure_main_model_part, fluid_main_model_part, project_parameters):
    return TrilinosPartitionedFSIDirichletNeumannSolver(structure_main_model_part, fluid_main_model_part, project_parameters)


class TrilinosPartitionedFSIDirichletNeumannSolver(trilinos_partitioned_fsi_base_solver.TrilinosPartitionedFSIBaseSolver):
    def __init__(self, structure_main_model_part, fluid_main_model_part, project_parameters):

        if (KratosMPI.mpi.rank == 0) : print("*** Trilinos partitioned Dirichlet-Neumann FSI solver construction starts...")

        super(TrilinosPartitionedFSIDirichletNeumannSolver, self).__init__(structure_main_model_part, fluid_main_model_part, project_parameters)

        if (KratosMPI.mpi.rank == 0) : print("*** Trilinos partitioned Dirichlet-Neumann FSI solver construction finished.")


    def Initialize(self):

        # Set the Trilinos space
        self.trilinos_space = KratosTrilinos.TrilinosSparseSpace()

        # Set the Epetra communicator
        self.epetra_communicator = KratosTrilinos.CreateCommunicator()

        # Construct the coupling partitioned strategy
        coupling_utility_parameters = self.settings["coupling_solver_settings"]["coupling_strategy"]
        self.coupling_utility = convergence_accelerator_factory.CreateTrilinosConvergenceAccelerator(self._GetFluidInterfaceSubmodelPart(),
                                                                                                     self.epetra_communicator, 
                                                                                                     coupling_utility_parameters)

        # Get the domain size
        self.domain_size = self._GetDomainSize()

        # Get the nodal update FSI utilities
        self.nodal_update_utilities = self._GetNodalUpdateUtilities()

        # Get the partitioned FSI utilities
        self.trilinos_partitioned_fsi_utilities = self._GetPartitionedFSIUtilities()

        # Python structure solver initialization
        self.structure_solver.Initialize()

        # Ensure that the fluid reaction fluxes are computed if D-N scheme is considered
        if self.fluid_solver.settings["compute_reactions"].GetBool() == False:
            self.fluid_solver.settings["compute_reactions"].SetBool(True)

        # Python fluid solver initialization
        self.fluid_solver.Initialize()

        # Python mesh solver initialization
        self.mesh_solver.Initialize()

        # Initialize the Dirichlet-Neumann interface
        self._InitializeDirichletNeumannInterface()

        # Construct the interface mapper
        self._SetUpMapper()

        # Set the Neumann B.C. in the structure interface
        self._SetStructureNeumannCondition() #TODO change when the interface is able to correctly transfer distributed forces

        # Initialize the iteration value vector
        self._InitializeIterationValueVector()


        # Compute the fluid domain NODAL_AREA values (required as weight in the residual norm computation)
        KratosMultiphysics.CalculateNodalAreaProcess(self.fluid_solver.GetComputingModelPart(), self.domain_size).Execute()

        # Strategies initialization
        super(TrilinosPartitionedFSIDirichletNeumannSolver, self).Initialize()


    def Solve(self):

        ## Solvers initialization
        self.InitializeSolutionStep()

        ## Solvers predict
        self.Predict()

        ## Compute mesh prediction ##
        if (self.double_faced_structure):
            self._ComputeMeshPredictionDoubleFaced()
        else:
            self._ComputeMeshPredictionSingleFaced()

        ## Non-Linear interface coupling iteration ##
        for nl_it in range(1,self.max_nl_it+1):

            if (KratosMPI.mpi.rank == 0) : print("     NL-ITERATION ",nl_it,"STARTS.")
            self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it
            self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it

            self.coupling_utility.InitializeNonLinearIteration()

            if (KratosMPI.mpi.rank == 0) : print("     Residual computation starts...")

            # Solve the mesh problem as well as the fluid problem
            self._SolveMeshAndFluid()

            # Solve the structure problem and computes the displacement residual
            if (self.double_faced_structure):
                self._SolveStructureDoubleFaced()
                dis_residual = self._ComputeDisplacementResidualDoubleFaced()
            else:
                self._SolveStructureSingleFaced()
                dis_residual = self._ComputeDisplacementResidualSingleFaced()

            # Residual computation
            nl_res_norm = self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM]
            interface_dofs = self.trilinos_partitioned_fsi_utilities.GetInterfaceResidualSize(self._GetFluidInterfaceSubmodelPart())

            # Check convergence
            if nl_res_norm/math.sqrt(interface_dofs) < self.nl_tol:
                if (KratosMPI.mpi.rank == 0):
                    print("     NON-LINEAR ITERATION CONVERGENCE ACHIEVED")
                    print("     Total non-linear iterations: ",nl_it," |res|/sqrt(Ndofs) = ",nl_res_norm/math.sqrt(interface_dofs))

                break
            else:
                # If convergence is not achieved, perform the correction of the prediction
                if (KratosMPI.mpi.rank == 0):
                    print("     Residual computation finished. |res|/sqrt(Ndofs) =", nl_res_norm/math.sqrt(interface_dofs))
                    print("     Performing non-linear iteration ",nl_it," correction.")

                self.coupling_utility.UpdateSolution(dis_residual.GetReference(), self.iteration_value.GetReference())

                if ((nl_it == self.max_nl_it) and (KratosMPI.mpi.rank == 0)):
                    print("***********************************************************")
                    print("***********************************************************")
                    print("       NON-LINEAR ITERATION CONVERGENCE NOT ACHIEVED       ")
                    print("***********************************************************")
                    print("***********************************************************")

            self.coupling_utility.FinalizeNonLinearIteration()

        ## Compute the mesh residual as final testing (it is expected to be 0)
        self.trilinos_partitioned_fsi_utilities.ComputeFluidInterfaceMeshVelocityResidualNorm(self._GetFluidInterfaceSubmodelPart())
        mesh_res_norm = self.fluid_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.FSI_INTERFACE_MESH_RESIDUAL_NORM)
        if (KratosMPI.mpi.rank == 0):
            print("     NL residual norm: ", nl_res_norm)
            print("     Mesh residual norm: ", mesh_res_norm)

        ## Finalize solution step
        self.fluid_solver.FinalizeSolutionStep()
        self.structure_solver.FinalizeSolutionStep()
        self.coupling_utility.FinalizeSolutionStep()


    def _InitializeIterationValueVector(self):
        # Note that the FSI problem is defined in terms of the fluid interface
        # Besides, one of the two interfaces is considered for the residual vector computation
        # in case a shell structure is analised. This is due to the fact that the same mesh_solver
        # displacement is map to both fluid sides.

        # Initialize the iteration value for the residual computation
        self.iteration_value = self.trilinos_space.CreateEmptyVectorPointer(self.epetra_communicator)
        self.trilinos_partitioned_fsi_utilities.SetUpInterfaceVector(self._GetFluidInterfaceSubmodelPart(),self.iteration_value)


    def _InitializeDirichletNeumannInterface(self):
        # Initialize Dirichlet fluid interface
        num_fl_interfaces = self.settings["coupling_solver_settings"]["fluid_interfaces_list"].size()
        for fl_interface_id in range(num_fl_interfaces):
            fl_interface_name = self.settings["coupling_solver_settings"]["fluid_interfaces_list"][fl_interface_id].GetString()
            fl_interface_submodelpart = self.fluid_solver.main_model_part.GetSubModelPart(fl_interface_name)

            # Fix the VELOCITY, MESH_DISPLACEMENT and MESH_VELOCITY variables in all the fluid interface submodelparts
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VELOCITY_X, True, fl_interface_submodelpart.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VELOCITY_Y, True, fl_interface_submodelpart.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.MESH_VELOCITY_X, True, fl_interface_submodelpart.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.MESH_VELOCITY_Y, True, fl_interface_submodelpart.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.MESH_DISPLACEMENT_X, True, fl_interface_submodelpart.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.MESH_DISPLACEMENT_Y, True, fl_interface_submodelpart.Nodes)
            if (self.domain_size == 3):
                KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VELOCITY_Z, True, fl_interface_submodelpart.Nodes)
                KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.MESH_VELOCITY_Z, True, fl_interface_submodelpart.Nodes)
                KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.MESH_DISPLACEMENT_Z, True, fl_interface_submodelpart.Nodes)

            # Set the interface flag
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, fl_interface_submodelpart.Nodes)

        # Initialize Neumann structure interface
        num_str_interfaces = self.settings["coupling_solver_settings"]["structure_interfaces_list"].size()
        for str_interface_id in range(num_str_interfaces):
            str_interface_name = self.settings["coupling_solver_settings"]["structure_interfaces_list"][str_interface_id].GetString()
            str_interface_submodelpart = self.structure_solver.main_model_part.GetSubModelPart(str_interface_name)

            # Set the interface flag
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, str_interface_submodelpart.Nodes)


    def _SolveMeshAndFluid(self):

        # Set the iteration_value displacement as MESH_DISPLACEMENT
        num_fl_interfaces = self.settings["coupling_solver_settings"]["fluid_interfaces_list"].size()
        for fl_interface_id in range(num_fl_interfaces):
            fl_interface_name = self.settings["coupling_solver_settings"]["fluid_interfaces_list"][fl_interface_id].GetString()
            fl_interface_submodelpart = self.fluid_solver.main_model_part.GetSubModelPart(fl_interface_name)
            self.trilinos_partitioned_fsi_utilities.UpdateInterfaceValues(fl_interface_submodelpart,
                                                                          KratosMultiphysics.MESH_DISPLACEMENT,
                                                                          self.iteration_value.GetReference())

        # Solve the mesh problem (or moves the interface nodes)
        if (self.solve_mesh_at_each_iteration == True):
            self.mesh_solver.Solve()
        else:
            self.mesh_solver.MoveMesh()

        # Update MESH_VELOCITY and MESH_ACCELERATION with Newmark formulas
        self.nodal_update_utilities.UpdateMeshTimeDerivatives(self.fluid_main_model_part, self.time_step)

        # Impose the structure MESH_VELOCITY and MESH_ACCELERATION in the fluid interface VELOCITY and ACCELERATION
        self.nodal_update_utilities.SetMeshTimeDerivativesOnInterface(self._GetFluidInterfaceSubmodelPart())

        # Solve fluid problem
        self.fluid_solver.SolveSolutionStep()


    def _SolveStructureSingleFaced(self):

        # Transfer fluid reaction to solid interface
        self.interface_mapper.Map(KratosMultiphysics.REACTION,
                                  KratosStructural.POINT_LOAD,
                                  KratosMapping.MapperFactory.SWAP_SIGN |
                                  KratosMapping.MapperFactory.CONSERVATIVE)

        # Solve the structure problem
        self.structure_solver.SolveSolutionStep()


    # TODO: UPDATE THIS ONCE THE DOUBLE SIDED SURFACES MAPPING IS AVAILABLE IN MAPPING APPLICATION
    def _SolveStructureDoubleFaced(self):

        # Transfer fluid reaction from both sides to the structure interface
        keep_sign = False
        distribute_load = True
        self.interface_mapper.PositiveFluidToStructure_VectorMap(KratosMultiphysics.REACTION,
                                                                 KratosFSI.POSITIVE_MAPPED_VECTOR_VARIABLE,
                                                                 keep_sign,
                                                                 distribute_load)
        self.interface_mapper.NegativeFluidToStructure_VectorMap(KratosMultiphysics.REACTION,
                                                                 KratosFSI.NEGATIVE_MAPPED_VECTOR_VARIABLE,
                                                                 keep_sign,
                                                                 distribute_load)

        # Add the two faces contributions to the POINT_LOAD variable
        # TODO: Add this to the variables utils
        for node in self._GetStructureInterfaceSubmodelPart().Nodes:
            pos_face_force = node.GetSolutionStepValue(KratosFSI.POSITIVE_MAPPED_VECTOR_VARIABLE)
            neg_face_force = node.GetSolutionStepValue(KratosFSI.NEGATIVE_MAPPED_VECTOR_VARIABLE)
            node.SetSolutionStepValue(KratosStructural.POINT_LOAD, 0, pos_face_force+neg_face_force)

        # Solve the structure problem
        self.structure_solver.SolveSolutionStep()


    def _ComputeDisplacementResidualSingleFaced(self):
        # Project the structure displacement onto the fluid interface
        self.interface_mapper.InverseMap(KratosMultiphysics.VECTOR_PROJECTED,
                                         KratosMultiphysics.DISPLACEMENT)

        # Compute the fluid interface residual vector by means of the VECTOR_PROJECTED variable
        # Besides, its norm is stored within the ProcessInfo.
        disp_residual = self.trilinos_space.CreateEmptyVectorPointer(self.epetra_communicator)
        self.trilinos_partitioned_fsi_utilities.SetUpInterfaceVector(self._GetFluidInterfaceSubmodelPart(), disp_residual)
        self.trilinos_partitioned_fsi_utilities.ComputeInterfaceVectorResidual(self._GetFluidInterfaceSubmodelPart(),
                                                                               KratosMultiphysics.MESH_DISPLACEMENT,
                                                                               KratosMultiphysics.VECTOR_PROJECTED,
                                                                               disp_residual.GetReference())

        return disp_residual


    # TODO: UPDATE THIS ONCE THE DOUBLE SIDED SURFACES MAPPING IS AVAILABLE IN MAPPING APPLICATION
    def _ComputeDisplacementResidualDoubleFaced(self):
        # Project the structure velocity onto the fluid interface
        keep_sign = True
        distribute_load = False
        #TODO: One of these mappings is not needed. At the moment both are performed to ensure that
        # the _GetFluidInterfaceSubmodelPart(), which is the one used to compute the interface residual,
        # gets the structure DISPLACEMENT. Think a way to properly identify the reference fluid interface.
        self.interface_mapper.StructureToPositiveFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                                 KratosMultiphysics.VECTOR_PROJECTED,
                                                                 keep_sign,
                                                                 distribute_load)
        self.interface_mapper.StructureToNegativeFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                                 KratosMultiphysics.VECTOR_PROJECTED,
                                                                 keep_sign,
                                                                 distribute_load)

        # Compute the fluid interface residual vector by means of the VECTOR_PROJECTED variable
        # Besides, its norm is stored within the ProcessInfo.
        disp_residual = self.trilinos_space.CreateEmptyVectorPointer(self.epetra_communicator)
        self.trilinos_partitioned_fsi_utilities.SetUpInterfaceVector(self._GetFluidInterfaceSubmodelPart(), disp_residual)
        self.trilinos_partitioned_fsi_utilities.ComputeInterfaceVectorResidual(self._GetFluidInterfaceSubmodelPart(),
                                                                               KratosMultiphysics.MESH_DISPLACEMENT,
                                                                               KratosMultiphysics.VECTOR_PROJECTED,
                                                                               disp_residual.GetReference())

        return disp_residual

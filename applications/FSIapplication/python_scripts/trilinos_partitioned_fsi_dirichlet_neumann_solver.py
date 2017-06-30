from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import utilities
# import NonConformant_OneSideMap                # Import non-conformant mapper #TODO: USE THE NEW MAPPER
import python_solvers_wrapper_fluid            # Import the fluid Python solvers wrapper

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.MetisApplication as KratosMetis
import KratosMultiphysics.TrilinosApplication as KratosTrilinos
import KratosMultiphysics.ALEApplication as KratosALE
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
import trilinos_partitioned_fsi_base_solver


def CreateSolver(structure_main_model_part, fluid_main_model_part, project_parameters):
    return TrilinosPartitionedFSIDirichletNeumannSolver(structure_main_model_part, fluid_main_model_part, project_parameters)


class TrilinosPartitionedFSIDirichletNeumannSolver(trilinos_partitioned_fsi_base_solver.TrilinosPartitionedFSIBaseSolver):
    def __init__(self, structure_main_model_part, fluid_main_model_part, project_parameters):

        print("*** Trilinos partitioned Dirichlet-Neumann FSI solver construction starts...")
        super(TrilinosPartitionedFSIDirichletNeumannSolver, self).__init__(structure_main_model_part, fluid_main_model_part, project_parameters)
        print("*** Trilinos partitioned Dirichlet-Neumann FSI solver construction finished.")


    def Initialize(self):

        # Set the Epetra communicator
        self.epetra_communicator = KratosTrilinos.CreateCommunicator()

        # Get the domain size
        self.domain_size = self._GetDomainSize()

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

        # Construct the interface mapper
        self._SetUpMapper()

        # Set the Neumann B.C. in the structure interface
        self._SetStructureNeumannCondition() #TODO change when the interface is able to correctly transfer distributed forces

        # Initialize the iteration value vector
        self._InitializeIterationValueVector()

        # Initialize the Dirichlet-Neumann interface
        self._InitializeDirichletNeumannInterface()

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

            print("     NL-ITERATION ",nl_it,"STARTS.")
            self.fluid_solver.main_model_part.ProcessInfo[KratosFSI.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it
            self.structure_solver.main_model_part.ProcessInfo[KratosFSI.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it

            self.coupling_utility.InitializeNonLinearIteration()

            print("     Residual computation starts...")
            # Sets self.iteration_value as MESH_DISPLACEMENT
            # Solves the mesh problem
            # Sets self.iteration_value derivatives as VELOCITY, ACCELERATION and MESH_VELOCITY
            # Solves the fluid problem and computes the nodal fluxes (REACTION)
            self._SolveMeshAndFluid()

            # Sets the structure POINT_LOAD as the fluid interface REACTION
            # Solves the structure problem

            if (self.double_faced_structure):
                self._SolveStructureDoubleFaced()
                dis_residual = self._ComputeDisplacementResidualDoubleFaced()
            else:
                self._SolveStructureSingleFaced()
                dis_residual = self._ComputeDisplacementResidualSingleFaced()

            # Residual computation
            nl_res_norm = self.fluid_solver.main_model_part.ProcessInfo[KratosFSI.FSI_INTERFACE_RESIDUAL_NORM]
            interface_dofs = self.trilinos_partitioned_fsi_utilities.GetInterfaceResidualSize(self._GetFluidInterfaceSubmodelPart())

            # Check convergence
            if nl_res_norm/math.sqrt(interface_dofs) < self.nl_tol:
                print("     NON-LINEAR ITERATION CONVERGENCE ACHIEVED")
                print("     Total non-linear iterations: ",nl_it," |res|/sqrt(Ndofs) = ",nl_res_norm/math.sqrt(interface_dofs))
                break
            else:
                # If convergence is not achieved, perform the correction of the prediction
                print("     Residual computation finished. |res|/sqrt(Ndofs) =", nl_res_norm/math.sqrt(interface_dofs))
                print("     Performing non-linear iteration ",nl_it," correction.")
                # self.coupling_utility.UpdateSolution(vel_residual, self.iteration_value)
                self.coupling_utility.UpdateSolution(dis_residual, self.iteration_value)
                self.coupling_utility.FinalizeNonLinearIteration()

        ## Compute the mesh residual as final testing (it is expected to be 0)
        self.trilinos_partitioned_fsi_utilities.ComputeFluidInterfaceMeshVelocityResidualNorm(self._GetFluidInterfaceSubmodelPart())
        mesh_res_norm = self.fluid_solver.main_model_part.ProcessInfo.GetValue(KratosFSI.FSI_INTERFACE_MESH_RESIDUAL_NORM)
        print("     NL residual norm: ", nl_res_norm)
        print("     Mesh residual norm: ", mesh_res_norm)

        ## Finalize solution step
        self.fluid_solver.SolverFinalizeSolutionStep()
        self.structure_solver.FinalizeSolutionStep()
        self.coupling_utility.FinalizeSolutionStep()


    def _InitializeIterationValueVector(self):
        # Note that the FSI problem is defined in terms of the fluid interface
        # Besides, one of the two interfaces is considered for the residual vector computation
        # in case a shell structure is analised. This is due to the fact that the same mesh_solver
        # displacement is map to both fluid sides.

        # Initialize the iteration value for the residual computation
        fluid_interface_residual_size = self.trilinos_partitioned_fsi_utilities.GetInterfaceResidualSize(self._GetFluidInterfaceSubmodelPart())
        self.iteration_value = KratosMultiphysics.Vector(fluid_interface_residual_size)
        for i in range(0,fluid_interface_residual_size):
            self.iteration_value[i] = 0.0


    def _InitializeDirichletNeumannInterface(self):
        # Fix the VELOCITY, MESH_DISPLACEMENT and MESH_VELOCITY variables in all the fluid interface submodelparts
        # TODO: This method allows to supress the set_interface_process
        num_fl_interfaces = self.settings["coupling_solver_settings"]["fluid_interfaces_list"].size()
        for fl_interface_id in range(num_fl_interfaces):
            fl_interface_name = self.settings["coupling_solver_settings"]["fluid_interfaces_list"][fl_interface_id].GetString()
            fl_interface_submodelpart = self.fluid_solver.main_model_part.GetSubModelPart(fl_interface_name)

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



    def _SolveMeshAndFluid(self):

        # Set the mesh displacement as the iteration_value displacement
        num_fl_interfaces = self.settings["coupling_solver_settings"]["fluid_interfaces_list"].size()
        for fl_interface_id in range(num_fl_interfaces):
            fl_interface_name = self.settings["coupling_solver_settings"]["fluid_interfaces_list"][fl_interface_id].GetString()
            fl_interface_submodelpart = self.fluid_solver.main_model_part.GetSubModelPart(fl_interface_name)
            self.trilinos_partitioned_fsi_utilities.SetInterfaceVectorVariable(fl_interface_submodelpart, KratosMultiphysics.MESH_DISPLACEMENT, self.iteration_value)

        # Solve the mesh problem (or moves the interface nodes)
        if (self.solve_mesh_at_each_iteration == True):
            self.mesh_solver.Solve()
        else:
            self.mesh_solver.MoveNodes()

        # Fluid domain velocities imposition
        # Compute the velocity associated to the iteration_value displacement and set it to VELOCITY and MESH_VELOCITY
        # Note that the VELOCITY and the MESH_VELOCITY values only coincide if the same time schemes are used
        # Currently, the mesh solver only includes the BDF2 so the MESH_VELOCITY values are forced to be the fluid ones
        self._ComputeCorrectedInterfaceDisplacementDerivatives()

        # Solve fluid problem
        self.fluid_solver.SolverSolveSolutionStep()


    def _SolveStructureSingleFaced(self):

        # Transfer fluid reaction to solid interface
        keep_sign = False
        distribute_load = True
        self.interface_mapper.FluidToStructure_VectorMap(KratosMultiphysics.REACTION,
                                                         KratosStructural.POINT_LOAD,
                                                         keep_sign,
                                                         distribute_load)

        # Solve the structure problem
        self.structure_solver.SolveSolutionStep()


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
        # Project the structure velocity onto the fluid interface
        keep_sign = True
        distribute_load = False
        self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                         KratosFSI.VECTOR_PROJECTED,
                                                         keep_sign,
                                                         distribute_load)

        # Compute the fluid interface residual vector by means of the VECTOR_PROJECTED variable
        # Besides, its norm is stored within the ProcessInfo.
        disp_residual = KratosMultiphysics.Vector(self.trilinos_partitioned_fsi_utilities.GetInterfaceResidualSize(self._GetFluidInterfaceSubmodelPart()))
        self.trilinos_partitioned_fsi_utilities.ComputeInterfaceVectorResidual(self._GetFluidInterfaceSubmodelPart(), KratosMultiphysics.MESH_DISPLACEMENT, KratosFSI.VECTOR_PROJECTED, disp_residual)

        return disp_residual


    def _ComputeDisplacementResidualDoubleFaced(self):
        # Project the structure velocity onto the fluid interface
        keep_sign = True
        distribute_load = False
        #TODO: One of these mappings is not needed. At the moment both are performed to ensure that
        # the _GetFluidInterfaceSubmodelPart(), which is the one used to compute the interface residual,
        # gets the structure DISPLACEMENT. Think a way to properly identify the reference fluid interface.
        self.interface_mapper.StructureToPositiveFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                                 KratosFSI.VECTOR_PROJECTED,
                                                                 keep_sign,
                                                                 distribute_load)
        self.interface_mapper.StructureToNegativeFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                                 KratosFSI.VECTOR_PROJECTED,
                                                                 keep_sign,
                                                                 distribute_load)

        # Compute the fluid interface residual vector by means of the VECTOR_PROJECTED variable
        # Besides, its norm is stored within the ProcessInfo.
        disp_residual = KratosMultiphysics.Vector(self.trilinos_partitioned_fsi_utilities.GetInterfaceResidualSize(self._GetFluidInterfaceSubmodelPart()))
        self.trilinos_partitioned_fsi_utilities.ComputeInterfaceVectorResidual(self._GetFluidInterfaceSubmodelPart(), KratosMultiphysics.MESH_DISPLACEMENT, KratosFSI.VECTOR_PROJECTED, disp_residual)

        return disp_residual


    ### INTERFACE MOVEMENT UTILITY ###
    # Function to update the velocity and acceleration according to the displacement in self.iteration_value.
    # Note that at the moment only the Bossak scheme is considered
    def _ComputeCorrectedInterfaceDisplacementDerivatives(self):

        # Bossak parameters
        alpha = self.settings["fluid_solver_settings"]["alpha"].GetDouble()
        gamma = 0.5*(1-2*alpha)
        beta = ((1-alpha)**2)/4

        i = 0
        if (self.domain_size == 2):
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                u_n = node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT,1)
                v_n = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,1)
                a_n = node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION,1)

                u_n1 = KratosMultiphysics.Vector(3)
                u_n1[0] = self.iteration_value[i]
                u_n1[1] = self.iteration_value[i+1]
                u_n1[2] = 0.0
                i+=2

                a_n1 = KratosMultiphysics.Vector(3)
                a_n1[0] = (u_n1[0] - u_n[0] - self.time_step*v_n[0] - (self.time_step**2)*(0.5-beta+beta*alpha)*a_n[0])/((self.time_step**2)*beta*(1-alpha))
                a_n1[1] = (u_n1[1] - u_n[1] - self.time_step*v_n[1] - (self.time_step**2)*(0.5-beta+beta*alpha)*a_n[1])/((self.time_step**2)*beta*(1-alpha))
                a_n1[2] = 0.0

                v_n1 = KratosMultiphysics.Vector(3)
                v_n1[0] = v_n[0] + self.time_step*(1-gamma)*a_n[0] + self.time_step*gamma*(1-alpha)*a_n1[0] + self.time_step*gamma*alpha*a_n[0]
                v_n1[1] = v_n[1] + self.time_step*(1-gamma)*a_n[1] + self.time_step*gamma*(1-alpha)*a_n1[1] + self.time_step*gamma*alpha*a_n[1]
                v_n1[2] = 0.0

                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0,v_n1)
                node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY,0,v_n1)
                node.SetSolutionStepValue(KratosMultiphysics.ACCELERATION,0,a_n1)

        else:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                u_n = node.GetSolutionStepValue(KratosMultiphysics.MESH_DISPLACEMENT,1)
                v_n = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,1)
                a_n = node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION,1)

                u_n1 = KratosMultiphysics.Vector(3)
                u_n1[0] = self.iteration_value[i]
                u_n1[1] = self.iteration_value[i+1]
                u_n1[2] = self.iteration_value[i+2]
                i+=3

                a_n1 = KratosMultiphysics.Vector(3)
                a_n1[0] = (u_n1[0] - u_n[0] - self.time_step*v_n[0] - (self.time_step**2)*(0.5-beta+beta*alpha)*a_n[0])/((self.time_step**2)*beta*(1-alpha))
                a_n1[1] = (u_n1[1] - u_n[1] - self.time_step*v_n[1] - (self.time_step**2)*(0.5-beta+beta*alpha)*a_n[1])/((self.time_step**2)*beta*(1-alpha))
                a_n1[2] = (u_n1[2] - u_n[2] - self.time_step*v_n[2] - (self.time_step**2)*(0.5-beta+beta*alpha)*a_n[2])/((self.time_step**2)*beta*(1-alpha))

                v_n1 = KratosMultiphysics.Vector(3)
                v_n1[0] = v_n[0] + self.time_step*(1-gamma)*a_n[0] + self.time_step*gamma*(1-alpha)*a_n1[0] + self.time_step*gamma*alpha*a_n[0]
                v_n1[1] = v_n[1] + self.time_step*(1-gamma)*a_n[1] + self.time_step*gamma*(1-alpha)*a_n1[1] + self.time_step*gamma*alpha*a_n[1]
                v_n1[2] = v_n[2] + self.time_step*(1-gamma)*a_n[2] + self.time_step*gamma*(1-alpha)*a_n1[2] + self.time_step*gamma*alpha*a_n[2]

                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0,v_n1)
                node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY,0,v_n1)
                node.SetSolutionStepValue(KratosMultiphysics.ACCELERATION,0,a_n1)

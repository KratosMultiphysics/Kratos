from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import utilities
import NonConformant_OneSideMap                # Import non-conformant mapper
import python_solvers_wrapper_fluid            # Import the fluid Python solvers wrapper

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as KratosALE
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
import partitioned_fsi_base_solver


def CreateSolver(structure_main_model_part, fluid_main_model_part, project_parameters):
    return PartitionedFSIDirichletNeumannSolver(structure_main_model_part, fluid_main_model_part, project_parameters)


class PartitionedFSIDirichletNeumannSolver(partitioned_fsi_base_solver.PartitionedFSIBaseSolver):
    def __init__(self, structure_main_model_part, fluid_main_model_part, project_parameters):

        print("*** Partitioned Dirichlet-Neumann FSI solver construction starts...")
        super(PartitionedFSIDirichletNeumannSolver, self).__init__(structure_main_model_part, fluid_main_model_part, project_parameters)
        print("*** Partitioned Dirichlet-Neumann FSI solver construction finished.")


    def Initialize(self):

        # Get the domain size
        self.domain_size = self._GetDomainSize()

        # Get the partitioned FSI utilities
        self.partitioned_fsi_utilities = self._GetPartitionedFSIUtilities()

        # Python structure solver initialization
        self.structure_solver.Initialize()

        # Ensure that the fluid reaction fluxes are computed if D-N scheme is considered
        if self.fluid_solver.settings["compute_reactions"].GetBool() == False:
            self.fluid_solver.settings["compute_reactions"].SetBool(True)

        # In the D-N scheme the interface correction is done over the velocity
        self.correction_over_velocity = True

        # Python fluid solver initialization
        self.fluid_solver.Initialize()

        # Python mesh solver initialization
        self.mesh_solver.Initialize()

        # Construct the interface mapper
        # Recall, to set the INTERFACE flag in both the fluid and solid interface before the mapper construction
        # Currently this is done with the FSI application Python process set_interface_process.py
        search_radius_factor = 2.0
        mapper_max_iterations = 100
        mapper_tolerance = 0.001*self.nl_tol
        self.interface_mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(self.fluid_solver.main_model_part,
                                                                                  self.structure_solver.main_model_part,
                                                                                  search_radius_factor,
                                                                                  mapper_max_iterations,
                                                                                  mapper_tolerance)

        # Set the Neumann B.C. in the structure interface
        self._SetStructureNeumannCondition() #TODO change when the interface is able to correctly transfer distributed forces

        # Note that the FSI problem is defined in terms of the fluid interface
        # Initialize the iteration value for the residual computation
        fluid_interface_residual_size = self.partitioned_fsi_utilities.GetFluidInterfaceResidualSize()
        self.iteration_value = KratosMultiphysics.Vector(fluid_interface_residual_size)     # Interface solution guess (it might be velocity or fluxes depending on the type of coupling)
        for i in range(0,fluid_interface_residual_size):
            self.iteration_value[i] = 0.0001

        # Compute the fluid domain NODAL_AREA values (required as weight in the residual norm computation)
        KratosMultiphysics.CalculateNodalAreaProcess(self.fluid_solver.GetComputingModelPart(), self.domain_size).Execute()

        # Strategies initialization
        super(PartitionedFSIDirichletNeumannSolver, self).Initialize()


    def Solve(self):

        ## Solvers initialization
        self.InitializeSolutionStep()

        ## Solvers predict
        self.Predict()

        ## Compute mesh prediction ##
        self._ComputeMeshPrediction()

        # Fills self.iteration_value with the structure velocity mapped to the fluid (VECTOR_PROJECTED)
        self._InitializeIterationValueVector()

        # Set the self.iteration value vector as initial fluid interface velocity and mesh velocity
        # self.partitioned_fsi_utilities.SetAndFixFluidInterfaceVectorVariable(KratosMultiphysics.VELOCITY, True, self.iteration_value)
        # self.partitioned_fsi_utilities.SetAndFixFluidInterfaceVectorVariable(KratosMultiphysics.MESH_VELOCITY, False, self.iteration_value)

        ## Non-Linear interface coupling iteration ##
        for nl_it in range(1,self.max_nl_it+1):

            print("     NL-ITERATION ",nl_it,"STARTS.")
            self.fluid_solver.main_model_part.ProcessInfo[KratosFSI.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it
            self.structure_solver.main_model_part.ProcessInfo[KratosFSI.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it

            self.coupling_utility.InitializeNonLinearIteration()

            print("     Residual computation starts...")
            # Sets self.iteration_value as fluid interface VELOCITY
            # Sets self.iteration_value as fluid interface MESH_VELOCITY
            # Solves the mesh problem
            # Solves the fluid problem and computes the nodal fluxes (REACTION)
            self._SolveMeshAndFluid()

            # Sets the structure POINT_LOAD as the fluid interface REACTION
            # Solves the structure problem
            self._SolveStructure()

            # Residual computation
            vel_residual = self._ComputeVelocityResidual()
            nl_res_norm = self.fluid_solver.main_model_part.ProcessInfo[KratosFSI.FSI_INTERFACE_RESIDUAL_NORM]
            FluidInterfaceArea = self.partitioned_fsi_utilities.GetFluidInterfaceArea()

            # Check convergence
            if nl_res_norm/FluidInterfaceArea < self.nl_tol:
                print("     NON-LINEAR ITERATION CONVERGENCE ACHIEVED")
                print("     Total non-linear iterations: ",nl_it," |res|/A = ",nl_res_norm/FluidInterfaceArea)
                break
            else:
                # If convergence is not achieved, perform the correction of the prediction
                print("     Residual computation finished. |res|/A =", nl_res_norm/FluidInterfaceArea)
                print("     Performing non-linear iteration ",nl_it," correction.")
                self.coupling_utility.UpdateSolution(vel_residual, self.iteration_value)
                self.coupling_utility.FinalizeNonLinearIteration()

            # Set corrected iteration value vector as fluid velocity and mesh velocity
            self.partitioned_fsi_utilities.SetAndFixFluidInterfaceVectorVariable(KratosMultiphysics.VELOCITY, True, self.iteration_value)


        ## Compute the mesh residual
        self.partitioned_fsi_utilities.ComputeFluidInterfaceMeshVelocityResidualNorm()
        mesh_res_norm = self.fluid_solver.main_model_part.ProcessInfo.GetValue(KratosFSI.FSI_INTERFACE_MESH_RESIDUAL_NORM)
        print("     NL residual norm: ", nl_res_norm)
        print("     Mesh residual norm: ", mesh_res_norm)

        ## Finalize solution step
        self.fluid_solver.SolverFinalizeSolutionStep()
        self.structure_solver.SolverFinalizeSolutionStep()
        self.coupling_utility.FinalizeSolutionStep()


    def _InitializeIterationValueVector(self):

        # Project the structure velocity onto the fluid interface
        keep_sign = True
        distribute_load = False
        self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.VELOCITY,
                                                         KratosFSI.VECTOR_PROJECTED,
                                                         keep_sign,
                                                         distribute_load)

        i = 0
        if (self.domain_size == 2):
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                    self.iteration_value[i*self.domain_size] = node.GetSolutionStepValue(KratosFSI.VECTOR_PROJECTED_X)
                    self.iteration_value[i*self.domain_size+1] = node.GetSolutionStepValue(KratosFSI.VECTOR_PROJECTED_Y)
                    i+=1
        else:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                self.iteration_value[i*self.domain_size] = node.GetSolutionStepValue(KratosFSI.VECTOR_PROJECTED_X)
                self.iteration_value[i*self.domain_size+1] = node.GetSolutionStepValue(KratosFSI.VECTOR_PROJECTED_Y)
                self.iteration_value[i*self.domain_size+2] = node.GetSolutionStepValue(KratosFSI.VECTOR_PROJECTED_Z)
                i+=1


    def _SolveMeshAndFluid(self):

        # Set the mesh displacement as the structure one
        keep_sign = True
        distribute_load = False
        self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                         KratosALE.MESH_DISPLACEMENT,
                                                         keep_sign,
                                                         distribute_load)

        # Solve the mesh problem (or moves the interface nodes)
        if (self.solve_mesh_at_each_iteration == True):
            self.mesh_solver.Solve()
        else:
            self.mesh_solver.MoveNodes()

        # Fluid domain velocities imposition
        self.partitioned_fsi_utilities.SetAndFixFluidInterfaceVectorVariable(KratosMultiphysics.VELOCITY, True, self.iteration_value)
        self.partitioned_fsi_utilities.SetAndFixFluidInterfaceVectorVariable(KratosMultiphysics.MESH_VELOCITY, False, self.iteration_value)
        # Note that the VELOCITY and the MESH_VELOCITY values only coincide if the same time schemes are used
        # Currently, the mesh solver only includes the BDF2 so the MESH_VELOCITY values are forced to be the fluid ones

        # Solve fluid problem
        self.fluid_solver.SolverSolveSolutionStep()


    def _SolveStructure(self):

        # Transfer fluid reaction to solid interface
        keep_sign = False
        distribute_load = True
        self.interface_mapper.FluidToStructure_VectorMap(KratosMultiphysics.REACTION,
                                                         KratosSolid.POINT_LOAD,
                                                         keep_sign,
                                                         distribute_load)

        # Solve the structure problem
        self.structure_solver.SolverSolveSolutionStep()


    def _ComputeVelocityResidual(self):
        # Project the structure velocity onto the fluid interface
        keep_sign = True
        distribute_load = False
        self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.VELOCITY,
                                                         KratosFSI.VECTOR_PROJECTED,
                                                         keep_sign,
                                                         distribute_load)

        # Compute the fluid interface residual vector by means of the VECTOR_PROJECTED variable
        # Besides, its norm is stored within the ProcessInfo.
        vel_residual = KratosMultiphysics.Vector(self.partitioned_fsi_utilities.GetFluidInterfaceResidualSize())
        self.partitioned_fsi_utilities.ComputeFluidInterfaceVelocityResidual(vel_residual)

        return vel_residual


    ### INTERFACE MOVEMENT UTILITY ###
    # Function to update the position of the interface during iterations
    def _ComputeCorrectedInterfacePosition(self):

        # Bossak parameters
        alpha = -1/3
        gamma = 0.5*(1-2*alpha)
        beta = ((1-alpha)**2)/4

        i = 0
        if self.domain_size == 2:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                u_n = node.GetSolutionStepValue(KratosALE.MESH_DISPLACEMENT,1)
                v_n = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,1)
                a_n = node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION,1)

                v_n1 = KratosMultiphysics.Vector(3)
                v_n1[0] = self.iteration_value[i]
                v_n1[1] = self.iteration_value[i+1]
                v_n1[2] = 0.0
                i+=2

                # Compute the current acceleration associated to the corrected interface velocity with the Bossak formulaes
                a_n1 = KratosMultiphysics.Vector(3)
                a_n1[0] = (v_n1[0] - v_n[0] - self.time_step*(1-gamma*(alpha-1))*a_n[0]) / ((1-alpha)*self.time_step*gamma)
                a_n1[1] = (v_n1[1] - v_n[1] - self.time_step*(1-gamma*(alpha-1))*a_n[1]) / ((1-alpha)*self.time_step*gamma)
                a_n1[2] = 0.0

                # Compute the current displacement associated to the corrected interface velocity with the Bossak formulaes
                u_n1 = KratosMultiphysics.Vector(3)
                u_n1[0] = u_n[0] + self.time_step*v_n[0] + (self.time_step**2)*(0.5-beta)*a_n[0] + (self.time_step**2)*beta*((1-alpha)*a_n1[0] + alpha*a_n[0])
                u_n1[1] = u_n[1] + self.time_step*v_n[1] + (self.time_step**2)*(0.5-beta)*a_n[1] + (self.time_step**2)*beta*((1-alpha)*a_n1[1] + alpha*a_n[1])
                u_n1[2] = 0.0

                # Set the obtained corrected interface displacement
                node.SetSolutionStepValue(KratosALE.MESH_DISPLACEMENT,0,u_n1)

        elif self.domain_size == 3:
            for node in self._GetFluidInterfaceSubmodelPart().Nodes:
                u_n = node.GetSolutionStepValue(KratosALE.MESH_DISPLACEMENT,1)
                v_n = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,1)
                a_n = node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION,1)

                v_n1 = KratosMultiphysics.Vector(3)
                v_n1[0] = self.iteration_value[i]
                v_n1[1] = self.iteration_value[i+1]
                v_n1[2] = self.iteration_value[i+2]
                i+=3

                # Compute the current acceleration associated to the corrected interface velocity with the Bossak formulaes
                a_n1 = KratosMultiphysics.Vector(3)
                a_n1[0] = (v_n1[0] - v_n[0] - self.time_step*(1-gamma*(alpha-1))*a_n[0]) / ((1-alpha)*self.time_step*gamma)
                a_n1[1] = (v_n1[1] - v_n[1] - self.time_step*(1-gamma*(alpha-1))*a_n[1]) / ((1-alpha)*self.time_step*gamma)
                a_n1[2] = (v_n1[2] - v_n[2] - self.time_step*(1-gamma*(alpha-1))*a_n[2]) / ((1-alpha)*self.time_step*gamma)

                # Compute the current displacement associated to the corrected interface velocity with the Bossak formulaes
                u_n1 = KratosMultiphysics.Vector(3)
                u_n1[0] = u_n[0] + self.time_step*v_n[0] + (self.time_step**2)*(0.5-beta)*a_n[0] + (self.time_step**2)*beta*((1-alpha)*a_n1[0] + alpha*a_n[0])
                u_n1[1] = u_n[1] + self.time_step*v_n[1] + (self.time_step**2)*(0.5-beta)*a_n[1] + (self.time_step**2)*beta*((1-alpha)*a_n1[1] + alpha*a_n[1])
                u_n1[2] = u_n[2] + self.time_step*v_n[2] + (self.time_step**2)*(0.5-beta)*a_n[2] + (self.time_step**2)*beta*((1-alpha)*a_n1[2] + alpha*a_n[2])

                # Set the obtained corrected interface displacement
                node.SetSolutionStepValue(KratosALE.MESH_DISPLACEMENT,0,u_n1)

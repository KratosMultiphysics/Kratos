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

        # Python fluid solver initialization
        self.fluid_solver.Initialize()

        # Python mesh solver initialization
        self.mesh_solver.Initialize()

        # Construct the interface mapper
        # Recall, to set the INTERFACE flag in both the fluid and solid interface before the mapper construction
        # Currently this is done with the FSI application Python process set_interface_process.py
        search_radius_factor = 2.0
        mapper_max_iterations = 200
        mapper_tolerance = 1e-12
        self.interface_mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(self.fluid_solver.main_model_part,
                                                                                  self.structure_solver.main_model_part,
                                                                                  search_radius_factor,
                                                                                  mapper_max_iterations,
                                                                                  mapper_tolerance)

        # Set the Neumann B.C. in the structure interface
        self._SetStructureNeumannCondition() #TODO change when the interface is able to correctly transfer distributed forces

        # Initialize the iteration value vector
        self._InitializeIterationValueVector()

        # Initialize the Dirichlet-Neumann interface
        self._InitializeDirichletNeumannInterface()

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
            self._SolveStructure()

            # Residual computation
            dis_residual = self._ComputeDisplacementResidual()
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
                # self.coupling_utility.UpdateSolution(vel_residual, self.iteration_value)
                self.coupling_utility.UpdateSolution(dis_residual, self.iteration_value)
                self.coupling_utility.FinalizeNonLinearIteration()

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
        # Note that the FSI problem is defined in terms of the fluid interface
        # Initialize the iteration value for the residual computation
        fluid_interface_residual_size = self.partitioned_fsi_utilities.GetFluidInterfaceResidualSize()
        self.iteration_value = KratosMultiphysics.Vector(fluid_interface_residual_size)
        for i in range(0,fluid_interface_residual_size):
            self.iteration_value[i] = 0.0


    def _InitializeDirichletNeumannInterface(self):
        # Fix the VELOCITY and MESH_VELOCITY variables
        # MESH_DISPLACEMENT variable is supposed to be already fixed by the set_interface_process
        for node in self._GetFluidInterfaceSubmodelPart().Nodes:
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.Fix(KratosMultiphysics.VELOCITY_Y)
            node.Fix(KratosMultiphysics.VELOCITY_Z)
            node.Fix(KratosMultiphysics.MESH_VELOCITY_X)
            node.Fix(KratosMultiphysics.MESH_VELOCITY_Y)
            node.Fix(KratosMultiphysics.MESH_VELOCITY_Z)


    def _SolveMeshAndFluid(self):

        # Set the mesh displacement as the iteration_value displacement
        self.partitioned_fsi_utilities.SetFluidInterfaceVectorVariable(KratosALE.MESH_DISPLACEMENT, self.iteration_value)

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


    def _ComputeDisplacementResidual(self):
        # Project the structure velocity onto the fluid interface
        keep_sign = True
        distribute_load = False
        self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                         KratosFSI.VECTOR_PROJECTED,
                                                         keep_sign,
                                                         distribute_load)

        # Compute the fluid interface residual vector by means of the VECTOR_PROJECTED variable
        # Besides, its norm is stored within the ProcessInfo.
        disp_residual = KratosMultiphysics.Vector(self.partitioned_fsi_utilities.GetFluidInterfaceResidualSize())
        self.partitioned_fsi_utilities.ComputeFluidInterfaceVectorResidual(KratosALE.MESH_DISPLACEMENT, KratosFSI.VECTOR_PROJECTED, disp_residual)

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
                u_n = node.GetSolutionStepValue(KratosALE.MESH_DISPLACEMENT,1)
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
                u_n = node.GetSolutionStepValue(KratosALE.MESH_DISPLACEMENT,1)
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

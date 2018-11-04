from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from math import sqrt

# Import utilities
import python_solvers_wrapper_fluid            # Import the fluid Python solvers wrapper
import python_solvers_wrapper_structural       # Import the structure Python solvers wrapper
import convergence_accelerator_factory         # Import the FSI convergence accelerator factory

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications(
    "MetisApplication",
    "TrilinosApplication",
    "MappingApplication",
    "FSIApplication",
    "MeshMovingApplication",
    "FluidDynamicsApplication",
    "StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.MetisApplication as KratosMetis
import KratosMultiphysics.TrilinosApplication as KratosTrilinos
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Import base class file
import trilinos_partitioned_fsi_base_solver     # Base class file

def CreateSolver(model, project_parameters):
    return TrilinosPartitionedFSIDirichletNeumannSolver(model, project_parameters)

class TrilinosPartitionedFSIDirichletNeumannSolver(trilinos_partitioned_fsi_base_solver.TrilinosPartitionedFSIBaseSolver):
    def __init__(self, model, project_parameters):
        super(TrilinosPartitionedFSIDirichletNeumannSolver, self).__init__(model, project_parameters)
        self._PrintInfoOnRankZero("::[TrilinosPartitionedFSIDirichletNeumannSolver]::", "Solver construction finished.")

    def Initialize(self):
        # Set the Trilinos space
        self.trilinos_space = KratosTrilinos.TrilinosSparseSpace()

        # Set the Epetra communicator
        self.epetra_communicator = KratosTrilinos.CreateCommunicator()

        # Construct the coupling partitioned strategy
        coupling_utility_parameters = self.settings["coupling_solver_settings"]["solver_settings"]["coupling_strategy"]
        self.coupling_utility = convergence_accelerator_factory.CreateTrilinosConvergenceAccelerator(self._GetFluidInterfaceSubmodelPart(),
                                                                                                     self.epetra_communicator,
                                                                                                     coupling_utility_parameters)

        # Get the domain size
        self.domain_size = self._GetDomainSize()

        # Get the nodal update FSI utilities
        self.nodal_update_utilities = self._GetNodalUpdateUtilities()

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

        # Initialize the Dirichlet-Neumann interface
        self._InitializeDirichletNeumannInterface()

        # Construct the interface mapper
        self._SetUpMapper()

        # Set the Neumann B.C. in the structure interface
        self._SetStructureNeumannCondition() #TODO change when the interface is able to correctly transfer distributed forces

        # Initialize the iteration value for the residual computation, which is defined in terms of the fluid interface
        # In case a shell structure is analised, only one of the two interfaces is considered for the residual vector computation
        # This is due to the fact that the same mesh_solver displacement is map to both fluid sides.
        self.iteration_value = self.partitioned_fsi_utilities.SetUpInterfaceVector(self._GetFluidInterfaceSubmodelPart())

        # Compute the fluid domain NODAL_AREA values (required as weight in the residual norm computation)
        KratosMultiphysics.CalculateNodalAreaProcess(self.fluid_solver.GetComputingModelPart(), self.domain_size).Execute()

    def _InitializeDirichletNeumannInterface(self):
        # Initialize Dirichlet fluid interface
        fluid_interfaces_list = self.settings["coupling_solver_settings"]["solver_settings"]["fluid_interfaces_list"]
        for fl_interface_id in range(fluid_interfaces_list.size()):
            fl_interface_name = fluid_interfaces_list[fl_interface_id].GetString()
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
        structure_interfaces_list = self.settings["coupling_solver_settings"]["solver_settings"]["structure_interfaces_list"]
        for str_interface_id in range(structure_interfaces_list.size()):
            str_interface_name = structure_interfaces_list[str_interface_id].GetString()
            str_interface_submodelpart = self.structure_solver.main_model_part.GetSubModelPart(str_interface_name)

            # Set the interface flag
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, str_interface_submodelpart.Nodes)

    def _SolveMeshAndFluid(self):
        # Set the iteration_value displacement as MESH_DISPLACEMENT
        coupling_solver_settings = self.settings["coupling_solver_settings"]["solver_settings"]
        num_fl_interfaces = coupling_solver_settings["fluid_interfaces_list"].size()
        for fl_interface_id in range(num_fl_interfaces):
            fl_interface_name = coupling_solver_settings["fluid_interfaces_list"][fl_interface_id].GetString()
            fl_interface_submodelpart = self.fluid_solver.main_model_part.GetSubModelPart(fl_interface_name)
            self.partitioned_fsi_utilities.UpdateInterfaceValues(fl_interface_submodelpart,
                                                                 KratosMultiphysics.MESH_DISPLACEMENT,
                                                                 self.iteration_value)

        # Solve the mesh problem (or moves the interface nodes)
        if self.solve_mesh_at_each_iteration:
            self.mesh_solver.Solve()
        else:
            self.mesh_solver.MoveMesh()

        # Update MESH_VELOCITY and MESH_ACCELERATION with Newmark formulas
        self.nodal_update_utilities.UpdateMeshTimeDerivatives(self.fluid_solver.GetComputingModelPart(),
                                                              self.time_step)

        # Impose the structure MESH_VELOCITY and MESH_ACCELERATION in the fluid interface VELOCITY and ACCELERATION
        self.nodal_update_utilities.SetMeshTimeDerivativesOnInterface(self._GetFluidInterfaceSubmodelPart())

        # Solve fluid problem
        self.fluid_solver.SolveSolutionStep()

    def _SolveStructureSingleFaced(self):
        # Set the redistribution settings
        redistribution_tolerance = 1e-8
        redistribution_max_iters = 50

        # Convert the nodal reaction to traction loads before transfering
        KratosMultiphysics.VariableRedistributionUtility.DistributePointValues(
            self._GetFluidInterfaceSubmodelPart(),
            KratosMultiphysics.REACTION,
            KratosMultiphysics.VAUX_EQ_TRACTION,
            redistribution_tolerance,
            redistribution_max_iters)

        # Transfer fluid tractions to the structure interface
        self.interface_mapper.Map(KratosMultiphysics.VAUX_EQ_TRACTION,
                                  KratosMultiphysics.VAUX_EQ_TRACTION,
                                  KratosMapping.Mapper.SWAP_SIGN)

        # Convert the transferred traction loads to point loads
        KratosMultiphysics.VariableRedistributionUtility.ConvertDistributedValuesToPoint(
            self._GetStructureInterfaceSubmodelPart(),
            KratosMultiphysics.VAUX_EQ_TRACTION,
            KratosStructural.POINT_LOAD)

        # Solve the structure problem
        is_converged = self.structure_solver.SolveSolutionStep()
        if not is_converged:
            self._PrintWarningOnRankZero("Structure solver did not converge.")

    def _SolveStructureDoubleFaced(self):
        # Set the redistribution settings
        redistribution_tolerance = 1e-8
        redistribution_max_iters = 50

        # Convert the nodal reaction to traction loads before transfering
        KratosMultiphysics.VariableRedistributionUtility.DistributePointValues(
            self._GetFluidPositiveInterfaceSubmodelPart(),
            KratosMultiphysics.REACTION,
            KratosMultiphysics.VAUX_EQ_TRACTION,
            redistribution_tolerance,
            redistribution_max_iters)

        KratosMultiphysics.VariableRedistributionUtility.DistributePointValues(
            self._GetFluidNegativeInterfaceSubmodelPart(),
            KratosMultiphysics.REACTION,
            KratosMultiphysics.VAUX_EQ_TRACTION,
            redistribution_tolerance,
            redistribution_max_iters)

        # Transfer fluid tractions to the structure interface
        # Note that the ADD_VALUES flag is only specified for the second mapper
        # since we want the first mapper to overwrite the existent values
        self.pos_interface_mapper.Map(KratosMultiphysics.VAUX_EQ_TRACTION,
                                      KratosMultiphysics.VAUX_EQ_TRACTION,
                                      KratosMapping.Mapper.SWAP_SIGN)

        # Transfer fluid tractions to the structure interface
        self.neg_interface_mapper.Map(KratosMultiphysics.VAUX_EQ_TRACTION,
                                      KratosMultiphysics.VAUX_EQ_TRACTION,
                                      KratosMapping.Mapper.SWAP_SIGN | KratosMapping.Mapper.ADD_VALUES)

        # Convert the transferred traction loads to point loads
        KratosMultiphysics.VariableRedistributionUtility.ConvertDistributedValuesToPoint(
            self._GetStructureInterfaceSubmodelPart(),
            KratosMultiphysics.VAUX_EQ_TRACTION,
            KratosStructural.POINT_LOAD)

        # Solve the structure problem
        is_converged = self.structure_solver.SolveSolutionStep()
        if not is_converged:
            self._PrintWarningOnRankZero("Structure solver did not converge.")

    def _ComputeDisplacementResidualSingleFaced(self):
        # Project the structure displacement onto the fluid interface
        self.interface_mapper.InverseMap(KratosMultiphysics.VECTOR_PROJECTED,
                                         KratosMultiphysics.DISPLACEMENT)

        # Compute the fluid interface residual vector by means of the VECTOR_PROJECTED variable
        # Besides, its norm is stored within the ProcessInfo.
        disp_residual = self.partitioned_fsi_utilities.SetUpInterfaceVector(self._GetFluidInterfaceSubmodelPart())
        self.partitioned_fsi_utilities.ComputeInterfaceResidualVector(self._GetFluidInterfaceSubmodelPart(),
                                                                      KratosMultiphysics.MESH_DISPLACEMENT,
                                                                      KratosMultiphysics.VECTOR_PROJECTED,
                                                                      disp_residual)

        return disp_residual

    def _ComputeDisplacementResidualDoubleFaced(self):
        # Project the structure displacement onto the fluid positive interface
        self.pos_interface_mapper.InverseMap(KratosMultiphysics.VECTOR_PROJECTED,
                                             KratosMultiphysics.DISPLACEMENT)

        self.neg_interface_mapper.InverseMap(KratosMultiphysics.VECTOR_PROJECTED,
                                             KratosMultiphysics.DISPLACEMENT)

        #TODO: One of these mappings is not needed. At the moment both are performed to ensure that
        # the _GetFluidInterfaceSubmodelPart(), which is the one used to compute the interface residual,
        # gets the structure DISPLACEMENT. Think a way to properly identify the reference fluid interface.
        disp_residual = self.partitioned_fsi_utilities.SetUpInterfaceVector(self._GetFluidInterfaceSubmodelPart())
        self.partitioned_fsi_utilities.ComputeInterfaceResidualVector(self._GetFluidInterfaceSubmodelPart(),
                                                                      KratosMultiphysics.MESH_DISPLACEMENT,
                                                                      KratosMultiphysics.VECTOR_PROJECTED,
                                                                      disp_residual)

        return disp_residual
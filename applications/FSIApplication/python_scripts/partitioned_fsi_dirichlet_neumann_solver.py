from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Import base class file
from KratosMultiphysics.FSIApplication import partitioned_fsi_base_solver

def CreateSolver(model, project_parameters):
    return PartitionedFSIDirichletNeumannSolver(model, project_parameters)

class PartitionedFSIDirichletNeumannSolver(partitioned_fsi_base_solver.PartitionedFSIBaseSolver):
    def __init__(self, model, project_parameters):
        super(PartitionedFSIDirichletNeumannSolver, self).__init__(model, project_parameters)
        KratosMultiphysics.Logger.PrintInfo("::[PartitionedFSIDirichletNeumannSolver]::", "Solver construction finished.")

    @classmethod
    def GetDefaultSettings(cls):
        """This function returns the default-settings used by this class
        """
        this_defaults = KratosMultiphysics.Parameters("""{
            "coupling_scheme": "dirichlet_neumann"
        }""")
        this_defaults.AddMissingParameters(super(PartitionedFSIDirichletNeumannSolver, cls).GetDefaultSettings())
        return this_defaults

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

        # Coupling utility initialization
        # The _GetConvergenceAccelerator is supposed to construct the convergence accelerator in here
        self._GetConvergenceAccelerator().Initialize()

    def _InitializeIterationValueVector(self):
        # Note that the FSI problem is defined in terms of the fluid interface
        # Besides, one of the two interfaces is considered for the residual vector computation
        # in case a shell structure is analised. This is due to the fact that the same mesh_solver
        # displacement is map to both fluid sides.

        # Initialize the iteration value for the residual computation
        fluid_interface_residual_size = self.partitioned_fsi_utilities.GetInterfaceResidualSize(self._GetFluidInterfaceSubmodelPart())
        self.iteration_value = KratosMultiphysics.Vector(fluid_interface_residual_size)
        for i in range(0,fluid_interface_residual_size):
            self.iteration_value[i] = 0.0

    def _InitializeDirichletNeumannInterface(self):
        # Initialize Dirichlet fluid interface
        coupling_settings = self.settings["coupling_settings"]
        num_fl_interfaces = coupling_settings["fluid_interfaces_list"].size()
        for fl_interface_id in range(num_fl_interfaces):
            fl_interface_name = coupling_settings["fluid_interfaces_list"][fl_interface_id].GetString()
            fl_interface_submodelpart = self.fluid_solver.main_model_part.GetSubModelPart(fl_interface_name)

            # Fix the VELOCITY, MESH_DISPLACEMENT and MESH_VELOCITY variables in all the fluid interface submodelparts
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VELOCITY_X, True, fl_interface_submodelpart.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VELOCITY_Y, True, fl_interface_submodelpart.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VELOCITY_Z, True, fl_interface_submodelpart.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.MESH_DISPLACEMENT_X, True, fl_interface_submodelpart.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.MESH_DISPLACEMENT_Y, True, fl_interface_submodelpart.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.MESH_DISPLACEMENT_Z, True, fl_interface_submodelpart.Nodes)

            # Set the interface flag
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, fl_interface_submodelpart.Nodes)

        # Initialize Neumann structure interface
        num_str_interfaces = coupling_settings["structure_interfaces_list"].size()
        for str_interface_id in range(num_str_interfaces):
            str_interface_name = coupling_settings["structure_interfaces_list"][str_interface_id].GetString()
            str_interface_submodelpart = self.structure_solver.main_model_part.GetSubModelPart(str_interface_name)

            # Set the interface flag
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, str_interface_submodelpart.Nodes)

    def _SolveMeshAndFluid(self):
        # Set the iteration_value displacement as MESH_DISPLACEMENT
        coupling_settings = self.settings["coupling_settings"]
        num_fl_interfaces = coupling_settings["fluid_interfaces_list"].size()
        for fl_interface_id in range(num_fl_interfaces):
            fl_interface_name = coupling_settings["fluid_interfaces_list"][fl_interface_id].GetString()
            fl_interface_submodelpart = self.fluid_solver.main_model_part.GetSubModelPart(fl_interface_name)
            self.partitioned_fsi_utilities.UpdateInterfaceValues(
                fl_interface_submodelpart,
                KratosMultiphysics.MESH_DISPLACEMENT,
                self.iteration_value)

        # Solve the mesh problem (or moves the interface nodes)
        if self.solve_mesh_at_each_iteration:
            self.mesh_solver.SolveSolutionStep()
        else:
            self.mesh_solver.MoveMesh()

        # Impose the structure MESH_VELOCITY and MESH_ACCELERATION in the fluid interface VELOCITY and ACCELERATION
        KratosMultiphysics.VariableUtils().CopyVectorVar(
                KratosMultiphysics.MESH_VELOCITY,
                KratosMultiphysics.VELOCITY,
                self._GetFluidInterfaceSubmodelPart().GetCommunicator().LocalMesh().Nodes)
        self._GetFluidInterfaceSubmodelPart().GetCommunicator().SynchronizeVariable(KratosMultiphysics.VELOCITY)

        # Solve fluid problem
        self.fluid_solver.SolveSolutionStep()

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
                                                                 KratosMultiphysics.POSITIVE_MAPPED_VECTOR_VARIABLE,
                                                                 keep_sign,
                                                                 distribute_load)
        self.interface_mapper.NegativeFluidToStructure_VectorMap(KratosMultiphysics.REACTION,
                                                                 KratosMultiphysics.NEGATIVE_MAPPED_VECTOR_VARIABLE,
                                                                 keep_sign,
                                                                 distribute_load)

        # Add the two faces contributions to the POINT_LOAD variable
        # TODO: Add this to the variables utils
        for node in self._GetStructureInterfaceSubmodelPart().Nodes:
            pos_face_force = node.GetSolutionStepValue(KratosMultiphysics.POSITIVE_MAPPED_VECTOR_VARIABLE)
            neg_face_force = node.GetSolutionStepValue(KratosMultiphysics.NEGATIVE_MAPPED_VECTOR_VARIABLE)
            node.SetSolutionStepValue(KratosStructural.POINT_LOAD, 0, pos_face_force+neg_face_force)

        # Solve the structure problem
        self.structure_solver.SolveSolutionStep()

    def _ComputeDisplacementResidualSingleFaced(self):
        # Project the structure velocity onto the fluid interface
        keep_sign = True
        distribute_load = False
        self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                         KratosMultiphysics.VECTOR_PROJECTED,
                                                         keep_sign,
                                                         distribute_load)

        # Compute the fluid interface residual vector by means of the VECTOR_PROJECTED variable
        # Besides, its norm is stored within the ProcessInfo.
        disp_residual = KratosMultiphysics.Vector(self.partitioned_fsi_utilities.GetInterfaceResidualSize(self._GetFluidInterfaceSubmodelPart()))
        self.partitioned_fsi_utilities.ComputeInterfaceResidualVector(
            self._GetFluidInterfaceSubmodelPart(),
            KratosMultiphysics.MESH_DISPLACEMENT,
            KratosMultiphysics.VECTOR_PROJECTED,
            KratosMultiphysics.FSI_INTERFACE_RESIDUAL,
            disp_residual,
            "nodal",
            KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM)

        return disp_residual

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
        disp_residual = KratosMultiphysics.Vector(self.partitioned_fsi_utilities.GetInterfaceResidualSize(self._GetFluidInterfaceSubmodelPart()))
        self.partitioned_fsi_utilities.ComputeInterfaceResidualVector(
            self._GetFluidInterfaceSubmodelPart(),
            KratosMultiphysics.MESH_DISPLACEMENT,
            KratosMultiphysics.VECTOR_PROJECTED,
            KratosMultiphysics.FSI_INTERFACE_MESH_RESIDUAL,
            disp_residual,
            "nodal",
            KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM)

        return disp_residual

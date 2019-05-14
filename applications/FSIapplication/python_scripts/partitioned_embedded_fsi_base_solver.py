from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from math import sqrt   # Import the square root from python library

# Import utilities
from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid            # Import the fluid Python solvers wrapper
from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural # Import the structure Python solvers wrapper
from KratosMultiphysics.MeshMovingApplication import python_solvers_wrapper_mesh_motion         # Import the mesh motion Python solvers wrapper
from KratosMultiphysics.FSIApplication import fsi_coupling_interface # Import the FSI coupling interface utility
from KratosMultiphysics.FSIApplication import convergence_accelerator_factory                   # Import the FSI convergence accelerator factory

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def CreateSolver(model, project_parameters):
    return PartitionedEmbeddedFSIBaseSolver(model, project_parameters)

class PartitionedEmbeddedFSIBaseSolver(PythonSolver):

    def _ValidateSettings(self, project_parameters):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "echo_level": 0,
            "parallel_type": "OpenMP",
            "solver_type": "partitioned_embedded",
            "coupling_scheme": "dirichlet_neumann",
            "structure_solver_settings": {
            },
            "fluid_solver_settings":{
            },
            "coupling_settings":{
            }
        }""")

        project_parameters.ValidateAndAssignDefaults(default_settings)

        if not project_parameters["structure_solver_settings"].Has("multi_point_constraints_used"):
            project_parameters["structure_solver_settings"].AddEmptyValue("multi_point_constraints_used")
            project_parameters["structure_solver_settings"]["multi_point_constraints_used"].SetBool(False)

        return project_parameters

    def __init__(self, model, project_parameters):
        # Validate settings
        project_parameters = self._ValidateSettings(project_parameters)

        # Call the base Python solver constructor
        super(PartitionedEmbeddedFSIBaseSolver,self).__init__(model, project_parameters)

        # Auxiliar variables
        self.parallel_type = self.settings["parallel_type"].GetString()
        coupling_settings = self.settings["coupling_settings"]
        self.max_nl_it = coupling_settings["nl_max_it"].GetInt()
        self.nl_tol = coupling_settings["nl_tol"].GetDouble()
        self.structure_interface_submodelpart_name = coupling_settings["structure_interfaces_list"][0].GetString()

        # Construct the structure solver
        self.structure_solver = python_solvers_wrapper_structural.CreateSolverByParameters(self.model, self.settings["structure_solver_settings"], self.parallel_type)
        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Structure solver construction finished')

        # Construct the fluid solver
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_solver_settings"], self.parallel_type)
        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Fluid solver construction finished')
        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Partitioned embedded FSI base solver construction finished')

    def GetMinimumBufferSize(self):
        buffer_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_structure = self.structure_solver.GetMinimumBufferSize()
        return max(buffer_structure,buffer_fluid)

    #TODO: CHECK WHICH VARIABLES ARE MISSING OR NOT NEEDED
    def AddVariables(self):
        # Structure solver variables addition
        self.structure_solver.AddVariables()
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)

        # Fluid solver variables addition
        self.fluid_solver.AddVariables()
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
        #TODO: REMOVE THIS! ONLY ADDED FOR DEBUGGING WITHOUT MODIFYING THE JSON WHEN THE FM-ALE ALGORITHM IS SWITCHED OFF.
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)
        #TODO: ADD FOR DEBUGGING
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LINEAR_MOMENTUM)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_MOMENTUM)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE_OLD_IT)

        # FSI coupling required variables addition
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.STRUCTURE_VELOCITY)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SCALAR_PROJECTED)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_PROJECTED)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SCALAR_INTERFACE_RESIDUAL)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FSI_INTERFACE_RESIDUAL)

    def ImportModelPart(self):
        # Fluid and structure solvers ImportModelPart() call
        self.fluid_solver.ImportModelPart()
        self.structure_solver.ImportModelPart()

    def PrepareModelPart(self):
        # Get the minimum buffer size between the mesh, fluid and structure solvers
        self._GetAndSetMinimumBufferSize()
        # Fluid and structure solvers PrepareModelPart() call
        self.fluid_solver.PrepareModelPart()
        self.structure_solver.PrepareModelPart()

        #TODO: REMOVE THIS CALL AFTER DEBUGGING
        self.domain_size = self._GetDomainSize()
        # self._get_embedded_skin_utility()

    def AddDofs(self):
        # Add DOFs structure
        self.structure_solver.AddDofs()
        # Add DOFs fluid
        self.fluid_solver.AddDofs()

    def Initialize(self):
        # Get the domain size
        self.domain_size = self._GetDomainSize()

        # Coupling utility initialization
        # The _get_convergence_accelerator is supposed to construct the convergence accelerator in here
        self._get_convergence_accelerator().Initialize()

        # FSI interface coupling interfaces initialization
        # The _get_fsi_coupling_interface_structure is supposed to construct the FSI coupling structure interface in here
        self._get_fsi_coupling_interface_fluid().GetInterfaceModelPart()
        self._get_fsi_coupling_interface_structure().GetInterfaceModelPart()

        # Python structure solver initialization
        self.structure_solver.Initialize()

        # Python fluid solver initialization
        self.fluid_solver.Initialize()

        # Initialize the Dirichlet-Neumann interface
        self._initialize_fsi_interfaces()

        # Initialize the iteration value vector
        self._initialize_iteration_value_vector()

        # Initialize the distance field
        update_distance_process = True
        self._get_distance_to_skin_process(update_distance_process).Execute()

        # Initialize the embedded intersections model part
        # self._get_embedded_skin_utility().GenerateSkin()

    def AdvanceInTime(self, current_time):
        fluid_new_time = self.fluid_solver.AdvanceInTime(current_time)
        structure_new_time = self.structure_solver.AdvanceInTime(current_time)

        self._get_fsi_coupling_interface_fluid().GetInterfaceModelPart().GetRootModelPart().CloneTimeStep(fluid_new_time)
        self._get_fsi_coupling_interface_structure().GetInterfaceModelPart().GetRootModelPart().CloneTimeStep(structure_new_time)

        # TODO: REMOVE AFTER DEBUGGING
        self.GetStructureIntersectionsModelPart().GetRootModelPart().CloneTimeStep(fluid_new_time)
        self.GetStructureIntersectionsModelPart().GetRootModelPart().ProcessInfo[KratosMultiphysics.STEP] += 1

        if abs(fluid_new_time - structure_new_time) > 1e-12:
            err_msg =  'Fluid new time is: ' + str(fluid_new_time) + '\n'
            err_msg += 'Structure new time is: ' + str(structure_new_time) + '\n'
            err_msg += 'No substepping has been implemented yet. Fluid and structure time must coincide.'
            raise Exception(err_msg)

        return fluid_new_time

    def InitializeSolutionStep(self):
        # Initialize solution step of fluid, structure and coupling solvers
        self.fluid_solver.InitializeSolutionStep()
        self.structure_solver.InitializeSolutionStep()
        self._get_convergence_accelerator().InitializeSolutionStep()

    def Predict(self):
        # Structure solver prediction. It is important to firstly perform the structure
        # prediction to update the current buffer position before the FM-ALE operations.
        # Otherwise position 0 and 1 of the buffer coincide since the advance in time
        # has been already performed but no update has been done yet. Besides, this will
        # give a better approximation of the level-set position at the end of step.
        self.structure_solver.Predict()

        # # Update the level set position
        self._update_level_set()

        # Fluid solver prediction
        self.fluid_solver.Predict()

    def GetComputingModelPart(self):
        err_msg =  'Calling GetComputingModelPart() method in a partitioned solver.\n'
        err_msg += 'Specify the domain of interest by calling:\n'
        err_msg += '\t- GetFluidComputingModelPart()\n'
        err_msg += '\t- GetStructureComputingModelPart()\n'
        raise Exception(err_msg)

    def GetFluidComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart()

    def GetStructureComputingModelPart(self):
        return self.structure_solver.GetComputingModelPart()

    def GetStructureSkinModelPart(self):
        return self.model.GetModelPart(self._get_structure_interface_model_part_name())

    def GetStructureSkinElementBasedModelPart(self):
        # Create an auxiliar model part to save the element based skin
        element_based_skin_model_part_name = self._get_structure_interface_model_part_name() + "ElementBased"
        if self.model.HasModelPart(element_based_skin_model_part_name):
            self.model.DeleteModelPart(element_based_skin_model_part_name)
        self.element_based_skin_model_part = self.model.CreateModelPart(element_based_skin_model_part_name)

        # Copy the skin model part conditions to an auxiliar model part elements.
        # This is required for the computation of the distance function, which
        # takes the elements of the second modelpart as skin. If this operation
        # is not performed, no elements are found, yielding a wrong level set.
        self._get_partitioned_fsi_utilities().CopySkinToElements(
            self.GetStructureSkinModelPart(),
            self.element_based_skin_model_part)

        return self.element_based_skin_model_part

    def GetStructureIntersectionsModelPart(self):
        if not hasattr(self, '_embedded_intersections_model_part'):
            embedded_intersections_root_part = self.model.CreateModelPart("EmbeddedIntersectionsModelPart")
            self._embedded_intersections_model_part = embedded_intersections_root_part.CreateSubModelPart("SkinEmbeddedIntersectionsModelPart")
            self._embedded_intersections_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
            self._embedded_intersections_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
            self._embedded_intersections_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
            self._embedded_intersections_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
            self._embedded_intersections_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_PROJECTED)
        return self._embedded_intersections_model_part

    def GetOutputVariables(self):
        pass

    def SaveRestart(self):
        pass

    def SolveSolutionStep(self):
        ## TODO: REMOVE AFTER DEBUGGING
        import gid_output_process

        ## Wall non-linear coupling iteration ##
        nl_it = 0
        while (nl_it < self.max_nl_it and self.fluid_solver._TimeBufferIsInitialized()):
            nl_it += 1
            KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', '\tFSI non-linear iteration = ' + str(nl_it))

            self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it
            self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it

            self._get_convergence_accelerator().InitializeNonLinearIteration()

            # Map the DISPLACEMENT from the structure FSI coupling interface to fluid FSI coupling interface
            # TODO: THIS SHOULD BE A MAPPER. NOW WE ASSUME THAT THE EMBEDDED FLUID INTERFACE COINCIDES WITH THE STRUCTURE ONE
            buffer_step = 0
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.DISPLACEMENT,
                self._get_fsi_coupling_interface_structure().GetInterfaceModelPart(),
                self._get_fsi_coupling_interface_fluid().GetInterfaceModelPart(),
                buffer_step)

            # TODO: REMOVE AFTER DEBUGGING (FAKE DISPLACEMENT FIELD)
            # t = self.fluid_solver.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
            # for node in self._get_fsi_coupling_interface_fluid().GetInterfaceModelPart().Nodes:
            #     aux_disp = 1e-3 * t * node.Y ** 2
            #     node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [aux_disp, 0.0, 0.0])

            # Update the fluid FSI coupling interface position
            self._get_fsi_coupling_interface_fluid().UpdatePosition()

            # Update the EMBEDDED_VELOCITY and solve the fluid problem
            self._solve_fluid()

            # Interpolate the pressure to the fluid FSI coupling interface
            self._get_partitioned_fsi_utilities().EmbeddedPressureToPositiveFacePressureInterpolator(
                self.GetFluidComputingModelPart(),
                self._get_fsi_coupling_interface_fluid().GetInterfaceModelPart())

            # Map PRESSURE from fluid FSI coupling interface to structure FSI coupling interface
            # TODO: THIS SHOULD BE A MAPPER. NOW WE ASSUME THAT THE EMBEDDED FLUID INTERFACE COINCIDES WITH THE STRUCTURE ONE
            for node_orig, node_dest in zip(
                self._get_fsi_coupling_interface_fluid().GetInterfaceModelPart().Nodes,
                self._get_fsi_coupling_interface_structure().GetInterfaceModelPart().Nodes):
                aux_pres = node_orig.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
                # t = self.fluid_solver.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
                # aux_pres = 0.5 * t * node_dest.Y
                node_dest.SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE, aux_pres)

            # Transfer PRESSURE to structure FSI coupling interface father model part
            self._get_fsi_coupling_interface_structure().TransferValuesToFatherModelPart(KratosMultiphysics.POSITIVE_FACE_PRESSURE)

            # Update the PRESSURE load and solve the structure problem
            self._solve_structure()

            # fluid_gid_output = gid_output_process.GiDOutputProcess(
            #     self.fluid_solver.GetComputingModelPart(),
            #     "fluid_step_" + str(self.fluid_solver.GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP]) + "_it_" + str(nl_it),
            #     KratosMultiphysics.Parameters("""{
            #         "result_file_configuration": {
            #             "gidpost_flags": {
            #                 "GiDPostMode": "GiD_PostBinary",
            #                 "WriteDeformedMeshFlag": "WriteDeformed",
            #                 "WriteConditionsFlag": "WriteConditions",
            #                 "MultiFileFlag": "SingleFile"
            #             },
            #             "file_label": "time",
            #             "output_control_type": "step",
            #             "output_frequency": 1.0,
            #             "body_output": true,
            #             "node_output": false,
            #             "skin_output": false,
            #             "plane_output": [],
            #             "nodal_results": ["VELOCITY","PRESSURE","MESH_VELOCITY","REACTION","DISTANCE","EMBEDDED_WET_VELOCITY","MESH_DISPLACEMENT","DISPLACEMENT","TEMPERATURE","TEMPERATURE_OLD_IT","LINEAR_MOMENTUM","ANGULAR_MOMENTUM"],
            #             "nodal_nonhistorical_results": ["EMBEDDED_IS_ACTIVE"],
            #             "elemental_conditional_flags_results": ["ACTIVE","INTERFACE","TO_SPLIT"],
            #             "gauss_point_results": ["EMBEDDED_VELOCITY"]
            #         }
            #     }""")
            # )
            # fluid_gid_output.ExecuteInitialize()
            # fluid_gid_output.ExecuteBeforeSolutionLoop()
            # fluid_gid_output.ExecuteInitializeSolutionStep()
            # fluid_gid_output.PrintOutput()
            # fluid_gid_output.ExecuteFinalizeSolutionStep()
            # fluid_gid_output.ExecuteFinalize()

            # virtual_fluid_gid_output = gid_output_process.GiDOutputProcess(
            #     self.model.GetModelPart("VirtualModelPart"),
            #     "virtual_fluid_step_" + str(self.model.GetModelPart("VirtualModelPart").ProcessInfo[KratosMultiphysics.STEP]) + "_it_" + str(nl_it),
            #     KratosMultiphysics.Parameters("""{
            #         "result_file_configuration": {
            #             "gidpost_flags": {
            #                 "GiDPostMode": "GiD_PostBinary",
            #                 "WriteDeformedMeshFlag": "WriteDeformed",
            #                 "WriteConditionsFlag": "WriteConditions",
            #                 "MultiFileFlag": "SingleFile"
            #             },
            #             "file_label": "time",
            #             "output_control_type": "step",
            #             "output_frequency": 1.0,
            #             "body_output": true,
            #             "node_output": false,
            #             "skin_output": false,
            #             "plane_output": [],
            #             "nodal_results": ["VELOCITY","PRESSURE","MESH_VELOCITY","MESH_DISPLACEMENT"]
            #         }
            #     }""")
            # )
            # virtual_fluid_gid_output.ExecuteInitialize()
            # virtual_fluid_gid_output.ExecuteBeforeSolutionLoop()
            # virtual_fluid_gid_output.ExecuteInitializeSolutionStep()
            # virtual_fluid_gid_output.PrintOutput()
            # virtual_fluid_gid_output.ExecuteFinalizeSolutionStep()
            # virtual_fluid_gid_output.ExecuteFinalize()

            # structure_gid_output = gid_output_process.GiDOutputProcess(
            #     self.structure_solver.GetComputingModelPart(),
            #     "structure_step_" + str(self.structure_solver.GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP]) + "_it_" + str(nl_it),
            #     KratosMultiphysics.Parameters("""{
            #         "result_file_configuration": {
            #             "gidpost_flags": {
            #                 "GiDPostMode": "GiD_PostBinary",
            #                 "WriteDeformedMeshFlag": "WriteDeformed",
            #                 "WriteConditionsFlag": "WriteConditions",
            #                 "MultiFileFlag": "SingleFile"
            #             },
            #             "file_label": "time",
            #             "output_control_type": "step",
            #             "output_frequency": 1.0,
            #             "body_output": true,
            #             "node_output": false,
            #             "skin_output": false,
            #             "plane_output": [],
            #             "nodal_results": ["DISPLACEMENT","VELOCITY","STRUCTURE_VELOCITY","ACCELERATION","REACTION","POSITIVE_FACE_PRESSURE","VECTOR_PROJECTED","FSI_INTERFACE_RESIDUAL"],
            #             "gauss_point_results": []
            #         }
            #     }""")
            # )
            # structure_gid_output.ExecuteInitialize()
            # structure_gid_output.ExecuteBeforeSolutionLoop()
            # structure_gid_output.ExecuteInitializeSolutionStep()
            # structure_gid_output.PrintOutput()
            # structure_gid_output.ExecuteFinalizeSolutionStep()
            # structure_gid_output.ExecuteFinalize()

            # Compute the residual and perform the update
            dis_residual_norm = self._get_fsi_coupling_interface_structure().Update()

            # End the FSI non-linear iteration
            self._get_convergence_accelerator().FinalizeNonLinearIteration()

            # Check convergence
            is_converged = self._check_FSI_convergence(dis_residual_norm)
            if (is_converged):
                return True

            # Compute the structure interface VELOCITY residual vector
            # dis_residual = self._compute_displacement_residual()

            # # Check convergence
            # nl_res_norm = self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM]
            # interface_dofs = self._get_partitioned_fsi_utilities().GetInterfaceResidualSize(self._GetStructureInterfaceSubmodelPart())
            # if nl_res_norm/sqrt(interface_dofs) < self.nl_tol:
            #     self._get_convergence_accelerator().FinalizeNonLinearIteration()
            #     KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', '\tNon-linear iteration convergence achieved')
            #     KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', '\tTotal non-linear iterations: ' + str(nl_it) + ' |res|/sqrt(nDOFS) = ' + str(nl_res_norm/sqrt(interface_dofs)))
            #     break
            # else:
            #     # If convergence is not achieved, perform the correction of the prediction
            #     KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', '\tResidual computation finished. |res|/sqrt(nDOFS) = ' + str(nl_res_norm/sqrt(interface_dofs)))
            #     self._get_convergence_accelerator().UpdateSolution(dis_residual, self.iteration_value)
            #     self._get_convergence_accelerator().FinalizeNonLinearIteration()

            #     if (nl_it == self.max_nl_it):
            #         KratosMultiphysics.Logger.PrintWarning('PartitionedEmbeddedFSIBaseSolver', '\tFSI NON-LINEAR ITERATION CONVERGENCE NOT ACHIEVED')

        return False

    def FinalizeSolutionStep(self):
        # Finalize solution step
        self.fluid_solver.FinalizeSolutionStep()
        self.structure_solver.FinalizeSolutionStep()
        self._get_convergence_accelerator().FinalizeSolutionStep()

        # TODO: REMOVE AFTER DEBUGGING
        for node in self.fluid_solver.GetComputingModelPart().Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.LINEAR_MOMENTUM, node.GetSolutionStepValue(KratosMultiphysics.VELOCITY, 1))
            node.SetSolutionStepValue(KratosMultiphysics.ANGULAR_MOMENTUM, node.GetSolutionStepValue(KratosMultiphysics.VELOCITY, 2))
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, node.GetSolutionStepValue(KratosMultiphysics.PRESSURE, 1))
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE_OLD_IT, node.GetSolutionStepValue(KratosMultiphysics.PRESSURE, 2))

    def SetEchoLevel(self, structure_echo_level, fluid_echo_level):
        self.fluid_solver.SetEchoLevel(self, fluid_echo_level)
        self.structure_solver.SetEchoLevel(self, structure_echo_level)

    def SetTimeStep(self, step):
        self.fluid_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
        self.structure_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)

    def Clear(self):
        self.fluid_solver.Clear()
        self.structure_solver.Clear()

    def Check(self):
        self.fluid_solver.Check()
        self.structure_solver.Check()

    #######################################################################
    ##############          PRIVATE METHODS SECTION          ##############
    #######################################################################

    def _get_distance_to_skin_process(self, update_distance_process = False):
        if update_distance_process:
            self._distance_to_skin_process = self._create_distance_to_skin_process()
        return self._distance_to_skin_process

    def _create_distance_to_skin_process(self):
        # Set the distance computation process
        if self.domain_size == 2:
            return KratosMultiphysics.CalculateDistanceToSkinProcess2D(
                self.GetFluidComputingModelPart(),
                self._get_fsi_coupling_interface_fluid().GetInterfaceModelPart())
        elif self.domain_size == 3:
            return KratosMultiphysics.CalculateDistanceToSkinProcess3D(
                self.GetFluidComputingModelPart(),
                self._get_fsi_coupling_interface_fluid().GetInterfaceModelPart())
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))

    # def _create_distance_to_skin_process(self):
    #     # Set the distance computation process
    #     if self.domain_size == 2:
    #         return KratosMultiphysics.CalculateDistanceToSkinProcess2D(
    #             self.GetFluidComputingModelPart(),
    #             self.GetStructureSkinElementBasedModelPart())
    #     elif self.domain_size == 3:
    #         return KratosMultiphysics.CalculateDistanceToSkinProcess3D(
    #             self.GetFluidComputingModelPart(),
    #             self.GetStructureSkinElementBasedModelPart())
    #     else:
    #         raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))

    def _get_variational_distance_process(self):
        if not hasattr(self, '_variational_distance_process'):
            self._variational_distance_process = self._create_variational_distance_process()
        return self._variational_distance_process

    def _create_variational_distance_process(self):
        # Set the variational distance calculation process
        # So far, we are using the same linear solver used in the CFD problem resolution
        maximum_iterations = 2 #TODO: Make this user-definable
        if self.domain_size == 2:
            return KratosMultiphysics.VariationalDistanceCalculationProcess2D(
                self.GetFluidComputingModelPart(),
                self.fluid_solver.linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess2D.NOT_CALCULATE_EXACT_DISTANCES_TO_PLANE)
        elif self.domain_size == 3:
            return KratosMultiphysics.VariationalDistanceCalculationProcess3D(
                self.GetFluidComputingModelPart(),
                self.fluid_solver.linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess3D.NOT_CALCULATE_EXACT_DISTANCES_TO_PLANE)
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))

    def _get_distance_modification_process(self):
        if not hasattr(self, '_distance_modification_process'):
            self._distance_modification_process = self._create_distance_modification_process()
        return self._distance_modification_process

    def _create_distance_modification_process(self):
        distance_modification_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "distance_factor": 2.0,
            "distance_threshold": 1e-5,
            "continuous_distance": true,
            "check_at_each_time_step": true,
            "avoid_almost_empty_elements": true,
            "deactivate_full_negative_elements": true
        }''')
        #TODO: MODIFY THIS TO RETRIEVE THE MODEL PART NAME FROM THE JSON
        distance_modification_settings["model_part_name"].SetString("Parts_Fluid")
        return KratosFluid.DistanceModificationProcess(self.model, distance_modification_settings)

    # def _get_embedded_skin_utility(self):
    #     if not hasattr(self, '_embedded_skin_utility'):
    #         self._embedded_skin_utility = self._create_embedded_skin_utility()
    #     return self._embedded_skin_utility

    # def _create_embedded_skin_utility(self):
    #     level_set_type = "continuous"
    #     if self.domain_size == 2:
    #         return KratosMultiphysics.EmbeddedSkinUtility2D(
    #             self.GetFluidComputingModelPart(),
    #             self.GetStructureIntersectionsModelPart(),
    #             level_set_type)
    #     elif self.domain_size == 3:
    #         return KratosMultiphysics.EmbeddedSkinUtility3D(
    #             self.GetFluidComputingModelPart(),
    #             self.GetStructureIntersectionsModelPart(),
    #             level_set_type)
    #     else:
    #         raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))

    def _get_structure_interface_model_part_name(self):
        str_int_list = self.settings["coupling_settings"]["structure_interfaces_list"]
        if (str_int_list.size() != 1):
            raise Exception("FSI embedded solver structure skin must be contained in a unique model part")
        return str_int_list[0].GetString()

    def _initialize_iteration_value_vector(self):
        # Note that the embedded FSI problem is defined in terms of the structure interface
        # Initialize the iteration value for the residual computation
        str_int_res_size = self._get_partitioned_fsi_utilities().GetInterfaceResidualSize(self._GetStructureInterfaceSubmodelPart())
        self.iteration_value = KratosMultiphysics.Vector(str_int_res_size)
        for i in range(0,str_int_res_size):
            self.iteration_value[i] = 0.0

    def _map_pressure_to_positive_face_pressure(self):
        # # Set fake pressure field
        # time = self._GetStructureInterfaceSubmodelPart().ProcessInfo[KratosMultiphysics.TIME]
        # for node in self._GetStructureInterfaceSubmodelPart().Nodes:
        #     node.SetSolutionStepValue(KratosMultiphysics.POSITIVE_FACE_PRESSURE, 5.0 * node.Y * time)

        self._get_partitioned_fsi_utilities().EmbeddedPressureToPositiveFacePressureInterpolator(
            self.GetFluidComputingModelPart(),
            self._GetStructureInterfaceSubmodelPart())

    def _initialize_fsi_interfaces(self):
        # Initialize Neumann structure interface
        str_interface_submodelpart = self.model.GetModelPart(self._get_structure_interface_model_part_name())

        # Set the INTERFACE flag to the structure skin
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, str_interface_submodelpart.Nodes)

    def _update_level_set(self):
        # Recompute the distance field with the obtained solution
        self._get_distance_to_skin_process().Execute()

        # Extend the distance field but keeping the zero level set
        self._get_variational_distance_process().Execute()

        # Correct the distance field
        self._get_distance_modification_process().Execute()

        # Recompute the new embedded intersections model part
        # self._get_embedded_skin_utility().GenerateSkin()

    def _solve_fluid(self):
        # Update the current iteration level-set position
        self._update_level_set()

        # # VELOCITY based EMBEDDED_VELOCITY.
        # self._get_distance_to_skin_process().CalculateEmbeddedVariableFromSkin(
        #     KratosMultiphysics.VELOCITY,
        #     KratosMultiphysics.EMBEDDED_VELOCITY)

        # # Linearised EMBEDDED_VELOCITY.
        # dt = self._GetStructureInterfaceSubmodelPart().ProcessInfo[KratosMultiphysics.DELTA_TIME]
        # for node in self._GetStructureInterfaceSubmodelPart().Nodes:
        #     d_0 = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0)
        #     d_1 = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 1)
        #     node.SetSolutionStepValue(KratosMultiphysics.STRUCTURE_VELOCITY, (d_0 - d_1) / dt)

        # self._get_distance_to_skin_process().CalculateEmbeddedVariableFromSkin(
        #     KratosMultiphysics.STRUCTURE_VELOCITY,
        #     KratosMultiphysics.EMBEDDED_VELOCITY)

        # Solve fluid problem
        self.fluid_solver.SolveSolutionStep() # This contains the FM-ALE operations

        # Undo the FM-ALE virtual mesh movement
        self.fluid_solver._get_mesh_moving_util().UndoMeshMovement()

        # Restore the fluid node fixity to its original status
        self._get_distance_modification_process().ExecuteFinalizeSolutionStep()

    def _solve_structure(self):
        # Interpolate the intersections model part PRESSURE values from the fluid mesh
        # self._get_embedded_skin_utility().InterpolateMeshVariableToSkin(
        #     KratosMultiphysics.PRESSURE,
        #     KratosMultiphysics.PRESSURE)

        # Map the intersections model part PRESSURE values to the structure skin
        # self._map_pressure_to_positive_face_pressure()

        # Solve the structure problem
        self.structure_solver.SolveSolutionStep()

    def _compute_displacement_residual(self):
        # Compute the structure interface residual vector. The residual vector is computed as the difference
        # between the obtained velocity value (VELOCITY) and the velocity in the beggining of the
        # iteration, which has been previously saved in the auxiliary variable VECTOR_PROJECTED. Besides, the
        # residual norm is stored in the ProcessInfo.
        dis_residual = KratosMultiphysics.Vector(
            self._get_partitioned_fsi_utilities().GetInterfaceResidualSize(self._GetStructureInterfaceSubmodelPart()))

        self._get_partitioned_fsi_utilities().ComputeInterfaceResidualVector(
            self._GetStructureInterfaceSubmodelPart(),
            KratosMultiphysics.DISPLACEMENT,
            KratosMultiphysics.VECTOR_PROJECTED,
            KratosMultiphysics.FSI_INTERFACE_RESIDUAL,
            dis_residual,
            "nodal",
            KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM)
        # self._get_partitioned_fsi_utilities().ComputeInterfaceResidualVector(
        #     self._GetStructureInterfaceSubmodelPart(),
        #     KratosMultiphysics.VECTOR_PROJECTED,
        #     KratosMultiphysics.DISPLACEMENT,
        #     KratosMultiphysics.FSI_INTERFACE_RESIDUAL,
        #     dis_residual,
        #     "nodal",
        #     KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM)

        return dis_residual

    def _check_FSI_convergence(self, residual_norm):
        interface_dofs = self._get_partitioned_fsi_utilities().GetInterfaceResidualSize(self._GetStructureInterfaceSubmodelPart())
        normalised_residual = residual_norm/sqrt(interface_dofs)
        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', '\t|res|/sqrt(nDOFS) = ' + str(normalised_residual))
        if normalised_residual < self.nl_tol:
            return True
        else:
            return False

    # This method returns the convergence accelerator.
    # If it is not created yet, it calls the _create_convergence_accelerator first
    def _get_convergence_accelerator(self):
        if not hasattr(self, '_convergence_accelerator'):
            self._convergence_accelerator = self._create_convergence_accelerator()
        return self._convergence_accelerator

    # This method constructs the convergence accelerator coupling utility
    def _create_convergence_accelerator(self):
        convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(self.settings["coupling_settings"]["coupling_strategy_settings"])
        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Coupling strategy construction finished')
        return convergence_accelerator

    def _get_fsi_coupling_interface_structure(self):
        if not hasattr(self, '_fsi_coupling_interface_structure'):
            self._fsi_coupling_interface_structure = self._create_fsi_coupling_interface_structure()
        return self._fsi_coupling_interface_structure

    def _create_fsi_coupling_interface_structure(self):
        # Set auxiliary settings
        aux_settings = KratosMultiphysics.Parameters(
        """{
            "model_part_name": "FSICouplingInterfaceStructure",
            "father_model_part_name": "",
            "input_variable_name": "POSITIVE_FACE_PRESSURE",
            "output_variable_name": "DISPLACEMENT"
        }""")
        aux_settings["father_model_part_name"].SetString(self.structure_interface_submodelpart_name)

        # Construct the FSI coupling interface
        fsi_coupling_interface_structure = fsi_coupling_interface.FSICouplingInterface(
            self.model,
            aux_settings,
            self._get_convergence_accelerator())

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Structure FSI coupling interface created')

        return fsi_coupling_interface_structure

    def _get_fsi_coupling_interface_fluid(self):
        if not hasattr(self, '_fsi_coupling_interface_fluid'):
            self._fsi_coupling_interface_fluid = self._create_fsi_coupling_interface_fluid()
        return self._fsi_coupling_interface_fluid

    def _create_fsi_coupling_interface_fluid(self):
        # Set auxiliary settings
        # Note that in the embedded case, the fluid interface is identical to the structure one
        # This is intentionally done, since this copy will be used in the level set computation
        aux_settings = KratosMultiphysics.Parameters(
        """{
            "model_part_name": "FSICouplingInterfaceFluid",
            "father_model_part_name": "",
            "input_variable_name": "DISPLACEMENT",
            "output_variable_name": "PRESSURE"
        }""")
        aux_settings["father_model_part_name"].SetString(self.structure_interface_submodelpart_name)

        # Construct the FSI coupling interface
        fsi_coupling_interface_fluid = fsi_coupling_interface.FSICouplingInterface(
            self.model,
            aux_settings)

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Fluid FSI coupling interface created')

        return fsi_coupling_interface_fluid

    # This method finds the maximum buffer size between mesh,
    # fluid and structure solvers and sets it to all the solvers.
    def _GetAndSetMinimumBufferSize(self):
        fluid_buffer_size = self.fluid_solver.min_buffer_size
        str_buffer_size = self.structure_solver.settings["buffer_size"].GetInt()

        buffer_size = max(fluid_buffer_size, str_buffer_size)

        self.fluid_solver.min_buffer_size = buffer_size
        self.structure_solver.settings["buffer_size"].SetInt(buffer_size)

    def _GetStructureInterfaceSubmodelPart(self):
        # Returns the structure interface submodelpart that will be used in the residual minimization
        return self.structure_solver.main_model_part.GetSubModelPart(self.structure_interface_submodelpart_name)

    def _GetDomainSize(self):
        fluid_domain_size = self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        structure_domain_size = self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        if fluid_domain_size !=structure_domain_size:
            raise("ERROR: Solid domain size and fluid domain size are not equal!")

        return fluid_domain_size

    def  _ComputeDeltaTime(self):
        fluid_time_step = self.fluid_solver._ComputeDeltaTime()
        structure_time_step = self.structure_solver.ComputeDeltaTime()

        if abs(fluid_time_step - structure_time_step) > 1e-12:
            err_msg =  'Fluid time step is: ' + str(fluid_time_step) + '\n'
            err_msg += 'Structure time step is: ' + str(structure_time_step) + '\n'
            err_msg += 'No substepping has been implemented yet. Fluid and structure time step must coincide.'
            raise Exception(err_msg)

        return fluid_time_step

    def _get_partitioned_fsi_utilities(self):
        if not hasattr(self, '_partitioned_fsi_utilities'):
            self._partitioned_fsi_utilities = self._create_partitioned_fsi_utilities()
        return self._partitioned_fsi_utilities

    def _create_partitioned_fsi_utilities(self):
        if self.domain_size == 2:
            return KratosFSI.PartitionedFSIUtilitiesArray2D()
        elif self.domain_size == 3:
            return KratosFSI.PartitionedFSIUtilitiesArray3D()
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))
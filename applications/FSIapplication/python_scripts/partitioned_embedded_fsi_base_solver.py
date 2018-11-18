from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from math import sqrt   # Import the square root from python library

# Import utilities
import NonConformant_OneSideMap                # Import non-conformant mapper
import python_solvers_wrapper_fluid            # Import the fluid Python solvers wrapper
import python_solvers_wrapper_structural       # Import the structure Python solvers wrapper
import python_solvers_wrapper_mesh_motion      # Import the mesh motion Python solvers wrapper
import convergence_accelerator_factory         # Import the FSI convergence accelerator factory

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Import base class file
from python_solver import PythonSolver

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
        self.fluid_interface_submodelpart_name = coupling_settings["fluid_interfaces_list"][0].GetString()
        self.structure_interface_submodelpart_name = coupling_settings["structure_interfaces_list"][0].GetString()

        # Construct the structure solver
        self.structure_solver = python_solvers_wrapper_structural.CreateSolverByParameters(self.model, self.settings["structure_solver_settings"], self.parallel_type)
        self._PrintInfoOnRankZero("::[PartitionedEmbeddedFSIBaseSolver]::", "Structure solver construction finished.")

        # Construct the fluid solver
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_solver_settings"], self.parallel_type)
        self._PrintInfoOnRankZero("::[PartitionedEmbeddedFSIBaseSolver]::", "Fluid solver construction finished.")
        self._PrintInfoOnRankZero("::[PartitionedEmbeddedFSIBaseSolver]::", "Partitioned embedded FSI base solver construction finished.")

    def GetMinimumBufferSize(self):
        buffer_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_structure = self.structure_solver.GetMinimumBufferSize()
        return max(buffer_structure,buffer_fluid)

    #TODO: CHECK WHICH VARIABLES ARE MISSING OR NOT NEEDED
    def AddVariables(self):
        ## Structure variables addition
        # Standard CSM variables addition
        self.structure_solver.AddVariables()

        ## Fluid variables addition
        # Standard CFD variables addition
        self.fluid_solver.AddVariables()
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_ACCELERATION) # TODO: This should be added in the mesh solvers

        ## FSIApplication variables addition
        NonConformant_OneSideMap.AddVariables(self.fluid_solver.main_model_part,self.structure_solver.main_model_part)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_PROJECTED)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FSI_INTERFACE_RESIDUAL)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FSI_INTERFACE_MESH_RESIDUAL)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_MAPPED_VECTOR_VARIABLE)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_MAPPED_VECTOR_VARIABLE)
        self.structure_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_PROJECTED)

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

    def AddDofs(self):
        # Add DOFs structure
        self.structure_solver.AddDofs()
        # Add DOFs fluid
        self.fluid_solver.AddDofs()

    def Initialize(self):
        # Get the domain size
        self.domain_size = self._GetDomainSize()

        # Python structure solver initialization
        self.structure_solver.Initialize()

        # Python fluid solver initialization
        self.fluid_solver.Initialize()

        # Initialize the Dirichlet-Neumann interface
        self._initialize_fsi_interfaces()

        # Set the Neumann B.C. in the structure interface
        #TODO: SHOULD THIS BE DONE BY THE INTERFACE?
        self._set_structure_neumann_condition()

        # Initialize the iteration value vector
        self._initialize_iteration_value_vector()

        # Compute the fluid domain NODAL_AREA values (required as weight in the residual norm computation)
        KratosMultiphysics.CalculateNodalAreaProcess(self.fluid_solver.GetComputingModelPart(), self.domain_size).Execute()

        # Coupling utility initialization
        # The _get_convergence_accelerator is supposed to construct the convergence accelerator in here
        self._get_convergence_accelerator().Initialize()

    def AdvanceInTime(self, current_time):
        fluid_new_time = self.fluid_solver.AdvanceInTime(current_time)
        structure_new_time = self.structure_solver.AdvanceInTime(current_time)

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

        # Recompute the distance field
        self._get_distance_to_skin_process().Execute()

        # Recompute the embedded intersections model part
        self._get_embedded_skin_utility().GenerateSkin()

    def Predict(self):
        # Perform fluid and structure solvers predictions
        self.fluid_solver.Predict()
        self.structure_solver.Predict()

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
        return self.model.GetSubModelPart(self._get_structure_interface_model_part_name())

    def GetStructureIntersectionsModelPart(self):
        if not hasattr(self, '_embedded_intersections_model_part'):
            self._embedded_intersections_model_part = self.model.CreateModelPart("EmbeddedIntersectionsModelPart")
        return self._embedded_intersections_model_part

    def GetOutputVariables(self):
        pass

    def SaveRestart(self):
        pass

    def SolveSolutionStep(self):
        ## Non-Linear interface coupling iteration ##
        for nl_it in range(1,self.max_nl_it+1):

            self._PrintInfoOnRankZero("","\tFSI non-linear iteration = ", nl_it)

            self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it
            self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it

            self._get_convergence_accelerator().InitializeNonLinearIteration()

            # Update the EMBEDDED_VELOCITY and solve the fluid problem
            self._solve_fluid()

            # Update the PRESSURE load and solve the structure problem 
            self._solve_structure()

            # Compute the structure interface VELOCITY residual vector
            vel_residual = self._compute_velocity_residual()

            # Check convergence
            nl_res_norm = self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.FSI_INTERFACE_RESIDUAL_NORM]
            interface_dofs = self._get_partitioned_fsi_utilities().GetInterfaceResidualSize(self._GetStructureInterfaceSubmodelPart())
            if nl_res_norm/sqrt(interface_dofs) < self.nl_tol:
                self._get_convergence_accelerator().FinalizeNonLinearIteration()
                self._PrintInfoOnRankZero("","\tNon-linear iteration convergence achieved")
                self._PrintInfoOnRankZero("","\tTotal non-linear iterations: ", nl_it, " |res|/sqrt(nDOFS) = ", nl_res_norm/sqrt(interface_dofs))
                break
            else:
                # If convergence is not achieved, perform the correction of the prediction
                self._PrintInfoOnRankZero("","\tResidual computation finished. |res|/sqrt(nDOFS) = ", nl_res_norm/sqrt(interface_dofs))
                self._get_convergence_accelerator().UpdateSolution(vel_residual, self.iteration_value)
                self._get_convergence_accelerator().FinalizeNonLinearIteration()

                if (nl_it == self.max_nl_it):
                    self._PrintWarningOnRankZero("","\tFSI NON-LINEAR ITERATION CONVERGENCE NOT ACHIEVED")

        ## Finalize solution step
        self.fluid_solver.FinalizeSolutionStep()
        self.structure_solver.FinalizeSolutionStep()
        self._get_convergence_accelerator().FinalizeSolutionStep()

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

    def _get_distance_to_skin_process(self):
        if not hasattr(self, '_distance_to_skin_process'):
            self._distance_to_skin_process = _create_distance_to_skin_process()
        return self._distance_to_skin_process

    def _create_distance_to_skin_process(self):
        if self.domain_size == 2:
            return KratosMultiphysics.CalculateDistanceToSkinProcess2D(
                self.GetFluidComputingModelPart(),
                self.GetStructureSkinModelPart())
        elif self.domain_size == 3:
            return KratosMultiphysics.CalculateDistanceToSkinProcess3D(
                self.GetFluidComputingModelPart(),
                self.GetStructureSkinModelPart())
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))

    def _get_embedded_skin_utility(self):
        if not hasattr(self, '_embedded_skin_utility'):
            self._embedded_skin_utility = self._create_embedded_skin_utility()
        return self._embedded_skin_utility

    def _create_embedded_skin_utility(self):
        level_set_type = "continuous"
        if self.domain_size == 2:
            return KratosMultiphysics.EmbeddedSkinUtility2D(
                self.GetFluidComputingModelPart(), 
                self.GetStructureIntersectionsModelPart(),
                level_set_type)
        elif self.domain_size == 3:
            return KratosMultiphysics.EmbeddedSkinUtility3D(
                self.GetFluidComputingModelPart(), 
                self.GetStructureIntersectionsModelPart(),
                level_set_type)
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))

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
        # Maps the PRESSURE value from the generated intersections skins to the structure skin.
        # Note that the SurfaceLoadCondition in the StructuralMechanicsApplication uses the 
        # variable POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE as nodal surface loads.
        map_parameters = KratosMultiphysics.Parameters("""{
            "echo_level"                       : 0,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2
        }""")

        if self.domain_size == 2:
            KratosMultiphysics.SimpleMortarMapperProcess2D2NDouble(
                self.GetStructureIntersectionsModelPart(),
                self.GetStructureSkinModelPart(),
                KratosMultiphysics.PRESSURE,
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                map_parameters).Execute()
        elif self.domain_size == 3:
            #TODO: INCLUDE THE MAPPING AGAINST TRIANGLE CONDITIONS!
            #TODO: WILL BE NEEDED TO MAP TO SHELLS (DOUBLE SIDED SURFACES)
            KratosMultiphysics.SimpleMortarMapperProcess3D3N4NDouble(
                self.GetStructureIntersectionsModelPart(),
                self.GetStructureSkinModelPart(),
                KratosMultiphysics.PRESSURE,
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                map_parameters).Execute()
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))

    def _map_velocity_to_vector_projected(self):
        # Maps the PRESSURE value from the generated intersections skins to the structure skin.
        # Note that the SurfaceLoadCondition in the StructuralMechanicsApplication uses the 
        # variable POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE as nodal surface loads.
        map_parameters = KratosMultiphysics.Parameters("""{
            "echo_level"                       : 0,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2
        }""")

        if self.domain_size == 2:
            KratosMultiphysics.SimpleMortarMapperProcess2D2NDouble(
                self.GetStructureIntersectionsModelPart(),
                self.GetStructureSkinModelPart(),
                KratosMultiphysics.VELOCITY,
                KratosMultiphysics.VECTOR_PROJECTED,
                map_parameters).Execute()
        elif self.domain_size == 3:
            #TODO: INCLUDE THE MAPPING AGAINST TRIANGLE CONDITIONS!
            #TODO: WILL BE NEEDED TO MAP TO SHELLS (DOUBLE SIDED SURFACES)
            KratosMultiphysics.SimpleMortarMapperProcess3D3N4NDouble(
                self.GetStructureIntersectionsModelPart(),
                self.GetStructureSkinModelPart(),
                KratosMultiphysics.VELOCITY,
                KratosMultiphysics.VECTOR_PROJECTED,
                map_parameters).Execute()
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))

    def _initialize_fsi_interfaces(self):
        # Initialize Neumann structure interface
        str_interface_submodelpart = self.model.GetModelPart(self._get_structure_interface_model_part_name())

        # Set the INTERFACE flag to the structure skin
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, str_interface_submodelpart.Nodes)

        # Initialize Dirichlet fluid interface (generate the intersections skin)
        #TODO
        self._get_embedded_skin_utility().GenerateSkin()

    def _solve_fluid(self):
        # Update de EMBEDDED_VELOCITY value
        self._get_distance_to_skin_process().CalculateEmbeddedVariableFromSkin(
            KratosMultiphysics.VELOCITY,
            KratosMultiphysics.EMBEDDED_VELOCITY)

        # Solve fluid problem
        self.fluid_solver.SolveSolutionStep()

        # First, interpolate the intersections model part VELOCITY values from the fluid mesh
        # Then, map the computed velocity to an auxiliar variable in the structure skin.
        # This will be used once the structure is solved to compute the interface residual.
        self._get_embedded_skin_utility().InterpolateMeshVariableToSkin(KratosMultiphysics.VELOCITY, KratosMultiphysics.VELOCITY)
        self._map_velocity_to_vector_projected()

    def _solve_structure(self):
        # Interpolate the intersections model part PRESSURE values from the fluid mesh
        # Then map the intersections model part PRESSURE values to the structure skin
        self._get_embedded_skin_utility().InterpolateMeshVariableToSkin(KratosMultiphysics.PRESSURE, KratosMultiphysics.PRESSURE)
        self._map_pressure_to_positive_face_pressure().Execute()

        # Solve the structure problem
        self.structure_solver.SolveSolutionStep()

    def _compute_velocity_residual(self):
        # Compute the structure interface residual vector. The residual vector is computed 
        # as the difference between the obtained velocity value (VELOCITY) and the velocity 
        # in the beggining of the iteration. This is mapped from the fluid and saved in 
        # VECTOR_PROJECTED. Besides, the residual norm is stored in the ProcessInfo.
        vel_residual = KratosMultiphysics.Vector(
            self._get_partitioned_fsi_utilities().GetInterfaceResidualSize(self._GetStructureInterfaceSubmodelPart()))

        self._get_partitioned_fsi_utilities().ComputeInterfaceResidualVector(
            self._GetStructureInterfaceSubmodelPart(),
            KratosMultiphysics.VECTOR_PROJECTED,
            KratosMultiphysics.VELOCITY,
            vel_residual)

        return vel_residual

    # This method is to be overwritten in the MPI solver
    def _PrintInfoOnRankZero(self, *args):
        KratosMultiphysics.Logger.PrintInfo(" ".join(map(str, args)))

    # This method is to be overwritten in the MPI solver
    def _PrintWarningOnRankZero(self, *args):
        KratosMultiphysics.Logger.PrintWarning(" ".join(map(str, args)))

    # This method returns the convergence accelerator.
    # If it is not created yet, it calls the _create_convergence_accelerator first
    def _get_convergence_accelerator(self):
        if not hasattr(self, '_convergence_accelerator'):
            self._convergence_accelerator = self._create_convergence_accelerator()
        return self._convergence_accelerator

    # This method constructs the convergence accelerator coupling utility
    def _create_convergence_accelerator(self):
        convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(self.settings["coupling_settings"]["coupling_strategy_settings"])
        self._PrintInfoOnRankZero("::[PartitionedEmbeddedFSIBaseSolver]::", "Coupling strategy construction finished.")
        return convergence_accelerator

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
            return KratosFSI.PartitionedFSIUtilities2D()
        elif self.domain_size == 3:
            return KratosFSI.PartitionedFSIUtilities3D()
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))

    #TODO: MUST BE SET IN THE INTERFACE?¿?¿
    def _set_structure_neumann_condition(self):
        print("NOT REQUIRED?¿¿")
        # structure_computational_submodelpart = self.structure_solver.GetComputingModelPart()

        # # Get the maximum condition id
        # max_cond_id = 0
        # for condition in self.structure_solver.main_model_part.Conditions:
        #     max_cond_id = max(max_cond_id, condition.Id)

        # max_cond_id = self.structure_solver.main_model_part.GetCommunicator().MaxAll(max_cond_id)

        # # Set up the point load condition in the structure interface
        # structure_interfaces_list = self.settings["coupling_settings"]["structure_interfaces_list"]
        # for i in range(structure_interfaces_list.size()):
        #     interface_submodelpart_name = structure_interfaces_list[i].GetString()
        #     interface_submodelpart_i = self.structure_solver.main_model_part.GetSubModelPart(interface_submodelpart_name)

        #     # Get the number of conditions to be set in each processor
        #     local_nodes_number_accumulated = -1
        #     local_nodes_number = len(interface_submodelpart_i.GetCommunicator().LocalMesh().Nodes)
        #     local_nodes_number_accumulated = interface_submodelpart_i.GetCommunicator().ScanSum(local_nodes_number, local_nodes_number_accumulated)

        #     # Create the point load condition
        #     aux_count = max_cond_id + local_nodes_number_accumulated
        #     if self.domain_size == 2:
        #         for node in interface_submodelpart_i.GetCommunicator().LocalMesh().Nodes:
        #             aux_count+=1
        #             structure_computational_submodelpart.CreateNewCondition("PointLoadCondition2D1N", 
        #                                                                     int(aux_count), 
        #                                                                     [node.Id], 
        #                                                                     self.structure_solver.main_model_part.Properties[0])
        #     elif self.domain_size == 3:
        #         for node in interface_submodelpart_i.GetCommunicator().LocalMesh().Nodes:
        #             aux_count+=1
        #             structure_computational_submodelpart.CreateNewCondition("PointLoadCondition3D1N",
        #                                                                     int(aux_count),
        #                                                                     [node.Id],
        #                                                                     self.structure_solver.main_model_part.Properties[0])
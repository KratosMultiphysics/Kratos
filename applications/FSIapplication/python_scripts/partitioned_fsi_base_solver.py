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
from python_solver import PythonSolver

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications(
    "FSIApplication",
    "MeshMovingApplication",
    "FluidDynamicsApplication",
    "StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

def CreateSolver(model, project_parameters):
    return PartitionedFSIBaseSolver(model, project_parameters)

class PartitionedFSIBaseSolver(PythonSolver):

    def _ValidateSettings(self, project_parameters):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "echo_level": 0,
            "parallel_type": "OpenMP",
            "solver_type": "partitioned",
            "coupling_scheme": "dirichlet_neumann",
            "structure_solver_settings": {
            },
            "fluid_solver_settings":{
            },
            "mesh_solver_settings":{
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
        super(PartitionedFSIBaseSolver,self).__init__(model, project_parameters)

        # Auxiliar variables
        self.parallel_type = self.settings["parallel_type"].GetString()
        coupling_settings = self.settings["coupling_settings"]
        self.max_nl_it = coupling_settings["nl_max_it"].GetInt()
        self.nl_tol = coupling_settings["nl_tol"].GetDouble()
        self.solve_mesh_at_each_iteration = coupling_settings["solve_mesh_at_each_iteration"].GetBool()
        self.fluid_interface_submodelpart_name = coupling_settings["fluid_interfaces_list"][0].GetString()
        self.structure_interface_submodelpart_name = coupling_settings["structure_interfaces_list"][0].GetString()

        # Construct the structure solver
        self.structure_solver = python_solvers_wrapper_structural.CreateSolverByParameters(self.model, self.settings["structure_solver_settings"], self.parallel_type)
        self._PrintInfoOnRankZero("::[PartitionedFSIBaseSolver]::", "Structure solver construction finished.")

        # Construct the fluid solver
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_solver_settings"], self.parallel_type)
        self._PrintInfoOnRankZero("::[PartitionedFSIBaseSolver]::", "Fluid solver construction finished.")

        # Construct the ALE mesh solver
        self.mesh_solver = python_solvers_wrapper_mesh_motion.CreateSolverByParameters(self.model, self.settings["mesh_solver_settings"], self.parallel_type)
        self._PrintInfoOnRankZero("::[PartitionedFSIBaseSolver]::", "ALE mesh solver construction finished.")
        self._PrintInfoOnRankZero("::[PartitionedFSIBaseSolver]::", "Partitioned FSI base solver construction finished.")

    def GetMinimumBufferSize(self):
        # Get structure buffer size
        buffer_structure = self.structure_solver.GetMinimumBufferSize()
        # Get fluid buffer size
        buffer_fluid = self.fluid_solver.GetMinimumBufferSize()

        return max(buffer_structure,buffer_fluid)

    def AddVariables(self):
        ## Structure variables addition
        # Standard CSM variables addition
        self.structure_solver.AddVariables()

        ## Fluid variables addition
        # Standard CFD variables addition
        self.fluid_solver.AddVariables()
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_ACCELERATION) # TODO: This should be added in the mesh solvers
        # Mesh solver variables addition
        self.mesh_solver.AddVariables()

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
        self.mesh_solver.AddDofs()

    def Initialize(self):
        err_msg =  'Calling the base partitioned FSI solver Initialize() method.\n'
        err_msg += 'Implement the custom Initialize() method in the derived solver.'
        raise Exception(err_msg)

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
        self._GetConvergenceAccelerator().InitializeSolutionStep()

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

    def GetOutputVariables(self):
        pass

    def SaveRestart(self):
        pass

    def SolveSolutionStep(self):
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

            self._PrintInfoOnRankZero("","\tFSI non-linear iteration = ", nl_it)

            self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it
            self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it

            self._GetConvergenceAccelerator().InitializeNonLinearIteration()

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
            interface_dofs = self.partitioned_fsi_utilities.GetInterfaceResidualSize(self._GetFluidInterfaceSubmodelPart())

            # Check convergence
            if nl_res_norm/sqrt(interface_dofs) < self.nl_tol:
                self._GetConvergenceAccelerator().FinalizeNonLinearIteration()
                self._PrintInfoOnRankZero("","\tNon-linear iteration convergence achieved")
                self._PrintInfoOnRankZero("","\tTotal non-linear iterations: ", nl_it, " |res|/sqrt(nDOFS) = ", nl_res_norm/sqrt(interface_dofs))
                break
            else:
                # If convergence is not achieved, perform the correction of the prediction
                self._PrintInfoOnRankZero("","\tResidual computation finished. |res|/sqrt(nDOFS) = ", nl_res_norm/sqrt(interface_dofs))
                self._GetConvergenceAccelerator().UpdateSolution(dis_residual, self.iteration_value)
                self._GetConvergenceAccelerator().FinalizeNonLinearIteration()

                if (nl_it == self.max_nl_it):
                    self._PrintWarningOnRankZero("","\tFSI NON-LINEAR ITERATION CONVERGENCE NOT ACHIEVED")

        ## Compute the mesh residual as final testing (it is expected to be 0)
        self.partitioned_fsi_utilities.ComputeFluidInterfaceMeshVelocityResidualNorm(self._GetFluidInterfaceSubmodelPart())
        mesh_res_norm = self.fluid_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.FSI_INTERFACE_MESH_RESIDUAL_NORM)
        self._PrintInfoOnRankZero("","\tNL residual norm: ", nl_res_norm)
        self._PrintInfoOnRankZero("","\tMesh residual norm: ", mesh_res_norm)

        ## Finalize solution step
        self.fluid_solver.FinalizeSolutionStep()
        self.structure_solver.FinalizeSolutionStep()
        self._GetConvergenceAccelerator().FinalizeSolutionStep()

    def SetEchoLevel(self, structure_echo_level, fluid_echo_level):
        self.structure_solver.SetEchoLevel(self, structure_echo_level)
        self.fluid_solver.SetEchoLevel(self, fluid_echo_level)

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

    # This method is to be overwritten in the MPI solver
    def _PrintInfoOnRankZero(self, *args):
        KratosMultiphysics.Logger.PrintInfo(" ".join(map(str, args)))

    # This method is to be overwritten in the MPI solver
    def _PrintWarningOnRankZero(self, *args):
        KratosMultiphysics.Logger.PrintWarning(" ".join(map(str, args)))

    # This method returns the convergence accelerator.
    # If it is not created yet, it calls the _CreateConvergenceAccelerator first
    def _GetConvergenceAccelerator(self):
        if not hasattr(self, '_convergence_accelerator'):
            self._convergence_accelerator = self._CreateConvergenceAccelerator()
        return self._convergence_accelerator

    # This method constructs the convergence accelerator coupling utility
    def _CreateConvergenceAccelerator(self):
        convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(self.settings["coupling_settings"]["coupling_strategy_settings"])
        self._PrintInfoOnRankZero("::[PartitionedFSIBaseSolver]::", "Coupling strategy construction finished.")
        return convergence_accelerator

    # This method finds the maximum buffer size between mesh, 
    # fluid and structure solvers and sets it to all the solvers.
    def _GetAndSetMinimumBufferSize(self):
        fluid_buffer_size = self.fluid_solver.min_buffer_size
        mesh_buffer_size = self.mesh_solver.settings["buffer_size"].GetInt()
        str_buffer_size = self.structure_solver.settings["buffer_size"].GetInt()

        buffer_size = max(fluid_buffer_size, mesh_buffer_size)
        buffer_size = max(buffer_size, str_buffer_size)

        self.fluid_solver.min_buffer_size = buffer_size
        self.mesh_solver.settings["buffer_size"].SetInt(buffer_size)
        self.structure_solver.settings["buffer_size"].SetInt(buffer_size)

    def _GetFluidInterfaceSubmodelPart(self):
        # Returns the fluid interface submodelpart that will be used in the residual minimization
        return self.fluid_solver.main_model_part.GetSubModelPart(self.fluid_interface_submodelpart_name)

    def _GetFluidPositiveInterfaceSubmodelPart(self):
        mapper_settings = self.settings["coupling_settings"]["mapper_settings"]

        # Get the fluid interface faces submodelpart names
        for mapper_id in range(2):
            if (mapper_settings[mapper_id]["mapper_face"].GetString() == "Positive"):
                pos_face_submodelpart_name = mapper_settings[mapper_id]["fluid_interface_submodelpart_name"].GetString()

        # Returns the fluid positive interface submodelpart
        return self.fluid_solver.main_model_part.GetSubModelPart(pos_face_submodelpart_name)

    def _GetFluidNegativeInterfaceSubmodelPart(self):
        mapper_settings = self.settings["coupling_settings"]["mapper_settings"]

        # Get the fluid interface faces submodelpart names
        for mapper_id in range(2):
            if (mapper_settings[mapper_id]["mapper_face"].GetString() == "Negative"):
                neg_face_submodelpart_name = mapper_settings[mapper_id]["fluid_interface_submodelpart_name"].GetString()

        # Returns the fluid negative interface submodelpart
        return self.fluid_solver.main_model_part.GetSubModelPart(neg_face_submodelpart_name)

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

    def _GetNodalUpdateUtilities(self):
        structure_time_scheme = self.structure_solver.dynamic_settings["scheme_type"].GetString()
        if (structure_time_scheme == "newmark"):
            damp_factor_m = 0.0
        elif (structure_time_scheme == "bossak"):
            damp_factor_m = -0.3
        else:
            err_msg =  "Requested structure time scheme type \"" + structure_time_scheme + "\" is not available!\n"
            err_msg += "Available options are: \"newmark\", \"bossak\", \"relaxation\""
            raise Exception(err_msg)

        if (self.domain_size == 2):
            return KratosFSI.NodalUpdateNewmark2D(damp_factor_m)
        else:
            return KratosFSI.NodalUpdateNewmark3D(damp_factor_m)

    def _GetPartitionedFSIUtilities(self):
        if (self.domain_size == 2):
            return KratosFSI.PartitionedFSIUtilities2D()
        else:
            return KratosFSI.PartitionedFSIUtilities3D()

    def _SetUpMapper(self):
        # Recall, to set the INTERFACE flag in both the fluid and solid interface before the mapper construction
        search_radius_factor = 2.0
        mapper_max_iterations = 200
        mapper_tolerance = 1e-12

        mappers_settings = self.settings["coupling_settings"]["mapper_settings"]

        if (mappers_settings.size() == 1):
            fluid_submodelpart_name = mappers_settings[0]["fluid_interface_submodelpart_name"].GetString()
            structure_submodelpart_name = mappers_settings[0]["structure_interface_submodelpart_name"].GetString()

            fluid_submodelpart = self.fluid_solver.main_model_part.GetSubModelPart(fluid_submodelpart_name)
            structure_submodelpart = self.structure_solver.main_model_part.GetSubModelPart(structure_submodelpart_name)

            self.interface_mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(fluid_submodelpart,
                                                                                      structure_submodelpart,
                                                                                      search_radius_factor,
                                                                                      mapper_max_iterations,
                                                                                      mapper_tolerance)

            self.double_faced_structure = False

        elif (mappers_settings.size() == 2):
            # Get the fluid interface faces submodelpart names
            for mapper_id in range(2):
                if (mappers_settings[mapper_id]["mapper_face"].GetString() == "Positive"):
                    pos_face_submodelpart_name = mappers_settings[mapper_id]["fluid_interface_submodelpart_name"].GetString()
                elif (mappers_settings[mapper_id]["mapper_face"].GetString() == "Negative"):
                    neg_face_submodelpart_name = mappers_settings[mapper_id]["fluid_interface_submodelpart_name"].GetString()
                else:
                    raise Exception("Unique mapper flag has been set but more than one mapper exist in mapper_settings.")
            # Get the structure submodelpart name
            structure_submodelpart_name = mappers_settings[0]["structure_interface_submodelpart_name"].GetString()

            # Grab the interface submodelparts
            pos_fluid_submodelpart = self.fluid_solver.main_model_part.GetSubModelPart(pos_face_submodelpart_name)
            neg_fluid_submodelpart = self.fluid_solver.main_model_part.GetSubModelPart(neg_face_submodelpart_name)
            structure_submodelpart = self.structure_solver.main_model_part.GetSubModelPart(structure_submodelpart_name)

            self.interface_mapper = NonConformant_OneSideMap.NonConformantTwoFaces_OneSideMap(pos_fluid_submodelpart,
                                                                                              neg_fluid_submodelpart,
                                                                                              structure_submodelpart,
                                                                                              search_radius_factor,
                                                                                              mapper_max_iterations,
                                                                                              mapper_tolerance)

            self.double_faced_structure = True

        else:
            raise Exception("Case with more than 2 mappers has not been implemented yet.\n \
                             Please, in case you are using single faced immersed bodies, set the skin entities in a unique submodelpart.\n \
                             In case you are considering double faced immersed bodies (shells or membranes), set all the positive faces \
                             in a unique submodelpart and all the negative ones in another submodelpart.")

    def _SetStructureNeumannCondition(self):

        structure_computational_submodelpart = self.structure_solver.GetComputingModelPart()

        # Get the maximum condition id
        max_cond_id = 0
        for condition in self.structure_solver.main_model_part.Conditions:
            max_cond_id = max(max_cond_id, condition.Id)

        max_cond_id = self.structure_solver.main_model_part.GetCommunicator().MaxAll(max_cond_id)

        # Set up the point load condition in the structure interface
        structure_interfaces_list = self.settings["coupling_settings"]["structure_interfaces_list"]
        for i in range(structure_interfaces_list.size()):
            interface_submodelpart_name = structure_interfaces_list[i].GetString()
            interface_submodelpart_i = self.structure_solver.main_model_part.GetSubModelPart(interface_submodelpart_name)

            # Get the number of conditions to be set in each processor
            local_nodes_number_accumulated = -1
            local_nodes_number = len(interface_submodelpart_i.GetCommunicator().LocalMesh().Nodes)
            local_nodes_number_accumulated = interface_submodelpart_i.GetCommunicator().ScanSum(local_nodes_number, local_nodes_number_accumulated)

            # Create the point load condition
            aux_count = max_cond_id + local_nodes_number_accumulated
            if self.domain_size == 2:
                for node in interface_submodelpart_i.GetCommunicator().LocalMesh().Nodes:
                    aux_count+=1
                    structure_computational_submodelpart.CreateNewCondition("PointLoadCondition2D1N", 
                                                                            int(aux_count), 
                                                                            [node.Id], 
                                                                            self.structure_solver.main_model_part.Properties[0])
            elif self.domain_size == 3:
                for node in interface_submodelpart_i.GetCommunicator().LocalMesh().Nodes:
                    aux_count+=1
                    structure_computational_submodelpart.CreateNewCondition("PointLoadCondition3D1N",
                                                                            int(aux_count),
                                                                            [node.Id],
                                                                            self.structure_solver.main_model_part.Properties[0])

    def _ComputeMeshPredictionSingleFaced(self):

            print("Computing time step ",self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]," prediction...")
            # Get the previous step fluid interface nodal fluxes
            keep_sign = False
            distribute_load = True
            self.interface_mapper.FluidToStructure_VectorMap(KratosMultiphysics.REACTION,
                                                             KratosStructural.POINT_LOAD,
                                                             keep_sign,
                                                             distribute_load)

            # Solve the current step structure problem with the previous step fluid interface nodal fluxes
            self.structure_solver.SolveSolutionStep()

            # Map the obtained structure displacement to the fluid interface
            keep_sign = True
            distribute_load = False
            self.interface_mapper.StructureToFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                             KratosMultiphysics.MESH_DISPLACEMENT,
                                                             keep_sign,
                                                             distribute_load)

            # Solve the mesh problem
            self.mesh_solver.InitializeSolutionStep()
            self.mesh_solver.Predict()
            self.mesh_solver.SolveSolutionStep()
            self.mesh_solver.FinalizeSolutionStep()

            print("Mesh prediction computed.")


    def _ComputeMeshPredictionDoubleFaced(self):

            print("Computing time step ",self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP],"double faced prediction...")
            # Get the previous step fluid interface nodal fluxes from both positive and negative faces
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

            # Solve the current step structure problem with the previous step fluid interface nodal fluxes
            self.structure_solver.SolveSolutionStep()

            # Map the obtained structure displacement to both positive and negative fluid interfaces
            keep_sign = True
            distribute_load = False
            self.interface_mapper.StructureToPositiveFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                                     KratosMultiphysics.MESH_DISPLACEMENT,
                                                                     keep_sign,
                                                                     distribute_load)
            self.interface_mapper.StructureToNegativeFluid_VectorMap(KratosMultiphysics.DISPLACEMENT,
                                                                     KratosMultiphysics.MESH_DISPLACEMENT,
                                                                     keep_sign,
                                                                     distribute_load)

            # Solve the mesh problem
            self.mesh_solver.InitializeSolutionStep()
            self.mesh_solver.Predict()
            self.mesh_solver.SolveSolutionStep()
            self.mesh_solver.FinalizeSolutionStep()

            print("Mesh prediction computed.")

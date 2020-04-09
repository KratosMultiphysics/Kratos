from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from math import sqrt   # Import the square root from python library

# Import utilities
from KratosMultiphysics.FSIApplication import NonConformant_OneSideMap                # Import non-conformant mapper
from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid            # Import the fluid Python solvers wrapper
from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural       # Import the structure Python solvers wrapper
from KratosMultiphysics.MeshMovingApplication import python_solvers_wrapper_mesh_motion      # Import the mesh motion Python solvers wrapper
from KratosMultiphysics.FSIApplication import convergence_accelerator_factory         # Import the FSI convergence accelerator factory

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

def CreateSolver(model, project_parameters):
    return PartitionedFSIBaseSolver(model, project_parameters)

class PartitionedFSIBaseSolver(PythonSolver):
    def __init__(self, model, project_parameters):
        # TODO: Remove this as soon as the MPCs are implemented in MPI
        # This has to be done prior to the defaults check to avoid the structural solver to throw an error in MPI
        if not project_parameters["structure_solver_settings"].Has("multi_point_constraints_used"):
            project_parameters["structure_solver_settings"].AddEmptyValue("multi_point_constraints_used")
            project_parameters["structure_solver_settings"]["multi_point_constraints_used"].SetBool(False)

        # Call the base Python solver constructor
        # Note that default settings in GetDefaultSettings() are validated in here
        self._validate_settings_in_baseclass = True
        super(PartitionedFSIBaseSolver,self).__init__(model, project_parameters)

        # Auxiliar variables
        self.parallel_type = self.settings["parallel_type"].GetString()
        coupling_settings = self.settings["coupling_settings"]
        self.max_nl_it = coupling_settings["nl_max_it"].GetInt()
        self.nl_tol = coupling_settings["nl_tol"].GetDouble()
        self.solve_mesh_at_each_iteration = coupling_settings["solve_mesh_at_each_iteration"].GetBool()
        self.fluid_interface_submodelpart_name = coupling_settings["fluid_interfaces_list"][0].GetString()
        self.structure_interface_submodelpart_name = coupling_settings["structure_interfaces_list"][0].GetString()

        ## Construct the structure solver
        self.structure_solver = python_solvers_wrapper_structural.CreateSolverByParameters(self.model, self.settings["structure_solver_settings"], self.parallel_type)
        KratosMultiphysics.Logger.PrintInfo("::[PartitionedFSIBaseSolver]::", "Structure solver construction finished.")

        ## Construct the fluid solver
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_solver_settings"], self.parallel_type)
        KratosMultiphysics.Logger.PrintInfo("::[PartitionedFSIBaseSolver]::", "Fluid solver construction finished.")

        ## Check the ALE settings before the mesh solver construction
        mesh_solver_settings = self.settings["mesh_solver_settings"]
        fluid_solver_settings = self.settings["fluid_solver_settings"]

        # Check that the ALE and the fluid are the same model part
        fluid_model_part_name =  fluid_solver_settings["model_part_name"].GetString()
        if mesh_solver_settings.Has("model_part_name"):
            if not fluid_model_part_name == mesh_solver_settings["model_part_name"].GetString():
                err_msg =  'Fluid and mesh solver have to use the same MainModelPart ("model_part_name")!\n'
                raise Exception(err_msg)
        else:
            mesh_solver_settings.AddValue("model_part_name", fluid_solver_settings["model_part_name"])

        # Ensure that the fluid domain is not imported twice
        if mesh_solver_settings.Has("model_import_settings"):
            if mesh_solver_settings["model_import_settings"].Has("input_type"):
                if not mesh_solver_settings["model_import_settings"]["input_type"].GetString() == "use_input_model_part":
                    mesh_solver_settings["model_import_settings"]["input_type"].SetString("use_input_model_part")
            else:
                mesh_solver_settings["model_import_settings"].AddEmptyValue("input_type").SetString("use_input_model_part")
        else:
            mesh_solver_settings.AddEmptyValue("model_import_settings").AddEmptyValue("input_type").SetString("use_input_model_part")

        # Check that the ALE and the fluid have the sime time scheme
        fluid_time_scheme =  fluid_solver_settings["time_scheme"].GetString()
        if mesh_solver_settings.Has("mesh_velocity_calculation"):
            if mesh_solver_settings["mesh_velocity_calculation"].Has("time_scheme"):
                if not fluid_time_scheme == mesh_solver_settings["mesh_velocity_calculation"]["time_scheme"].GetString():
                    err_msg = 'Fluid and mesh solver require to use the same time scheme ("time_scheme") for consistency!\n'
                    raise Exception(err_msg)
            else:
                mesh_solver_settings["mesh_velocity_calculation"].AddValue("time_scheme", fluid_solver_settings["time_scheme"])
        else:
            mesh_solver_settings.AddEmptyValue("mesh_velocity_calculation")
            mesh_solver_settings["mesh_velocity_calculation"].AddValue("time_scheme", fluid_solver_settings["time_scheme"])

        # Check domain size
        fluid_domain_size = fluid_solver_settings["domain_size"].GetInt()
        if mesh_solver_settings.Has("domain_size"):
            if not fluid_domain_size == mesh_solver_settings["domain_size"].GetInt():
                raise Exception('Fluid and mesh solver have different "domain_size"!')
        else:
            mesh_solver_settings.AddValue("domain_size", fluid_solver_settings["domain_size"])

        # Ensure that the MESH_VELOCITY is computed
        if mesh_solver_settings.Has("calculate_mesh_velocity"):
            if not mesh_solver_settings["calculate_mesh_velocity"].GetBool():
                mesh_solver_settings.SetValue("calculate_mesh_velocity", True)
                KratosMultiphysics.Logger.PrintWarning("","Mesh velocity calculation was deactivated. Switching \"calculate_mesh_velocity\" on")
        else:
            mesh_solver_settings.AddEmptyValue("calculate_mesh_velocity").SetBool(True)

        ## Construct the ALE mesh solver
        self.mesh_solver = python_solvers_wrapper_mesh_motion.CreateSolverByParameters(self.model, self.settings["mesh_solver_settings"], self.parallel_type)
        KratosMultiphysics.Logger.PrintInfo("::[PartitionedFSIBaseSolver]::", "ALE mesh solver construction finished.")
        KratosMultiphysics.Logger.PrintInfo("::[PartitionedFSIBaseSolver]::", "Partitioned FSI base solver construction finished.")

    @classmethod
    def GetDefaultSettings(cls):
        """This function returns the default-settings used by this class
        """
        this_defaults = KratosMultiphysics.Parameters("""{
            "echo_level": 0,
            "parallel_type": "OpenMP",
            "solver_type": "partitioned",
            "structure_solver_settings": {
            },
            "fluid_solver_settings":{
            },
            "mesh_solver_settings":{
            },
            "coupling_settings":{
            }
        }""")
        this_defaults.AddMissingParameters(super(PartitionedFSIBaseSolver, cls).GetDefaultSettings())
        return this_defaults

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
        self.mesh_solver.ImportModelPart()

    def PrepareModelPart(self):
        # Get the minimum buffer size between the mesh, fluid and structure solvers
        self._GetAndSetMinimumBufferSize()
        # Fluid and structure solvers PrepareModelPart() call
        self.structure_solver.PrepareModelPart()
        self.fluid_solver.PrepareModelPart()
        self.mesh_solver.PrepareModelPart()

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
        self.structure_solver.InitializeSolutionStep()
        self.fluid_solver.InitializeSolutionStep()
        self.mesh_solver.InitializeSolutionStep()
        self._GetConvergenceAccelerator().InitializeSolutionStep()

    def Predict(self):
        # Perform fluid and structure solvers predictions
        self.structure_solver.Predict()
        self.fluid_solver.Predict()
        self.mesh_solver.Predict()

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
        KratosMultiphysics.Logger.PrintInfo("","\tMesh prediction computed.")

        ## Non-Linear interface coupling iteration ##
        nl_it = 0
        is_converged = False
        while not is_converged:
            nl_it += 1
            KratosMultiphysics.Logger.PrintInfo("","\tFSI non-linear iteration = ", nl_it)

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
                is_converged = True
                self._GetConvergenceAccelerator().FinalizeNonLinearIteration()
                KratosMultiphysics.Logger.PrintInfo("","\tNon-linear iteration convergence achieved")
                KratosMultiphysics.Logger.PrintInfo("","\tTotal non-linear iterations: ", nl_it, " |res|/sqrt(nDOFS) = ", nl_res_norm/sqrt(interface_dofs))
                break
            else:
                # If convergence is not achieved, perform the correction of the prediction
                KratosMultiphysics.Logger.PrintInfo("","\tResidual computation finished. |res|/sqrt(nDOFS) = ", nl_res_norm/sqrt(interface_dofs))
                self._GetConvergenceAccelerator().UpdateSolution(dis_residual, self.iteration_value)
                self._GetConvergenceAccelerator().FinalizeNonLinearIteration()

                if (nl_it == self.max_nl_it):
                    KratosMultiphysics.Logger.PrintWarning("","\tFSI NON-LINEAR ITERATION CONVERGENCE NOT ACHIEVED")
                    break

        ## Compute the mesh residual as final testing (it is expected to be 0)
        mesh_res_norm = self.partitioned_fsi_utilities.ComputeInterfaceResidualNorm(
            self._GetFluidInterfaceSubmodelPart(),
            KratosMultiphysics.VELOCITY,
            KratosMultiphysics.MESH_VELOCITY,
            KratosMultiphysics.FSI_INTERFACE_MESH_RESIDUAL,
            "nodal")
        KratosMultiphysics.Logger.PrintInfo("","\tNL residual norm: ", nl_res_norm)
        KratosMultiphysics.Logger.PrintInfo("","\tMesh residual norm: ", mesh_res_norm)

        return is_converged

    def FinalizeSolutionStep(self):
        self.structure_solver.FinalizeSolutionStep()
        self.fluid_solver.FinalizeSolutionStep()
        self.mesh_solver.FinalizeSolutionStep()
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

    # This method returns the convergence accelerator.
    # If it is not created yet, it calls the _CreateConvergenceAccelerator first
    def _GetConvergenceAccelerator(self):
        if not hasattr(self, '_convergence_accelerator'):
            self._convergence_accelerator = self._CreateConvergenceAccelerator()
        return self._convergence_accelerator

    # This method constructs the convergence accelerator coupling utility
    def _CreateConvergenceAccelerator(self):
        convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(self.settings["coupling_settings"]["coupling_strategy_settings"])
        KratosMultiphysics.Logger.PrintInfo("::[PartitionedFSIBaseSolver]::", "Coupling strategy construction finished.")
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

    def _GetPartitionedFSIUtilities(self):
        if (self.domain_size == 2):
            return KratosFSI.PartitionedFSIUtilitiesArray2D()
        else:
            return KratosFSI.PartitionedFSIUtilitiesArray3D()

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

        max_cond_id = self.structure_solver.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(max_cond_id)

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

    def _ComputeMeshPredictionDoubleFaced(self):
            # Get the previous step fluid interface nodal fluxes from both positive and negative faces
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

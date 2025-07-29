from math import sqrt   # Import the square root from python library

# Import utilities
from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid            # Import the fluid Python solvers wrapper
from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural       # Import the structure Python solvers wrapper
from KratosMultiphysics.MeshMovingApplication import python_solvers_wrapper_mesh_motion      # Import the mesh motion Python solvers wrapper
from KratosMultiphysics.FSIApplication import fsi_coupling_interface                            # Import the FSI coupling interface utility
from KratosMultiphysics.FSIApplication import convergence_accelerator_factory         # Import the FSI convergence accelerator factory

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.MappingApplication as KratosMapping
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
        # Note that default settings in GetDefaultParameters() are validated in here
        super().__init__(model, project_parameters)

        # This method encapsulates all the operations required to be done in the constructor
        # It is intended to be overwritten in the derived classes
        self._AuxiliaryInitOperations()

        KratosMultiphysics.Logger.PrintInfo("PartitionedFSIBaseSolver", "Partitioned FSI base solver construction finished.")

    @classmethod
    def GetDefaultParameters(cls):

        # Note that only the coupling settings are validated
        # The subdomain solver settings will be validated while instantiating these
        this_defaults = KratosMultiphysics.Parameters("""{
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
                "coupling_strategy_settings": {
                    "solver_type": "MVQN",
                    "w_0": 0.825
                },
                "mapper_settings": [{
                    "fluid_interface_submodelpart_name": "",
                    "mapper_face": "unique",
                    "structure_interface_submodelpart_name": ""
                }],
                "nl_max_it": 25,
                "nl_tol": 1e-8,
                "solve_mesh_at_each_iteration": true,
                "fluid_interfaces_list": [],
                "structure_interfaces_list": []
            }
        }""")

        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def ValidateSettings(self):
        default_settings = self.GetDefaultParameters()

        ## Base class settings validation
        super().ValidateSettings()

        ## Validate coupling settings
        self.settings["coupling_settings"].ValidateAndAssignDefaults(default_settings["coupling_settings"])

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
        # Mesh solver variables addition
        self.mesh_solver.AddVariables()

        ## Fluid traction variables addition
        # These are required for the reaction variable redistribution in the fluid model part skin
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_MAPPED_VECTOR_VARIABLE)
        self.fluid_solver.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_MAPPED_VECTOR_VARIABLE)

    def ImportModelPart(self):
        # Fluid and structure solvers ImportModelPart() call
        self.fluid_solver.ImportModelPart()
        self.structure_solver.ImportModelPart()
        self.mesh_solver.ImportModelPart()

    def PrepareModelPart(self):
        # Fluid and structure solvers PrepareModelPart() call
        self.structure_solver.PrepareModelPart()
        self.fluid_solver.PrepareModelPart()
        self.mesh_solver.PrepareModelPart()

        # FSI interface coupling interfaces initialization
        # The getter methods are to construct the FSI coupling interfaces in here
        self._GetFSICouplingInterfaceFluidPositive().GetInterfaceModelPart()
        if self.double_faced_structure:
            self._GetFSICouplingInterfaceFluidNegative().GetInterfaceModelPart()
        self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart()

    def AddDofs(self):
        # Add DOFs structure
        self.structure_solver.AddDofs()

        # Add DOFs fluid
        self.fluid_solver.AddDofs()
        self.mesh_solver.AddDofs()

    def Initialize(self):
        # Create the required mappers
        self._GetStructureToFluidPositiveInterfaceMapper()
        if self.double_faced_structure:
            self._GetStructureToFluidNegativeInterfaceMapper()

        # Coupling utility initialization
        # The _GetConvergenceAccelerator is supposed to construct the convergence accelerator in here
        self._GetConvergenceAccelerator().Initialize()

        # Python structure solver initialization
        self.structure_solver.Initialize()

        # Ensure that the fluid reaction fluxes are computed before initializing the fluid strategy
        if self.fluid_solver.settings["compute_reactions"].GetBool() == False:
            self.fluid_solver.settings["compute_reactions"].SetBool(True)

        # Python fluid solver initialization
        self.fluid_solver.Initialize()

        # Python mesh solver initialization
        self.mesh_solver.Initialize()

        # Perform all the operations required to set up the coupling interfaces
        self._InitializeCouplingInterfaces()

    def AdvanceInTime(self, current_time):
        # Subdomains time advance
        fluid_new_time = self.fluid_solver.AdvanceInTime(current_time)
        structure_new_time = self.structure_solver.AdvanceInTime(current_time)

        if abs(fluid_new_time - structure_new_time) > 1e-12:
            err_msg =  'Fluid new time is: ' + str(fluid_new_time) + '\n'
            err_msg += 'Structure new time is: ' + str(structure_new_time) + '\n'
            err_msg += 'No substepping has been implemented yet. Fluid and structure time must coincide.'
            raise Exception(err_msg)

        # Update the auxiliary coupling interface model parts database in case these are to be output
        # Note that there are also situations (e.g. embedded) in which this is mandatory for the proper performance of the solver
        self._AdvanceInTimeCouplingInterfaces(fluid_new_time)

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
        ## Initialize residual
        dis_residual_norm = self._ComputeInitialResidual()

        ## FSI nonlinear iteration loop
        nl_it = 1
        is_converged = False
        while (nl_it < self.max_nl_it):
            # Check convergence
            if self._CheckFSIConvergence(nl_it, dis_residual_norm):
                KratosMultiphysics.Logger.PrintInfo('PartitionedFSIBaseSolver', 'FSI non-linear iteration convergence achieved in {0} iterations.'.format(nl_it))
                is_converged = True
                break
            else:
                nl_it += 1

            # Perform the displacement convergence accelerator update
            self._GetConvergenceAccelerator().InitializeNonLinearIteration()
            if not self._GetConvergenceAccelerator().IsBlockNewton():
                self._GetFSICouplingInterfaceStructure().Update()
            else:
                self._GetFSICouplingInterfaceStructure().UpdateDisplacement()

            # Update the structure interface position
            self._GetFSICouplingInterfaceStructure().UpdatePosition(KratosMultiphysics.RELAXED_DISPLACEMENT)

            # Map the relaxed displacement from the structure coupling interface to the fluid one
            self._MapStructureInterfaceDisplacement()

            # Solve the mesh and fluid problem
            self._SolveFluid()

            # Calculate the fluid interface traction from the interface reaction
            #TODO: THINK ABOUT COMPUTING THE SUM OF NEGATIVE AND POSITIVE DISTRIBUTED LOAD AND TRANSFER THIS TO THE STRUCTURE (INSTEAD OF PASSING POSITIVE AND NEGATIVE)
            self._CalculateFluidInterfaceTraction()

            # Transfer the fluid traction to the structure interface
            self._MapFluidInterfaceTraction()

            # Compute the current iteration traction
            if not self._GetConvergenceAccelerator().IsBlockNewton():
                # Directly send the map load from the structure FSI coupling interface to the parent one
                self._GetFSICouplingInterfaceStructure().TransferValuesToFatherModelPart(self._GetTractionVariable())
            else:
                # Perform the traction convergence accelerator update
                self._GetFSICouplingInterfaceStructure().UpdateTraction()
                KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                    KratosMultiphysics.RELAXED_TRACTION,
                    self._GetTractionVariable(),
                    self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                    self._GetFSICouplingInterfaceStructure().GetFatherModelPart(),
                    0)

            self._GetConvergenceAccelerator().FinalizeNonLinearIteration()

            # Solve the structure problem
            self._SolveStructure()

            # Compute the residual vector
            dis_residual_norm = self._GetFSICouplingInterfaceStructure().ComputeResidualVector()

            # Check if maximum iterations are reached
            if nl_it == self.max_nl_it:
                KratosMultiphysics.Logger.PrintWarning('PartitionedFSIBaseSolver', 'FSI non-linear maximum iterations {0} reached'.format(self.max_nl_it))

        return is_converged

    def FinalizeSolutionStep(self):
        self.structure_solver.FinalizeSolutionStep()
        self.fluid_solver.FinalizeSolutionStep()
        self.mesh_solver.FinalizeSolutionStep()
        self._GetConvergenceAccelerator().FinalizeSolutionStep()

    def Finalize(self):
        self.structure_solver.Finalize()
        self.fluid_solver.Finalize()
        self.mesh_solver.Finalize()
        self._GetConvergenceAccelerator().Finalize()

    def SetEchoLevel(self, structure_echo_level, fluid_echo_level):
        self.structure_solver.SetEchoLevel(self, structure_echo_level)
        self.fluid_solver.SetEchoLevel(self, fluid_echo_level)

    def Clear(self):
        self.fluid_solver.Clear()
        self.structure_solver.Clear()

    def Check(self):
        self.fluid_solver.Check()
        self.structure_solver.Check()

    #######################################################################
    ##############          PRIVATE METHODS SECTION          ##############
    #######################################################################

    def _AuxiliaryInitOperations(self):
        # Auxiliar variables
        self.parallel_type = self.settings["parallel_type"].GetString()
        coupling_settings = self.settings["coupling_settings"]
        self.max_nl_it = coupling_settings["nl_max_it"].GetInt()
        self.nl_tol = coupling_settings["nl_tol"].GetDouble()
        self.solve_mesh_at_each_iteration = coupling_settings["solve_mesh_at_each_iteration"].GetBool()

        ## Save the FSI interfaces information in a dictionary
        mappers_settings = self.settings["coupling_settings"]["mapper_settings"]
        self.interfaces_dict = {}
        self.interfaces_dict['structure'] = coupling_settings["structure_interfaces_list"][0].GetString()
        if mappers_settings.size() == 1:
            self.double_faced_structure = False
            self.interfaces_dict['fluid_positive'] = mappers_settings[0]["fluid_interface_submodelpart_name"].GetString()
            self.interfaces_dict['fluid_negative'] = None
        elif mappers_settings.size() == 2:
            self.double_faced_structure = True
            for mapper_id in range(2):
                if (mappers_settings[mapper_id]["mapper_face"].GetString() == "Positive"):
                    self.interfaces_dict['fluid_positive'] = mappers_settings[mapper_id]["fluid_interface_submodelpart_name"].GetString()
                elif (mappers_settings[mapper_id]["mapper_face"].GetString() == "Negative"):
                    self.interfaces_dict['fluid_negative'] = mappers_settings[mapper_id]["fluid_interface_submodelpart_name"].GetString()
                else:
                    err_msg = "The mapper face is not \'Positve\' nor \'Negative\'."
                    raise Exception(err_msg)
        else:
            err_msg = "Case with more than 2 mappers has not been implemented yet.\n"
            err_msg += "Please, in case you are using single faced immersed bodies, set the skin entities in a unique submodelpart.\n"
            err_msg += "In case you are considering double faced immersed bodies (shells or membranes), set all the positive faces in a unique submodelpart and all the negative ones in another submodelpart."
            raise Exception(err_msg)

        ## Construct the structure solver
        self.structure_solver = python_solvers_wrapper_structural.CreateSolverByParameters(self.model, self.settings["structure_solver_settings"], self.parallel_type)
        KratosMultiphysics.Logger.PrintInfo("PartitionedFSIBaseSolver", "Structure solver construction finished.")

        ## Construct the fluid solver
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_solver_settings"], self.parallel_type)
        KratosMultiphysics.Logger.PrintInfo("PartitionedFSIBaseSolver", "Fluid solver construction finished.")

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
        if fluid_solver_settings.Has("time_scheme"):
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
        KratosMultiphysics.Logger.PrintInfo("PartitionedFSIBaseSolver", "ALE mesh solver construction finished.")

    def _AdvanceInTimeCouplingInterfaces(self, new_time):
        self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart().CloneTimeStep(new_time)
        self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart().ProcessInfo[KratosMultiphysics.STEP] = self._GetFSICouplingInterfaceStructure().GetFatherModelPart().ProcessInfo[KratosMultiphysics.STEP]
        self._GetFSICouplingInterfaceFluidPositive().GetInterfaceModelPart().CloneTimeStep(new_time)
        self._GetFSICouplingInterfaceFluidPositive().GetInterfaceModelPart().ProcessInfo[KratosMultiphysics.STEP] = self._GetFSICouplingInterfaceFluidPositive().GetFatherModelPart().ProcessInfo[KratosMultiphysics.STEP]
        if self.double_faced_structure:
            self._GetFSICouplingInterfaceFluidNegative().GetInterfaceModelPart().CloneTimeStep(new_time)
            self._GetFSICouplingInterfaceFluidNegative().GetInterfaceModelPart().ProcessInfo[KratosMultiphysics.STEP] = self._GetFSICouplingInterfaceFluidNegative().GetFatherModelPart().ProcessInfo[KratosMultiphysics.STEP]

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

    #TODO: SHOULDN'T BE NEEDED AS WE ALWAYS USE THE POSITIVE ONE
    # def _GetFluidInterfaceSubmodelPart(self):
    #     # Returns the fluid interface submodelpart that will be used in the residual minimization
    #     return self.fluid_solver.main_model_part.GetSubModelPart(self.fluid_interface_submodelpart_name)

    def _InitializeCouplingInterfaces(self):
        # Prepare the Dirichlet fluid interfaces
        self.__SetFluidInterfaceFixity(self._GetFluidPositiveInterfaceSubmodelPart())
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, self._GetFluidPositiveInterfaceSubmodelPart().Nodes)
        if self.double_faced_structure:
            self.__SetFluidInterfaceFixity(self._GetFluidNegativeInterfaceSubmodelPart())
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, self._GetFluidNegativeInterfaceSubmodelPart().Nodes)

        # Prepare the Neumann structure interface
        self._SetStructureNeumannCondition()
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, self._GetStructureInterfaceSubmodelPart().Nodes)

    def __SetFluidInterfaceFixity(self, fluid_interface_submodelpart):
        # Fix the VELOCITY, MESH_DISPLACEMENT and MESH_VELOCITY variables in all the fluid interface submodelparts
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VELOCITY_X, True, fluid_interface_submodelpart.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VELOCITY_Y, True, fluid_interface_submodelpart.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.MESH_DISPLACEMENT_X, True, fluid_interface_submodelpart.Nodes)
        KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.MESH_DISPLACEMENT_Y, True, fluid_interface_submodelpart.Nodes)
        if self._GetDomainSize() == 3:
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.VELOCITY_Z, True, fluid_interface_submodelpart.Nodes)
            KratosMultiphysics.VariableUtils().ApplyFixity(KratosMultiphysics.MESH_DISPLACEMENT_Z, True, fluid_interface_submodelpart.Nodes)

    def _GetFluidPositiveInterfaceSubmodelPart(self):
        # Returns the fluid positive interface submodelpart
        positive_interface_name = self.interfaces_dict['fluid_positive']
        return self.model.GetModelPart(positive_interface_name)

    def _GetFluidNegativeInterfaceSubmodelPart(self):
        negative_interface_name = self.interfaces_dict['fluid_negative']
        return self.model.GetModelPart(negative_interface_name)

    def _GetStructureInterfaceSubmodelPart(self):
        # Returns the structure interface submodelpart that will be used in the residual minimization
        return self.model.GetModelPart(self.__GetStructureInterfaceModelPartName())

    def __GetStructureInterfaceModelPartName(self):
        if not hasattr(self, '_structure_interface_submodelpart_name'):
            str_int_list = self.settings["coupling_settings"]["structure_interfaces_list"]
            if (str_int_list.size() != 1):
                raise Exception("FSI structure skin must be contained in a unique model part")
            self._structure_interface_submodelpart_name = str_int_list[0].GetString()
        return self._structure_interface_submodelpart_name

    def _GetDomainSize(self):
        if not hasattr(self, '_domain_size'):
            fluid_domain_size = self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            structure_domain_size = self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
            if fluid_domain_size !=structure_domain_size:
                raise("ERROR: Solid domain size and fluid domain size are not equal!")
            self._domain_size = fluid_domain_size
        return self._domain_size

    def  _ComputeDeltaTime(self):
        fluid_time_step = self.fluid_solver._ComputeDeltaTime()
        structure_time_step = self.structure_solver.ComputeDeltaTime()

        if abs(fluid_time_step - structure_time_step) > 1e-12:
            err_msg =  'Fluid time step is: ' + str(fluid_time_step) + '\n'
            err_msg += 'Structure time step is: ' + str(structure_time_step) + '\n'
            err_msg += 'No substepping has been implemented yet. Fluid and structure time step must coincide.'
            raise Exception(err_msg)

        return fluid_time_step

    #TODO: THIS MUST BE DONE IN C++ (FSIPartitionedUtils)
    #TODO: THINK ABOUT ASSUMING THESE TO BE ALREADY IN THE MPDA (AS THE EMBEDDED DOES) OR TO ALWAYS CREATE THEM (ALSO IN THE EMBEDDED)
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
            interface_submodelpart_i = self.model.GetModelPart(interface_submodelpart_name)

            # Get the number of conditions to be set in each processor
            local_nodes_number_accumulated = -1
            local_nodes_number = len(interface_submodelpart_i.GetCommunicator().LocalMesh().Nodes)
            local_nodes_number_accumulated = interface_submodelpart_i.GetCommunicator().GetDataCommunicator().ScanSum(local_nodes_number)

            # Create the point load condition
            aux_count = max_cond_id + local_nodes_number_accumulated
            if self._GetDomainSize() == 2:
                for node in interface_submodelpart_i.GetCommunicator().LocalMesh().Nodes:
                    aux_count+=1
                    structure_computational_submodelpart.CreateNewCondition("PointLoadCondition2D1N",
                                                                            int(aux_count),
                                                                            [node.Id],
                                                                            self.structure_solver.main_model_part.Properties[0])
            elif self._GetDomainSize() == 3:
                for node in interface_submodelpart_i.GetCommunicator().LocalMesh().Nodes:
                    aux_count+=1
                    structure_computational_submodelpart.CreateNewCondition("PointLoadCondition3D1N",
                                                                            int(aux_count),
                                                                            [node.Id],
                                                                            self.structure_solver.main_model_part.Properties[0])

    def _GetStructureToFluidPositiveInterfaceMapper(self):
        if not hasattr(self, '_structure_to_fluid_positive_interface_mapper'):
            self._structure_to_fluid_positive_interface_mapper = self._CreateStructureToFluidInterfaceMapper(
                self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                self._GetFSICouplingInterfaceFluidPositive().GetInterfaceModelPart())
        return self._structure_to_fluid_positive_interface_mapper

    def _GetStructureToFluidNegativeInterfaceMapper(self):
        if not hasattr(self, '_structure_to_fluid_negative_interface_mapper'):
            self._structure_to_fluid_negative_interface_mapper = self._CreateStructureToFluidInterfaceMapper(
                self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                self._GetFSICouplingInterfaceFluidNegative().GetInterfaceModelPart())
        return self._structure_to_fluid_negative_interface_mapper

    @classmethod
    def _CreateStructureToFluidInterfaceMapper(self, structure_interface, fluid_interface):
        mapper_params = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "echo_level" : 0
        }""")
        structure_to_fluid_interface_mapper = KratosMapping.MapperFactory.CreateMapper(
            structure_interface,
            fluid_interface,
            mapper_params)

        return structure_to_fluid_interface_mapper

    def _GetFSICouplingInterfaceStructure(self):
        if not hasattr(self, '_fsi_coupling_interface_structure'):
            self._fsi_coupling_interface_structure = self._CreateFSICouplingInterfaceStructure()
        return self._fsi_coupling_interface_structure

    def _CreateFSICouplingInterfaceStructure(self):
        # Set auxiliary settings
        aux_settings = KratosMultiphysics.Parameters(
        """{
            "model_part_name": "FSICouplingInterfaceStructure",
            "parent_model_part_name": "",
            "input_variable_list": ["SURFACE_LOAD"],
            "output_variable_list": ["DISPLACEMENT"]
        }""")
        aux_settings["parent_model_part_name"].SetString(self.interfaces_dict['structure'])

        # Get the convergence accelerator instance (note that it is created in this call)
        convergence_accelerator = self._GetConvergenceAccelerator()

        # Construct the FSI coupling interface
        fsi_coupling_interface_structure = fsi_coupling_interface.FSICouplingInterface(
            self.model,
            aux_settings,
            convergence_accelerator)

        KratosMultiphysics.Logger.PrintInfo('PartitionedFSIBaseSolver', 'Structure FSI coupling interface created')

        return fsi_coupling_interface_structure

    @classmethod
    def _GetFSICouplingInterfaceFluid(self):
        err_msg = 'Body-fitted partitioned FSI solver does not implement the \'_GetFSICouplingInterfaceFluid\' method.\n'
        err_msg += 'Use \'_GetFSICouplingInterfaceFluidPositive\' as well as \'_GetFSICouplingInterfaceFluidNegative\' if the structure is thin-walled.'
        raise Exception(err_msg)

    def _GetFSICouplingInterfaceFluidPositive(self):
        if not hasattr(self, '_fsi_coupling_interface_fluid_positive'):
            fluid_interface_side = 'fluid_positive'
            self._fsi_coupling_interface_fluid_positive = self.__CreateFSICouplingInterfaceFluid(fluid_interface_side)
        return self._fsi_coupling_interface_fluid_positive

    def _GetFSICouplingInterfaceFluidNegative(self):
        if not hasattr(self, '_fsi_coupling_interface_fluid_negative'):
            fluid_interface_side = 'fluid_negative'
            self._fsi_coupling_interface_fluid_negative = self.__CreateFSICouplingInterfaceFluid(fluid_interface_side)
        return self._fsi_coupling_interface_fluid_negative

    def __CreateFSICouplingInterfaceFluid(self, fluid_interface_side):
        # Set auxiliary settings
        aux_settings = KratosMultiphysics.Parameters(
        """{
            "model_part_name": "",
            "parent_model_part_name": "",
            "input_variable_list": ["MESH_DISPLACEMENT"],
            "output_variable_list": []
        }""")
        aux_settings["model_part_name"].SetString("FSICouplingInterface{}".format(fluid_interface_side.capitalize()))
        aux_settings["parent_model_part_name"].SetString(self.interfaces_dict[fluid_interface_side])
        aux_settings["output_variable_list"].Append("{}_MAPPED_VECTOR_VARIABLE".format("POSITIVE" if fluid_interface_side == 'fluid_positive' else "NEGATIVE"))

        # Construct the FSI coupling interface
        # Note that no convergence accelerator is created in the fluid side
        fsi_coupling_interface_fluid = fsi_coupling_interface.FSICouplingInterface(
            self.model,
            aux_settings)

        KratosMultiphysics.Logger.PrintInfo('PartitionedFSIBaseSolver', 'Fluid {} FSI coupling interface created'.format(fluid_interface_side))

        return fsi_coupling_interface_fluid

    # This method returns the convergence accelerator.
    # If it is not created yet, it calls the _CreateConvergenceAccelerator first
    def _GetConvergenceAccelerator(self):
        if not hasattr(self, '_convergence_accelerator'):
            self._convergence_accelerator = self._CreateConvergenceAccelerator()
        return self._convergence_accelerator

    # This method constructs the convergence accelerator coupling utility
    def _CreateConvergenceAccelerator(self):
        convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(self.settings["coupling_settings"]["coupling_strategy_settings"])
        KratosMultiphysics.Logger.PrintInfo('PartitionedFSIBaseSolver', 'Coupling strategy construction finished')
        return convergence_accelerator

    def _MapStructureInterfaceDisplacement(self):
        # Map the RELAXED_DISPLACEMENT from the structure FSI coupling interface to fluid FSI coupling interface
        self._GetStructureToFluidPositiveInterfaceMapper().Map(
            KratosMultiphysics.RELAXED_DISPLACEMENT,
            KratosMultiphysics.MESH_DISPLACEMENT)

        # In case the structure is double faced, we take advantage of the equal positive and negative sides meshes
        if self.double_faced_structure:
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.MESH_DISPLACEMENT,
                KratosMultiphysics.MESH_DISPLACEMENT,
                self._GetFSICouplingInterfaceFluidPositive.GetInterfaceModelPart(),
                self._GetFSICouplingInterfaceFluidNegative.GetInterfaceModelPart())

        # Update the fluid interface position
        self._GetFSICouplingInterfaceFluidPositive().UpdatePosition(KratosMultiphysics.MESH_DISPLACEMENT)
        if self.double_faced_structure:
            self._GetFSICouplingInterfaceFluidNegative().UpdatePosition(KratosMultiphysics.MESH_DISPLACEMENT)

    def _SolveFluid(self):
        # Transfer values from fluid FSI coupling interface to fluid skin
        self._GetFSICouplingInterfaceFluidPositive().TransferValuesToFatherModelPart(KratosMultiphysics.MESH_DISPLACEMENT)
        if self.double_faced_structure:
            self._GetFSICouplingInterfaceFluidNegative().TransferValuesToFatherModelPart(KratosMultiphysics.MESH_DISPLACEMENT)

        # Solve the mesh problem (or moves the interface nodes)
        if self.solve_mesh_at_each_iteration:
            self.mesh_solver.SolveSolutionStep()
        else:
            self.mesh_solver.MoveMesh()

        # Impose the structure MESH_VELOCITY in the fluid interface VELOCITY (Lagrangian interface)
        KratosMultiphysics.VariableUtils().CopyVariable(
                KratosMultiphysics.MESH_VELOCITY,
                KratosMultiphysics.VELOCITY,
                self._GetFluidPositiveInterfaceSubmodelPart().GetCommunicator().LocalMesh().Nodes)
        self._GetFluidPositiveInterfaceSubmodelPart().GetCommunicator().SynchronizeVariable(KratosMultiphysics.VELOCITY)
        if self.double_faced_structure:
            KratosMultiphysics.VariableUtils().CopyVariable(
                    KratosMultiphysics.MESH_VELOCITY,
                    KratosMultiphysics.VELOCITY,
                    self._GetFluidNegativeInterfaceSubmodelPart().GetCommunicator().LocalMesh().Nodes)
            self._GetFluidNegativeInterfaceSubmodelPart().GetCommunicator().SynchronizeVariable(KratosMultiphysics.VELOCITY)

        # Solve fluid problem
        self.fluid_solver.SolveSolutionStep()

    def _CalculateFluidInterfaceTraction(self):
        # Distribute the REACTION point load
        distribution_tolerance = 1.0e-12
        distribution_max_iterations = 500

        KratosMultiphysics.VariableRedistributionUtility.DistributePointValues(
            self._GetFluidPositiveInterfaceSubmodelPart(),
            self._GetFluidPositiveInterfaceSubmodelPart().Conditions,
            KratosMultiphysics.REACTION,
            KratosMultiphysics.POSITIVE_MAPPED_VECTOR_VARIABLE,
            distribution_tolerance,
            distribution_max_iterations)
        if self.double_faced_structure:
            KratosMultiphysics.VariableRedistributionUtility.DistributePointValues(
                self._GetFluidNegativeInterfaceSubmodelPart(),
                self._GetFluidNegativeInterfaceSubmodelPart().Conditions,
                KratosMultiphysics.REACTION,
                KratosMultiphysics.NEGATIVE_MAPPED_VECTOR_VARIABLE,
                distribution_tolerance,
                distribution_max_iterations)

        # Transfer traction values from fluid interface to fluid FSI coupling interface
        # NOTE: POSITIVE_MAPPED_VECTOR_VARIABLE is the positive side traction
        # NOTE: NEGATIVE_MAPPED_VECTOR_VARIABLE is the negative side traction
        self._GetFSICouplingInterfaceFluidPositive().GetValuesFromFatherModelPart(KratosMultiphysics.POSITIVE_MAPPED_VECTOR_VARIABLE)
        if self.double_faced_structure:
            self._GetFSICouplingInterfaceFluidNegative().GetValuesFromFatherModelPart(KratosMultiphysics.NEGATIVE_MAPPED_VECTOR_VARIABLE)

    def _MapFluidInterfaceTraction(self):
        # Map the distributed REACTION from the fluid FSI coupling interface to structure FSI coupling interface
        self._GetStructureToFluidPositiveInterfaceMapper().InverseMap(
            self._GetTractionVariable(),
            KratosMultiphysics.POSITIVE_MAPPED_VECTOR_VARIABLE,
            KratosMultiphysics.Mapper.SWAP_SIGN)
        if self.double_faced_structure:
            self._GetStructureToFluidNegativeInterfaceMapper().InverseMap(
                self._GetTractionVariable(),
                KratosMultiphysics.NEGATIVE_MAPPED_VECTOR_VARIABLE,
                KratosMultiphysics.Mapper.ADD_VALUES | KratosMultiphysics.Mapper.SWAP_SIGN)

        # Send the mapped load from the structure FSI coupling interface to the parent one
        self._GetFSICouplingInterfaceStructure().TransferValuesToFatherModelPart(self._GetTractionVariable())

    def _SolveStructure(self):
        #TODO: Once we use the distributed (relaxed) traction in the body fitted we can remove this as well as the _SolveStructure in the embedded solver
        # Convert distributed traction to point values
        KratosMultiphysics.VariableRedistributionUtility.ConvertDistributedValuesToPoint(
            self._GetStructureInterfaceSubmodelPart(),
            self._GetStructureInterfaceSubmodelPart().Conditions,
            self._GetTractionVariable(),
            KratosStructural.POINT_LOAD)

        # Solve the structure problem
        self.structure_solver.SolveSolutionStep()

    def _ComputeInitialResidual(self):
        # Save as RELAXED_DISPLACEMENT the DISPLACEMENT coming from the structure Predict()
        # Note that this would be required in order to set the first observation matrices
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
            KratosMultiphysics.DISPLACEMENT,
            KratosMultiphysics.RELAXED_DISPLACEMENT,
            self._GetStructureInterfaceSubmodelPart(),
            self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
            0)

        # Update the structure interface position with the DISPLACEMENT values from the predict
        self._GetFSICouplingInterfaceStructure().UpdatePosition(KratosMultiphysics.RELAXED_DISPLACEMENT)

        # Map the relaxed displacement from the structure coupling interface to the fluid one
        self._MapStructureInterfaceDisplacement()

        # Solve the mesh and fluid problem
        self._SolveFluid()

        # Calculate the fluid interface traction from the interface reaction
        #TODO: THINK ABOUT COMPUTING THE SUM OF NEGATIVE AND POSITIVE DISTRIBUTED LOAD AND TRANSFER THIS TO THE STRUCTURE (INSTEAD OF PASSING POSITIVE AND NEGATIVE)
        self._CalculateFluidInterfaceTraction()

        # Transfer the fluid traction to the structure interface
        self._MapFluidInterfaceTraction()

        # Save as RELAXED_TRATION the TRACTION coming from the fluid
        # Note that this would be required in order to set the first observation matrices
        if self._GetConvergenceAccelerator().IsBlockNewton():
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                self._GetTractionVariable(),
                KratosMultiphysics.RELAXED_TRACTION,
                self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                0)

        # Directly send the map load from the structure FSI coupling interface to the parent one
        self._GetFSICouplingInterfaceStructure().TransferValuesToFatherModelPart(self._GetTractionVariable())

        # Solve the structure problem
        self._SolveStructure()

        # Compute the residual vector
        dis_residual_norm = self._GetFSICouplingInterfaceStructure().ComputeResidualVector()

        return dis_residual_norm

    def _CheckFSIConvergence(self, nl_it, residual_norm):
        interface_dofs = self._GetPartitionedFSIUtilities().GetInterfaceResidualSize(self._GetStructureInterfaceSubmodelPart())
        normalised_residual = residual_norm/sqrt(interface_dofs)
        KratosMultiphysics.Logger.PrintInfo('PartitionedFSIBaseSolver', 'FSI non-linear iteration = {0} |res|/sqrt(nDOFS) = {1}'.format(nl_it, normalised_residual))
        return normalised_residual < self.nl_tol

    def _GetPartitionedFSIUtilities(self):
        if not hasattr(self, '_partitioned_fsi_utilities'):
            self._partitioned_fsi_utilities = self._CreatePartitionedFSIUtilities()
        return self._partitioned_fsi_utilities

    def _CreatePartitionedFSIUtilities(self):
        if self._GetDomainSize() == 2:
            return KratosFSI.PartitionedFSIUtilitiesArray2D()
        elif self._GetDomainSize() == 3:
            return KratosFSI.PartitionedFSIUtilitiesArray3D()
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self._GetDomainSize()))

    @classmethod
    def _GetTractionVariable(self):
        #NOTE: No distinction is made between LINE_LOAD and SURFACE_LOAD as this is only used for the redistribution
        #NOTE: Hence, this has only an auxiliary use as the traction will be converted to point load prior to the structure solve
        #TODO: Use LINE_LOAD in 2D and SURFACE_LOAD in 3D when using distributed loads conditions
        return KratosStructural.SURFACE_LOAD

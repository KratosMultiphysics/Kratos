from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from math import sqrt   # Import the square root from python library

# Import utilities
from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid            # Import the fluid Python solvers wrapper
from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural # Import the structure Python solvers wrapper
from KratosMultiphysics.FSIApplication import fsi_coupling_interface                            # Import the FSI coupling interface utility
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
        self.level_set_type = self.settings["fluid_solver_settings"]["formulation"]["level_set_type"].GetString()

        # First call to create the embedded intersections model part
        self.__GetEmbedddedSkinUtilityModelPart()

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Fluid solver construction finished')
        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Partitioned embedded FSI base solver construction finished')

    def GetMinimumBufferSize(self):
        buffer_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_structure = self.structure_solver.GetMinimumBufferSize()
        return max(buffer_structure,buffer_fluid)

    def AddVariables(self):
        # Fluid and structure solvers variables addition
        self.fluid_solver.AddVariables()
        self.structure_solver.AddVariables()

    def ImportModelPart(self):
        # Fluid and structure solvers ImportModelPart() call
        self.fluid_solver.ImportModelPart()
        self.structure_solver.ImportModelPart()

    def PrepareModelPart(self):
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
        self.domain_size = self.__GetDomainSize()

        # Coupling utility initialization
        # The __GetConvergenceAccelerator is supposed to construct the convergence accelerator in here
        self.__GetConvergenceAccelerator().Initialize()

        # FSI interface coupling interfaces initialization
        # The __GetFSICouplingInterfaceStructure is supposed to construct the FSI coupling structure interface in here
        self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart()
        self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart()

        # Python structure solver initialization
        self.structure_solver.Initialize()

        # Python fluid solver initialization
        self.fluid_solver.Initialize()

        # Compute the fluid domain NODAL_AREA values
        # Required by the parallel distance calculator if the distance has to be extended
        if (self.level_set_type == "continuous"):
            KratosMultiphysics.CalculateNodalAreaProcess(self.GetFluidComputingModelPart(), self.domain_size).Execute()

        # Initialize the Dirichlet-Neumann interface
        self.__InitializeFSIInterfaces()

        # Initialize the iteration value vector
        self.__InitializeIterationValueVector()

        # Initialize the distance field
        update_distance_process = True
        self.__GetDistanceToSkinProcess(update_distance_process).Execute()
        if (self.level_set_type == "continuous"):
            self.__ExtendLevelSet()

        # Initialize the embedded skin utility
        self.__GetEmbeddedSkinUtility()

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', "Finished initialization.")

    def AdvanceInTime(self, current_time):
        fluid_new_time = self.fluid_solver.AdvanceInTime(current_time)
        structure_new_time = self.structure_solver.AdvanceInTime(current_time)

        self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart().GetRootModelPart().CloneTimeStep(fluid_new_time)
        self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart().GetRootModelPart().CloneTimeStep(structure_new_time)

        if abs(fluid_new_time - structure_new_time) > 1e-12:
            err_msg =  'Fluid new time is: ' + str(fluid_new_time) + '\n'
            err_msg += 'Structure new time is: ' + str(structure_new_time) + '\n'
            err_msg += 'No substepping has been implemented yet. Fluid and structure time step must coincide.'
            raise Exception(err_msg)

        return fluid_new_time

    def InitializeSolutionStep(self):
        # Initialize solution step of fluid, structure and coupling solvers
        self.fluid_solver.InitializeSolutionStep()
        self.structure_solver.InitializeSolutionStep()
        self.__GetConvergenceAccelerator().InitializeSolutionStep()

    def Predict(self):
        # Structure solver prediction. It is important to firstly perform the structure
        # prediction to update the current buffer position before the FM-ALE operations.
        # Otherwise position 0 and 1 of the buffer coincide since the advance in time
        # has been already performed but no update has been done yet. Besides, this will
        # give a better approximation of the level-set position at the end of step.
        self.structure_solver.Predict()

        # Update the level set position with the structure prediction
        self.__UpdateLevelSet()

        # Correct the updated level set
        self.fluid_solver._GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

        # Fluid solver prediction
        self.fluid_solver.Predict()

        # Restore the fluid node fixity to its original status
        self.fluid_solver._GetDistanceModificationProcess().ExecuteFinalizeSolutionStep()

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
        return self.model.GetModelPart(self.__GetStructureInterfaceModelPartName())

    def GetStructureSkinElementBasedModelPart(self):
        # Create an auxiliar model part to save the element based skin
        element_based_skin_model_part_name = self.__GetStructureInterfaceModelPartName() + "ElementBased"
        if self.model.HasModelPart(element_based_skin_model_part_name):
            self.model.DeleteModelPart(element_based_skin_model_part_name)
        self.element_based_skin_model_part = self.model.CreateModelPart(element_based_skin_model_part_name)

        # Copy the skin model part conditions to an auxiliar model part elements.
        # This is required for the computation of the distance function, which
        # takes the elements of the second modelpart as skin. If this operation
        # is not performed, no elements are found, yielding a wrong level set.
        self.__GetPartitionedFSIUtilities().CopySkinToElements(
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

    def SolveSolutionStep(self):
        ## Non-linear coupling iteration ##
        nl_it = 0
        while (nl_it < self.max_nl_it and self.fluid_solver._TimeBufferIsInitialized()):
            nl_it += 1
            KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', '\tFSI non-linear iteration = ' + str(nl_it))

            self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it
            self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.CONVERGENCE_ACCELERATOR_ITERATION] = nl_it

            self.__GetConvergenceAccelerator().InitializeNonLinearIteration()

            # Map the RELAXED_DISP from the structure FSI coupling interface to fluid FSI coupling interface
            # This RELAXED_DISP is intended to be used to update the skin position before the level set update
            # Note that we take advance of the fact that the coupling interfaces coincide in the embedded case
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.RELAXED_DISPLACEMENT,
                KratosMultiphysics.DISPLACEMENT,
                self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                0)

            # Update the fluid FSI coupling interface position
            self.__GetFSICouplingInterfaceFluid().UpdatePosition()

            # Update the EMBEDDED_VELOCITY and solve the fluid problem
            self.__SolveFluid()

            # Transfer the fluid load to the structure FSI coupling
            self.__TransferFluidLoad()

            # Transfer fluid from the structure FSI coupling interface to father model part
            if (self.level_set_type == "continuous"):
                self.__GetFSICouplingInterfaceStructure().TransferValuesToFatherModelPart(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
            elif (self.level_set_type == "discontinuous"):
                self.__GetFSICouplingInterfaceStructure().TransferValuesToFatherModelPart(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
                self.__GetFSICouplingInterfaceStructure().TransferValuesToFatherModelPart(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
            else:
                err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
                raise Exception(err_msg)

            # Solve the structure problem
            self.__SolveStructure()

            # Compute the residual and perform the update
            dis_residual_norm = self.__GetFSICouplingInterfaceStructure().Update()

            # End the FSI non-linear iteration
            self.__GetConvergenceAccelerator().FinalizeNonLinearIteration()

            # Check convergence
            is_converged = self.__CheckFSIConvergence(dis_residual_norm)
            if (is_converged):
                return True

        return False

    def FinalizeSolutionStep(self):
        # Finalize solution step
        self.fluid_solver.FinalizeSolutionStep()
        self.structure_solver.FinalizeSolutionStep()
        self.__GetConvergenceAccelerator().FinalizeSolutionStep()

    def SetEchoLevel(self, structure_echo_level, fluid_echo_level):
        self.fluid_solver.SetEchoLevel(self, fluid_echo_level)
        self.structure_solver.SetEchoLevel(self, structure_echo_level)

    def Clear(self):
        self.fluid_solver.Clear()
        self.structure_solver.Clear()

    def Check(self):
        self.fluid_solver.Check()
        self.structure_solver.Check()

    #######################################################################
    ##############          PRIVATE METHODS SECTION          ##############
    #######################################################################

    def __GetDistanceToSkinProcess(self, update_distance_process = False):
        if update_distance_process:
            self._distance_to_skin_process = self.__CreateDistanceToSkinProcess()
        return self._distance_to_skin_process

    def __CreateDistanceToSkinProcess(self):
        # Set the distance computation process
        if (self.level_set_type == "continuous"):
            raycasting_relative_tolerance = 1.0e-10
            if self.domain_size == 2:
                return KratosMultiphysics.CalculateDistanceToSkinProcess2D(
                    self.GetFluidComputingModelPart(),
                    self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                    raycasting_relative_tolerance)
            elif self.domain_size == 3:
                return KratosMultiphysics.CalculateDistanceToSkinProcess3D(
                    self.GetFluidComputingModelPart(),
                    self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                    raycasting_relative_tolerance)
            else:
                raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))
        elif (self.level_set_type == "discontinuous"):
            if self.domain_size == 2:
                return KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess2D(
                    self.GetFluidComputingModelPart(),
                    self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart())
            elif self.domain_size == 3:
                return KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess3D(
                    self.GetFluidComputingModelPart(),
                    self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart())
            else:
                raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))
        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

    def __GetParallelDistanceCalculator(self):
        if not hasattr(self, '_parallel_distance_calculator'):
            self._parallel_distance_calculator = self.__CreateParallelDistanceCalculator()
        return self._parallel_distance_calculator

    def __CreateParallelDistanceCalculator(self):
        if self.domain_size == 2:
            return KratosMultiphysics.ParallelDistanceCalculator2D()
        elif self.domain_size == 3:
            return KratosMultiphysics.ParallelDistanceCalculator3D()
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))

    def __GetEmbeddedSkinUtility(self):
        if not hasattr(self, '_embedded_skin_utility'):
            self._embedded_skin_utility = self.__CreateEmbeddedSkinUtility()
        return self._embedded_skin_utility

    def __CreateEmbeddedSkinUtility(self):
        if self.domain_size == 2:
            return KratosMultiphysics.EmbeddedSkinUtility2D(
                self.GetFluidComputingModelPart(),
                self.__GetEmbedddedSkinUtilityModelPart(),
                self.level_set_type)
        elif self.domain_size == 3:
            return KratosMultiphysics.EmbeddedSkinUtility3D(
                self.GetFluidComputingModelPart(),
                self.__GetEmbedddedSkinUtilityModelPart(),
                self.level_set_type)
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))

    def __GetEmbedddedSkinUtilityModelPart(self):
        if not hasattr(self, '_embedded_skin_utility_model_part'):
            self._embedded_skin_utility_model_part = self.__CreateEmbeddedSkinUtilityModelPart()
        return self._embedded_skin_utility_model_part

    def __CreateEmbeddedSkinUtilityModelPart(self):
        embedded_skin_utility_skin_model_part_name = "EmbeddedSkinUtilityModelPart"
        embedded_skin_utility_skin_model_part =self.model.CreateModelPart(embedded_skin_utility_skin_model_part_name)
        embedded_skin_utility_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        if (self.level_set_type == "continuous"):
            embedded_skin_utility_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        elif (self.level_set_type == "discontinuous"):
            embedded_skin_utility_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
            embedded_skin_utility_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)
        return embedded_skin_utility_skin_model_part

    def __GetStructureInterfaceModelPartName(self):
        str_int_list = self.settings["coupling_settings"]["structure_interfaces_list"]
        if (str_int_list.size() != 1):
            raise Exception("FSI embedded solver structure skin must be contained in a unique model part")
        return str_int_list[0].GetString()

    def __InitializeIterationValueVector(self):
        # Note that the embedded FSI problem is defined in terms of the structure interface
        # Initialize the iteration value for the residual computation
        str_int_res_size = self.__GetPartitionedFSIUtilities().GetInterfaceResidualSize(self.__GetStructureInterfaceSubmodelPart())
        self.iteration_value = KratosMultiphysics.Vector(str_int_res_size)
        for i in range(0,str_int_res_size):
            self.iteration_value[i] = 0.0

    def __InitializeFSIInterfaces(self):
        # Initialize Neumann structure interface
        str_interface_submodelpart = self.model.GetModelPart(self.__GetStructureInterfaceModelPartName())

        # Set the INTERFACE flag to the structure skin
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, str_interface_submodelpart.Nodes)

    def __UpdateLevelSet(self):
        # Recompute the distance field with the obtained solution
        self.__GetDistanceToSkinProcess().Execute()

        # Extend the level set to the first layer of non-intersected elements
        # This is required in case the distance modification process moves the level set
        # to a non-intersected element to prevent almost empty fluid elements.
        # Note that non-intersected elements have a large default distance value, which might
        # alter the zero isosurface when the distance modification avoids almost empty elements.
        if (self.level_set_type == "continuous"):
            self.__ExtendLevelSet()

    def __ExtendLevelSet(self):
        max_layers = 2
        max_distance = 1.0e+12
        self.__GetParallelDistanceCalculator().CalculateDistances(
            self.GetFluidComputingModelPart(),
            KratosMultiphysics.DISTANCE,
            KratosMultiphysics.NODAL_AREA,
            max_layers,
            max_distance)

    def __SolveFluid(self):
        # Update the current iteration level-set position
        self.__UpdateLevelSet()

        # Solve fluid problem
        self.fluid_solver.SolveSolutionStep() # This contains the FM-ALE operations and the level set correction

    def __TransferFluidLoad(self):
        if (self.level_set_type == "continuous"):
            # Interpolate the pressure to the fluid FSI coupling interface
            self.__GetPartitionedFSIUtilities().EmbeddedPressureToPositiveFacePressureInterpolator(
                self.GetFluidComputingModelPart(),
                self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart())

            # Map PRESSURE from fluid FSI coupling interface to structure FSI coupling interface
            # Note that in here we take advantage of the fact that the coupling interfaces coincide in the embedded case
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.PRESSURE,
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                0)

        elif (self.level_set_type == "discontinuous"):
            # Generate the intersections skin to map from
            self.__GetEmbeddedSkinUtility().GenerateSkin()

            # Interpolate POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE from the background mesh
            self.__GetEmbeddedSkinUtility().InterpolateDiscontinuousMeshVariableToSkin(
                KratosMultiphysics.PRESSURE, KratosMultiphysics.POSITIVE_FACE_PRESSURE, "positive")
            self.__GetEmbeddedSkinUtility().InterpolateDiscontinuousMeshVariableToSkin(
                KratosMultiphysics.PRESSURE, KratosMultiphysics.NEGATIVE_FACE_PRESSURE, "negative")

            mapper_params = KratosMultiphysics.Parameters("""{
                "mapper_type": "nearest_element",
                "echo_level" : 0
            }""")
            mapper = KratosMapping.MapperFactory.CreateMapper(
                self.__GetEmbedddedSkinUtilityModelPart(),
                self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                mapper_params)
            mapper.Map(KratosMultiphysics.POSITIVE_FACE_PRESSURE, KratosMultiphysics.POSITIVE_FACE_PRESSURE)
            mapper.Map(KratosMultiphysics.NEGATIVE_FACE_PRESSURE, KratosMultiphysics.NEGATIVE_FACE_PRESSURE)

            # Transfer POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE from the fluid coupling interface to the structure one
            # Note that in here we take advantage of the fact that the coupling interfaces coincide in the embedded case
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                0)
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.NEGATIVE_FACE_PRESSURE,
                self.__GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                self.__GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                0)
        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

    def __SolveStructure(self):
        # Solve the structure problem
        self.structure_solver.SolveSolutionStep()

    def __CheckFSIConvergence(self, residual_norm):
        interface_dofs = self.__GetPartitionedFSIUtilities().GetInterfaceResidualSize(self.__GetStructureInterfaceSubmodelPart())
        normalised_residual = residual_norm/sqrt(interface_dofs)
        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', '\t|res|/sqrt(nDOFS) = ' + str(normalised_residual))
        return normalised_residual < self.nl_tol

    # This method returns the convergence accelerator.
    # If it is not created yet, it calls the __CreateConvergenceAccelerator first
    def __GetConvergenceAccelerator(self):
        if not hasattr(self, '_convergence_accelerator'):
            self._convergence_accelerator = self.__CreateConvergenceAccelerator()
        return self._convergence_accelerator

    # This method constructs the convergence accelerator coupling utility
    def __CreateConvergenceAccelerator(self):
        convergence_accelerator = convergence_accelerator_factory.CreateConvergenceAccelerator(self.settings["coupling_settings"]["coupling_strategy_settings"])
        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Coupling strategy construction finished')
        return convergence_accelerator

    def __GetFSICouplingInterfaceStructure(self):
        if not hasattr(self, '_fsi_coupling_interface_structure'):
            self._fsi_coupling_interface_structure = self.__CreateFSICouplingInterfaceStructure()
        return self._fsi_coupling_interface_structure

    def __CreateFSICouplingInterfaceStructure(self):
        # Set auxiliary settings
        if (self.level_set_type == "continuous"):
            aux_settings = KratosMultiphysics.Parameters(
            """{
                "model_part_name": "FSICouplingInterfaceStructure",
                "parent_model_part_name": "",
                "input_variable_list": ["POSITIVE_FACE_PRESSURE"],
                "output_variable_list": ["DISPLACEMENT"]
            }""")
        elif (self.level_set_type == "discontinuous"):
            aux_settings = KratosMultiphysics.Parameters(
            """{
                "model_part_name": "FSICouplingInterfaceStructure",
                "parent_model_part_name": "",
                "input_variable_list": ["POSITIVE_FACE_PRESSURE","NEGATIVE_FACE_PRESSURE"],
                "output_variable_list": ["DISPLACEMENT"]
            }""")
        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

        aux_settings["parent_model_part_name"].SetString(self.structure_interface_submodelpart_name)

        # Construct the FSI coupling interface
        fsi_coupling_interface_structure = fsi_coupling_interface.FSICouplingInterface(
            self.model,
            aux_settings,
            self.__GetConvergenceAccelerator())

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Structure FSI coupling interface created')

        return fsi_coupling_interface_structure

    def __GetFSICouplingInterfaceFluid(self):
        if not hasattr(self, '_fsi_coupling_interface_fluid'):
            self._fsi_coupling_interface_fluid = self.__CreateFSICouplingInterfaceFluid()
        return self._fsi_coupling_interface_fluid

    def __CreateFSICouplingInterfaceFluid(self):
        # Set auxiliary settings
        # Note that in the embedded case, the fluid interface is identical to the structure one
        # This is intentionally done, since this copy will be used in the level set computation
        if (self.level_set_type == "continuous"):
            aux_settings = KratosMultiphysics.Parameters(
            """{
                "model_part_name": "FSICouplingInterfaceFluid",
                "parent_model_part_name": "",
                "input_variable_list": ["DISPLACEMENT"],
                "output_variable_list": ["PRESSURE"]
            }""")
        elif (self.level_set_type == "discontinuous"):
            aux_settings = KratosMultiphysics.Parameters(
            """{
                "model_part_name": "FSICouplingInterfaceFluid",
                "parent_model_part_name": "",
                "input_variable_list": ["DISPLACEMENT"],
                "output_variable_list": ["POSITIVE_FACE_PRESSURE","NEGATIVE_FACE_PRESSURE"]
            }""")
        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

        aux_settings["parent_model_part_name"].SetString(self.structure_interface_submodelpart_name)

        # Construct the FSI coupling interface
        fsi_coupling_interface_fluid = fsi_coupling_interface.FSICouplingInterface(
            self.model,
            aux_settings)

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Fluid FSI coupling interface created')

        return fsi_coupling_interface_fluid

    def __GetStructureInterfaceSubmodelPart(self):
        # Returns the structure interface submodelpart that will be used in the residual minimization
        return self.model.GetModelPart(self.structure_interface_submodelpart_name)

    def __GetDomainSize(self):
        fluid_domain_size = self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        structure_domain_size = self.structure_solver.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        if fluid_domain_size !=structure_domain_size:
            raise("ERROR: Solid domain size and fluid domain size are not equal!")

        return fluid_domain_size

    def __GetPartitionedFSIUtilities(self):
        if not hasattr(self, '_partitioned_fsi_utilities'):
            self._partitioned_fsi_utilities = self.__CreatePartitionedFSIUtilities()
        return self._partitioned_fsi_utilities

    def __CreatePartitionedFSIUtilities(self):
        if self.domain_size == 2:
            return KratosFSI.PartitionedFSIUtilitiesArray2D()
        elif self.domain_size == 3:
            return KratosFSI.PartitionedFSIUtilitiesArray3D()
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))
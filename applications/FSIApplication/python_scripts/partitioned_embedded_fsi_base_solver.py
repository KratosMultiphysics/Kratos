from math import sqrt   # Import the square root from python library

# Import utilities
from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid            # Import the fluid Python solvers wrapper
from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural # Import the structure Python solvers wrapper
from KratosMultiphysics.FSIApplication import fsi_coupling_interface                            # Import the FSI coupling interface utility
from KratosMultiphysics.FSIApplication import convergence_accelerator_factory                   # Import the FSI convergence accelerator factory

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

# Import base class
from KratosMultiphysics.FSIApplication.partitioned_fsi_base_solver import PartitionedFSIBaseSolver

def CreateSolver(model, project_parameters):
    return PartitionedEmbeddedFSIBaseSolver(model, project_parameters)

class PartitionedEmbeddedFSIBaseSolver(PartitionedFSIBaseSolver):

    def __init__(self, model, project_parameters):
        # Call the base solver constructor
        super().__init__(model, project_parameters)

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Partitioned embedded FSI base solver construction finished')

    @classmethod
    def GetDefaultParameters(cls):

        # Note that only the coupling settings are validated
        # The subdomain solver settings will be validated while instantiating these
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
                "coupling_strategy_settings": {
                    "abs_cut_off_tol": 1e-06,
                    "solver_type": "MVQN",
                    "w_0": 0.5
                },
                "nl_max_it": 30,
                "nl_tol": 1e-07,
                "structure_interfaces_list": []
            }
        }""")

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    #TODO: Use the base solver one once we use the fluid ALE solver
    def AddVariables(self):
        # Fluid and structure solvers variables addition
        self.fluid_solver.AddVariables()
        self.structure_solver.AddVariables()

    #TODO: Use the base solver one once we use the fluid ALE solver
    def ImportModelPart(self):
        # Fluid and structure solvers ImportModelPart() call
        self.fluid_solver.ImportModelPart()
        self.structure_solver.ImportModelPart()

    #TODO: Use the base solver one once we use the fluid ALE solver
    def PrepareModelPart(self):
        # Fluid and structure solvers PrepareModelPart() call
        self.fluid_solver.PrepareModelPart()
        self.structure_solver.PrepareModelPart()

        # Perform all the operations required to set up the coupling interfaces
        self._InitializeCouplingInterfaces()

    #TODO: Use the base solver one once we use the fluid ALE solver
    def AddDofs(self):
        # Add DOFs structure
        self.structure_solver.AddDofs()
        # Add DOFs fluid
        self.fluid_solver.AddDofs()

    def Initialize(self):
        # Coupling utility initialization
        # The _GetConvergenceAccelerator is supposed to construct the convergence accelerator in here
        self._GetConvergenceAccelerator().Initialize()

        # Python structure solver initialization
        self.structure_solver.Initialize()

        # Python fluid solver initialization
        self.fluid_solver.Initialize()

        # Compute the fluid domain NODAL_AREA values
        # Required by the parallel distance calculator if the distance has to be extended
        if (self.level_set_type == "continuous"):
            KratosMultiphysics.CalculateNodalAreaProcess(self.GetFluidComputingModelPart(), self._GetDomainSize()).Execute()

        # Initialize the distance field
        update_distance_process = True
        self.__GetDistanceToSkinProcess(update_distance_process).Execute()
        if (self.level_set_type == "continuous"):
            self.__ExtendLevelSet()

        # Initialize the embedded skin utility
        self.__GetEmbeddedSkinUtility()

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', "Finished initialization.")

    #TODO: Use the base one once the body fitted uses the fluid ALE solver
    def InitializeSolutionStep(self):
        # Initialize solution step of fluid, structure and coupling solvers
        self.fluid_solver.InitializeSolutionStep()
        self.structure_solver.InitializeSolutionStep()
        self._GetConvergenceAccelerator().InitializeSolutionStep()

    def Predict(self):
        # Structure solver prediction. It is important to firstly perform the structure
        # prediction to update the current buffer position before the FM-ALE operations.
        # Otherwise position 0 and 1 of the buffer coincide since the advance in time
        # has been already performed but no update has been done yet. Besides, this will
        # give a better approximation of the level-set position at the end of step.
        self.structure_solver.Predict()

        # Update the level set position with the structure prediction
        self.__UpdateLevelSet()

        # Fluid solver prediction
        self.fluid_solver.Predict()

    #TODO: Use the base solver one once we use the fluid ALE solver for the fluid
    def FinalizeSolutionStep(self):
        # Finalize solution step
        self.fluid_solver.FinalizeSolutionStep()
        self.structure_solver.FinalizeSolutionStep()
        self._GetConvergenceAccelerator().FinalizeSolutionStep()

    #TODO: Use the base solver one once we use the fluid ALE solver for the fluid
    def Finalize(self):
        self.fluid_solver.Finalize()
        self.structure_solver.Finalize()
        self._GetConvergenceAccelerator().Finalize()

    #######################################################################
    ##############          PRIVATE METHODS SECTION          ##############
    #######################################################################

    def _AuxiliaryInitOperations(self):
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
        self.level_set_type = self.settings["fluid_solver_settings"]["formulation"]["level_set_type"].GetString()

        # First call to create the embedded intersections model part
        self.__GetEmbedddedSkinUtilityModelPart()

    def _AdvanceInTimeCouplingInterfaces(self, new_time):
        # Even though these are auxiliary model parts, this is mandatory to be done to properly set up the database
        # Note that if this operations are removed, some auxiliary utils (e.g. FM-ALE algorithm in embedded) will perform wrong
        self._GetFSICouplingInterfaceFluid().GetInterfaceModelPart().GetRootModelPart().CloneTimeStep(new_time)
        self._GetFSICouplingInterfaceFluid().GetInterfaceModelPart().ProcessInfo[KratosMultiphysics.STEP] = self._GetFSICouplingInterfaceFluid().GetFatherModelPart().ProcessInfo[KratosMultiphysics.STEP]
        self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart().GetRootModelPart().CloneTimeStep(new_time)
        self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart().ProcessInfo[KratosMultiphysics.STEP] = self._GetFSICouplingInterfaceStructure().GetFatherModelPart().ProcessInfo[KratosMultiphysics.STEP]

    def _InitializeCouplingInterfaces(self):
        # FSI interface coupling interfaces initialization
        # The getter methods are to construct the FSI coupling structure interface in here
        self._GetFSICouplingInterfaceFluid().GetInterfaceModelPart()
        self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart()

        # Set the INTERFACE flag to the structure skin
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.INTERFACE, True, self._GetStructureInterfaceSubmodelPart().Nodes)

    def __GetDistanceToSkinProcess(self, update_distance_process = False):
        if update_distance_process:
            self._distance_to_skin_process = self.__CreateDistanceToSkinProcess()
        return self._distance_to_skin_process

    def __CreateDistanceToSkinProcess(self):
        # Set the distance computation process
        if (self.level_set_type == "continuous"):
            raycasting_relative_tolerance = 1.0e-10
            if self._GetDomainSize() == 2:
                return KratosMultiphysics.CalculateDistanceToSkinProcess2D(
                    self.GetFluidComputingModelPart(),
                    self._GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                    raycasting_relative_tolerance)
            elif self._GetDomainSize() == 3:
                return KratosMultiphysics.CalculateDistanceToSkinProcess3D(
                    self.GetFluidComputingModelPart(),
                    self._GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                    raycasting_relative_tolerance)
            else:
                raise Exception("Domain size expected to be 2 or 3. Got " + str(self._GetDomainSize()))
        elif (self.level_set_type == "discontinuous"):
            discontinuous_distance_settings = KratosMultiphysics.Parameters("""{
                "calculate_elemental_edge_distances" : true,
                "calculate_elemental_edge_distances_extrapolated" : true
            }""")
            if self._GetDomainSize() == 2:
                return KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess2D(
                    self.GetFluidComputingModelPart(),
                    self._GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                    discontinuous_distance_settings)
            elif self._GetDomainSize() == 3:
                return KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess3D(
                    self.GetFluidComputingModelPart(),
                    self._GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                    discontinuous_distance_settings)
            else:
                raise Exception("Domain size expected to be 2 or 3. Got " + str(self._GetDomainSize()))
        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

    def __GetParallelDistanceCalculator(self):
        if not hasattr(self, '_parallel_distance_calculator'):
            self._parallel_distance_calculator = self.__CreateParallelDistanceCalculator()
        return self._parallel_distance_calculator

    def __CreateParallelDistanceCalculator(self):
        parallel_redistance_settings = KratosMultiphysics.Parameters("""{
            "max_levels" : 2,
            "max_distance": 1e12
        }""")
        if self._GetDomainSize() == 2:
            return KratosMultiphysics.ParallelDistanceCalculationProcess2D(
                self.GetFluidComputingModelPart(),
                parallel_redistance_settings)
        elif self._GetDomainSize() == 3:
            return KratosMultiphysics.ParallelDistanceCalculationProcess3D(
                self.GetFluidComputingModelPart(),
                parallel_redistance_settings)
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self._GetDomainSize()))

    def __GetEmbeddedSkinUtility(self):
        if not hasattr(self, '_embedded_skin_utility'):
            self._embedded_skin_utility = self.__CreateEmbeddedSkinUtility()
        return self._embedded_skin_utility

    def __CreateEmbeddedSkinUtility(self):
        if self._GetDomainSize() == 2:
            return KratosMultiphysics.EmbeddedSkinUtility2D(
                self.GetFluidComputingModelPart(),
                self.__GetEmbedddedSkinUtilityModelPart(),
                self.level_set_type)
        elif self._GetDomainSize() == 3:
            return KratosMultiphysics.EmbeddedSkinUtility3D(
                self.GetFluidComputingModelPart(),
                self.__GetEmbedddedSkinUtilityModelPart(),
                self.level_set_type)
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self._GetDomainSize()))

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
        self.__GetParallelDistanceCalculator().Execute()

    def _MapStructureInterfaceDisplacement(self):
        # Map the RELAXED_DISP from the structure FSI coupling interface to fluid FSI coupling interface
        # Note that we take advance of the fact that the coupling interfaces coincide in the embedded case
        # Then update the fluid FSI coupling interface position
        KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
            KratosMultiphysics.RELAXED_DISPLACEMENT,
            KratosMultiphysics.DISPLACEMENT,
            self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
            self._GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
            0)
        self._GetFSICouplingInterfaceFluid().UpdatePosition(KratosMultiphysics.DISPLACEMENT)

    def _SolveFluid(self):
        # Update the current iteration level-set position
        self.__UpdateLevelSet()

        # Solve fluid problem
        self.fluid_solver.SolveSolutionStep() # This contains the FM-ALE operations and the level set correction

    def _CalculateFluidInterfaceTraction(self):
        if (self.level_set_type == "continuous"):
            # Interpolate the pressure to the fluid FSI coupling interface
            self._GetPartitionedFSIUtilities().EmbeddedPressureToPositiveFacePressureInterpolator(
                self.GetFluidComputingModelPart(),
                self._GetFSICouplingInterfaceFluid().GetInterfaceModelPart())

        elif (self.level_set_type == "discontinuous"):
            # Generate the intersections skin to map from
            self.__GetEmbeddedSkinUtility().GenerateSkin()

            # Interpolate POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE from the background mesh
            self.__GetEmbeddedSkinUtility().InterpolateDiscontinuousMeshVariableToSkin(
                KratosMultiphysics.PRESSURE, KratosMultiphysics.POSITIVE_FACE_PRESSURE, "positive")
            self.__GetEmbeddedSkinUtility().InterpolateDiscontinuousMeshVariableToSkin(
                KratosMultiphysics.PRESSURE, KratosMultiphysics.NEGATIVE_FACE_PRESSURE, "negative")

        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

    def _MapFluidInterfaceTraction(self):
        if (self.level_set_type == "continuous"):
            # Map PRESSURE from fluid FSI coupling interface to structure FSI coupling interface
            # Note that in here we take advantage of the fact that the coupling interfaces coincide in the embedded case
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                self._GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                0)

            # Convert the pressure scalar load to a traction vector one
            swap_traction_sign = True
            self._GetPartitionedFSIUtilities().CalculateTractionFromPressureValues(
                self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                self._GetTractionVariable(),
                swap_traction_sign)

        elif (self.level_set_type == "discontinuous"):
            # Map the POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE from the auxiliary embedded skin model part,
            # which is created from the elemental level set intersections, to the fluid FSI coupling interface
            # Note that the mapper instance is created each time as the embedded skin mesh potentially changes at each iteration
            mapper_params = KratosMultiphysics.Parameters("""{
                "mapper_type": "nearest_element",
                "echo_level" : 0
            }""")
            mapper = KratosMultiphysics.MapperFactory.CreateMapper(
                self.__GetEmbedddedSkinUtilityModelPart(),
                self._GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                mapper_params)
            mapper.Map(KratosMultiphysics.POSITIVE_FACE_PRESSURE, KratosMultiphysics.POSITIVE_FACE_PRESSURE)
            mapper.Map(KratosMultiphysics.NEGATIVE_FACE_PRESSURE, KratosMultiphysics.NEGATIVE_FACE_PRESSURE)

            # Transfer POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE from the fluid coupling interface to the structure one
            # Note that in here we take advantage of the fact that the coupling interfaces coincide in the embedded case
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                self._GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                0)
            KratosMultiphysics.VariableUtils().CopyModelPartNodalVar(
                KratosMultiphysics.NEGATIVE_FACE_PRESSURE,
                self._GetFSICouplingInterfaceFluid().GetInterfaceModelPart(),
                self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                0)

            # Convert the pressure scalar load to a traction vector one
            swap_traction_sign = True
            self._GetPartitionedFSIUtilities().CalculateTractionFromPressureValues(
                self._GetFSICouplingInterfaceStructure().GetInterfaceModelPart(),
                KratosMultiphysics.POSITIVE_FACE_PRESSURE,
                KratosMultiphysics.NEGATIVE_FACE_PRESSURE,
                self._GetTractionVariable(),
                swap_traction_sign)

        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

    #TODO: Use the base solver one once the body-fitted does not require to do the traction to point load conversion
    def _SolveStructure(self):
        # Solve the structure problem
        self.structure_solver.SolveSolutionStep()

    def _CreateFSICouplingInterfaceStructure(self):
        # Set auxiliary settings
        if (self.level_set_type == "continuous"):
            aux_settings = KratosMultiphysics.Parameters(
            """{
                "model_part_name": "FSICouplingInterfaceStructure",
                "parent_model_part_name": "",
                "input_variable_list": [],
                "output_variable_list": ["DISPLACEMENT"],
                "auxiliary_variable_list": ["POSITIVE_FACE_PRESSURE","NORMAL"]
            }""")
        elif (self.level_set_type == "discontinuous"):
            aux_settings = KratosMultiphysics.Parameters(
            """{
                "model_part_name": "FSICouplingInterfaceStructure",
                "parent_model_part_name": "",
                "input_variable_list": [],
                "output_variable_list": ["DISPLACEMENT"],
                "auxiliary_variable_list": ["POSITIVE_FACE_PRESSURE","NEGATIVE_FACE_PRESSURE","NORMAL"]
            }""")
        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

        aux_settings["parent_model_part_name"].SetString(self.structure_interface_submodelpart_name)
        aux_settings["input_variable_list"].Append(self._GetTractionVariable().Name())

        # Construct the FSI coupling interface
        fsi_coupling_interface_structure = fsi_coupling_interface.FSICouplingInterface(
            self.model,
            aux_settings,
            self._GetConvergenceAccelerator())

        KratosMultiphysics.Logger.PrintInfo('PartitionedEmbeddedFSIBaseSolver', 'Structure FSI coupling interface created')

        return fsi_coupling_interface_structure

    def _GetFSICouplingInterfaceFluid(self):
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
                "output_variable_list": ["POSITIVE_FACE_PRESSURE"]
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

    @classmethod
    def _GetFSICouplingInterfaceFluidPositive(self):
        err_msg = 'Embedded partitioned FSI solver does not implement the \'_GetFSICouplingInterfaceFluidPositive\' method.\n'
        err_msg += 'Use \'_GetFSICouplingInterfaceFluid\' method as the fluid coupling interface is unique.'
        raise Exception(err_msg)

    @classmethod
    def _GetFSICouplingInterfaceFluidNegative(self):
        err_msg = 'Embedded partitioned FSI solver does not implement the \'_GetFSICouplingInterfaceFluidNegative\' method.\n'
        err_msg += 'Use \'_GetFSICouplingInterfaceFluid\' method as the fluid coupling interface is unique.'
        raise Exception(err_msg)

    def _GetTractionVariable(self):
        if self._GetDomainSize() == 2:
            return KratosStructural.LINE_LOAD
        elif self._GetDomainSize() == 3:
            return KratosStructural.SURFACE_LOAD
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self._GetDomainSize()))

        # if self._GetConvergenceAccelerator().IsBlockNewton():
        #     if self._GetDomainSize() == 2:
        #         return KratosStructural.LINE_LOAD
        #     elif self._GetDomainSize() == 3:
        #         return KratosStructural.SURFACE_LOAD
        #     else:
        #         raise Exception("Domain size expected to be 2 or 3. Got " + str(self._GetDomainSize()))
        # else:
        #     if self.level_set_type == "continuous":
        #         return KratosMultiphysics.POSITIVE_FACE_PRESSURE
        #     elif self.level_set_type == "discontinuous":
        #         if self._GetDomainSize() == 2:
        #             return KratosStructural.LINE_LOAD
        #         elif self._GetDomainSize() == 3:
        #             return KratosStructural.SURFACE_LOAD
        #         else:
        #             raise Exception(
        #                 "Domain size expected to be 2 or 3. Got " + str(self._GetDomainSize()))
        #     else:
        #         raise Exception("Wrong level set type '{}'".format(self.level_set_type))

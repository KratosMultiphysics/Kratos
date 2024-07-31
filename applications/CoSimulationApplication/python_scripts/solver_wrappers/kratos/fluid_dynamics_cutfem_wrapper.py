# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing FluidDynamics
if not CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    raise ImportError("The FluidDynamicsApplication is not available!")
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

# Importing MeshMoving
#import KratosMultiphysics.MeshMovingApplication as KMM

# Importing FSI
if not CheckIfApplicationsAvailable("FSIApplication"):
    raise ImportError("The FSIApplication is not available!")
from KratosMultiphysics.FSIApplication import PartitionedFSIUtilitiesArray2D, PartitionedFSIUtilitiesArray3D

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import data_communicator_utilities


def Create(settings, model, solver_name):
    return FluidDynamicsCutFEMWrapper(settings, model, solver_name)


class FluidDynamicsCutFEMWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the FluidDynamicsApplication of Kratos when using the navier_stokes_embedded_solver.
    TODO Not yet MPI parallelized (--> trilinos for FSI utilities and FM-ALE??).
    TODO Test for continuous, discontinuous, MPI!
    TODO Are predictor, mapper and convergence accelerator working as expected?
    """

    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)
        self.interface_name = settings["solver_wrapper_settings"]["fluid_interface_name"].GetString()
        self.structure_interface_name = settings["solver_wrapper_settings"]["structure_interface_name"].GetString()
        self.level_set_type = settings["solver_wrapper_settings"]["level_set_type"].GetString()

        # Get and check domain size
        self.domain_size = self.project_parameters["solver_settings"]["domain_size"].GetInt()
        if self.domain_size != 2 and self.domain_size != 3:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self.domain_size))

        # For moving the coupling interface model part inside the wrapper
        #self.previous_displacement = [0.0, 0.0, 0.0]
        #self.previous_time = 0.0

        #TODO - for discontinuous only?
        # First call to create the embedded intersections model part
        self.__GetEmbedddedSkinUtilityModelPart()

    def Initialize(self):
        # Create interface model part for embedded fluid solver (distance calculation and FM_ALE)
        #NOTE this is done before the initialize of the solvers
        #NOTE the interface entities (nodes and elements) are taken from the structural interface on CoSim lever after solver initialization
        self.__CreateCutFEMInterfaceModelPart()

        super().Initialize()

        # Get domain size from computing model part (after solver initialization)
        # self.domain_size = self.__GetFluidComputingModelPart().ProcessInfo[KM.DOMAIN_SIZE]

        # Compute the fluid domain NODAL_AREA values
        # Required by the parallel distance calculator if the distance has to be extended
        if (self.level_set_type == "continuous"):
            KM.CalculateNodalAreaProcess(self.__GetFluidComputingModelPart(), self.domain_size).Execute()
        self.__UpdateLevelSet()

        # Create fsi utility for pressure interpolation from fluid background mesh to interface
        #self.__GetPartitionedFSIUtility()

        #TODO Initialize the embedded skin utility - for discontinuous only??
        self.__GetEmbeddedSkinUtility()

    def AdvanceInTime(self, current_time):
        new_time = super().AdvanceInTime(current_time)

        #self.__MoveInterfaceModelPart()  # for rigid-body solver
        self.__AdvanceCutFEMInterfaceInTime(new_time)

        return new_time

    def Predict(self):
        # Update level set before predicting the embedded fluid
        # NOTE should be done after the structure has been predicted and the interface has been moved because of FM-ALE operations!
        # Otherwise position 0 and 1 of the buffer coincide for FM-ALE operations since the advance in time has been already performed but no update has been done yet.
        #TODO Is there actually an update of interface position after Predict() happening?!
        self.__UpdateLevelSet()
        super().Predict()

    def SolveSolutionStep(self):
        # Update the level set from the interface model part before solving the CutFEM fluid
        # NOTE should be done after the structure has been solved and its movement has been transferred to the interface model part
        self.__UpdateLevelSet()

        # Solve fluid for current coupling iteration
        super().SolveSolutionStep()

        # Calculate resultants of the fluid for the interface model part
        self.__CalculateResultantsOfCutFEMInterfaceModelPart()

    def _CreateAnalysisStage(self):
        return FluidDynamicsAnalysis(self.model, self.project_parameters)

    def _GetDataCommunicator(self):
        if not KM.IsDistributedRun():
            return KM.ParallelEnvironment.GetDataCommunicator("Serial")

        #TODO does not work for MPI ?! --> TrilinosPartitionedFSIUtilitiesArray3D --> FM-ALE operations??

        # now we know that Kratos runs in MPI
        parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()

        # first check if the solver uses MPI
        if parallel_type != "MPI":
            return data_communicator_utilities.GetRankZeroDataCommunicator()

        # now we know that the solver uses MPI, only question left is whether to use all ranks or a subset
        if self.project_parameters["solver_settings"]["solver_type"].GetString() == "ale_fluid":
            model_import_settings = self.project_parameters["solver_settings"]["fluid_solver_settings"]["model_import_settings"]
        else:
            model_import_settings = self.project_parameters["solver_settings"]["model_import_settings"]
        self._CheckDataCommunicatorIsConsistentlyDefined(model_import_settings, self.settings["mpi_settings"])

        return super()._GetDataCommunicator()

    def __UpdateLevelSet(self):
        # Update level set distances
        if not hasattr(self, '_distance_to_skin_process'):
            self._distance_to_skin_process = self.__CreateDistanceToSkinProcess()
            self._distance_to_skin_process.Execute()
        else:
            self._distance_to_skin_process.Execute()

        # Extend level set
        if self.level_set_type == "continuous":
            if not hasattr(self, '_parallel_distance_calculator'):
                self._parallel_distance_calculator = self.__CreateParallelDistanceCalculator()
                self._parallel_distance_calculator.Execute()
            else:
                self._parallel_distance_calculator.Execute()

    def __CreateDistanceToSkinProcess(self):
        if self.level_set_type == "continuous":
            raycasting_relative_tolerance = 1.0e-10
            if self.domain_size == 2:
                distance_to_skin_process = KM.CalculateDistanceToSkinProcess2D(
                    self.__GetFluidComputingModelPart(),
                    self._interface_model_part,
                    raycasting_relative_tolerance)
            else:
                distance_to_skin_process = KM.CalculateDistanceToSkinProcess3D(
                    self.__GetFluidComputingModelPart(),
                    self._interface_model_part,
                    raycasting_relative_tolerance)

        elif self.level_set_type == "discontinuous":
            discontinuous_distance_settings = KM.Parameters("""{
                "calculate_elemental_edge_distances" : true,
                "calculate_elemental_edge_distances_extrapolated" : true
            }""")
            if self.domain_size == 2:
                return KM.CalculateDiscontinuousDistanceToSkinProcess2D(
                    self.__GetFluidComputingModelPart(),
                    self._interface_model_part,
                    discontinuous_distance_settings)
            else:
                return KM.CalculateDiscontinuousDistanceToSkinProcess3D(
                    self.__GetFluidComputingModelPart(),
                    self._interface_model_part,
                    discontinuous_distance_settings)
        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

        return distance_to_skin_process

    def __CreateParallelDistanceCalculator(self):
        parallel_redistance_settings = KM.Parameters("""{
            "max_levels" : 2,
            "max_distance": 1e12
        }""")
        if self.domain_size == 2:
            return KM.ParallelDistanceCalculationProcess2D(
                self.__GetFluidComputingModelPart(),
                parallel_redistance_settings)
        else:
            return KM.ParallelDistanceCalculationProcess3D(
                self.__GetFluidComputingModelPart(),
                parallel_redistance_settings)

    def __GetPartitionedFSIUtility(self):
        if not hasattr(self, '_partitioned_fsi_utility'):
            self.__CreatePartitionedFSIUtility()
        return self._partitioned_fsi_utility

    def __CreatePartitionedFSIUtility(self):
        #TODO TrilinosPartitionedFSIUtilities for MPI (of structure?)
        if self.domain_size == 2:
            self._partitioned_fsi_utility = PartitionedFSIUtilitiesArray2D()
        else:
            self._partitioned_fsi_utility = PartitionedFSIUtilitiesArray3D()

    def __GetFluidComputingModelPart(self):
        return self._analysis_stage._GetSolver().GetComputingModelPart()

    def __CreateCutFEMInterfaceModelPart(self):
        # Create interface model part for calculating the embedded distances, the FM-ALE utility and data mapping
        self._interface_model_part = self.model.CreateModelPart(self.interface_name)
        self._interface_model_part.ProcessInfo[KM.DOMAIN_SIZE] = self.domain_size

        # Add required variables
        self._interface_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self._interface_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self._interface_model_part.AddNodalSolutionStepVariable(KM.POSITIVE_FACE_PRESSURE)
        if self.level_set_type == "discontinuous":
            self._interface_model_part.AddNodalSolutionStepVariable(KM.NEGATIVE_FACE_PRESSURE)
        self._interface_model_part.AddNodalSolutionStepVariable(KM.REACTION)
        self._interface_model_part.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)

        # # Read the square geometry for rigid-body solver
        # KM.ModelPartIO('square', KM.ModelPartIO.READ | KM.ModelPartIO.SKIP_TIMER).ReadModelPart(self._interface_model_part)

    def __MoveInterfaceModelPart(self):
        # # Get global displacements
        # displacement = self._interface_model_part[KMC.GLOBAL_DISPLACEMENT]

        # # Rotation is not considered
        # angle = 0
        # axis = [1,0,0]

        # # Reference point
        # reference_point = [0.0, 0.0, 0.0]

        # # Apply displacement to mesh --> calculates MESH_DISPLACEMENT
        # KMM.MoveModelPart(
        #     self._interface_model_part,
        #     axis,                 # rotation axis
        #     angle,                # rotation angle
        #     reference_point, # one point of the rotation axis
        #     displacement)         # translation

        # Advance interface model part in time
        t = self.__GetFluidComputingModelPart().ProcessInfo[KM.TIME]
        self._interface_model_part.CloneTimeStep(t)
        self._interface_model_part.ProcessInfo[KM.STEP] = self.__GetFluidComputingModelPart().ProcessInfo[KM.STEP]

        # Move the interface model part
        # Rotate model part as well?
        #displacement = [0.25*t, 0.0, 0.0]
        displacement = self._interface_model_part[KMC.GLOBAL_DISPLACEMENT]
        dt = t - self.previous_time
        velocity = [(u_n-u_nn)/dt for u_n, u_nn in zip(displacement, self.previous_displacement)]
        self.previous_displacement = displacement
        self.previous_time = t
        for node in self._interface_model_part.Nodes:
            node.X = node.X0 + displacement[0]
            node.Y = node.Y0 + displacement[1]
            node.SetSolutionStepValue(KM.DISPLACEMENT, displacement)
            #node.SetSolutionStepValue(KM.VELOCITY, velocity)

    def __AdvanceCutFEMInterfaceInTime(self, new_time):
        # NOTE it is mandatory to advance this auxiliary model part in time to properly set up the database,
        # otherwise e.g. FM-ALE operations will perform wrong
        self._interface_model_part.GetRootModelPart().CloneTimeStep(new_time)
        self._interface_model_part.ProcessInfo[KM.STEP] = self.__GetFluidComputingModelPart().ProcessInfo[KM.STEP]

    def __CalculateResultantsOfCutFEMInterfaceModelPart(self):
        if self.level_set_type == "continuous":
            # Interpolate the pressure from the fluid background mesh to the fluid coupling interface
            self.__GetPartitionedFSIUtility().EmbeddedPressureToPositiveFacePressureInterpolator(
                self.__GetFluidComputingModelPart(),
                self._interface_model_part)

            # Calculate REACTION (traction vector) and REACTION_MOMENT from POSITIVE_FACE_PRESSURE (scalar)
            # NOTE shear is neglected here?!
            # NOTE traction sign is not swapped as the normal considered is the positive interface outwards one
            # which already points to the structure (unlike standard body-fitted solver)
            self.__GetPartitionedFSIUtility().CalculateTractionFromPressureValues(
                self._interface_model_part,
                KM.POSITIVE_FACE_PRESSURE,
                KM.REACTION,
                False  # SwapTractionSign
            )

        elif self.level_set_type == "discontinuous":
            # Generate the intersections skin to map from
            self.__GetEmbeddedSkinUtility().GenerateSkin()

            # Interpolate POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE from the background mesh
            self.__GetEmbeddedSkinUtility().InterpolateDiscontinuousMeshVariableToSkin(
                KM.PRESSURE, KM.POSITIVE_FACE_PRESSURE, "positive")
            self.__GetEmbeddedSkinUtility().InterpolateDiscontinuousMeshVariableToSkin(
                KM.PRESSURE, KM.NEGATIVE_FACE_PRESSURE, "negative")

            # Map the POSITIVE_FACE_PRESSURE and NEGATIVE_FACE_PRESSURE from the auxiliary embedded skin model part,
            # which is created from the elemental level set intersections, to the fluid FSI coupling interface
            # Note that the mapper instance is created each time as the embedded skin mesh potentially changes at each iteration
            mapper_params = KM.Parameters("""{
                "mapper_type": "nearest_element",
                "echo_level" : 0
            }""")
            mapper = KM.MapperFactory.CreateMapper(
                self.__GetEmbedddedSkinUtilityModelPart(),
                self._interface_model_part,
                mapper_params)
            mapper.Map(KM.POSITIVE_FACE_PRESSURE, KM.POSITIVE_FACE_PRESSURE)
            mapper.Map(KM.NEGATIVE_FACE_PRESSURE, KM.NEGATIVE_FACE_PRESSURE)

            # Calculate REACTION (traction vector) and REACTION_MOMENT from POSITIVE_FACE_PRESSURE (scalar)
            # NOTE traction sign is not swapped as the normal considered is the positive interface outwards one
            # which already points to the structure (unlike standard body-fitted solver)
            swap_traction_sign = False
            self.__GetPartitionedFSIUtility().CalculateTractionFromPressureValues(
                self._interface_model_part,
                KM.POSITIVE_FACE_PRESSURE,
                KM.NEGATIVE_FACE_PRESSURE,
                KM.REACTION,
                swap_traction_sign)

        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

        # # Calculate the total reaction for rigid-body solver
        #reference_point = [0.0, 0.0, 0.0]
        # reaction_force, reaction_moment = KM.ForceAndTorqueUtils.ComputeEquivalentForceAndTorque(
        #     self._interface_model_part,
        #     reference_point,
        #     KM.REACTION,
        #     KM.REACTION_MOMENT
        # )

        # # Save the forces in the model part for rigid-body solver
        # self._interface_model_part[KMC.RESULTANT_FORCE] = reaction_force
        # self._interface_model_part[KMC.RESULTANT_MOMENT] = reaction_moment

    def __GetEmbeddedSkinUtility(self):
        if not hasattr(self, '_embedded_skin_utility'):
            self.__CreateEmbeddedSkinUtility()
        return self._embedded_skin_utility

    def __CreateEmbeddedSkinUtility(self):
        if self.domain_size == 2:
            self._embedded_skin_utility = KM.EmbeddedSkinUtility2D(
                self.__GetFluidComputingModelPart(),
                self.__GetEmbedddedSkinUtilityModelPart(),
                self.level_set_type)
        elif self.domain_size == 3:
             self._embedded_skin_utility = KM.EmbeddedSkinUtility3D(
                self.__GetFluidComputingModelPart(),
                self.__GetEmbedddedSkinUtilityModelPart(),
                self.level_set_type)
        else:
            raise Exception("Domain size expected to be 2 or 3. Got " + str(self._GetDomainSize()))

    def __GetEmbedddedSkinUtilityModelPart(self):
        if not hasattr(self, '_embedded_skin_utility_model_part'):
            self.__CreateEmbeddedSkinUtilityModelPart()
        return self._embedded_skin_utility_model_part

    def __CreateEmbeddedSkinUtilityModelPart(self):
        skin_model_part =self.model.CreateModelPart("EmbeddedSkinUtilityModelPart")
        skin_model_part.AddNodalSolutionStepVariable(KM.NORMAL)
        if (self.level_set_type == "continuous"):
            skin_model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        elif (self.level_set_type == "discontinuous"):
            skin_model_part.AddNodalSolutionStepVariable(KM.POSITIVE_FACE_PRESSURE)
            skin_model_part.AddNodalSolutionStepVariable(KM.NEGATIVE_FACE_PRESSURE)
        else:
            err_msg = 'Level set type is: \'' + self.level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)
        self._embedded_skin_utility_model_part = skin_model_part


# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

import sys
from pathlib import Path

# Test files live in subdirectories mirroring python_scripts/ where a package
# earned one. Discovering them beats listing them: a new subdirectory whose
# tests were never imported would look like a passing run.
for _directory in sorted(Path(__file__).parent.rglob("*")):
    if _directory.is_dir() and _directory.name != "__pycache__" \
            and any(_directory.glob("test_*.py")):
        sys.path.insert(0, str(_directory))

from test_torch_bridge import TestTorchBridge
from test_torch_bridge import TestTorchBridgeWithoutTorch
from test_dataset_export_process import TestDatasetExportProcess
from test_tessellation import TestTessellationTables
from test_tessellation import TestTessellationVolumeConservation
from test_tessellation import TestDompierreTessellation
from test_tessellation import TestHigherOrderSubdivision
from test_tessellation import TestVectorizedProvenanceFastPath
from test_curved_tessellation import TestCurvedShapeFunctions
from test_curved_tessellation import TestCurvedLevelOneEquivalence
from test_curved_tessellation import TestCurvedMeasureConvergence
from test_curved_tessellation import TestCurvedWatertightness
from test_curved_tessellation import TestCurvedProvenanceAndGather
from test_curved_tessellation import TestCurvedValidation
from test_provenance_roundtrip import TestProvenanceRoundTrip
from test_physicsnemo_mesh_construction import TestPhysicsNemoMeshConstruction
from test_domain_mesh import TestDomainMesh
from test_model_registry import TestModelRegistry
from test_model_registry import TestModelCards
from test_model_registry import TestOutputNormalization
from test_model_registry import TestOutputNormalizationTensorPath
from test_model_registry import TestInputNormalization
from test_domain_parallel_utils import TestDomainParallelUtilsSerial
from test_model_registry import TestLoadOutputNormalization
from test_training_utils import TestTrainModel
from test_training_utils import TestSaveTrainedModel
from test_rollout_utils import TestEvaluateRollout
from test_inference_process import TestInferenceProcess
from test_inference_process import TestOutputNormalizationThroughProcess
from test_hybrid_initialization_process import TestHybridInitializationProcess
from test_torch_dataset import TestTorchDataset
from test_mesh_curation import TestMeshAugmentations
from test_mesh_curation import TestMeshDatasetCuration
from test_import_contract import TestImportContract
from test_validation_metrics_process import TestValidationMetricsProcess
from test_mesh_export_process import TestMeshExportProcess
from test_distributed_utils import TestConsistencyCheck
from test_distributed_utils import TestInitializeDistributedManager
from test_grid_bridge import TestGridBridge
from test_grid_bridge import TestGridDerivatives
from test_calculus_bridge import TestCalculusBridge
from test_adaptive_remeshing import TestTargetSizeField
from test_adaptive_remeshing import TestMmgAdaptation
from test_adaptive_remeshing import TestAdaptiveRemeshProcess
from test_adaptive_remeshing import TestSurfacePartition
from test_adaptive_remeshing import TestRemeshSurface
from test_graph_bridge import TestGraphBridge
from test_graph_inference_process import TestGraphInferenceProcess
from test_graph_inference_process import TestGraphBridgePyGConversion
from test_graph_scalable_gnns import TestBistrideHierarchy
from test_graph_scalable_gnns import TestWorldEdges
from test_graph_scalable_gnns import TestScalableGnnsThroughProcess
from test_graph_scalable_gnns import TestScalableGnnOnRealFluidSolve
from test_time_series_inference_process import TestTimeSeriesInferenceProcess
from test_real_solver_integration import TestRealSolverDatasetExport
from test_real_solver_integration import TestRealSolverActiveLearning
from test_real_solver_integration import TestRealSolverSurrogateValidation
from test_superresolution_process import TestSuperResolutionProcess
from test_superresolution_process import TestSRResNetThroughProcess
from test_sample_io import TestApplyParameterOverrides
from test_sample_io import TestLoadFieldsFromNpzDirectory
from test_execution_backends import TestInProcessBackend
from test_execution_backends import TestSubprocessBackend
from test_kratos_label_strategy import TestKratosLabelStrategy
from test_active_learning_driver_integration import TestActiveLearningDriverIntegration
from test_query_strategies import TestQueryStrategySelection
from test_query_strategies import TestEnsembleDisagreementStrategy
from test_query_strategies import TestEntropyStrategy
from test_query_strategies import TestSolverResidualStrategy
from test_metrology import TestValidationMetricsMetrology
from test_parallel_labeling import TestRunCasesBatch
from test_parallel_labeling import TestParallelSubprocessBackend
from test_solver_residuals import TestResidualEvaluator
from test_training_utils import TestTrainModelCallbacks
from test_training_utils import TestSaveShardedModel
from test_streaming_dataset import TestLiveSampleQueue
from test_streaming_dataset import TestStreamingDataset
from test_streaming_dataset import TestStreamingMatchesTheFileRoundTrip
from test_streaming_dataset import TestWarmRestart
from test_grid_dataset_export_process import TestGridDatasetExportProcess
from test_sequence_inference_process import TestGridSequenceDataset
from test_sequence_inference_process import TestSequenceInferenceProcess
from test_sequence_inference_process import TestOne2ManyRNNThroughProcess
from test_point_cloud_inference_process import TestPointCloudInferenceProcess
from test_point_cloud_inference_process import TestTransolverThroughProcess
from test_diffusion_bridge import TestGridPairDataset
from test_diffusion_bridge import TestTrainDiffusionModel
from test_diffusion_bridge import TestDiffusionInferenceProcess
from test_structural_solver_integration import TestStructuralSolverCase
from test_structural_solver_integration import TestStructuralSolverSurrogate
from test_mapping_bridge import TestBackgroundGrid
from test_mapping_bridge import TestMappingBridge
from test_distributed_utils import TestGeometryNameMap
from test_distributed_utils import TestSerialGatherModelPartPassThrough
from test_distributed_utils import TestProcessGroupHelpers
from test_cae_dataset_export_process import TestCaeDatasetExportProcess
from test_cae_dataset_export_process import TestCaeDatasetThroughDatapipes
from test_rom_bridge import TestRomBridge
from test_rom_bridge import TestRowsOfIds
from test_rom_surrogate_process import TestRomSurrogateProcess
from test_rom_temporal import TestRomTrajectoryDataset
from test_rom_temporal import TestSequenceModelTraining
from test_point_cloud_inference_process import TestFIGConvUNetThroughProcess
from test_rom_real_integration import TestRomApplicationBasisInterop
from test_rom_real_integration import TestRomSurrogateOnRealSolves
from test_cfd_bridge import TestCfdBridgePolyData
from test_cfd_bridge import TestCfdMetrics
from test_cfd_bridge import TestHybridInitializationDelegation
from test_cfd_bridge import TestValidationProcessCfdMetrics
from test_cfd_bridge import TestCfdEvaluationWrappers
from test_curator_bridge import TestCuratorSource
from test_curator_bridge import TestCuratorPipeline
from test_curator_bridge import TestCuratorExportProcess
from test_usd_export import TestUsdExportCore
from test_usd_export import TestUsdPointsKind
from test_usd_export import TestUsdExportProcess
from test_nim_client import TestNimClient
from test_nim_client import TestNimInferenceProcess
from test_nim_client import TestNimLive
from test_model_registry import TestTorchCompile
from test_model_registry import TestNvtxUtils
from test_point_cloud_inference_process import TestGeoTransolverThroughProcess
from test_point_cloud_inference_process import TestCoordinateNormalizationOffTheUnitCube
from test_onnx_inference_process import TestExportOnnxModelValidation
from test_onnx_inference_process import TestOnnxExportAndSession
from test_onnx_inference_process import TestOnnxInferenceProcess
from test_onnx_inference_process import TestOnnxDeviceParsing
from test_onnx_inference_process import TestOnnxCudaFallbackIsReported
from test_onnx_inference_process import TestOnnxGpuInference
from test_triton_export import TestTritonConfigGeneration
from test_triton_export import TestTritonRepositoryExport
from test_triton_export import TestTritonInferenceProcess
from test_cosim_surrogate_wrapper import TestCoSimSurrogateWrapper
from test_cosim_surrogate_wrapper import TestCoSimSurrogateCoupledLoop
from test_ood_guard_utils import TestOODGuardUtils
from test_ood_guard_utils import TestOODGuardEndToEnd
from test_uncertainty_utils import TestMonteCarloPredict
from test_uncertainty_utils import TestConcreteDropoutTraining
from test_uncertainty_utils import TestProcessUncertainty
from test_gp_uq_head import TestGpHead
from test_gp_uq_head import TestCalibrationMetrics
from test_gp_uq_head import TestCalibrationThroughValidationProcess
from test_gp_uq_head import TestGpUncertaintyThroughProcess
from test_validation_metrics_process import TestEnsembleMetrics
from test_ensemble_metrics import TestEnsembleMetricValues
from test_ensemble_metrics import TestEnsembleComparisonsThroughProcess
from test_ensemble_metrics import TestRetainEnsemble
from test_superresolution_process import TestGridOperators2D
from test_superresolution_process import TestGridOperatorZoo
from test_sequence_inference_process import TestWindowAsTimeAxis
from test_domino_inference_process import TestDominoSurfaceScatter
from test_domino_inference_process import TestDominoThroughProcess
from test_domino_inference_process import TestDominoDenormalization
from test_domino_inference_process import TestDominoRealCheckpointDenormalization
from test_domino_finetune import TestDominoCorrector
from test_domino_finetune import TestDominoLora
from test_domino_finetune import TestDominoLoraOnRealCheckpoint
from test_diffusion_bridge import TestDitDenoiser
from test_diffusion_bridge import TestUnet3dDenoiser
from test_particle_bridge import TestParticleBridge
from test_array_backend_utils import TestArrayBackendSelection
from test_array_backend_utils import TestArrayBackendCuPy
from test_cupy_backend_parity import TestCuPyBackendParity
from test_cupy_backend_parity import TestCuPyBackendStaysOptIn
from test_particle_bridge import TestParticleTrajectoryDataset
from test_particle_inference_process import TestParticleInferenceProcess
from test_particle_inference_process import TestMeshGraphNetParticles
from test_vfgn_bridge import TestNormalizationStats
from test_vfgn_bridge import TestVfgnSimulator
from test_vfgn_bridge import TestUpstreamForwardIsBroken
from test_sintering_surrogate import TestSinteringCase
from test_sintering_surrogate import TestVfgnOnRealSinteringTrajectories
from test_graph_partition_utils import TestHaloSubgraphSerial
from test_physics_informed import TestPhysicsInformed
from test_pinn_solve_process import TestPinnSolveProcess
from test_pinn_solve_process import TestPinnOnANonUnitDomain
from test_differentiable_residual import TestDifferentiableResidual
from test_differentiable_residual import TestDifferentiableResidualStructural
from test_differentiable_residual import TestTransientResidualElementIntegrated
from test_differentiable_residual import TestTransientResidualDynamicScheme
from test_temporal_training import TestTrajectoryWindowDataset
from test_temporal_training import TestAutoregressiveTraining
from test_temporal_training import TestTransientStructuralSurrogate
from test_sensitivity_utils import TestSensitivityUtils
from test_sensitivity_utils import TestShapeSensitivityFieldOnAnotherPhysics
from test_mesh_deformation import TestDeformPoints
from test_mesh_deformation import TestRegularizationEnergy
from test_mesh_deformation import TestCoordinateWriteBack
from test_mesh_deformation import TestShapeSensitivities
from test_mesh_deformation import TestMeshMovingInteriorSmoothing
from test_mesh_deformation import TestExactControlSensitivities
from test_vertex_morphing_comparison import TestVertexMorphingComparison
from test_adjoint_cross_validation import TestAdjointCrossValidation
from test_adjoint_bridge import TestAdjointBridgeConversion
from test_adjoint_bridge import TestAdjointBridgeFactory
from test_adjoint_bridge import TestAdjointBridgeFields
from test_adjoint_bridge import TestAdjointObjectiveWeights
from test_adjoint_integration import TestAdjointBridgeStructuralMechanics
from test_adjoint_integration import TestAdjointBridgeConvectionDiffusion
from test_adjoint_integration import TestAdjointSensitivityProcess
from test_sobolev_training import TestSensitivityGradient
from test_sobolev_training import TestSensitivityLossTerm
from test_sobolev_training import TestSobolevTrainingImprovesGradients
from test_surrogate_response_function import TestSurrogateResponseFunction
from test_surrogate_response_function import TestSurrogateResponseFunctionExactMode
from test_provenance_cache import TestProvenanceCache
from test_suite_registration import TestSuiteRegistration
from test_suite_registration import TestDocumentedIdentifiersExist
from test_suite_registration import TestDocumentedImportPathsResolve
from test_suite_registration import TestCrossModuleAttributesExist
from test_suite_registration import TestBenchmarkStillRuns
from test_diffusion_protocol import TestWrapDiffusionModel
from test_diffusion_protocol import TestTrainDiffusionModelProtocol
from test_diffusion_protocol import TestGenerateEnsembleProtocol
from test_diffusion_protocol import TestDpsGuidance
from test_diffusion_protocol import TestMultiDiffusionPatching
from test_diffusion_protocol import TestKratosResidualGuidance
from test_diffusion_protocol import TestDiffusionInferenceProcessGuidance
from test_diffusion_protocol import TestResidualGuidanceThroughTheProcess
from test_grid_bridge import TestFlatBoundingBoxInterpolation
from test_grid_bridge import TestTorchInterpolation
from test_point_subsampling import TestSelectPointSubset
from test_point_subsampling import TestExpandToAllPoints
from test_point_subsampling import TestSubsamplingThroughTheProcess
from test_geometry_guard import TestGeometryGuardUtils
from test_geometry_guard import TestGeometryGuardCheck
from test_operator_interfaces import TestDeepONetContract
from test_operator_interfaces import TestDeepONetThroughTheProcess
from test_operator_interfaces import TestGlobeBridge
from test_operator_interfaces import TestGlobeOnAKratosModelPart
from test_mesh_operations import TestMeshOperations
from test_mesh_operations import TestCurvatureOnAModelPart
from test_mesh_operations import TestExtrusion
from test_mesh_operations import TestMeshDeformers
from test_mesh_operations import TestBvhSampling
from test_mesh_io_paths import TestZarrMeshIo
from test_mesh_io_paths import TestCuratorExportProcessZarrBackend
from test_mesh_io_paths import TestVtkBridge
from test_validation_distribution_metrics import TestDistributionMetrics
from test_validation_distribution_metrics import TestSpectralMetrics
from test_validation_distribution_metrics import TestSpectralComparisonsThroughTheProcess
from test_training_performance import TestTrainingPerformanceLayer
from test_training_performance import TestTrainingSchedulers
from test_operator_model_recipes import TestDpotInterface
from test_operator_model_recipes import TestTopoDiffDenoiser
from test_generative_design import TestComplianceCase
from test_generative_design import TestTopoDiffOnComplianceData
from test_aerojepa_pretraining import TestAeroJepaPretraining
from test_pinn_cavity import TestPinnCavity
from test_vortex_shedding import TestCylinderCase
from test_vortex_shedding import TestVortexSheddingSurrogate
from test_contact_surrogate import TestContactCase
from test_contact_surrogate import TestContactFieldsThroughTheBridge
from test_airfoil_surrogate import TestAirfoilCase
from test_airfoil_surrogate import TestAirfoilSurrogate
from test_airfoil_surrogate import TestAirfoilAdjointDispatch
from test_geomechanics_surrogate import TestConsolidationCase
from test_geomechanics_surrogate import TestConsolidationSurrogate
from test_iga_gather import TestIgaAnalysisGather
from test_rans_closure import TestRansCase
from test_rans_closure import TestLearnedClosure
from test_docs_links import TestDocumentationLinks
from test_docs_figures import TestDocsFigures
from test_notebooks import TestNotebooks
from test_mesh_spatial import TestSignedDistance
from test_mesh_spatial import TestGridVectorOperators
from test_mesh_generate import TestImplicitDomainGeneration
from test_mesh_generate import TestRefitToImplicit
from test_mesh_generate import TestSurfaceFromLevelSet
from test_mesh_generate import TestFillBoundaryLoop
from test_mesh_generate import TestFillSurfaceWithTetrahedra
from test_mesh_generate import TestTetgenFill
from test_mesh_generate import TestTetgenMissingPackage
from test_nurbs_sampling import TestNurbsSampling
from test_nurbs_sampling import TestNurbsSamplingWithIgaModeler
from test_nurbs_sampling import TestNurbsMeshBuild
from test_mesh_generate import TestPopulateModelPart
from test_mesh_generate import TestFieldTransferToGeneratedMesh
from test_mesh_generate import TestGeneratedMeshFeedsMmg
from test_graphcast_recipe import TestShallowWaterReference
from test_graphcast_recipe import TestGraphCastShallowWater
from test_lennard_jones_recipe import TestLennardJonesReference
from test_lennard_jones_recipe import TestPeriodicParticleGraph
from test_lennard_jones_recipe import TestLennardJonesMeshGraphNet
from test_corrdiff_recipe import TestCorrDiffTwoStage
from test_corrdiff_recipe import TestRegressionSettingsThroughProcess
from test_corrdiff_recipe import TestFwiInversionRecipe


def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTorchBridge]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestArrayBackendSelection]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestArrayBackendCuPy]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCuPyBackendParity]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCuPyBackendStaysOptIn]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTorchBridgeWithoutTorch]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDatasetExportProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTessellationTables]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTessellationVolumeConservation]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDompierreTessellation]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHigherOrderSubdivision]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestVectorizedProvenanceFastPath]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCurvedShapeFunctions]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCurvedLevelOneEquivalence]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCurvedMeasureConvergence]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCurvedWatertightness]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCurvedProvenanceAndGather]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCurvedValidation]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestProvenanceRoundTrip]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPhysicsNemoMeshConstruction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestApplyParameterOverrides]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLoadFieldsFromNpzDirectory]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestInProcessBackend]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSubprocessBackend]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestKratosLabelStrategy]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestActiveLearningDriverIntegration]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestQueryStrategySelection]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestEnsembleDisagreementStrategy]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestEntropyStrategy]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSolverResidualStrategy]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestValidationMetricsMetrology]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRunCasesBatch]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestParallelSubprocessBackend]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestResidualEvaluator]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTrainModelCallbacks]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSaveShardedModel]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLiveSampleQueue]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestStreamingDataset]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestStreamingMatchesTheFileRoundTrip]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestWarmRestart]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGridDatasetExportProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGridSequenceDataset]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSequenceInferenceProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOne2ManyRNNThroughProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPointCloudInferenceProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTransolverThroughProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGridPairDataset]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTrainDiffusionModel]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDiffusionInferenceProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestStructuralSolverCase]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestStructuralSolverSurrogate]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestBackgroundGrid]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMappingBridge]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGeometryNameMap]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSerialGatherModelPartPassThrough]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestProcessGroupHelpers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCaeDatasetExportProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCaeDatasetThroughDatapipes]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRomBridge]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRowsOfIds]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRomSurrogateProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRomTrajectoryDataset]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSequenceModelTraining]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestFIGConvUNetThroughProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRomApplicationBasisInterop]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRomSurrogateOnRealSolves]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDomainMesh]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestModelRegistry]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestModelCards]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOutputNormalization]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOutputNormalizationTensorPath]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestInputNormalization]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDomainParallelUtilsSerial]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLoadOutputNormalization]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTrainModel]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSaveTrainedModel]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestEvaluateRollout]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestInferenceProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOutputNormalizationThroughProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHybridInitializationProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTorchDataset]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMeshAugmentations]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMeshDatasetCuration]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestImportContract]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestValidationMetricsProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCfdBridgePolyData]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCfdMetrics]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHybridInitializationDelegation]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestValidationProcessCfdMetrics]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCfdEvaluationWrappers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCuratorSource]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCuratorPipeline]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCuratorExportProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestUsdExportCore]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestUsdPointsKind]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestUsdExportProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestNimClient]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestNimInferenceProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestNimLive]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTorchCompile]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestNvtxUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGeoTransolverThroughProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoordinateNormalizationOffTheUnitCube]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestExportOnnxModelValidation]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOnnxExportAndSession]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOnnxInferenceProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOnnxDeviceParsing]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOnnxCudaFallbackIsReported]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOnnxGpuInference]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTritonConfigGeneration]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTritonRepositoryExport]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTritonInferenceProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoSimSurrogateWrapper]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoSimSurrogateCoupledLoop]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOODGuardUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestOODGuardEndToEnd]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMonteCarloPredict]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConcreteDropoutTraining]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestProcessUncertainty]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGpHead]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCalibrationMetrics]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCalibrationThroughValidationProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGpUncertaintyThroughProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestEnsembleMetrics]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestEnsembleMetricValues]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestEnsembleComparisonsThroughProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRetainEnsemble]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGridOperators2D]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGridOperatorZoo]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestWindowAsTimeAxis]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDominoSurfaceScatter]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDominoThroughProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDominoDenormalization]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDominoRealCheckpointDenormalization]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDominoCorrector]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDominoLora]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDominoLoraOnRealCheckpoint]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDitDenoiser]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestUnet3dDenoiser]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestParticleBridge]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestParticleTrajectoryDataset]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestParticleInferenceProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMeshGraphNetParticles]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestNormalizationStats]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestVfgnSimulator]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestUpstreamForwardIsBroken]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSinteringCase]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestVfgnOnRealSinteringTrajectories]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestHaloSubgraphSerial]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPhysicsInformed]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPinnSolveProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPinnOnANonUnitDomain]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDifferentiableResidual]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDifferentiableResidualStructural]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTransientResidualElementIntegrated]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTransientResidualDynamicScheme]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTrajectoryWindowDataset]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAutoregressiveTraining]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTransientStructuralSurrogate]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSensitivityUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestShapeSensitivityFieldOnAnotherPhysics]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDeformPoints]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRegularizationEnergy]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCoordinateWriteBack]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestShapeSensitivities]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMeshMovingInteriorSmoothing]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestExactControlSensitivities]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestVertexMorphingComparison]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointCrossValidation]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointBridgeConversion]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointBridgeFactory]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointBridgeFields]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointObjectiveWeights]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointBridgeStructuralMechanics]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointBridgeConvectionDiffusion]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdjointSensitivityProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSensitivityGradient]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSensitivityLossTerm]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSobolevTrainingImprovesGradients]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSurrogateResponseFunction]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSurrogateResponseFunctionExactMode]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSignedDistance]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGridVectorOperators]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestImplicitDomainGeneration]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRefitToImplicit]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSurfaceFromLevelSet]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestFillBoundaryLoop]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestFillSurfaceWithTetrahedra]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTetgenFill]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTetgenMissingPackage]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestNurbsSampling]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestNurbsSamplingWithIgaModeler]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestNurbsMeshBuild]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPopulateModelPart]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestFieldTransferToGeneratedMesh]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGeneratedMeshFeedsMmg]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestShallowWaterReference]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGraphCastShallowWater]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLennardJonesReference]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPeriodicParticleGraph]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLennardJonesMeshGraphNet]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCorrDiffTwoStage]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRegressionSettingsThroughProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestFwiInversionRecipe]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMeshExportProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConsistencyCheck]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestInitializeDistributedManager]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGridBridge]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGridDerivatives]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCalculusBridge]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTargetSizeField]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMmgAdaptation]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAdaptiveRemeshProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSurfacePartition]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRemeshSurface]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSuperResolutionProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSRResNetThroughProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGraphBridge]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGraphInferenceProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGraphBridgePyGConversion]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestBistrideHierarchy]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestWorldEdges]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestScalableGnnsThroughProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestScalableGnnOnRealFluidSolve]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTimeSeriesInferenceProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRealSolverDatasetExport]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRealSolverActiveLearning]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRealSolverSurrogateValidation]))

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSuiteRegistration]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestProvenanceCache]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDocumentedIdentifiersExist]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDocumentedImportPathsResolve]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCrossModuleAttributesExist]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestBenchmarkStillRuns]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestWrapDiffusionModel]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTrainDiffusionModelProtocol]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGenerateEnsembleProtocol]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDpsGuidance]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMultiDiffusionPatching]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestKratosResidualGuidance]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDiffusionInferenceProcessGuidance]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestResidualGuidanceThroughTheProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestFlatBoundingBoxInterpolation]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTorchInterpolation]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSelectPointSubset]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestExpandToAllPoints]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSubsamplingThroughTheProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGeometryGuardUtils]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGeometryGuardCheck]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDeepONetContract]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDeepONetThroughTheProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGlobeBridge]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestGlobeOnAKratosModelPart]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMeshOperations]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCurvatureOnAModelPart]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestExtrusion]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestMeshDeformers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestBvhSampling]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestZarrMeshIo]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCuratorExportProcessZarrBackend]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestVtkBridge]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDistributionMetrics]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSpectralMetrics]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestSpectralComparisonsThroughTheProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTrainingPerformanceLayer]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTrainingSchedulers]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDpotInterface]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTopoDiffDenoiser]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestComplianceCase]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestTopoDiffOnComplianceData]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAeroJepaPretraining]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestPinnCavity]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestCylinderCase]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestVortexSheddingSurrogate]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestContactCase]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestContactFieldsThroughTheBridge]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAirfoilCase]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAirfoilSurrogate]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestAirfoilAdjointDispatch]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConsolidationCase]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestConsolidationSurrogate]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestIgaAnalysisGather]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestRansCase]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestLearnedClosure]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDocumentationLinks]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestDocsFigures]))

    nightSuite = suites['nightly'] # These tests are executed in the nightly build
    nightSuite.addTests(smallSuite)

    validationSuite = suites['validation'] # These tests are very long and should not be in nightly, for validation
    validationSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TestNotebooks]))

    allSuite = suites['all'] # Create a test suite that contains all the tests
    allSuite.addTests(nightSuite) # Already contains the smallSuite
    allSuite.addTests(validationSuite)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())

# Application dependent names and paths
import KratosMultiphysics
from KratosMultiphysics import _ImportApplication, python_registry_utilities
import KratosMultiphysics.StructuralMechanicsApplication
import KratosGeoMechanicsApplication as _KratosGeoMechanicsApplication

# Explicitly export required bindings by name.
_required_bindings = (
    "KratosGeoMechanicsApplication",
    "KratosExecuteCallBackFunctions",
    "KratosExecuteCriticalHeadInfo",

    "BackwardEulerQuasistaticPwScheme",
    "GeneralizedNewmarkTScheme",
    "GeoLoadSteppingScheme",
    "GeoStaticScheme",
    "GeoMechanicsNewtonRaphsonStrategy",
    "GeoMechanicsNewtonRaphsonErosionProcessStrategy",
    "NewmarkDynamicUPwScheme",
    "NewmarkQuasistaticDampedUPwScheme",
    "NewmarkQuasistaticPwScheme",
    "NewmarkQuasistaticUPwScheme",
    "CustomWorkflowFactory",
    "DeactivateConditionsOnInactiveElements",
    "ProcessUtilities",
    "ResidualBasedBlockBuilderAndSolverWithMassAndDamping",
    "NodeUtilities",

    "ApplyCPhiReductionProcess",
    "ApplyFinalStressesOfPreviousStageToInitialState",
    "ApplyExcavationProcess",
    "ApplyInitialUniformStressField",
    "ApplyK0ProcedureProcess",
    "ApplyNormalLoadTableProcess",
    "ApplyScalarConstraintTableProcess",
    "ApplyVectorConstraintTableProcess",
    "BackwardEulerQuasistaticUPwScheme",
    "BackwardEulerTScheme",
    "CalculateIncrementalMotionProcess",
    "CalculateTotalMotionProcess",
    "FindNeighbourElementsOfConditionsProcess",
    "FindNeighboursOfInterfacesProcess",
    "GeoExtrapolateIntegrationPointValuesToNodesProcess",
    "SetAbsorbingBoundaryParametersProcess",
    "SetMultipleMovingLoadsProcess",
    "SetParameterFieldProcess",

    "AIR_HUMIDITY",
    "AIR_TEMPERATURE",
    "DT_PRESSURE_COEFFICIENT",
    "DT_TEMPERATURE",
    "DT_TEMPERATURE_COEFFICIENT",
    "DT_WATER_PRESSURE",
    "FLUID_FLUX_VECTOR",
    "GEO_PLASTICITY_STATUS",
    "GEO_RELATIVE_DISPLACEMENT_VECTOR",
    "GEO_EFFECTIVE_TRACTION_VECTOR",
    "HYDRAULIC_DISCHARGE",
    "HYDRAULIC_HEAD",
    "INCREMENTAL_DISPLACEMENT",
    "INCREMENTAL_ROTATION",
    "IS_CONVERGED",
    "LOCAL_STRESS_VECTOR",
    "NORMAL_HEAT_FLUX",
    "NORMAL_FLUID_FLUX",
    "PIPE_ACTIVE",
    "PIPE_ELEMENT_LENGTH",
    "PRECIPITATION",
    "SOLAR_RADIATION",
    "TIME_UNIT_CONVERTER",
    "TOTAL_DISPLACEMENT",
    "TOTAL_DISPLACEMENT_X",
    "TOTAL_DISPLACEMENT_Y",
    "TOTAL_DISPLACEMENT_Z",
    "TOTAL_ROTATION",
    "TOTAL_ROTATION_X",
    "TOTAL_ROTATION_Y",
    "TOTAL_ROTATION_Z",
    "TOTAL_STRESS_TENSOR",
    "VELOCITY_COEFFICIENT",
    "WIND_SPEED",

    "UMAT_PARAMETERS",
)
for _name in _required_bindings:
    globals()[_name] = getattr(_KratosGeoMechanicsApplication, _name)

KratosGeoMechanicsApplication = _KratosGeoMechanicsApplication.KratosGeoMechanicsApplication

application = KratosGeoMechanicsApplication()
application_name = "KratosGeoMechanicsApplication"

_ImportApplication(application, application_name)

if not KratosMultiphysics.Registry.HasItem("Stages.KratosMultiphysics.GeoMechanicsApplication.GeoMechanicsAnalysis"):
    from . import python_registry_lists
    python_registry_utilities.RegisterAll("KratosMultiphysics.GeoMechanicsApplication", python_registry_lists)
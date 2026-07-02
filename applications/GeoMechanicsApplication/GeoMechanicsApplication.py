# Application dependent names and paths
import KratosMultiphysics
from KratosMultiphysics import _ImportApplication, python_registry_utilities
import KratosMultiphysics.StructuralMechanicsApplication
import KratosGeoMechanicsApplication as KratosGeo

KratosGeoMechanicsApplication = KratosGeo.KratosGeoMechanicsApplication

# Re-export specific items from KratosGeoMechanicsApplication
ProcessUtilities = KratosGeo.ProcessUtilities
CustomWorkflowFactory = KratosGeo.CustomWorkflowFactory
NodeUtilities = KratosGeo.NodeUtilities
DeactivateConditionsOnInactiveElements = KratosGeo.DeactivateConditionsOnInactiveElements
ResidualBasedBlockBuilderAndSolverWithMassAndDamping = KratosGeo.ResidualBasedBlockBuilderAndSolverWithMassAndDamping
FindNeighbourElementsOfConditionsProcess = KratosGeo.FindNeighbourElementsOfConditionsProcess

# Schemes
BackwardEulerQuasistaticPwScheme = KratosGeo.BackwardEulerQuasistaticPwScheme
BackwardEulerQuasistaticUPwScheme = KratosGeo.BackwardEulerQuasistaticUPwScheme
BackwardEulerTScheme = KratosGeo.BackwardEulerTScheme
GeneralizedNewmarkTScheme = KratosGeo.GeneralizedNewmarkTScheme
GeoLoadSteppingScheme = KratosGeo.GeoLoadSteppingScheme
GeoStaticScheme = KratosGeo.GeoStaticScheme
NewmarkDynamicUPwScheme = KratosGeo.NewmarkDynamicUPwScheme
NewmarkQuasistaticDampedUPwScheme = KratosGeo.NewmarkQuasistaticDampedUPwScheme
NewmarkQuasistaticPwScheme = KratosGeo.NewmarkQuasistaticPwScheme
NewmarkQuasistaticUPwScheme = KratosGeo.NewmarkQuasistaticUPwScheme

# Strategies
GeoMechanicsNewtonRaphsonStrategy = KratosGeo.GeoMechanicsNewtonRaphsonStrategy
GeoMechanicsNewtonRaphsonErosionProcessStrategy = KratosGeo.GeoMechanicsNewtonRaphsonErosionProcessStrategy

# Processes
ApplyCPhiReductionProcess = KratosGeo.ApplyCPhiReductionProcess
ApplyExcavationProcess = KratosGeo.ApplyExcavationProcess
ApplyFinalStressesOfPreviousStageToInitialState = KratosGeo.ApplyFinalStressesOfPreviousStageToInitialState
ApplyInitialUniformStressField = KratosGeo.ApplyInitialUniformStressField
ApplyK0ProcedureProcess = KratosGeo.ApplyK0ProcedureProcess
ApplyNormalLoadTableProcess = KratosGeo.ApplyNormalLoadTableProcess
ApplyScalarConstraintTableProcess = KratosGeo.ApplyScalarConstraintTableProcess
ApplyVectorConstraintTableProcess = KratosGeo.ApplyVectorConstraintTableProcess
CalculateIncrementalMotionProcess = KratosGeo.CalculateIncrementalMotionProcess
CalculateTotalMotionProcess = KratosGeo.CalculateTotalMotionProcess
FindNeighboursOfInterfacesProcess = KratosGeo.FindNeighboursOfInterfacesProcess
GeoExtrapolateIntegrationPointValuesToNodesProcess = KratosGeo.GeoExtrapolateIntegrationPointValuesToNodesProcess
SetAbsorbingBoundaryParametersProcess = KratosGeo.SetAbsorbingBoundaryParametersProcess
SetMultipleMovingLoadsProcess = KratosGeo.SetMultipleMovingLoadsProcess
SetParameterFieldProcess = KratosGeo.SetParameterFieldProcess

# Callback functions
KratosExecuteCallBackFunctions = KratosGeo.KratosExecuteCallBackFunctions
KratosExecuteCriticalHeadInfo = KratosGeo.KratosExecuteCriticalHeadInfo

# Constants - Variables
VELOCITY_COEFFICIENT = KratosGeo.VELOCITY_COEFFICIENT
TOTAL_DISPLACEMENT = KratosGeo.TOTAL_DISPLACEMENT
TOTAL_ROTATION = KratosGeo.TOTAL_ROTATION
INCREMENTAL_DISPLACEMENT = KratosGeo.INCREMENTAL_DISPLACEMENT
INCREMENTAL_ROTATION = KratosGeo.INCREMENTAL_ROTATION
GEO_RELATIVE_DISPLACEMENT_VECTOR = KratosGeo.GEO_RELATIVE_DISPLACEMENT_VECTOR
GEO_EFFECTIVE_TRACTION_VECTOR = KratosGeo.GEO_EFFECTIVE_TRACTION_VECTOR
TIME_UNIT_CONVERTER = KratosGeo.TIME_UNIT_CONVERTER
AIR_HUMIDITY = KratosGeo.AIR_HUMIDITY
AIR_TEMPERATURE = KratosGeo.AIR_TEMPERATURE
DT_PRESSURE_COEFFICIENT = KratosGeo.DT_PRESSURE_COEFFICIENT
DT_TEMPERATURE = KratosGeo.DT_TEMPERATURE
DT_TEMPERATURE_COEFFICIENT = KratosGeo.DT_TEMPERATURE_COEFFICIENT
DT_WATER_PRESSURE = KratosGeo.DT_WATER_PRESSURE
FLUID_FLUX_VECTOR = KratosGeo.FLUID_FLUX_VECTOR
GEO_PLASTICITY_STATUS = KratosGeo.GEO_PLASTICITY_STATUS
HYDRAULIC_DISCHARGE = KratosGeo.HYDRAULIC_DISCHARGE
HYDRAULIC_HEAD = KratosGeo.HYDRAULIC_HEAD
LOCAL_STRESS_VECTOR = KratosGeo.LOCAL_STRESS_VECTOR
NORMAL_FLUID_FLUX = KratosGeo.NORMAL_FLUID_FLUX
NORMAL_HEAT_FLUX = KratosGeo.NORMAL_HEAT_FLUX
PERMEABILITY_MATRIX = KratosGeo.PERMEABILITY_MATRIX
PRECIPITATION = KratosGeo.PRECIPITATION
PIPE_ELEMENT_LENGTH = KratosGeo.PIPE_ELEMENT_LENGTH
PIPE_ACTIVE = KratosGeo.PIPE_ACTIVE
IS_CONVERGED = KratosGeo.IS_CONVERGED
SOLAR_RADIATION = KratosGeo.SOLAR_RADIATION
TOTAL_STRESS_TENSOR = KratosGeo.TOTAL_STRESS_TENSOR
TOTAL_DISPLACEMENT_Y = KratosGeo.TOTAL_DISPLACEMENT_Y
UMAT_PARAMETERS = KratosGeo.UMAT_PARAMETERS
WIND_SPEED = KratosGeo.WIND_SPEED

# Try to import component variables from the base C++ module
_component_vars = ['TOTAL_DISPLACEMENT_X', 'TOTAL_DISPLACEMENT_Z',
                   'TOTAL_ROTATION_X', 'TOTAL_ROTATION_Y', 'TOTAL_ROTATION_Z',
                   'INCREMENTAL_DISPLACEMENT_X', 'INCREMENTAL_DISPLACEMENT_Z',
                   'INCREMENTAL_ROTATION_X', 'INCREMENTAL_ROTATION_Y', 'INCREMENTAL_ROTATION_Z']

for _var in _component_vars:
    if hasattr(KratosGeo, _var):
        globals()[_var] = getattr(KratosGeo, _var)

application = KratosGeoMechanicsApplication()
application_name = "KratosGeoMechanicsApplication"

_ImportApplication(application, application_name)

if not KratosMultiphysics.Registry.HasItem("Stages.KratosMultiphysics.GeoMechanicsApplication.GeoMechanicsAnalysis"):
    from . import python_registry_lists
    python_registry_utilities.RegisterAll("KratosMultiphysics.GeoMechanicsApplication", python_registry_lists)
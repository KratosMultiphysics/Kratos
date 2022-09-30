//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//


// System includes
#include "pybind11/stl.h"

// External includes
#ifdef KRATOS_USE_AMATRIX
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#endif // KRATOS_USE_AMATRIX

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "processes/process.h"
#include "includes/model_part.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "custom_utilities/fluid_auxiliary_utilities.h"
#include "custom_utilities/drag_utilities.h"
#include "custom_utilities/dynamic_smagorinsky_utilities.h"
#include "custom_utilities/estimate_dt_utilities.h"
#include "custom_utilities/fluid_characteristic_numbers_utilities.h"
#include "custom_utilities/fractional_step_settings_periodic.h"
#include "custom_utilities/fractional_step_settings.h"
#include "custom_utilities/integration_point_to_node_transformation_utility.h"
#include "custom_utilities/periodic_condition_utilities.h"
#include "custom_utilities/compressible_element_rotation_utility.h"
#include "custom_utilities/acceleration_limitation_utilities.h"
#include "custom_utilities/fluid_test_utilities.h"
#include "custom_utilities/fluid_adjoint_utilities.h"
#include "custom_utilities/fluid_fft_utilities.h"
#include "custom_utilities/fluid_lss_variable_utilities.h"
#include "custom_utilities/fluid_lss_sensitivity.h"
#include "custom_utilities/fluid_lss_shape_sensitivity.h"
#include "custom_utilities/fluid_model_part_preprocessing_utilities.h"

#include "utilities/split_tetrahedra.h"


namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double> > SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    // Dynamic Smagorinsky utilitites
    py::class_<DynamicSmagorinskyUtils>(m,"DynamicSmagorinskyUtils")
        .def(py::init<ModelPart&,unsigned int>())
        .def("StoreCoarseMesh",&DynamicSmagorinskyUtils::StoreCoarseMesh)
        .def("CalculateC",&DynamicSmagorinskyUtils::CalculateC)
        .def("CorrectFlagValues",&DynamicSmagorinskyUtils::CorrectFlagValues)
        ;

    // Estimate time step utilities
    py::class_<EstimateDtUtility>(m,"EstimateDtUtility")
        .def(py::init< ModelPart&, const double, const double, const double >())
        .def(py::init< ModelPart&, Parameters& >())
        .def("SetCFL",&EstimateDtUtility::SetCFL)
        .def("SetDtMax",&EstimateDtUtility::SetDtMin)
        .def("SetDtMax",&EstimateDtUtility::SetDtMax)
        .def("EstimateDt",&EstimateDtUtility::EstimateDt)
        ;

    // Fluid characteristic numbers utilities
    py::class_<FluidCharacteristicNumbersUtilities>(m,"FluidCharacteristicNumbersUtilities")
        .def_static("CalculateLocalCFL",(void (*)(ModelPart&)) &FluidCharacteristicNumbersUtilities::CalculateLocalCFL)
        ;

    // Periodic boundary conditions utilities
    typedef void (PeriodicConditionUtilities::*AddDoubleVariableType)(Properties&,Variable<double>&);

    AddDoubleVariableType AddDoubleVariable = &PeriodicConditionUtilities::AddPeriodicVariable;

    py::class_<PeriodicConditionUtilities>(m,"PeriodicConditionUtilities")
        .def(py::init<ModelPart&,unsigned int>())
        .def("SetUpSearchStructure",&PeriodicConditionUtilities::SetUpSearchStructure)
        .def("DefinePeriodicBoundary",&PeriodicConditionUtilities::DefinePeriodicBoundary)
        .def("AddPeriodicVariable",AddDoubleVariable)
    ;

    // Base settings
    typedef SolverSettings<SparseSpaceType,LocalSpaceType,LinearSolverType> BaseSettingsType;

    py::class_ < BaseSettingsType >(m, "BaseSettingsType" );

    // Fractional step settings
    py::enum_<FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel>(m,"StrategyLabel")
        .value("Velocity",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::Velocity)
        .value("Pressure",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::Pressure)
        //.value("EddyViscosity",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::EddyViscosity)
    ;

    py::enum_<FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::TurbulenceModelLabel>(m,"TurbulenceModelLabel")
        .value("SpalartAllmaras",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SpalartAllmaras)
    ;

    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*SetStrategyByParamsType)(FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel const&,LinearSolverType::Pointer,const double,const unsigned int);
    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*BuildTurbModelType)(FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::TurbulenceModelLabel const&, LinearSolverType::Pointer, const double, const unsigned int);
    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*PassTurbModelType)(Process::Pointer);
    SetStrategyByParamsType ThisSetStrategyOverload = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetStrategy;
    BuildTurbModelType SetTurbModel_Build = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetTurbulenceModel;
    PassTurbModelType SetTurbModel_Pass = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetTurbulenceModel;

    py::class_< FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>,BaseSettingsType>
        (m,"FractionalStepSettings")
        .def(py::init<ModelPart&,unsigned int,unsigned int,bool,bool,bool>())
        .def("SetStrategy",ThisSetStrategyOverload)
        .def("SetTurbulenceModel",SetTurbModel_Build)
        .def("SetTurbulenceModel",SetTurbModel_Pass)
        .def("GetStrategy",&FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::pGetStrategy)
        .def("SetEchoLevel",&FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetEchoLevel)
    ;

    py::class_< FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>,BaseSettingsType>
        (m,"FractionalStepSettingsPeriodic")
        .def(py::init<ModelPart&,unsigned int,unsigned int,bool,bool,bool,const Kratos::Variable<int>&>())
        .def("SetStrategy",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetStrategy)
        .def("GetStrategy",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::pGetStrategy)
        .def("SetEchoLevel",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetEchoLevel)
    ;

    // Transform from integration point to nodes utilities
    typedef IntegrationPointToNodeTransformationUtility<2,3> IntegrationPointToNodeTransformationUtility2DType;
    typedef IntegrationPointToNodeTransformationUtility<3,4> IntegrationPointToNodeTransformationUtility3DType;
    py::class_<IntegrationPointToNodeTransformationUtility2DType>(m,"IntegrationPointToNodeTransformationUtility2D")
        .def(py::init<>())
        .def("TransformFromIntegrationPointsToNodes",&IntegrationPointToNodeTransformationUtility2DType::TransformFromIntegrationPointsToNodes<double>)
        ;
    py::class_<IntegrationPointToNodeTransformationUtility3DType>(m,"IntegrationPointToNodeTransformationUtility3D")
        .def(py::init<>())
        .def("TransformFromIntegrationPointsToNodes",&IntegrationPointToNodeTransformationUtility3DType::TransformFromIntegrationPointsToNodes<double>)
        ;

    // Calculate embedded drag utilities
    py::class_< DragUtilities> (m,"DragUtilities")
        .def(py::init<>())
        .def("CalculateBodyFittedDrag", &DragUtilities::CalculateBodyFittedDrag)
        .def("CalculateEmbeddedDrag", &DragUtilities::CalculateEmbeddedDrag)
        ;

    py::class_<
        CompressibleElementRotationUtility<LocalSpaceType::MatrixType,LocalSpaceType::VectorType>,
        CompressibleElementRotationUtility<LocalSpaceType::MatrixType,LocalSpaceType::VectorType>::Pointer,
        CoordinateTransformationUtils<LocalSpaceType::MatrixType,LocalSpaceType::VectorType,double> >
        (m,"CompressibleElementRotationUtility")
        .def(py::init<const unsigned int,const Kratos::Flags&>())
        ;

    // Limit the maximal accelearation inside a time step to a physically possible value
    py::class_<
        AccelerationLimitationUtilities>
        (m,"AccelerationLimitationUtilities")
        .def(py::init< ModelPart&, const double >())
        .def("SetLimitAsMultipleOfGravitionalAcceleration", &AccelerationLimitationUtilities::SetLimitAsMultipleOfGravitionalAcceleration)
        .def("Execute", &AccelerationLimitationUtilities::Execute)
        ;

    // Auxiliary utilities
    py::class_<FluidAuxiliaryUtilities>(m, "FluidAuxiliaryUtilities")
        .def_static("CalculateFlowRate", &FluidAuxiliaryUtilities::CalculateFlowRate)
        .def_static("CalculateFlowRatePositiveSkin", [](const ModelPart& rModelPart){return FluidAuxiliaryUtilities::CalculateFlowRatePositiveSkin(rModelPart);})
        .def_static("CalculateFlowRatePositiveSkin", [](const ModelPart& rModelPart, const Flags& rSkinFlag){return FluidAuxiliaryUtilities::CalculateFlowRatePositiveSkin(rModelPart, rSkinFlag);})
        .def_static("CalculateFlowRateNegativeSkin", [](const ModelPart& rModelPart){return FluidAuxiliaryUtilities::CalculateFlowRateNegativeSkin(rModelPart);})
        .def_static("CalculateFlowRateNegativeSkin", [](const ModelPart& rModelPart, const Flags& rSkinFlag){return FluidAuxiliaryUtilities::CalculateFlowRateNegativeSkin(rModelPart, rSkinFlag);})
        .def_static("CalculateFluidVolume", &FluidAuxiliaryUtilities::CalculateFluidVolume)
        .def_static("CalculateFluidPositiveVolume", &FluidAuxiliaryUtilities::CalculateFluidPositiveVolume)
        .def_static("CalculateFluidNegativeVolume", &FluidAuxiliaryUtilities::CalculateFluidNegativeVolume)
        .def_static("MapVelocityFromSkinToVolumeRBF", &FluidAuxiliaryUtilities::MapVelocityFromSkinToVolumeRBF)
        ;

    py::class_<FluidTestUtilities>(m, "FluidTestUtilities")
        .def_static("RandomFillHistoricalVariable", &FluidTestUtilities::RandomFillHistoricalVariable<double>)
        .def_static("RandomFillHistoricalVariable", &FluidTestUtilities::RandomFillHistoricalVariable<array_1d<double, 3>>)
        .def_static("RandomFillNonHistoricalVariable", [](ModelPart::NodesContainerType& rNodesContainer, const Variable<double>& rVariable, const IndexType DomainSize, const double MinValue, const double MaxValue) { FluidTestUtilities::RandomFillNonHistoricalVariable(rNodesContainer, rVariable, DomainSize, MinValue, MaxValue);})
        .def_static("RandomFillNonHistoricalVariable", [](ModelPart::NodesContainerType& rNodesContainer, const Variable<array_1d<double, 3>>& rVariable, const IndexType DomainSize, const double MinValue, const double MaxValue) { FluidTestUtilities::RandomFillNonHistoricalVariable(rNodesContainer, rVariable, DomainSize, MinValue, MaxValue);})
        .def_static("RandomFillNonHistoricalVariable", [](ModelPart::ConditionsContainerType& rConditionsContainer, const Variable<double>& rVariable, const IndexType DomainSize, const double MinValue, const double MaxValue) { FluidTestUtilities::RandomFillNonHistoricalVariable(rConditionsContainer, rVariable, DomainSize, MinValue, MaxValue);})
        .def_static("RandomFillNonHistoricalVariable", [](ModelPart::ConditionsContainerType& rConditionsContainer, const Variable<array_1d<double, 3>>& rVariable, const IndexType DomainSize, const double MinValue, const double MaxValue) { FluidTestUtilities::RandomFillNonHistoricalVariable(rConditionsContainer, rVariable, DomainSize, MinValue, MaxValue);})
        .def_static("RandomFillNonHistoricalVariable", [](ModelPart::ElementsContainerType& rElementsContainer, const Variable<double>& rVariable, const IndexType DomainSize, const double MinValue, const double MaxValue) { FluidTestUtilities::RandomFillNonHistoricalVariable(rElementsContainer, rVariable, DomainSize, MinValue, MaxValue);})
        .def_static("RandomFillNonHistoricalVariable", [](ModelPart::ElementsContainerType& rElementsContainer, const Variable<array_1d<double, 3>>& rVariable, const IndexType DomainSize, const double MinValue, const double MaxValue) { FluidTestUtilities::RandomFillNonHistoricalVariable(rElementsContainer, rVariable, DomainSize, MinValue, MaxValue);})
        ;

    py::class_<FluidAdjointUtilities<2>>(m, "FluidAdjointUtilities2D")
        .def_static("CalculateTriangleAreaDerivative", &FluidAdjointUtilities<2>::CalculateTriangleAreaDerivative)
        ;

    py::class_<FluidAdjointUtilities<3>>(m, "FluidAdjointUtilities3D")
        .def_static("CalculateTriangleAreaDerivative", &FluidAdjointUtilities<3>::CalculateTriangleAreaDerivative)
        ;

    py::class_<FluidFFTUtilities>(m,"FluidFFTUtilities")
        .def(py::init<const double, const double, const double>())
        .def("CalculateFFTFrequencyDistribution",&FluidFFTUtilities::CalculateFFTFrequencyDistribution)
        .def("IsWithinWindowingRange",&FluidFFTUtilities::IsWithinWindowingRange)
        .def("CalculateHannWindowCoefficient",&FluidFFTUtilities::CalculateHannWindowCoefficient)
        .def("CalculateFFTRealCoefficient",&FluidFFTUtilities::CalculateFFTRealCoefficient)
        .def("CalculateFFTImagCoefficient",&FluidFFTUtilities::CalculateFFTImagCoefficient)
        .def("CalculateFFTAmplitudeSquare",&FluidFFTUtilities::CalculateFFTAmplitudeSquare)
        .def("CalculateFFTAmplitudeSquareDerivative",&FluidFFTUtilities::CalculateFFTAmplitudeSquareDerivative)
        .def("GetFrequencyResolution",&FluidFFTUtilities::GetFrequencyResolution)
        .def("GetFrequency",&FluidFFTUtilities::GetFrequency)
        .def("GetMaximumFrequency",&FluidFFTUtilities::GetMaximumFrequency)
        .def("GetTotalNumberOfSteps",&FluidFFTUtilities::GetTotalNumberOfSteps)
        .def("GetNumberOfWindowingSteps",&FluidFFTUtilities::GetNumberOfWindowingSteps)
        ;

    py::class_<FluidLSSVariableUtilities, FluidLSSVariableUtilities::Pointer>(m, "FluidLSSVariableUtilities")
        .def(py::init<const std::vector<const Variable<double>*>&, const std::vector<const Variable<double>*>&, const std::vector<const Variable<double>*>&, const std::vector<const Variable<double>*>&, const std::vector<const Variable<double>*>&, const std::vector<const Variable<double>*>&>())
        .def("GetPrimalValues", &FluidLSSVariableUtilities::GetPrimalValues<ModelPart::ConditionType>)
        .def("GetPrimalFirstDerivativeValues", &FluidLSSVariableUtilities::GetPrimalFirstDerivativeValues<ModelPart::ConditionType>)
        .def("GetAdjointValues", &FluidLSSVariableUtilities::GetAdjointValues<ModelPart::ConditionType>)
        .def("GetAdjointFirstDerivativeValues", &FluidLSSVariableUtilities::GetAdjointFirstDerivativeValues<ModelPart::ConditionType>)
        .def("GetLSSValues", &FluidLSSVariableUtilities::GetLSSValues<ModelPart::ConditionType>)
        .def("GetLSSFirstDerivativeValues", &FluidLSSVariableUtilities::GetLSSFirstDerivativeValues<ModelPart::ConditionType>)
        .def("GetPrimalValues", &FluidLSSVariableUtilities::GetPrimalValues<ModelPart::ElementType>)
        .def("GetPrimalFirstDerivativeValues", &FluidLSSVariableUtilities::GetPrimalFirstDerivativeValues<ModelPart::ElementType>)
        .def("GetAdjointValues", &FluidLSSVariableUtilities::GetAdjointValues<ModelPart::ElementType>)
        .def("GetAdjointFirstDerivativeValues", &FluidLSSVariableUtilities::GetAdjointFirstDerivativeValues<ModelPart::ElementType>)
        .def("GetLSSValues", &FluidLSSVariableUtilities::GetLSSValues<ModelPart::ElementType>)
        .def("GetLSSFirstDerivativeValues", &FluidLSSVariableUtilities::GetLSSFirstDerivativeValues<ModelPart::ElementType>)
        ;

    py::class_<FluidLSSSensitivity, FluidLSSSensitivity::Pointer>(m, "FluidLSSSensitivity")
        .def(py::init<>())
        .def("GetDerivativeVariable", &FluidLSSSensitivity::GetDerivativeVariable)
        .def("CalculateResidualSensitivity", (void(FluidLSSSensitivity::*)(Vector&, ModelPart::ConditionType&, const FluidAdjointSlipUtilities&, const ProcessInfo&))(&FluidLSSSensitivity::CalculateResidualSensitivity))
        .def("CalculateResidualSensitivity", (void(FluidLSSSensitivity::*)(Vector&, ModelPart::ElementType&, const FluidAdjointSlipUtilities&, const ProcessInfo&))(&FluidLSSSensitivity::CalculateResidualSensitivity))
        .def("CalculateResponseSensitivity", (double(FluidLSSSensitivity::*)(ModelPart::ConditionType&, AdjointResponseFunction&,const FluidAdjointSlipUtilities&, const ProcessInfo&))(&FluidLSSSensitivity::CalculateResponseSensitivity))
        .def("CalculateResponseSensitivity", (double(FluidLSSSensitivity::*)(ModelPart::ElementType&, AdjointResponseFunction&, const FluidAdjointSlipUtilities&, const ProcessInfo&))(&FluidLSSSensitivity::CalculateResponseSensitivity))
        ;

    py::class_<FluidLSSShapeSensitivity, FluidLSSShapeSensitivity::Pointer, FluidLSSSensitivity>(m, "FluidLSSShapeSensitivity")
        .def(py::init<Parameters, const std::size_t>())
        ;

    py::class_<FluidModelPartPreProcessingUtilities>(m, "FluidModelPartPreProcessingUtilities")
        .def_static("CreateModelPartForCommenInterface", &FluidModelPartPreProcessingUtilities::CreateModelPartForCommenInterface)
        .def_static("GetElementIdsWithAllNodesOnBoundaries", &FluidModelPartPreProcessingUtilities::GetElementIdsWithAllNodesOnBoundaries)
        .def_static("BreakElements", &FluidModelPartPreProcessingUtilities::BreakElements)
        ;


}

}  // namespace Python.

} // Namespace Kratos

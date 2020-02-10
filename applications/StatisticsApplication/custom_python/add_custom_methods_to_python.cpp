//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "custom_python/add_custom_methods_to_python.h"
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "custom_methods/spatial_methods.h"
#include "custom_methods/temporal_methods.h"
#include "custom_utilities/method_utilities.h"

namespace Kratos
{
namespace Python
{
void AddCustomMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Adding spatial methods
    auto spatial_method_module = m.def_submodule("SpatialMethods");

    using HistoricalSpatialMethods = SpatialMethods::HistoricalSpatialMethods;
    auto spatial_historical_method = spatial_method_module.def_submodule("Historical");
    spatial_historical_method.def_submodule("ValueMethods")
        .def("Sum", &HistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &HistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &HistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &HistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &HistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &HistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        ;
    spatial_historical_method.def_submodule("NormMethods")
        .def("Sum", &HistoricalSpatialMethods::CalculateNormSum<double>)
        .def("Sum", &HistoricalSpatialMethods::CalculateNormSum<array_1d<double, 3>>)
        .def("Sum", &HistoricalSpatialMethods::CalculateNormSum<Vector>)
        .def("Sum", &HistoricalSpatialMethods::CalculateNormSum<Matrix>)
        .def("Mean", &HistoricalSpatialMethods::CalculateNormMean<double>)
        .def("Mean", &HistoricalSpatialMethods::CalculateNormMean<array_1d<double, 3>>)
        .def("Mean", &HistoricalSpatialMethods::CalculateNormMean<Vector>)
        .def("Mean", &HistoricalSpatialMethods::CalculateNormMean<Matrix>)
        .def("Variance", &HistoricalSpatialMethods::CalculateNormVariance<double>)
        .def("Variance", &HistoricalSpatialMethods::CalculateNormVariance<array_1d<double, 3>>)
        .def("Variance", &HistoricalSpatialMethods::CalculateNormVariance<Vector>)
        .def("Variance", &HistoricalSpatialMethods::CalculateNormVariance<Matrix>)
        .def("Min", &HistoricalSpatialMethods::GetMin<double>)
        .def("Min", &HistoricalSpatialMethods::GetMin<array_1d<double, 3>>)
        .def("Min", &HistoricalSpatialMethods::GetMin<Vector>)
        .def("Min", &HistoricalSpatialMethods::GetMin<Matrix>)
        .def("Max", &HistoricalSpatialMethods::GetMax<double>)
        .def("Max", &HistoricalSpatialMethods::GetMax<array_1d<double, 3>>)
        .def("Max", &HistoricalSpatialMethods::GetMax<Vector>)
        .def("Max", &HistoricalSpatialMethods::GetMax<Matrix>)
        ;

    using NodalNonHistoricalSpatialMethods = SpatialMethods::NodalNonHistoricalSpatialMethods;
    auto non_historical_spatial_method_module = spatial_method_module.def_submodule("NonHistorical");
    auto spatial_non_historical_nodal_method = non_historical_spatial_method_module.def_submodule("Nodes");
    spatial_non_historical_nodal_method.def_submodule("ValueMethods")
        .def("Sum", &NodalNonHistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &NodalNonHistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        ;
    spatial_non_historical_nodal_method.def_submodule("NormMethods")
        .def("Sum", &NodalNonHistoricalSpatialMethods::CalculateNormSum<double>)
        .def("Sum", &NodalNonHistoricalSpatialMethods::CalculateNormSum<array_1d<double, 3>>)
        .def("Sum", &NodalNonHistoricalSpatialMethods::CalculateNormSum<Vector>)
        .def("Sum", &NodalNonHistoricalSpatialMethods::CalculateNormSum<Matrix>)
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateNormMean<double>)
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateNormMean<array_1d<double, 3>>)
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateNormMean<Vector>)
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateNormMean<Matrix>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateNormVariance<double>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateNormVariance<array_1d<double, 3>>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateNormVariance<Vector>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateNormVariance<Matrix>)
        .def("Min", &NodalNonHistoricalSpatialMethods::GetMin<double>)
        .def("Min", &NodalNonHistoricalSpatialMethods::GetMin<array_1d<double, 3>>)
        .def("Min", &NodalNonHistoricalSpatialMethods::GetMin<Vector>)
        .def("Min", &NodalNonHistoricalSpatialMethods::GetMin<Matrix>)
        .def("Max", &NodalNonHistoricalSpatialMethods::GetMax<double>)
        .def("Max", &NodalNonHistoricalSpatialMethods::GetMax<array_1d<double, 3>>)
        .def("Max", &NodalNonHistoricalSpatialMethods::GetMax<Vector>)
        .def("Max", &NodalNonHistoricalSpatialMethods::GetMax<Matrix>)
        ;

    using ConditionNonHistoricalSpatialMethods = SpatialMethods::ConditionNonHistoricalSpatialMethods;
    auto spatial_non_historical_condition_method = non_historical_spatial_method_module.def_submodule("Conditions");
    spatial_non_historical_condition_method.def_submodule("ValueMethods")
        .def("Sum", &ConditionNonHistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &ConditionNonHistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        ;
    spatial_non_historical_condition_method.def_submodule("NormMethods")
        .def("Sum", &ConditionNonHistoricalSpatialMethods::CalculateNormSum<double>)
        .def("Sum", &ConditionNonHistoricalSpatialMethods::CalculateNormSum<array_1d<double, 3>>)
        .def("Sum", &ConditionNonHistoricalSpatialMethods::CalculateNormSum<Vector>)
        .def("Sum", &ConditionNonHistoricalSpatialMethods::CalculateNormSum<Matrix>)
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateNormMean<double>)
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateNormMean<array_1d<double, 3>>)
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateNormMean<Vector>)
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateNormMean<Matrix>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateNormVariance<double>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateNormVariance<array_1d<double, 3>>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateNormVariance<Vector>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateNormVariance<Matrix>)
        .def("Min", &ConditionNonHistoricalSpatialMethods::GetMin<double>)
        .def("Min", &ConditionNonHistoricalSpatialMethods::GetMin<array_1d<double, 3>>)
        .def("Min", &ConditionNonHistoricalSpatialMethods::GetMin<Vector>)
        .def("Min", &ConditionNonHistoricalSpatialMethods::GetMin<Matrix>)
        .def("Max", &ConditionNonHistoricalSpatialMethods::GetMax<double>)
        .def("Max", &ConditionNonHistoricalSpatialMethods::GetMax<array_1d<double, 3>>)
        .def("Max", &ConditionNonHistoricalSpatialMethods::GetMax<Vector>)
        .def("Max", &ConditionNonHistoricalSpatialMethods::GetMax<Matrix>)
        ;

    using ElementNonHistoricalSpatialMethods = SpatialMethods::ElementNonHistoricalSpatialMethods;
    auto spatial_non_historical_element_method = non_historical_spatial_method_module.def_submodule("Elements");
    spatial_non_historical_element_method.def_submodule("ValueMethods")
        .def("Sum", &ElementNonHistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &ElementNonHistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        ;
    spatial_non_historical_element_method.def_submodule("NormMethods")
        .def("Sum", &ElementNonHistoricalSpatialMethods::CalculateNormSum<double>)
        .def("Sum", &ElementNonHistoricalSpatialMethods::CalculateNormSum<array_1d<double, 3>>)
        .def("Sum", &ElementNonHistoricalSpatialMethods::CalculateNormSum<Vector>)
        .def("Sum", &ElementNonHistoricalSpatialMethods::CalculateNormSum<Matrix>)
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateNormMean<double>)
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateNormMean<array_1d<double, 3>>)
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateNormMean<Vector>)
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateNormMean<Matrix>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateNormVariance<double>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateNormVariance<array_1d<double, 3>>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateNormVariance<Vector>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateNormVariance<Matrix>)
        .def("Min", &ElementNonHistoricalSpatialMethods::GetMin<double>)
        .def("Min", &ElementNonHistoricalSpatialMethods::GetMin<array_1d<double, 3>>)
        .def("Min", &ElementNonHistoricalSpatialMethods::GetMin<Vector>)
        .def("Min", &ElementNonHistoricalSpatialMethods::GetMin<Matrix>)
        .def("Max", &ElementNonHistoricalSpatialMethods::GetMax<double>)
        .def("Max", &ElementNonHistoricalSpatialMethods::GetMax<array_1d<double, 3>>)
        .def("Max", &ElementNonHistoricalSpatialMethods::GetMax<Vector>)
        .def("Max", &ElementNonHistoricalSpatialMethods::GetMax<Matrix>)
        ;

    // Adding temporal methods

    // adding base temporal method
    auto temporal_method_module = m.def_submodule("TemporalMethods");

    py::class_<TemporalMethods::TemporalMethod, TemporalMethods::TemporalMethod::Pointer>(temporal_method_module,"TemporalMethod")
    .def(py::init<ModelPart&>())
    .def("GetModelPart", &TemporalMethods::TemporalMethod::GetModelPart)
    .def("GetTotalTime", &TemporalMethods::TemporalMethod::GetTotalTime)
    .def("InitializeStatisticsMethod", &TemporalMethods::TemporalMethod::InitializeStatisticsMethod)
    .def("FinalizeStatisticsTimeStep", &TemporalMethods::TemporalMethod::FinalizeStatisticsTimeStep)
    ;

    auto temporal_historical_method = temporal_method_module.def_submodule("HistoricalInput");
    auto temporal_historical_historical_output_method = temporal_historical_method.def_submodule("HistoricalOutput");

    auto temporal_historical_historical_output_value_method = temporal_historical_historical_output_method.def_submodule("ValueMethods");
    using HistoricalInputHistoricalOutputTemporalMethods = TemporalMethods::HistoricalInputHistoricalOutputTemporalMethods;
    py::class_<HistoricalInputHistoricalOutputTemporalMethods::VarianceMethod, HistoricalInputHistoricalOutputTemporalMethods::VarianceMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_value_method, "Variance")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::VarianceMethod::CalculateValueStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::VarianceMethod::CalculateValueStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::VarianceMethod::CalculateValueStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::VarianceMethod::CalculateValueStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::VarianceMethod::InitializeValueStatisticsVariables<double>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::VarianceMethod::InitializeValueStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::VarianceMethod::InitializeValueStatisticsVariables<Vector>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::VarianceMethod::InitializeValueStatisticsVariables<Matrix>)
            ;
    py::class_<HistoricalInputHistoricalOutputTemporalMethods::MeanMethod, HistoricalInputHistoricalOutputTemporalMethods::MeanMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_value_method, "Mean")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MeanMethod::CalculateValueStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MeanMethod::CalculateValueStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MeanMethod::CalculateValueStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MeanMethod::CalculateValueStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::MeanMethod::InitializeValueStatisticsVariables<double>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::MeanMethod::InitializeValueStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::MeanMethod::InitializeValueStatisticsVariables<Vector>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::MeanMethod::InitializeValueStatisticsVariables<Matrix>)
            ;
    py::class_<HistoricalInputHistoricalOutputTemporalMethods::SumMethod, HistoricalInputHistoricalOutputTemporalMethods::SumMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_value_method, "Sum")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::SumMethod::CalculateValueStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::SumMethod::CalculateValueStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::SumMethod::CalculateValueStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::SumMethod::CalculateValueStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::SumMethod::InitializeValueStatisticsVariables<double>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::SumMethod::InitializeValueStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::SumMethod::InitializeValueStatisticsVariables<Vector>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::SumMethod::InitializeValueStatisticsVariables<Matrix>)
            ;

    auto temporal_historical_historical_output_norm_method = temporal_historical_historical_output_method.def_submodule("NormMethods");
    using HistoricalInputNonHistoricalOutputTemporalMethods = TemporalMethods::HistoricalInputNonHistoricalOutputTemporalMethods;
    py::class_<HistoricalInputNonHistoricalOutputTemporalMethods::VarianceMethod, HistoricalInputNonHistoricalOutputTemporalMethods::VarianceMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_norm_method, "Variance")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceMethod::CalculateNormStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceMethod::CalculateNormStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceMethod::CalculateNormStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceMethod::CalculateNormStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceMethod::InitializeNormStatisticsVariables)
            ;
    py::class_<HistoricalInputNonHistoricalOutputTemporalMethods::MeanMethod, HistoricalInputNonHistoricalOutputTemporalMethods::MeanMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_norm_method, "Mean")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanMethod::CalculateNormStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanMethod::CalculateNormStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanMethod::CalculateNormStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanMethod::CalculateNormStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanMethod::InitializeNormStatisticsVariables)
            ;
    py::class_<HistoricalInputNonHistoricalOutputTemporalMethods::SumMethod, HistoricalInputNonHistoricalOutputTemporalMethods::SumMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_norm_method, "Sum")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::SumMethod::CalculateNormStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::SumMethod::CalculateNormStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::SumMethod::CalculateNormStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::SumMethod::CalculateNormStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::SumMethod::InitializeNormStatisticsVariables)
            ;
    py::class_<HistoricalInputNonHistoricalOutputTemporalMethods::MinMethod, HistoricalInputNonHistoricalOutputTemporalMethods::MinMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_norm_method, "Min")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MinMethod::CalculateNormStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MinMethod::CalculateNormStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MinMethod::CalculateNormStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MinMethod::CalculateNormStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::MinMethod::InitializeNormStatisticsVariables)
            ;

    py::class_<HistoricalInputNonHistoricalOutputTemporalMethods::MaxMethod, HistoricalInputNonHistoricalOutputTemporalMethods::MaxMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_norm_method, "Max")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MaxMethod::CalculateNormStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MaxMethod::CalculateNormStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MaxMethod::CalculateNormStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MaxMethod::CalculateNormStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::MaxMethod::InitializeNormStatisticsVariables)
            ;
}

} // namespace Python.
} // Namespace Kratos

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

// Application includes
#include "custom_methods/temporal_methods.h"

// Include base h
#include "custom_python/add_custom_temporal_non_historical_condition_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomTemporalNonHistoricalConditionMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto method_module = m.def_submodule("Conditions");

    auto value_method_module = method_module.def_submodule("ValueMethods");
    using ConditionNonHistoricalTemporalMethods = TemporalMethods::ConditionNonHistoricalTemporalMethods;
    py::class_<ConditionNonHistoricalTemporalMethods::VarianceValueMethod, ConditionNonHistoricalTemporalMethods::VarianceValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            value_method_module, "Variance")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::VarianceValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::VarianceValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::VarianceValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::VarianceValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<Matrix>)
            ;
    py::class_<ConditionNonHistoricalTemporalMethods::MeanValueMethod, ConditionNonHistoricalTemporalMethods::MeanValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            value_method_module, "Mean")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MeanValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MeanValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MeanValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MeanValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<Matrix>)
            ;
    py::class_<ConditionNonHistoricalTemporalMethods::SumValueMethod, ConditionNonHistoricalTemporalMethods::SumValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            value_method_module, "Sum")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::SumValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::SumValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::SumValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::SumValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::SumValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::SumValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::SumValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::SumValueMethod::InitializeStatisticsVariables<Matrix>)
            ;

    auto norm_method_module = method_module.def_submodule("NormMethods");
    using ConditionNonHistoricalTemporalMethods = TemporalMethods::ConditionNonHistoricalTemporalMethods;
    py::class_<ConditionNonHistoricalTemporalMethods::VarianceNormMethod, ConditionNonHistoricalTemporalMethods::VarianceNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Variance")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::VarianceNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::VarianceNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::VarianceNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::VarianceNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::VarianceNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<ConditionNonHistoricalTemporalMethods::MeanNormMethod, ConditionNonHistoricalTemporalMethods::MeanNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Mean")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MeanNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MeanNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MeanNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MeanNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::MeanNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<ConditionNonHistoricalTemporalMethods::SumNormMethod, ConditionNonHistoricalTemporalMethods::SumNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Sum")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::SumNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::SumNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::SumNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::SumNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::SumNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<ConditionNonHistoricalTemporalMethods::MinNormMethod, ConditionNonHistoricalTemporalMethods::MinNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Min")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MinNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MinNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MinNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MinNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::MinNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<ConditionNonHistoricalTemporalMethods::MaxNormMethod, ConditionNonHistoricalTemporalMethods::MaxNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Max")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MaxNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MaxNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MaxNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ConditionNonHistoricalTemporalMethods::MaxNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ConditionNonHistoricalTemporalMethods::MaxNormMethod::InitializeStatisticsVariables)
            ;
}

} // namespace Python.
} // Namespace Kratos

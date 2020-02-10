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
#include "custom_python/add_custom_temporal_non_historical_element_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomTemporalNonHistoricalElementMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto method_module = m.def_submodule("Elements");

    auto value_method_module = method_module.def_submodule("ValueMethods");
    using ElementNonHistoricalTemporalMethods = TemporalMethods::ElementNonHistoricalTemporalMethods;
    py::class_<ElementNonHistoricalTemporalMethods::VarianceValueMethod, ElementNonHistoricalTemporalMethods::VarianceValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            value_method_module, "Variance")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::VarianceValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::VarianceValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::VarianceValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::VarianceValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<Matrix>)
            ;
    py::class_<ElementNonHistoricalTemporalMethods::MeanValueMethod, ElementNonHistoricalTemporalMethods::MeanValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            value_method_module, "Mean")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MeanValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MeanValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MeanValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MeanValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<Matrix>)
            ;
    py::class_<ElementNonHistoricalTemporalMethods::SumValueMethod, ElementNonHistoricalTemporalMethods::SumValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            value_method_module, "Sum")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::SumValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::SumValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::SumValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::SumValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::SumValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::SumValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::SumValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::SumValueMethod::InitializeStatisticsVariables<Matrix>)
            ;

    auto norm_method_module = method_module.def_submodule("NormMethods");
    using ElementNonHistoricalTemporalMethods = TemporalMethods::ElementNonHistoricalTemporalMethods;
    py::class_<ElementNonHistoricalTemporalMethods::VarianceNormMethod, ElementNonHistoricalTemporalMethods::VarianceNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Variance")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::VarianceNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::VarianceNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::VarianceNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::VarianceNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::VarianceNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<ElementNonHistoricalTemporalMethods::MeanNormMethod, ElementNonHistoricalTemporalMethods::MeanNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Mean")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MeanNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MeanNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MeanNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MeanNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::MeanNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<ElementNonHistoricalTemporalMethods::SumNormMethod, ElementNonHistoricalTemporalMethods::SumNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Sum")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::SumNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::SumNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::SumNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::SumNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::SumNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<ElementNonHistoricalTemporalMethods::MinNormMethod, ElementNonHistoricalTemporalMethods::MinNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Min")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MinNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MinNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MinNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MinNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::MinNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<ElementNonHistoricalTemporalMethods::MaxNormMethod, ElementNonHistoricalTemporalMethods::MaxNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Max")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MaxNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MaxNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MaxNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &ElementNonHistoricalTemporalMethods::MaxNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &ElementNonHistoricalTemporalMethods::MaxNormMethod::InitializeStatisticsVariables)
            ;
}

} // namespace Python.
} // Namespace Kratos

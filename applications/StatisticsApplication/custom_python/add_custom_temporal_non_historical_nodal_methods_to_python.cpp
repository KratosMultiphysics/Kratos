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
#include "custom_python/add_custom_temporal_non_historical_nodal_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomTemporalNonHistoricalNodalMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto method_module = m.def_submodule("Nodes");

    auto value_method_module = method_module.def_submodule("ValueMethods");
    using NodalNonHistoricalTemporalMethods = TemporalMethods::NodalNonHistoricalTemporalMethods;
    py::class_<NodalNonHistoricalTemporalMethods::VarianceValueMethod, NodalNonHistoricalTemporalMethods::VarianceValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            value_method_module, "Variance")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::VarianceValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::VarianceValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::VarianceValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::VarianceValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<Matrix>)
            ;
    py::class_<NodalNonHistoricalTemporalMethods::MeanValueMethod, NodalNonHistoricalTemporalMethods::MeanValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            value_method_module, "Mean")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MeanValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MeanValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MeanValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MeanValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<Matrix>)
            ;
    py::class_<NodalNonHistoricalTemporalMethods::SumValueMethod, NodalNonHistoricalTemporalMethods::SumValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            value_method_module, "Sum")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::SumValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::SumValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::SumValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::SumValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::SumValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::SumValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::SumValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::SumValueMethod::InitializeStatisticsVariables<Matrix>)
            ;

    auto norm_method_module = method_module.def_submodule("NormMethods");
    using NodalNonHistoricalTemporalMethods = TemporalMethods::NodalNonHistoricalTemporalMethods;
    py::class_<NodalNonHistoricalTemporalMethods::VarianceNormMethod, NodalNonHistoricalTemporalMethods::VarianceNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Variance")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::VarianceNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::VarianceNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::VarianceNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::VarianceNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::VarianceNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<NodalNonHistoricalTemporalMethods::MeanNormMethod, NodalNonHistoricalTemporalMethods::MeanNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Mean")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MeanNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MeanNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MeanNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MeanNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::MeanNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<NodalNonHistoricalTemporalMethods::SumNormMethod, NodalNonHistoricalTemporalMethods::SumNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Sum")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::SumNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::SumNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::SumNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::SumNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::SumNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<NodalNonHistoricalTemporalMethods::MinNormMethod, NodalNonHistoricalTemporalMethods::MinNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Min")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MinNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MinNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MinNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MinNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::MinNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<NodalNonHistoricalTemporalMethods::MaxNormMethod, NodalNonHistoricalTemporalMethods::MaxNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            norm_method_module, "Max")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MaxNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MaxNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MaxNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &NodalNonHistoricalTemporalMethods::MaxNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &NodalNonHistoricalTemporalMethods::MaxNormMethod::InitializeStatisticsVariables)
            ;
}

} // namespace Python.
} // Namespace Kratos

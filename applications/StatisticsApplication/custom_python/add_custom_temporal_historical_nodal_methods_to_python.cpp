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
#include "custom_python/add_custom_temporal_historical_nodal_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomTemporalHistoricalNodalMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto temporal_historical_historical_output_method = m.def_submodule("HistoricalOutput");

    auto temporal_historical_historical_output_value_method = temporal_historical_historical_output_method.def_submodule("ValueMethods");
    using HistoricalInputHistoricalOutputTemporalMethods = TemporalMethods::HistoricalInputHistoricalOutputTemporalMethods;
    py::class_<HistoricalInputHistoricalOutputTemporalMethods::VarianceValueMethod, HistoricalInputHistoricalOutputTemporalMethods::VarianceValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_value_method, "Variance")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::VarianceValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::VarianceValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::VarianceValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::VarianceValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<Matrix>)
            ;
    py::class_<HistoricalInputHistoricalOutputTemporalMethods::MeanValueMethod, HistoricalInputHistoricalOutputTemporalMethods::MeanValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_value_method, "Mean")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MeanValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MeanValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MeanValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MeanValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<Matrix>)
            ;
    py::class_<HistoricalInputHistoricalOutputTemporalMethods::SumValueMethod, HistoricalInputHistoricalOutputTemporalMethods::SumValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_value_method, "Sum")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::SumValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::SumValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::SumValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::SumValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::SumValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::SumValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::SumValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::SumValueMethod::InitializeStatisticsVariables<Matrix>)
            ;

    auto temporal_historical_historical_output_norm_method = temporal_historical_historical_output_method.def_submodule("NormMethods");
    using HistoricalInputHistoricalOutputTemporalMethods = TemporalMethods::HistoricalInputHistoricalOutputTemporalMethods;
    py::class_<HistoricalInputHistoricalOutputTemporalMethods::VarianceNormMethod, HistoricalInputHistoricalOutputTemporalMethods::VarianceNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_norm_method, "Variance")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::VarianceNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::VarianceNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::VarianceNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::VarianceNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::VarianceNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<HistoricalInputHistoricalOutputTemporalMethods::MeanNormMethod, HistoricalInputHistoricalOutputTemporalMethods::MeanNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_norm_method, "Mean")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MeanNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MeanNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MeanNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MeanNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::MeanNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<HistoricalInputHistoricalOutputTemporalMethods::SumNormMethod, HistoricalInputHistoricalOutputTemporalMethods::SumNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_norm_method, "Sum")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::SumNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::SumNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::SumNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::SumNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::SumNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<HistoricalInputHistoricalOutputTemporalMethods::MinNormMethod, HistoricalInputHistoricalOutputTemporalMethods::MinNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_norm_method, "Min")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MinNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MinNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MinNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MinNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::MinNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<HistoricalInputHistoricalOutputTemporalMethods::MaxNormMethod, HistoricalInputHistoricalOutputTemporalMethods::MaxNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_historical_output_norm_method, "Max")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MaxNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MaxNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MaxNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputHistoricalOutputTemporalMethods::MaxNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputHistoricalOutputTemporalMethods::MaxNormMethod::InitializeStatisticsVariables)
            ;

    auto temporal_historical_non_historical_output_method = m.def_submodule("NonHistoricalOutput");

    auto temporal_historical_non_historical_output_value_method = temporal_historical_non_historical_output_method.def_submodule("ValueMethods");
    using HistoricalInputNonHistoricalOutputTemporalMethods = TemporalMethods::HistoricalInputNonHistoricalOutputTemporalMethods;
    py::class_<HistoricalInputNonHistoricalOutputTemporalMethods::VarianceValueMethod, HistoricalInputNonHistoricalOutputTemporalMethods::VarianceValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_non_historical_output_value_method, "Variance")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceValueMethod::InitializeStatisticsVariables<Matrix>)
            ;
    py::class_<HistoricalInputNonHistoricalOutputTemporalMethods::MeanValueMethod, HistoricalInputNonHistoricalOutputTemporalMethods::MeanValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_non_historical_output_value_method, "Mean")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanValueMethod::InitializeStatisticsVariables<Matrix>)
            ;
    py::class_<HistoricalInputNonHistoricalOutputTemporalMethods::SumValueMethod, HistoricalInputNonHistoricalOutputTemporalMethods::SumValueMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_non_historical_output_value_method, "Sum")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::SumValueMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::SumValueMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::SumValueMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::SumValueMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::SumValueMethod::InitializeStatisticsVariables<double>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::SumValueMethod::InitializeStatisticsVariables<array_1d<double, 3>>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::SumValueMethod::InitializeStatisticsVariables<Vector>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::SumValueMethod::InitializeStatisticsVariables<Matrix>)
            ;

    auto temporal_historical_non_historical_output_norm_method = temporal_historical_non_historical_output_method.def_submodule("NormMethods");
    using HistoricalInputNonHistoricalOutputTemporalMethods = TemporalMethods::HistoricalInputNonHistoricalOutputTemporalMethods;
    py::class_<HistoricalInputNonHistoricalOutputTemporalMethods::VarianceNormMethod, HistoricalInputNonHistoricalOutputTemporalMethods::VarianceNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_non_historical_output_norm_method, "Variance")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::VarianceNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<HistoricalInputNonHistoricalOutputTemporalMethods::MeanNormMethod, HistoricalInputNonHistoricalOutputTemporalMethods::MeanNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_non_historical_output_norm_method, "Mean")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::MeanNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<HistoricalInputNonHistoricalOutputTemporalMethods::SumNormMethod, HistoricalInputNonHistoricalOutputTemporalMethods::SumNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_non_historical_output_norm_method, "Sum")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::SumNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::SumNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::SumNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::SumNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::SumNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<HistoricalInputNonHistoricalOutputTemporalMethods::MinNormMethod, HistoricalInputNonHistoricalOutputTemporalMethods::MinNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_non_historical_output_norm_method, "Min")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MinNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MinNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MinNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MinNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::MinNormMethod::InitializeStatisticsVariables)
            ;
    py::class_<HistoricalInputNonHistoricalOutputTemporalMethods::MaxNormMethod, HistoricalInputNonHistoricalOutputTemporalMethods::MaxNormMethod::Pointer, TemporalMethods::TemporalMethod>(
            temporal_historical_non_historical_output_norm_method, "Max")
            .def(py::init<ModelPart&>())
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MaxNormMethod::CalculateStatistics<double>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MaxNormMethod::CalculateStatistics<array_1d<double, 3>>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MaxNormMethod::CalculateStatistics<Vector>)
            .def("CalculateStatistics", &HistoricalInputNonHistoricalOutputTemporalMethods::MaxNormMethod::CalculateStatistics<Matrix>)
            .def("InitializeVariables", &HistoricalInputNonHistoricalOutputTemporalMethods::MaxNormMethod::InitializeStatisticsVariables)
            ;
}

} // namespace Python.
} // Namespace Kratos

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
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "custom_methods/temporal_methods.h"

// Include base h
#include "custom_python/add_custom_temporal_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomTemporalMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Adding temporal methods
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

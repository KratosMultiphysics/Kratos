//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes

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

    auto temporal_method_module = m.def_submodule("TemporalMethods");

    py::class_<TemporalMethods::TemporalMethod, TemporalMethods::TemporalMethod::Pointer>(temporal_method_module,"TemporalMethod")
        .def(py::init<ModelPart&, const int>())
        .def("GetModelPart", &TemporalMethods::TemporalMethod::GetModelPart)
        .def("GetTotalTime", &TemporalMethods::TemporalMethod::GetTotalTime)
        .def("GetEchoLevel", &TemporalMethods::TemporalMethod::GetEchoLevel)
        .def("InitializeStatisticsMethod", &TemporalMethods::TemporalMethod::InitializeStatisticsMethod)
        .def("CalculateStatistics", &TemporalMethods::TemporalMethod::CalculateStatistics)
        .def("FinalizeStatisticsTimeStep", &TemporalMethods::TemporalMethod::FinalizeStatisticsTimeStep)
    ;

    auto temporal_historical_method = temporal_method_module.def_submodule("Historical");
    auto temporal_historical_historical_method = temporal_historical_method.def_submodule("HistoricalOutput");
    temporal_historical_historical_method.def_submodule("ValueMethods");
    temporal_historical_historical_method.def_submodule("NormMethods");
    auto temporal_historical_non_historical_method = temporal_historical_method.def_submodule("NonHistoricalOutput");
    temporal_historical_non_historical_method.def_submodule("ValueMethods");
    temporal_historical_non_historical_method.def_submodule("NormMethods");
    auto temporal_non_historical_method = temporal_method_module.def_submodule("NonHistorical");
    auto temporal_non_historical_nodal_method = temporal_non_historical_method.def_submodule("Nodes");
    temporal_non_historical_nodal_method.def_submodule("ValueMethods");
    temporal_non_historical_nodal_method.def_submodule("NormMethods");
    auto temporal_non_historical_condition_method = temporal_non_historical_method.def_submodule("Conditions");
    temporal_non_historical_condition_method.def_submodule("ValueMethods");
    temporal_non_historical_condition_method.def_submodule("NormMethods");
    auto temporal_non_historical_element_method = temporal_non_historical_method.def_submodule("Elements");
    temporal_non_historical_element_method.def_submodule("ValueMethods");
    temporal_non_historical_element_method.def_submodule("NormMethods");

    ADD_KRATOS_STATISTICS_TEMPORAL_VALUE_METHOD_PYTHON_INTERFACE(SumMethod, "Sum", m, KRATOS_STATISTICS_ONE_OUTPUT_VALUE_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_VALUE_METHOD_PYTHON_INTERFACE(MeanMethod, "Mean", m, KRATOS_STATISTICS_ONE_OUTPUT_VALUE_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_VALUE_METHOD_PYTHON_INTERFACE(VarianceMethod, "Variance", m, KRATOS_STATISTICS_TWO_OUTPUT_VALUE_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_VALUE_METHOD_PYTHON_INTERFACE(RootMeanSquareMethod, "RootMeanSquare", m, KRATOS_STATISTICS_ONE_OUTPUT_VALUE_TYPE)

    ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(SumMethod, "Sum", m, KRATOS_STATISTICS_ONE_OUTPUT_NORM_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(MeanMethod, "Mean", m, KRATOS_STATISTICS_ONE_OUTPUT_NORM_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(VarianceMethod, "Variance", m, KRATOS_STATISTICS_TWO_OUTPUT_NORM_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(RootMeanSquareMethod, "RootMeanSquare", m, KRATOS_STATISTICS_ONE_OUTPUT_NORM_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(MinMethod, "Min", m, KRATOS_STATISTICS_TWO_OUTPUT_NORM_TYPE)
    ADD_KRATOS_STATISTICS_TEMPORAL_NORM_METHOD_PYTHON_INTERFACE(MaxMethod, "Max", m, KRATOS_STATISTICS_TWO_OUTPUT_NORM_TYPE)

    ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(SumMethod, "Sum", m)
    ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(MeanMethod, "Mean", m)
    ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(VarianceMethod, "Variance", m)
    ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(RootMeanSquareMethod, "RootMeanSquare", m)
    ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(MinMethod, "Min", m)
    ADD_KRATOS_STATISTICS_TEMPORAL_CREATE_METHOD_PYTHON_INTERFACE(MaxMethod, "Max", m)
}

} // namespace Python.
} // Namespace Kratos

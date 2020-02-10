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
#include "custom_methods/spatial_methods.h"

// Include base h
#include "custom_python/add_custom_spatial_non_historical_condition_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomSpatialNonHistoricalConditionMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Adding spatial methods
    auto method_module = m.def_submodule("Conditions");

    using ConditionNonHistoricalSpatialMethods = SpatialMethods::ConditionNonHistoricalSpatialMethods;

    method_module.def_submodule("ValueMethods")
        .def("Sum", &ConditionNonHistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &ConditionNonHistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        ;
    method_module.def_submodule("NormMethods")
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
        .def("Min", &ConditionNonHistoricalSpatialMethods::GetNormMin<double>)
        .def("Min", &ConditionNonHistoricalSpatialMethods::GetNormMin<array_1d<double, 3>>)
        .def("Min", &ConditionNonHistoricalSpatialMethods::GetNormMin<Vector>)
        .def("Min", &ConditionNonHistoricalSpatialMethods::GetNormMin<Matrix>)
        .def("Max", &ConditionNonHistoricalSpatialMethods::GetNormMax<double>)
        .def("Max", &ConditionNonHistoricalSpatialMethods::GetNormMax<array_1d<double, 3>>)
        .def("Max", &ConditionNonHistoricalSpatialMethods::GetNormMax<Vector>)
        .def("Max", &ConditionNonHistoricalSpatialMethods::GetNormMax<Matrix>)
        .def("Median", &ConditionNonHistoricalSpatialMethods::GetNormMedian<double>)
        .def("Median", &ConditionNonHistoricalSpatialMethods::GetNormMedian<array_1d<double, 3>>)
        .def("Median", &ConditionNonHistoricalSpatialMethods::GetNormMedian<Vector>)
        .def("Median", &ConditionNonHistoricalSpatialMethods::GetNormMedian<Matrix>)
        .def("Distribution", &ConditionNonHistoricalSpatialMethods::GetNormDistribution<double>)
        .def("Distribution", &ConditionNonHistoricalSpatialMethods::GetNormDistribution<array_1d<double, 3>>)
        .def("Distribution", &ConditionNonHistoricalSpatialMethods::GetNormDistribution<Vector>)
        .def("Distribution", &ConditionNonHistoricalSpatialMethods::GetNormDistribution<Matrix>)
        ;
}

} // namespace Python.
} // Namespace Kratos

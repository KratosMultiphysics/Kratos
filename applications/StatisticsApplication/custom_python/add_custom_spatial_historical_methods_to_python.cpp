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
#include "custom_methods/spatial_methods.h"

// Include base h
#include "add_custom_spatial_historical_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomSpatialHistoricalMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Adding spatial methods
    auto spatial_historical_method = m.def_submodule("Historical");

    using HistoricalSpatialMethods = SpatialMethods::HistoricalSpatialMethods;
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
        .def("Min", &HistoricalSpatialMethods::GetNormMin<double>)
        .def("Min", &HistoricalSpatialMethods::GetNormMin<array_1d<double, 3>>)
        .def("Min", &HistoricalSpatialMethods::GetNormMin<Vector>)
        .def("Min", &HistoricalSpatialMethods::GetNormMin<Matrix>)
        .def("Max", &HistoricalSpatialMethods::GetNormMax<double>)
        .def("Max", &HistoricalSpatialMethods::GetNormMax<array_1d<double, 3>>)
        .def("Max", &HistoricalSpatialMethods::GetNormMax<Vector>)
        .def("Max", &HistoricalSpatialMethods::GetNormMax<Matrix>)
        .def("Median", &HistoricalSpatialMethods::GetNormMedian<double>)
        .def("Median", &HistoricalSpatialMethods::GetNormMedian<array_1d<double, 3>>)
        .def("Median", &HistoricalSpatialMethods::GetNormMedian<Vector>)
        .def("Median", &HistoricalSpatialMethods::GetNormMedian<Matrix>)
        .def("Distribution", &HistoricalSpatialMethods::GetNormDistribution<double>)
        .def("Distribution", &HistoricalSpatialMethods::GetNormDistribution<array_1d<double, 3>>)
        .def("Distribution", &HistoricalSpatialMethods::GetNormDistribution<Vector>)
        .def("Distribution", &HistoricalSpatialMethods::GetNormDistribution<Matrix>)
        ;
}

} // namespace Python.
} // Namespace Kratos

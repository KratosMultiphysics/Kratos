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
#include "custom_python/add_custom_spatial_non_historical_element_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomSpatialNonHistoricalElementMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Adding spatial methods
    auto method_module = m.def_submodule("Elements");

    using ElementNonHistoricalSpatialMethods = SpatialMethods::ElementNonHistoricalSpatialMethods;

    method_module.def_submodule("ValueMethods")
        .def("Sum", &ElementNonHistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &ElementNonHistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        ;
    method_module.def_submodule("NormMethods")
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
        .def("Min", &ElementNonHistoricalSpatialMethods::GetNormMin<double>)
        .def("Min", &ElementNonHistoricalSpatialMethods::GetNormMin<array_1d<double, 3>>)
        .def("Min", &ElementNonHistoricalSpatialMethods::GetNormMin<Vector>)
        .def("Min", &ElementNonHistoricalSpatialMethods::GetNormMin<Matrix>)
        .def("Max", &ElementNonHistoricalSpatialMethods::GetNormMax<double>)
        .def("Max", &ElementNonHistoricalSpatialMethods::GetNormMax<array_1d<double, 3>>)
        .def("Max", &ElementNonHistoricalSpatialMethods::GetNormMax<Vector>)
        .def("Max", &ElementNonHistoricalSpatialMethods::GetNormMax<Matrix>)
        .def("Median", &ElementNonHistoricalSpatialMethods::GetNormMedian<double>)
        .def("Median", &ElementNonHistoricalSpatialMethods::GetNormMedian<array_1d<double, 3>>)
        .def("Median", &ElementNonHistoricalSpatialMethods::GetNormMedian<Vector>)
        .def("Median", &ElementNonHistoricalSpatialMethods::GetNormMedian<Matrix>)
        .def("Distribution", &ElementNonHistoricalSpatialMethods::GetNormDistribution<double>)
        .def("Distribution", &ElementNonHistoricalSpatialMethods::GetNormDistribution<array_1d<double, 3>>)
        .def("Distribution", &ElementNonHistoricalSpatialMethods::GetNormDistribution<Vector>)
        .def("Distribution", &ElementNonHistoricalSpatialMethods::GetNormDistribution<Matrix>)
        ;

}

} // namespace Python.
} // Namespace Kratos

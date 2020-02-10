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
#include "custom_python/add_custom_spatial_non_historical_nodal_methods_to_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomSpatialNonHistoricalNodalMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Adding spatial methods
    auto method_module = m.def_submodule("Nodes");

    using NodalNonHistoricalSpatialMethods = SpatialMethods::NodalNonHistoricalSpatialMethods;
    method_module.def_submodule("ValueMethods")
        .def("Sum", &NodalNonHistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &NodalNonHistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        ;
    method_module.def_submodule("NormMethods")
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
        .def("Min", &NodalNonHistoricalSpatialMethods::GetNormMin<double>)
        .def("Min", &NodalNonHistoricalSpatialMethods::GetNormMin<array_1d<double, 3>>)
        .def("Min", &NodalNonHistoricalSpatialMethods::GetNormMin<Vector>)
        .def("Min", &NodalNonHistoricalSpatialMethods::GetNormMin<Matrix>)
        .def("Max", &NodalNonHistoricalSpatialMethods::GetNormMax<double>)
        .def("Max", &NodalNonHistoricalSpatialMethods::GetNormMax<array_1d<double, 3>>)
        .def("Max", &NodalNonHistoricalSpatialMethods::GetNormMax<Vector>)
        .def("Max", &NodalNonHistoricalSpatialMethods::GetNormMax<Matrix>)
        .def("Median", &NodalNonHistoricalSpatialMethods::GetNormMedian<double>)
        .def("Median", &NodalNonHistoricalSpatialMethods::GetNormMedian<array_1d<double, 3>>)
        .def("Median", &NodalNonHistoricalSpatialMethods::GetNormMedian<Vector>)
        .def("Median", &NodalNonHistoricalSpatialMethods::GetNormMedian<Matrix>)
        .def("Distribution", &NodalNonHistoricalSpatialMethods::GetNormDistribution<double>)
        .def("Distribution", &NodalNonHistoricalSpatialMethods::GetNormDistribution<array_1d<double, 3>>)
        .def("Distribution", &NodalNonHistoricalSpatialMethods::GetNormDistribution<Vector>)
        .def("Distribution", &NodalNonHistoricalSpatialMethods::GetNormDistribution<Matrix>)
        ;
}

} // namespace Python.
} // Namespace Kratos

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
#include "custom_python/add_custom_methods_to_python.h"
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "custom_methods/spatial_methods.h"
#include "custom_methods/temporal_methods.h"
#include "custom_utilities/method_utilities.h"

namespace Kratos
{
namespace Python
{
void AddCustomSpatialNonHistoricalMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Adding spatial methods
    auto non_historical_spatial_method_module = m.def_submodule("NonHistorical");

    using NodalNonHistoricalSpatialMethods = SpatialMethods::NodalNonHistoricalSpatialMethods;
    auto spatial_non_historical_nodal_method = non_historical_spatial_method_module.def_submodule("Nodes");
    spatial_non_historical_nodal_method.def_submodule("ValueMethods")
        .def("Sum", &NodalNonHistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &NodalNonHistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        ;
    spatial_non_historical_nodal_method.def_submodule("NormMethods")
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

    using ConditionNonHistoricalSpatialMethods = SpatialMethods::ConditionNonHistoricalSpatialMethods;
    auto spatial_non_historical_condition_method = non_historical_spatial_method_module.def_submodule("Conditions");
    spatial_non_historical_condition_method.def_submodule("ValueMethods")
        .def("Sum", &ConditionNonHistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &ConditionNonHistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        ;
    spatial_non_historical_condition_method.def_submodule("NormMethods")
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

    using ElementNonHistoricalSpatialMethods = SpatialMethods::ElementNonHistoricalSpatialMethods;
    auto spatial_non_historical_element_method = non_historical_spatial_method_module.def_submodule("Elements");
    spatial_non_historical_element_method.def_submodule("ValueMethods")
        .def("Sum", &ElementNonHistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &ElementNonHistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        ;
    spatial_non_historical_element_method.def_submodule("NormMethods")
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

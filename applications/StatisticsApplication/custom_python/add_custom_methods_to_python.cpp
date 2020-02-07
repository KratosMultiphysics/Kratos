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
#include "custom_utilities/method_utilities.h"

namespace Kratos
{
namespace Python
{
void AddCustomMethodsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using NodeType = ModelPart::NodeType;
    using ElementType = ModelPart::ElementType;
    using ConditionType = ModelPart::ConditionType;

    using NodesContainerType = ModelPart::NodesContainerType;
    using ElementsContainerType = ModelPart::ElementsContainerType;
    using ConditionsContainerType = ModelPart::ConditionsContainerType;

    auto spatial_method_module = m.def_submodule("SpatialMethods");

    using HistoricalSpatialMethods = SpatialMethods::ContainerSpatialMethods<NodesContainerType, NodeType, MethodsUtilities::HistoricalDataValueRetrievalFunctor>;
    auto spatial_historical_method = spatial_method_module.def_submodule("Historical");
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
        .def("Mean", &HistoricalSpatialMethods::CalculateNormMean<double>)
        .def("Mean", &HistoricalSpatialMethods::CalculateNormMean<array_1d<double, 3>>)
        .def("Variance", &HistoricalSpatialMethods::CalculateNormVariance<double>)
        .def("Variance", &HistoricalSpatialMethods::CalculateNormVariance<array_1d<double, 3>>)
        .def("Min", &HistoricalSpatialMethods::GetMin<double>)
        .def("Min", &HistoricalSpatialMethods::GetMin<array_1d<double, 3>>)
        .def("Max", &HistoricalSpatialMethods::GetMax<double>)
        .def("Max", &HistoricalSpatialMethods::GetMax<array_1d<double, 3>>)
        ;

    using NodalNonHistoricalSpatialMethods = SpatialMethods::ContainerSpatialMethods<NodesContainerType, NodeType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    auto non_historical_spatial_method_module = spatial_method_module.def_submodule("NonHistorical");
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
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateNormMean<double>)
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateNormMean<array_1d<double, 3>>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateNormVariance<double>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateNormVariance<array_1d<double, 3>>)
        .def("Min", &NodalNonHistoricalSpatialMethods::GetMin<double>)
        .def("Min", &NodalNonHistoricalSpatialMethods::GetMin<array_1d<double, 3>>)
        .def("Max", &NodalNonHistoricalSpatialMethods::GetMax<double>)
        .def("Max", &NodalNonHistoricalSpatialMethods::GetMax<array_1d<double, 3>>)
        ;

    using ConditionNonHistoricalSpatialMethods = SpatialMethods::ContainerSpatialMethods<ConditionsContainerType, ConditionType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
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
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateNormMean<double>)
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateNormMean<array_1d<double, 3>>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateNormVariance<double>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateNormVariance<array_1d<double, 3>>)
        .def("Min", &ConditionNonHistoricalSpatialMethods::GetMin<double>)
        .def("Min", &ConditionNonHistoricalSpatialMethods::GetMin<array_1d<double, 3>>)
        .def("Max", &ConditionNonHistoricalSpatialMethods::GetMax<double>)
        .def("Max", &ConditionNonHistoricalSpatialMethods::GetMax<array_1d<double, 3>>)
        ;

    using ElementNonHistoricalSpatialMethods = SpatialMethods::ContainerSpatialMethods<ElementsContainerType, ElementType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
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
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateNormMean<double>)
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateNormMean<array_1d<double, 3>>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateNormVariance<double>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateNormVariance<array_1d<double, 3>>)
        .def("Min", &ElementNonHistoricalSpatialMethods::GetMin<double>)
        .def("Min", &ElementNonHistoricalSpatialMethods::GetMin<array_1d<double, 3>>)
        .def("Max", &ElementNonHistoricalSpatialMethods::GetMax<double>)
        .def("Max", &ElementNonHistoricalSpatialMethods::GetMax<array_1d<double, 3>>)
        ;
}

} // namespace Python.
} // Namespace Kratos

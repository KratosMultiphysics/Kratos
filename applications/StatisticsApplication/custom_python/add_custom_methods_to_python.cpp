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
    spatial_method_module.def_submodule("Historical")
        .def("Sum", &HistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &HistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &HistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &HistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &HistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &HistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        .def("Min", &HistoricalSpatialMethods::GetMin<double>)
        .def("Min", &HistoricalSpatialMethods::GetMin<array_1d<double, 3>>)
        .def("Max", &HistoricalSpatialMethods::GetMax<double>)
        .def("Max", &HistoricalSpatialMethods::GetMax<array_1d<double, 3>>)
        ;

    using NodalNonHistoricalSpatialMethods = SpatialMethods::ContainerSpatialMethods<NodesContainerType, NodeType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    auto non_historical_spatial_method_module = spatial_method_module.def_submodule("NonHistorical");
    non_historical_spatial_method_module.def_submodule("Nodes")
        .def("Sum", &NodalNonHistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &NodalNonHistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &NodalNonHistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &NodalNonHistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        .def("Min", &NodalNonHistoricalSpatialMethods::GetMin<double>)
        .def("Min", &NodalNonHistoricalSpatialMethods::GetMin<array_1d<double, 3>>)
        .def("Max", &NodalNonHistoricalSpatialMethods::GetMax<double>)
        .def("Max", &NodalNonHistoricalSpatialMethods::GetMax<array_1d<double, 3>>)
        ;

    using ConditionNonHistoricalSpatialMethods = SpatialMethods::ContainerSpatialMethods<ConditionsContainerType, ConditionType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    non_historical_spatial_method_module.def_submodule("Conditions")
        .def("Sum", &ConditionNonHistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &ConditionNonHistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &ConditionNonHistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &ConditionNonHistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        .def("Min", &ConditionNonHistoricalSpatialMethods::GetMin<double>)
        .def("Min", &ConditionNonHistoricalSpatialMethods::GetMin<array_1d<double, 3>>)
        .def("Max", &ConditionNonHistoricalSpatialMethods::GetMax<double>)
        .def("Max", &ConditionNonHistoricalSpatialMethods::GetMax<array_1d<double, 3>>)
        ;

    using ElementNonHistoricalSpatialMethods = SpatialMethods::ContainerSpatialMethods<ElementsContainerType, ElementType, MethodsUtilities::NonHistoricalDataValueRetrievalFunctor>;
    non_historical_spatial_method_module.def_submodule("Elements")
        .def("Sum", &ElementNonHistoricalSpatialMethods::CalculateSum<double>)
        .def("Sum", &ElementNonHistoricalSpatialMethods::CalculateSum<array_1d<double, 3>>)
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateMean<double>)
        .def("Mean", &ElementNonHistoricalSpatialMethods::CalculateMean<array_1d<double, 3>>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateVariance<double>)
        .def("Variance", &ElementNonHistoricalSpatialMethods::CalculateVariance<array_1d<double, 3>>)
        .def("Min", &ElementNonHistoricalSpatialMethods::GetMin<double>)
        .def("Min", &ElementNonHistoricalSpatialMethods::GetMin<array_1d<double, 3>>)
        .def("Max", &ElementNonHistoricalSpatialMethods::GetMax<double>)
        .def("Max", &ElementNonHistoricalSpatialMethods::GetMax<array_1d<double, 3>>)
        ;
}

} // namespace Python.
} // Namespace Kratos

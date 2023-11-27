//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//                   Ihar Antonau
//                   Fabian Meister
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes
#include "includes/define.h"
#include "utilities/string_utilities.h"

// Application includes
#include "custom_utilities/control_utils.h"
#include "custom_utilities/sensor_utils.h"
#include "custom_utilities/domain_sensor_view_cluster_data.h"
#include "custom_utilities/sensor_view_cluster.h"

// Include base h
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos::Python {

template<class TContainerType>
void AddSensorClusterUtilsToPython(
    pybind11::module& m)
{
    namespace py = pybind11;

    std::string upper_prefix, lower_prefix;
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        upper_prefix = "Nodal";
        lower_prefix = "nodal";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        upper_prefix = "Condition";
        lower_prefix = "condition";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        upper_prefix = "Element";
        lower_prefix = "element";
    }

    const std::string& cluster_data_name = upper_prefix + "SensorViewClusterData";
    using sensor_cluster_data = DomainSensorViewClusterData<TContainerType>;
    py::class_<sensor_cluster_data, typename sensor_cluster_data::Pointer>(m, cluster_data_name.c_str())
        .def(py::init<const typename sensor_cluster_data::SensorViewVectorType&>(), py::arg("sensor_views_list"))
        .def("AddDistances", &sensor_cluster_data::AddDistances, py::arg("distances_type"), py::arg("compressed_distances_matrix"))
        .def("GetContainer", &sensor_cluster_data::GetContainer, py::return_value_policy::reference)
        .def("GetSensorViews", &sensor_cluster_data::GetSensorViews, py::return_value_policy::reference)
        .def("GetModelPart", &sensor_cluster_data::GetModelPart, py::return_value_policy::reference_internal)
        ;

    const std::string& cluster_name = upper_prefix + "SensorViewCluster";
    const std::string& sensor_views_py_hint = lower_prefix + "_sensor_views";
    using cluster = SensorViewCluster<TContainerType>;
    py::class_<cluster, typename cluster::Pointer, IndexedObject>(m, cluster_name.c_str())
        .def(py::init<const IndexType, typename DomainSensorViewClusterData<TContainerType>::Pointer>(), py::arg("cluster_id"), py::arg(StringUtilities::ConvertCamelCaseToSnakeCase(cluster_data_name).c_str()))
        .def("Clone", &cluster::Clone)
        .def("Clear", &cluster::Clear)
        .def("SetSensorViews", &cluster::SetSensorViews, py::arg(sensor_views_py_hint.c_str()))
        .def("GetSensorViews", &cluster::GetSensorViews)
        .def("GetDistances", &cluster::GetDistances, py::arg("distances_type"))
        .def("SetEntities", &cluster::SetEntities, py::arg((lower_prefix + "_array").c_str()))
        .def("GetEntities", py::overload_cast<>(&cluster::GetEntities, py::const_), py::return_value_policy::reference)
        .def("GetDataContainer", &cluster::GetDataContainer)
        ;
}

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto control_utils = m.def_submodule("ControlUtils");
    control_utils.def("AssignEquivalentProperties", &ControlUtils::AssignEquivalentProperties<ModelPart::ConditionsContainerType>, py::arg("source_conditions"), py::arg("destination_conditions"));
    control_utils.def("AssignEquivalentProperties", &ControlUtils::AssignEquivalentProperties<ModelPart::ElementsContainerType>, py::arg("source_elements"), py::arg("destination_elements"));
    control_utils.def("ClipContainerExpression", &ControlUtils::ClipContainerExpression<ModelPart::NodesContainerType>, py::arg("nodal_expression"), py::arg("min"), py::arg("max"));
    control_utils.def("ClipContainerExpression", &ControlUtils::ClipContainerExpression<ModelPart::ConditionsContainerType>, py::arg("condition_expression"), py::arg("min"), py::arg("max"));
    control_utils.def("ClipContainerExpression", &ControlUtils::ClipContainerExpression<ModelPart::ElementsContainerType>, py::arg("element_expression"), py::arg("min"), py::arg("max"));

    auto sensor_utils = m.def_submodule("SensorUtils");
    sensor_utils.def("GetBestSensorViewForEveryEntity", [](const std::vector<SensorView<ModelPart::NodesContainerType>::Pointer>& rNormalizedSensorViews) { std::vector<SensorView<ModelPart::NodesContainerType>::Pointer> output;  SensorUtils::IdentifyBestSensorViewForEveryEntity<ModelPart::NodesContainerType>(output, rNormalizedSensorViews); return output;}, py::arg("normalized_nodal_sensor_views_list"));
    sensor_utils.def("GetBestSensorViewForEveryEntity", [](const std::vector<SensorView<ModelPart::ConditionsContainerType>::Pointer>& rNormalizedSensorViews) { std::vector<SensorView<ModelPart::ConditionsContainerType>::Pointer> output;  SensorUtils::IdentifyBestSensorViewForEveryEntity<ModelPart::ConditionsContainerType>(output, rNormalizedSensorViews); return output;}, py::arg("normalized_condition_sensor_views_list"));
    sensor_utils.def("GetBestSensorViewForEveryEntity", [](const std::vector<SensorView<ModelPart::ElementsContainerType>::Pointer>& rNormalizedSensorViews) { std::vector<SensorView<ModelPart::ElementsContainerType>::Pointer> output;  SensorUtils::IdentifyBestSensorViewForEveryEntity<ModelPart::ElementsContainerType>(output, rNormalizedSensorViews); return output;}, py::arg("normalized_element_sensor_views_list"));
    sensor_utils.def("GetSensorViewsModelPart", &SensorUtils::GetSensorViewsModelPart<ModelPart::NodesContainerType>, py::arg("nodal_sensor_views_list"), py::return_value_policy::reference_internal);
    sensor_utils.def("GetSensorViewsModelPart", &SensorUtils::GetSensorViewsModelPart<ModelPart::ConditionsContainerType>, py::arg("condition_sensor_views_list"), py::return_value_policy::reference_internal);
    sensor_utils.def("GetSensorViewsModelPart", &SensorUtils::GetSensorViewsModelPart<ModelPart::ElementsContainerType>, py::arg("element_sensor_views_list"), py::return_value_policy::reference_internal);
    sensor_utils.def("GetDomainSize", &SensorUtils::GetDomainSize<ModelPart::NodesContainerType>, py::arg("nodes_array"), py::arg("data_communicator"));
    sensor_utils.def("GetDomainSize", &SensorUtils::GetDomainSize<ModelPart::ConditionsContainerType>, py::arg("conditions_array"), py::arg("data_communicator"));
    sensor_utils.def("GetDomainSize", &SensorUtils::GetDomainSize<ModelPart::ElementsContainerType>, py::arg("elements_array"), py::arg("data_communicator"));
    sensor_utils.def("AssignEntitiesToClusters", &SensorUtils::AssignEntitiesToClusters<ModelPart::NodesContainerType>, py::arg("nodal_sensor_view_clusters"), py::arg("nodal_expressions_list"));
    sensor_utils.def("AssignEntitiesToClusters", &SensorUtils::AssignEntitiesToClusters<ModelPart::ConditionsContainerType>, py::arg("condition_sensor_view_clusters"), py::arg("condition_expression_list"));
    sensor_utils.def("AssignEntitiesToClusters", &SensorUtils::AssignEntitiesToClusters<ModelPart::ElementsContainerType>, py::arg("element_sensor_view_clusters"), py::arg("elements_expression_list"));

    auto cluster_utils = m.def_submodule("ClusterUtils");
    AddSensorClusterUtilsToPython<ModelPart::NodesContainerType>(cluster_utils);
    AddSensorClusterUtilsToPython<ModelPart::ConditionsContainerType>(cluster_utils);
    AddSensorClusterUtilsToPython<ModelPart::ElementsContainerType>(cluster_utils);

}

} // namespace Kratos::Python

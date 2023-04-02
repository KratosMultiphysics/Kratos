//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl,
//                   Suneth Warnakulasuriya
//

// System includes
#include "sstream"

// External includes
#include <pybind11/stl.h>

// Project includes
#include "includes/model_part.h"
#include "containers/container_expression/container_data_io.h"
#include "containers/container_expression/specialized_container_expression.h"

// Application includes
#include "custom_utilities/geometrical/symmetry_utility.h"
#include "custom_utilities/geometrical/model_part_utils.h"
#include "custom_utilities/optimization_utils.h"
#include "custom_utilities/optimization_info.h"

// Include base h
#include "add_custom_response_utilities_to_python.h"

namespace Kratos {
namespace Python {

template<class TDataType>
std::string GetDataTypeName(const std::string& rPrefix);

template<> std::string GetDataTypeName<bool>(const std::string& rPrefix) { return std::string(rPrefix).append("Bool"); }
template<> std::string GetDataTypeName<int>(const std::string& rPrefix) { return std::string(rPrefix).append("Int"); }
template<> std::string GetDataTypeName<double>(const std::string& rPrefix) { return std::string(rPrefix).append("Double"); }
template<> std::string GetDataTypeName<std::string>(const std::string& rPrefix) { return std::string(rPrefix).append("String"); }
template<> std::string GetDataTypeName<SpecializedContainerExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::Historical>>::Pointer>(const std::string& rPrefix) { return std::string(rPrefix).append("HistoricalExpression"); }
template<> std::string GetDataTypeName<SpecializedContainerExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>::Pointer>(const std::string& rPrefix) { return std::string(rPrefix).append("NodalNonHistoricalExpression"); }
template<> std::string GetDataTypeName<SpecializedContainerExpression<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>::Pointer>(const std::string& rPrefix) { return std::string(rPrefix).append("ConditionNonHistoricalExpression"); }
template<> std::string GetDataTypeName<SpecializedContainerExpression<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>::Pointer>(const std::string& rPrefix) { return std::string(rPrefix).append("ElementNonHistoricalExpression"); }

template<class... TArgs>
void AddOptimizationInfoToPython(pybind11::module m)
{
    namespace py = pybind11;

    using OptimizationInfoType = OptimizationInfo<TArgs...>;
    auto optimization_info_class = py::class_<OptimizationInfoType, typename OptimizationInfoType::Pointer>(m, "OptimizationInfo")
        .def(py::init<>())
        .def(py::init<std::size_t>(), py::arg("buffer_size"))
        .def("AdvanceStep", &OptimizationInfoType::AdvanceStep)
        .def("SetBufferSize", &OptimizationInfoType::SetBufferSize, py::arg("buffer_size"), py::arg("resize_sub_items") = false)
        .def("GetBufferSize", &OptimizationInfoType::GetBufferSize)
        .def("HasValue", &OptimizationInfoType::HasValue, py::arg("name"), py::arg("step_index") = 0)
        .def("SetValue", &OptimizationInfoType::SetValue, py::arg("name"), py::arg("value"), py::arg("step_index") = 0, py::arg("overwrite") = false)
        .def("IsSubItem", &OptimizationInfoType::template IsValue<typename OptimizationInfoType::Pointer>, py::arg("name"), py::arg("step_index") = 0)
        .def("GetSubItem", &OptimizationInfoType::template GetValue<typename OptimizationInfoType::Pointer>, py::arg("name"), py::arg("step_index") = 0)
        .def("__str__", [](const OptimizationInfoType& rSelf) { return rSelf.Info(); })
        ;

    (optimization_info_class.def(GetDataTypeName<TArgs>("Is").c_str(), &OptimizationInfoType::template IsValue<TArgs>, py::arg("name"), py::arg("step_index") = 0), ...);
    (optimization_info_class.def(GetDataTypeName<TArgs>("Get").c_str(), &OptimizationInfoType::template GetValue<TArgs>, py::arg("name"), py::arg("step_index") = 0), ...);
}

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using namespace pybind11::literals;

    py::class_<SymmetryUtility >(m, "SymmetryUtility")
        .def(py::init<std::string, ModelPart&, Parameters>())
        .def("Initialize", &SymmetryUtility::Initialize)
        .def("Update", &SymmetryUtility::Update)
        .def("ApplyOnVectorField", &SymmetryUtility::ApplyOnVectorField)
        .def("ApplyOnScalarField", &SymmetryUtility::ApplyOnScalarField)
        ;

    m.def_submodule("ModelPartUtils")
        .def("GetModelPartsWithCommonReferenceEntities", [](
            const std::vector<ModelPart*>& rEvaluatedModelPartsList,
            const std::vector<ModelPart*>& rReferenceModelPartsList,
            const bool AreNodesConsidered,
            const bool AreConditionsConsidered,
            const bool AreElementsConsidered,
            const bool AreParentsConsidered,
            const IndexType EchoLevel){
                const auto& r_model_parts = ModelPartUtils::GetModelPartsWithCommonReferenceEntities(
                    rEvaluatedModelPartsList,
                    rReferenceModelPartsList,
                    AreNodesConsidered,
                    AreConditionsConsidered,
                    AreElementsConsidered,
                    AreParentsConsidered,
                    EchoLevel);

                auto pylist = py::list();
                for (auto p_model_part : r_model_parts) {
                    auto pyobj = py::cast(*p_model_part, py::return_value_policy::reference);
                    pylist.append(pyobj);
                }
                return pylist;
            },
            "examined_model_parts_list"_a,
            "reference_model_parts"_a,
            "are_nodes_considered"_a,
            "are_conditions_considered"_a,
            "are_elements_considered"_a,
            "are_parents_considered"_a,
            "echo_level"_a = 0)
        .def("RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList", &ModelPartUtils::RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList,
            "model_parts_list"_a)
        ;

    py::class_<OptimizationUtils >(m, "OptimizationUtils")
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ConditionsContainerType, double>)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ElementsContainerType,double>)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ElementsContainerType,array_1d<double, 3>>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ConditionsContainerType, double>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ElementsContainerType,double>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ElementsContainerType,array_1d<double, 3>>)
        .def_static("AreAllEntitiesOfSameGeometryType", [](ModelPart::ConditionsContainerType& rContainer, const DataCommunicator& rDataCommunicator) { return OptimizationUtils::GetContainerEntityGeometryType(rContainer, rDataCommunicator) != GeometryData::KratosGeometryType::Kratos_generic_type; } )
        .def_static("AreAllEntitiesOfSameGeometryType", [](ModelPart::ElementsContainerType& rContainer, const DataCommunicator& rDataCommunicator) { return OptimizationUtils::GetContainerEntityGeometryType(rContainer, rDataCommunicator) != GeometryData::KratosGeometryType::Kratos_generic_type; } )
        .def_static("CreateEntitySpecificPropertiesForContainer", &OptimizationUtils::CreateEntitySpecificPropertiesForContainer<ModelPart::ConditionsContainerType>)
        .def_static("CreateEntitySpecificPropertiesForContainer", &OptimizationUtils::CreateEntitySpecificPropertiesForContainer<ModelPart::ElementsContainerType>)
        .def_static("GetVariableDimension", &OptimizationUtils::GetVariableDimension<double>)
        .def_static("GetVariableDimension", &OptimizationUtils::GetVariableDimension<array_1d<double, 3>>)
        .def_static("CopySolutionStepVariablesList", &OptimizationUtils::CopySolutionStepVariablesList)
        ;

    AddOptimizationInfoToPython<bool,
                                int,
                                double,
                                std::string,
                                SpecializedContainerExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::Historical>>::Pointer,
                                SpecializedContainerExpression<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>::Pointer,
                                SpecializedContainerExpression<ModelPart::ConditionsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>::Pointer,
                                SpecializedContainerExpression<ModelPart::ElementsContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>>::Pointer
                                >(m);
}

}  // namespace Python.
} // Namespace Kratos


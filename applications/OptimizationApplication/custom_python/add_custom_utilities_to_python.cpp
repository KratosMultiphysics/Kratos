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

// External includes
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

// Project includes
#include "python/numpy_utils.h"

// Application includes
#include "custom_utilities/geometrical/opt_app_model_part_utils.h"
#include "custom_utilities/geometrical/symmetry_utility.h"
#include "custom_utilities/implicit_filter_utils.h"
#include "custom_utilities/optimization_utils.h"

// Include base h
#include "add_custom_response_utilities_to_python.h"

namespace Kratos {
namespace Python {

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

    m.def_submodule("OptAppModelPartUtils")
        .def("LogModelPartStatus", &OptAppModelPartUtils::LogModelPartStatus, py::arg("model_part"), py::arg("status_to_log"))
        .def("GetModelPartStatusLog", &OptAppModelPartUtils::GetModelPartStatusLog, py::arg("model_part"))
        .def("CheckModelPartStatus", &OptAppModelPartUtils::CheckModelPartStatus, py::arg("model_part"), py::arg("status_to_check"))
        .def("GetModelPartsWithCommonReferenceEntities", [](
            const std::vector<ModelPart*>& rEvaluatedModelPartsList,
            const std::vector<ModelPart*>& rReferenceModelPartsList,
            const bool AreNodesConsidered,
            const bool AreConditionsConsidered,
            const bool AreElementsConsidered,
            const bool AreParentsConsidered,
            const IndexType EchoLevel){
                const auto& r_model_parts = OptAppModelPartUtils::GetModelPartsWithCommonReferenceEntities(
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
        .def("RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList", &OptAppModelPartUtils::RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList,
            "model_parts_list"_a)
        .def("GenerateModelPart",
            [](ModelPart::ConditionsContainerType& rOriginConditions, ModelPart& rDestinationModelPart, const std::string& rElementName) {
                OptAppModelPartUtils::GenerateModelPart(rOriginConditions, rDestinationModelPart, KratosComponents<Element>::Get(rElementName));
            },
            "conditions_container"_a,
            "destination_model_part"_a,
            "element_name"_a)
        ;

    m.def_submodule("OptimizationUtils")
        .def("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ConditionsContainerType, double>)
        .def("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ElementsContainerType,double>)
        .def("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ElementsContainerType,array_1d<double, 3>>)
        .def("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ConditionsContainerType, double>)
        .def("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ElementsContainerType,double>)
        .def("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ElementsContainerType,array_1d<double, 3>>)
        .def("AreAllEntitiesOfSameGeometryType", [](ModelPart::ConditionsContainerType& rContainer, const DataCommunicator& rDataCommunicator) { return OptimizationUtils::GetContainerEntityGeometryType(rContainer, rDataCommunicator) != GeometryData::KratosGeometryType::Kratos_generic_type; } )
        .def("AreAllEntitiesOfSameGeometryType", [](ModelPart::ElementsContainerType& rContainer, const DataCommunicator& rDataCommunicator) { return OptimizationUtils::GetContainerEntityGeometryType(rContainer, rDataCommunicator) != GeometryData::KratosGeometryType::Kratos_generic_type; } )
        .def("CreateEntitySpecificPropertiesForContainer", &OptimizationUtils::CreateEntitySpecificPropertiesForContainer<ModelPart::ConditionsContainerType>, py::arg("model_part"), py::arg("container"), py::arg("is_recursive"))
        .def("CreateEntitySpecificPropertiesForContainer", &OptimizationUtils::CreateEntitySpecificPropertiesForContainer<ModelPart::ElementsContainerType>, py::arg("model_part"), py::arg("container"), py::arg("is_recursive"))
        .def("UpdatePropertiesVariableWithRootValueRecursively", &OptimizationUtils::UpdatePropertiesVariableWithRootValueRecursively<ModelPart::ConditionsContainerType, double>, py::arg("container"), py::arg("variable"))
        .def("UpdatePropertiesVariableWithRootValueRecursively", &OptimizationUtils::UpdatePropertiesVariableWithRootValueRecursively<ModelPart::ElementsContainerType, double>, py::arg("container"), py::arg("variable"))
        .def("GetVariableDimension", &OptimizationUtils::GetVariableDimension<double>)
        .def("GetVariableDimension", &OptimizationUtils::GetVariableDimension<array_1d<double, 3>>)
        .def("SetSolutionStepVariablesList", &OptimizationUtils::SetSolutionStepVariablesList, py::arg("destination_model_part"), py::arg("origin_model_part"))
        .def("IsSolutionStepVariablesListASubSet", &OptimizationUtils::IsSolutionStepVariablesListASubSet, py::arg("main_set_model_part"), py::arg("sub_set_model_part"))
        .def("GetSolutionStepVariableNamesList", &OptimizationUtils::GetSolutionStepVariableNamesList, py::arg("model_part"))
        .def("GetComponentWiseModelParts", &OptimizationUtils::GetComponentWiseModelParts,
            py::arg("model"),
            py::arg("parameters"),
            py::return_value_policy::reference)
        .def("MapContainerDataToNodalData", OptimizationUtils::MapContainerDataToNodalData,
            py::arg("input_tensor_adaptor"),
            py::arg("nodes"))
        .def("MapNodalDataToContainerData", OptimizationUtils::MapNodalDataToContainerData<ModelPart::ConditionsContainerType::Pointer>,
            py::arg("input_tensor_adaptor"),
            py::arg("conditions"),
            py::arg("neighbour_count_tensor_adaptor"))
        .def("MapNodalDataToContainerData", OptimizationUtils::MapNodalDataToContainerData<ModelPart::ElementsContainerType::Pointer>,
            py::arg("input_tensor_adaptor"),
            py::arg("conditions"),
            py::arg("neighbour_count_tensor_adaptor"))
        ;

    m.def_submodule("ImplicitFilterUtils")
        .def("SetBulkRadiusForShapeFiltering", &ImplicitFilterUtils::SetBulkRadiusForShapeFiltering, py::arg("input_model_part"))
        ;
}

}  // namespace Python.
} // Namespace Kratos


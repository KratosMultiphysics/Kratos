//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes
#include "includes/define_python.h"
#include "containers/model.h"
#include "includes/kratos_components.h"
#include "utilities/quaternion.h"
#include "python/add_model_part_to_python.h"
#include "python/containers_interface.h"

namespace Kratos::Python
{

template<class TDataType>
void AddNodalSolutionStepVariable(ModelPart& rModelPart, Variable<TDataType> const& rThisVariable)
{
    rModelPart.AddNodalSolutionStepVariable(rThisVariable);
}

template<class TDataType>
bool HasNodalSolutionStepVariable(ModelPart& rModelPart, Variable<TDataType> const& rThisVariable)
{
    return rModelPart.HasNodalSolutionStepVariable(rThisVariable);
}

/** Retrieve the variable names of the entities in the given container.
 *
 * Retrieve the variable names of the entities in `rContainer`. If the
 * `doFullSearch` is enabled, it will iterate and check all the entities
 * in the container. If not enabled it will be assumed that first entity of
 * the container is representative of the list of variables in every intenty
 */
template<class TContainerType>
const std::unordered_set<std::string> GetNonHistoricalVariablesNames(ModelPart& rModelPart, TContainerType& rContainer, bool doFullSearch=false) {

    std::unordered_set<std::string> variable_names;

    if(doFullSearch) {
        if(rContainer.size() == 0) {
            KRATOS_WARNING("DEBUG") << "Checking and empty container" << std::endl;
        } else {
            for(auto & variable: rContainer.begin()->GetData()) {
                variable_names.insert(variable.first->Name());
            }
        }
    } else {
        if(rContainer.size() == 0) {
            KRATOS_WARNING("DEBUG") << "Checking and empty container" << std::endl;
        }
        for(auto & entity : rContainer) {
            for(auto & variable: entity.GetData()) {
                variable_names.insert(variable.first->Name());
            }
        }
    }

    return variable_names;
}

class SubModelPartView {
public:
    using iterator = ModelPart::SubModelPartsContainerType::const_iterator;

    SubModelPartView(const ModelPart::SubModelPartsContainerType& rSubModelParts)
        : mrSubModelParts(rSubModelParts) {}

    iterator begin() const { return mrSubModelParts.begin(); }
    iterator end() const { return mrSubModelParts.end(); }

    std::size_t size() const { return mrSubModelParts.size(); }
private:
    const ModelPart::SubModelPartsContainerType& mrSubModelParts;
};

void AddModelPartToPython(pybind11::module& m)
{

    ModelPart::IndexType(ModelPart::*pointer_to_clone_time_step_1)(void) = &ModelPart::CloneTimeStep;
    ModelPart::IndexType(ModelPart::*pointer_to_clone_time_step_2)(double) = &ModelPart::CloneTimeStep;
    ProcessInfo::Pointer(ModelPart::*pointer_to_get_process_info)(void) = &ModelPart::pGetProcessInfo;
    void (ModelPart::*pointer_to_set_process_info)(ProcessInfo::Pointer) = &ModelPart::SetProcessInfo;
    // ModelPart::MeshType::Pointer (ModelPart::*pointer_to_get_mesh)() = &ModelPart::pGetMesh;
    //      std::string& (ModelPart::*pointer_to_name)(void) = &ModelPart::Name;


    namespace py = pybind11;

    py::class_<typename ModelPart::SubModelPartsContainerType >(m, "SubModelPartsContainerType")
        .def("__iter__", [](typename ModelPart::SubModelPartsContainerType& self){ return py::make_iterator(self.begin(), self.end());},  py::keep_alive<0,1>())
        ;

    MapInterface<ModelPart::GeometryContainerType>().CreateInterface(m,"GeometryContainerType");
    PointerVectorSetPythonInterface<ModelPart::MasterSlaveConstraintContainerType>().CreateInterface(m,"MasterSlaveConstraintsArray");

    py::class_<Kratos::Python::SubModelPartView>(m, "SubModelPartView")
        .def("__iter__", [](const Kratos::Python::SubModelPartView &self) {
            return py::make_iterator<py::return_value_policy::reference_internal,
                                    SubModelPartView::iterator, SubModelPartView::iterator>(
                self.begin(), self.end()); }, py::keep_alive<0, 1>())
        .def("__len__", &Kratos::Python::SubModelPartView::size)
        .def("__getitem__", [](const Kratos::Python::SubModelPartView &self, std::size_t i) -> const ModelPart& {
            if (i >= self.size())
                throw py::index_error();
            return *std::next(self.begin(), i); }, py::return_value_policy::reference_internal)
        ;

    py::class_<ModelPart, Kratos::shared_ptr<ModelPart>, DataValueContainer, Flags>(m, "ModelPart")
        .def("__copy__", [](const ModelPart &) {
        throw py::type_error("ModelPart cannot be copied");})
        .def("__deepcopy__", [](const ModelPart &, py::dict) {
            throw py::type_error("ModelPart cannot be deep-copied");})
        .def_property("Name",
            [](ModelPart const& rModelPart) { return rModelPart.Name(); },
            [](ModelPart& rModelPart, std::string const& NewName) { rModelPart.Name() = NewName; }
        )
        .def("FullName", &ModelPart::FullName)
        //  .def_property("ProcessInfo", GetProcessInfo, SetProcessInfo)
        .def_property("ProcessInfo", pointer_to_get_process_info, pointer_to_set_process_info)
        .def("Clear", &ModelPart::Clear)
        .def("CreateSolutionStep", &ModelPart::CreateSolutionStep)
        .def("CloneSolutionStep", &ModelPart::CloneSolutionStep)
        .def("CreateTimeStep", &ModelPart::CreateTimeStep)
        .def("ReduceTimeStep", &ModelPart::ReduceTimeStep)
        .def("CloneTimeStep", pointer_to_clone_time_step_1)
        .def("CloneTimeStep", pointer_to_clone_time_step_2)
        //       .def("CopySolutionStepData",&ModelPart::CopySolutionStepData)
        .def("NumberOfNodes", &ModelPart::NumberOfNodes)
        .def("NumberOfNodes", [](ModelPart& rModelPart){ return rModelPart.NumberOfNodes();})
        .def("SetBufferSize", &ModelPart::SetBufferSize)
        .def("GetBufferSize", &ModelPart::GetBufferSize)
        .def("NumberOfElements", &ModelPart::NumberOfElements)
        .def("NumberOfElements", [](ModelPart& rModelPart){ return rModelPart.NumberOfElements();})
        .def("NumberOfConditions", &ModelPart::NumberOfConditions)
        .def("NumberOfConditions", [](ModelPart& rModelPart){ return rModelPart.NumberOfConditions();})
        .def("NumberOfGeometries", &ModelPart::NumberOfGeometries)
        .def("NumberOfMasterSlaveConstraints", &ModelPart::NumberOfMasterSlaveConstraints)
        .def("NumberOfMasterSlaveConstraints", [](ModelPart& rModelPart){ return rModelPart.NumberOfMasterSlaveConstraints();})
        .def("NumberOfMeshes", &ModelPart::NumberOfMeshes)
        .def("NumberOfProperties", &ModelPart::NumberOfProperties, py::arg("ThisIndex") = 0)
        .def("GetMesh", [](ModelPart& rModelPart) { return rModelPart.pGetMesh(); })
        .def("GetMesh", [](ModelPart& rModelPart, ModelPart::IndexType MeshIndex) {
            ModelPart::IndexType number_of_meshes = rModelPart.NumberOfMeshes();
            // adding necessary meshes to the model part.
            ModelPart::MeshType empty_mesh;
            for(ModelPart::IndexType i = number_of_meshes ; i < MeshIndex + 1 ; i++)
                rModelPart.GetMeshes().push_back(Kratos::make_shared<ModelPart::MeshType>(empty_mesh.Clone()));
            return rModelPart.pGetMesh(MeshIndex);
        })
        .def_property("Nodes",
            [](ModelPart& rModelPart){ return rModelPart.pNodes(); },
            [](ModelPart& rModelPart, ModelPart::NodesContainerType::Pointer pOtherNodes){ rModelPart.SetNodes(pOtherNodes); }
        )
        .def("HasNode", [](ModelPart& rModelPart, ModelPart::IndexType NodeId){ return rModelPart.HasNode(NodeId); })
        .def("HasNode", [](ModelPart& rModelPart, ModelPart::IndexType NodeId, ModelPart::IndexType ThisIndex){ return rModelPart.HasNode(NodeId, ThisIndex); })
        .def("GetNode", [](ModelPart& rModelPart, ModelPart::IndexType NodeId){ return rModelPart.pGetNode(NodeId); }, py::return_value_policy::reference_internal)
        .def("GetNode", [](ModelPart& rModelPart, ModelPart::IndexType NodeId, ModelPart::IndexType ThisIndex){ return rModelPart.pGetNode(NodeId, ThisIndex); }, py::return_value_policy::reference_internal)
        .def("GetNodes", [](ModelPart& rModelPart){ return rModelPart.pNodes(); }, py::return_value_policy::reference_internal)
        .def("SetNodes", [](ModelPart& rModelPart, ModelPart::NodesContainerType::Pointer pOtherNodes){ rModelPart.SetNodes(pOtherNodes); })
        .def("GetNodes", [](ModelPart& rModelPart, ModelPart::IndexType ThisIndex){ return rModelPart.pNodes(ThisIndex); }, py::return_value_policy::reference_internal)
        .def("SetNodes", [](ModelPart& rModelPart, ModelPart::NodesContainerType::Pointer pOtherNodes, ModelPart::IndexType ThisIndex){ rModelPart.SetNodes(pOtherNodes, ThisIndex); })
        .def("AddNode", [](ModelPart& rModelPart, ModelPart::NodeType::Pointer pNode){ rModelPart.AddNode(pNode); })
        .def("AddNode", [](ModelPart& rModelPart, ModelPart::NodeType::Pointer pNode, ModelPart::IndexType MeshId){ rModelPart.AddNode(pNode, MeshId); })
        .def("RemoveNode", [](ModelPart& rModelPart, ModelPart::IndexType NodeId){ rModelPart.RemoveNode(NodeId); })
        .def("RemoveNode", [](ModelPart& rModelPart, ModelPart::IndexType NodeId, ModelPart::IndexType ThisIndex){ rModelPart.RemoveNode(NodeId, ThisIndex); })
        .def("RemoveNode", [](ModelPart& rModelPart, ModelPart::NodeType::Pointer pThisNode){ rModelPart.RemoveNode(pThisNode); })
        .def("RemoveNode", [](ModelPart& rModelPart, ModelPart::NodeType::Pointer pThisNode, ModelPart::IndexType ThisIndex){ rModelPart.RemoveNode(pThisNode, ThisIndex); })
        .def("RemoveNodes", &ModelPart::RemoveNodes) // This likely takes a flag or similar, keep as is or verify specific overload if it was a helper
        .def("RemoveNodeFromAllLevels", [](ModelPart& rModelPart, ModelPart::IndexType NodeId){ rModelPart.RemoveNodeFromAllLevels(NodeId); })
        .def("RemoveNodeFromAllLevels", [](ModelPart& rModelPart, ModelPart::IndexType NodeId, ModelPart::IndexType ThisIndex){ rModelPart.RemoveNodeFromAllLevels(NodeId, ThisIndex); })
        .def("RemoveNodeFromAllLevels", [](ModelPart& rModelPart, ModelPart::NodeType::Pointer pThisNode){ rModelPart.RemoveNodeFromAllLevels(pThisNode); })
        .def("RemoveNodeFromAllLevels", [](ModelPart& rModelPart, ModelPart::NodeType::Pointer pThisNode, ModelPart::IndexType ThisIndex){ rModelPart.RemoveNodeFromAllLevels(pThisNode, ThisIndex); })
        .def("RemoveNodesFromAllLevels", [](ModelPart& rModelPart, Flags identifier_flag){ rModelPart.RemoveNodesFromAllLevels(identifier_flag); })
        .def("NodesArray", &ModelPart::NodesArray, py::return_value_policy::reference_internal)
        .def("NumberOfTables", &ModelPart::NumberOfTables)
        .def("AddTable", &ModelPart::AddTable)
        .def("GetTable", &ModelPart::pGetTable)
        .def("HasProperties", [](ModelPart &rSelf, int Id)
             { return rSelf.HasProperties(Id); })
        .def("HasProperties", [](ModelPart &rSelf, const std::string &rAddress)
             { return rSelf.HasProperties(rAddress); })
        .def("RecursivelyHasProperties", [](ModelPart &rSelf, int Id)
             { return rSelf.RecursivelyHasProperties(Id); })
        .def("CreateNewProperties", [](ModelPart &rSelf, int Id)
             { return rSelf.CreateNewProperties(Id); })
        .def("GetProperties", [](ModelPart &rSelf, int Id)
             { return rSelf.pGetProperties(Id); })
        .def("GetProperties", [](ModelPart &rSelf, const std::string &rAddress)
             { return rSelf.pGetProperties(rAddress); })
        .def("GetProperties", [](ModelPart &rSelf)
             { return rSelf.pProperties(); })
        .def("AddProperties", [](ModelPart &rSelf, Properties::Pointer pProperties)
             { rSelf.AddProperties(pProperties); })
        .def("RemoveProperties", [](ModelPart &rSelf, int Id)
             { return rSelf.RemoveProperties(Id); })
        .def("RemoveProperties", [](ModelPart &rSelf, Properties::Pointer pProperties)
             { return rSelf.RemoveProperties(pProperties); })
        .def("RemovePropertiesFromAllLevels", [](ModelPart &rSelf, int Id)
             { return rSelf.RemovePropertiesFromAllLevels(Id); })
        .def("RemovePropertiesFromAllLevels", [](ModelPart &rSelf, Properties::Pointer pProperties)
             { return rSelf.RemovePropertiesFromAllLevels(pProperties); })
        .def_property("Properties",
            [](ModelPart& rModelPart){ return rModelPart.pProperties(); },
            [](ModelPart& rModelPart, ModelPart::PropertiesContainerType::Pointer pOtherProperties){ rModelPart.SetProperties(pOtherProperties); }
        )
        .def("SetProperties", // This is now redundant due to def_property, but kept if used elsewhere explicitly
            [](ModelPart& rModelPart, ModelPart::PropertiesContainerType::Pointer pOtherProperties){ rModelPart.SetProperties(pOtherProperties); }
        )
        .def("PropertiesArray", &ModelPart::PropertiesArray, py::return_value_policy::reference_internal)
        .def_property("Elements",
            [](ModelPart& rModelPart){ return rModelPart.pElements(); },
            [](ModelPart& rModelPart, ModelPart::ElementsContainerType::Pointer pOtherElements){ rModelPart.SetElements(pOtherElements); }
        )
        .def("HasElement", [](ModelPart& rModelPart, ModelPart::IndexType ElementId){ return rModelPart.HasElement(ElementId); })
        .def("HasElement", [](ModelPart& rModelPart, ModelPart::IndexType ElementId, ModelPart::IndexType ThisIndex){ return rModelPart.HasElement(ElementId, ThisIndex); })
        .def("GetElement", [](ModelPart& rModelPart, ModelPart::IndexType ElementId){ return rModelPart.pGetElement(ElementId); }, py::return_value_policy::reference_internal)
        .def("GetElement", [](ModelPart& rModelPart, ModelPart::IndexType ElementId, ModelPart::IndexType ThisIndex){ return rModelPart.pGetElement(ElementId, ThisIndex); }, py::return_value_policy::reference_internal)
        .def("GetElements", [](ModelPart& rModelPart){ return rModelPart.pElements(); }, py::return_value_policy::reference_internal)
        .def("SetElements", [](ModelPart& rModelPart, ModelPart::ElementsContainerType::Pointer pOtherElements){ rModelPart.SetElements(pOtherElements); })
        .def("GetElements", [](ModelPart& rModelPart, ModelPart::IndexType ThisIndex){ return rModelPart.pElements(ThisIndex); }, py::return_value_policy::reference_internal)
        .def("SetElements", [](ModelPart& rModelPart, ModelPart::ElementsContainerType::Pointer pOtherElements, ModelPart::IndexType ThisIndex){ rModelPart.SetElements(pOtherElements, ThisIndex); })
        .def("AddElement", [](ModelPart& rModelPart, ModelPart::ElementType::Pointer pElement){ rModelPart.AddElement(pElement); })
        .def("AddElement", [](ModelPart& rModelPart, ModelPart::ElementType::Pointer pElement, ModelPart::IndexType MeshId){ rModelPart.AddElement(pElement, MeshId); })
        .def("RemoveElement", [](ModelPart& rModelPart, ModelPart::IndexType ElementId){ rModelPart.RemoveElement(ElementId); })
        .def("RemoveElement", [](ModelPart& rModelPart, ModelPart::IndexType ElementId, ModelPart::IndexType ThisIndex){ rModelPart.RemoveElement(ElementId, ThisIndex); })
        .def("RemoveElement", [](ModelPart& rModelPart, ModelPart::ElementType::Pointer pThisElement){ rModelPart.RemoveElement(pThisElement); })
        .def("RemoveElement", [](ModelPart& rModelPart, ModelPart::ElementType::Pointer pThisElement, ModelPart::IndexType ThisIndex){ rModelPart.RemoveElement(pThisElement, ThisIndex); })
        .def("RemoveElements", &ModelPart::RemoveElements) // Keep direct member call
        .def("RemoveElementFromAllLevels", [](ModelPart& rModelPart, ModelPart::IndexType ElementId){ rModelPart.RemoveElementFromAllLevels(ElementId); })
        .def("RemoveElementFromAllLevels", [](ModelPart& rModelPart, ModelPart::IndexType ElementId, ModelPart::IndexType ThisIndex){ rModelPart.RemoveElementFromAllLevels(ElementId, ThisIndex); })
        .def("RemoveElementFromAllLevels", [](ModelPart& rModelPart, ModelPart::ElementType::Pointer pThisElement){ rModelPart.RemoveElementFromAllLevels(pThisElement); })
        .def("RemoveElementFromAllLevels", [](ModelPart& rModelPart, ModelPart::ElementType::Pointer pThisElement, ModelPart::IndexType ThisIndex){ rModelPart.RemoveElementFromAllLevels(pThisElement, ThisIndex); })
        .def("RemoveElementsFromAllLevels", [](ModelPart& rModelPart, Flags identifier_flag){ rModelPart.RemoveElementsFromAllLevels(identifier_flag); })
        .def("ElementsArray", &ModelPart::ElementsArray, py::return_value_policy::reference_internal)
        .def_property("Conditions",
            [](ModelPart& rModelPart){ return rModelPart.pConditions(); },
            [](ModelPart& rModelPart, ModelPart::ConditionsContainerType::Pointer pOtherConditions){ rModelPart.SetConditions(pOtherConditions); }
        )
        .def("HasCondition", [](ModelPart& rModelPart, ModelPart::IndexType ConditionId){ return rModelPart.HasCondition(ConditionId); })
        .def("HasCondition", [](ModelPart& rModelPart, ModelPart::IndexType ConditionId, ModelPart::IndexType ThisIndex){ return rModelPart.HasCondition(ConditionId, ThisIndex); })
        .def("GetCondition", [](ModelPart& rModelPart, ModelPart::IndexType ConditionId){ return rModelPart.pGetCondition(ConditionId); }, py::return_value_policy::reference_internal)
        .def("GetCondition", [](ModelPart& rModelPart, ModelPart::IndexType ConditionId, ModelPart::IndexType ThisIndex){ return rModelPart.pGetCondition(ConditionId, ThisIndex); }, py::return_value_policy::reference_internal)
        .def("GetConditions", [](ModelPart& rModelPart){ return rModelPart.pConditions(); }, py::return_value_policy::reference_internal)
        .def("SetConditions", [](ModelPart& rModelPart, ModelPart::ConditionsContainerType::Pointer pOtherConditions){ rModelPart.SetConditions(pOtherConditions); })
        .def("GetConditions", [](ModelPart& rModelPart, ModelPart::IndexType ThisIndex){ return rModelPart.pConditions(ThisIndex); }, py::return_value_policy::reference_internal)
        .def("SetConditions", [](ModelPart& rModelPart, ModelPart::ConditionsContainerType::Pointer pOtherConditions, ModelPart::IndexType ThisIndex){ rModelPart.SetConditions(pOtherConditions, ThisIndex); })
        .def("AddCondition", [](ModelPart& rModelPart, Condition::Pointer newCondition){ rModelPart.AddCondition(newCondition); })
        .def("AddCondition", [](ModelPart& rModelPart, Condition::Pointer newCondition, unsigned int ThisIndex){ rModelPart.AddCondition(newCondition, ThisIndex); })
        .def("RemoveCondition", [](ModelPart& rModelPart, ModelPart::IndexType ConditionId){ rModelPart.RemoveCondition(ConditionId); })
        .def("RemoveCondition", [](ModelPart& rModelPart, ModelPart::IndexType ConditionId, ModelPart::IndexType ThisIndex){ rModelPart.RemoveCondition(ConditionId, ThisIndex); })
        .def("RemoveCondition", [](ModelPart& rModelPart, ModelPart::ConditionType::Pointer pThisCondition){ rModelPart.RemoveCondition(pThisCondition); })
        .def("RemoveCondition", [](ModelPart& rModelPart, ModelPart::ConditionType::Pointer pThisCondition, ModelPart::IndexType ThisIndex){ rModelPart.RemoveCondition(pThisCondition, ThisIndex); })
        .def("RemoveConditions", &ModelPart::RemoveConditions) // Keep direct member call
        .def("RemoveConditionFromAllLevels", [](ModelPart& rModelPart, ModelPart::IndexType ConditionId){ rModelPart.RemoveConditionFromAllLevels(ConditionId); })
        .def("RemoveConditionFromAllLevels", [](ModelPart& rModelPart, ModelPart::IndexType ConditionId, ModelPart::IndexType ThisIndex){ rModelPart.RemoveConditionFromAllLevels(ConditionId, ThisIndex); })
        .def("RemoveConditionFromAllLevels", [](ModelPart& rModelPart, ModelPart::ConditionType::Pointer pThisCondition){ rModelPart.RemoveConditionFromAllLevels(pThisCondition); })
        .def("RemoveConditionFromAllLevels", [](ModelPart& rModelPart, ModelPart::ConditionType::Pointer pThisCondition, ModelPart::IndexType ThisIndex){ rModelPart.RemoveConditionFromAllLevels(pThisCondition, ThisIndex); })
        .def("RemoveConditionsFromAllLevels", [](ModelPart& rModelPart, Flags identifier_flag){ rModelPart.RemoveConditionsFromAllLevels(identifier_flag); })
        .def("AddGeometry", [](ModelPart& rModelPart, ModelPart::GeometryType::Pointer pNewGeometry){ rModelPart.AddGeometry(pNewGeometry); })
        .def("GetGeometry", [](ModelPart& rModelPart, ModelPart::IndexType GeometryId){ return rModelPart.pGetGeometry(GeometryId); }, py::return_value_policy::reference_internal)
        .def("GetGeometry", [](ModelPart& rModelPart, const std::string& GeometryName){ return rModelPart.pGetGeometry(GeometryName); }, py::return_value_policy::reference_internal)
        .def("HasGeometry", [](ModelPart& rModelPart, ModelPart::IndexType GeometryId){ return rModelPart.HasGeometry(GeometryId); })
        .def("HasGeometry", [](ModelPart& rModelPart, const std::string& GeometryName){ return rModelPart.HasGeometry(GeometryName); })
        .def("RemoveGeometry", [](ModelPart& rModelPart, ModelPart::IndexType GeometryId){ rModelPart.RemoveGeometry(GeometryId); })
        .def("RemoveGeometry", [](ModelPart& rModelPart, const std::string& GeometryName){ rModelPart.RemoveGeometry(GeometryName); })
        .def("RemoveGeometryFromAllLevels", [](ModelPart& rModelPart, ModelPart::IndexType GeometryId){ rModelPart.RemoveGeometryFromAllLevels(GeometryId); })
        .def("RemoveGeometryFromAllLevels", [](ModelPart& rModelPart, const std::string& GeometryName){ rModelPart.RemoveGeometryFromAllLevels(GeometryName); })
        .def_property("Geometries", [](ModelPart &self)
                      { return self.Geometries(); }, [](ModelPart &self, ModelPart::GeometryContainerType &geometries)
                      { KRATOS_ERROR << "Setting geometries is not allowed! Trying to set value of ModelPart::Geometries."; })
        .def("CreateSubModelPart", &ModelPart::CreateSubModelPart, py::return_value_policy::reference_internal)
        .def("NumberOfSubModelParts", &ModelPart::NumberOfSubModelParts)
        .def("GetSubModelPart", py::overload_cast<const std::string &>(&ModelPart::GetSubModelPart), py::return_value_policy::reference_internal) // non-const version
        .def("RemoveSubModelPart", [](ModelPart& rModelPart, std::string const& ThisSubModelPartName){ rModelPart.RemoveSubModelPart(ThisSubModelPartName); })
        .def("RemoveSubModelPart", [](ModelPart& rModelPart, ModelPart& ThisSubModelPart){ rModelPart.RemoveSubModelPart(ThisSubModelPart); })
        .def("HasSubModelPart", &ModelPart::HasSubModelPart)
        .def("GetSubModelPartNames", &ModelPart::GetSubModelPartNames)
        .def("ConditionsArray", &ModelPart::ConditionsArray, py::return_value_policy::reference_internal)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<bool>)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<int>)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<double>)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<array_1d<double, 3>>)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<Vector>)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<Matrix>)
        .def("AddNodalSolutionStepVariable", AddNodalSolutionStepVariable<Quaternion<double>>)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<bool>)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<int>)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<double>)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<array_1d<double, 3>>)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<Vector>)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<Matrix>)
        .def("HasNodalSolutionStepVariable", HasNodalSolutionStepVariable<Quaternion<double>>)
        .def("GetNodalSolutionStepDataSize", &ModelPart::GetNodalSolutionStepDataSize)
        .def("GetNodalSolutionStepTotalDataSize", &ModelPart::GetNodalSolutionStepTotalDataSize)
        .def("OverwriteSolutionStepData", &ModelPart::OverwriteSolutionStepData)
        .def("CreateNewNode", [](ModelPart& rModelPart, int Id, double x, double y, double z){ return rModelPart.CreateNewNode(Id, x, y, z); }, py::return_value_policy::reference_internal)
        .def("CreateNewGeometry", [](ModelPart& rModelPart, const std::string& GeometryTypeName, std::vector<ModelPart::IndexType>& NodeIdList) {
            Geometry<Node>::PointsArrayType pGeometryNodeList;
            for (std::size_t i = 0; i < NodeIdList.size(); i++) {
                pGeometryNodeList.push_back(rModelPart.pGetNode(NodeIdList[i]));
            }
            return rModelPart.CreateNewGeometry(GeometryTypeName, pGeometryNodeList);
        }, py::return_value_policy::reference_internal)
        .def("CreateNewGeometry", [](ModelPart& rModelPart, const std::string& GeometryTypeName, ModelPart::IndexType GeometryId, std::vector<ModelPart::IndexType>& NodeIdList) {
            Geometry<Node>::PointsArrayType pGeometryNodeList;
            for(std::size_t i = 0; i < NodeIdList.size(); i++) {
                pGeometryNodeList.push_back(rModelPart.pGetNode(NodeIdList[i]));
            }
            return rModelPart.CreateNewGeometry(GeometryTypeName, GeometryId, pGeometryNodeList);
        }, py::return_value_policy::reference_internal)
        .def("CreateNewGeometry", [](ModelPart& rModelPart, const std::string& GeometryTypeName, const std::string& GeometryIdentifierName, std::vector<ModelPart::IndexType>& NodeIdList) {
            Geometry<Node>::PointsArrayType pGeometryNodeList;
            for (std::size_t i = 0; i < NodeIdList.size(); i++) {
                pGeometryNodeList.push_back(rModelPart.pGetNode(NodeIdList[i]));
            }
            return rModelPart.CreateNewGeometry(GeometryTypeName, GeometryIdentifierName, pGeometryNodeList);
        }, py::return_value_policy::reference_internal)
        .def("CreateNewGeometry", [](ModelPart& rModelPart, const std::string& GeometryTypeName, ModelPart::GeometryType::Pointer pGeometry){
            return rModelPart.CreateNewGeometry(GeometryTypeName, pGeometry);
        }, py::return_value_policy::reference_internal)
        .def("CreateNewGeometry", [](ModelPart& rModelPart, const std::string& GeometryTypeName, ModelPart::IndexType GeometryId, ModelPart::GeometryType::Pointer pGeometry){
            return rModelPart.CreateNewGeometry(GeometryTypeName, GeometryId, pGeometry);
        }, py::return_value_policy::reference_internal)
        .def("CreateNewGeometry", [](ModelPart& rModelPart, const std::string& GeometryTypeName, const std::string& GeometryIdentifierName, ModelPart::GeometryType::Pointer pGeometry){
            return rModelPart.CreateNewGeometry(GeometryTypeName, GeometryIdentifierName, pGeometry);
        }, py::return_value_policy::reference_internal)
        .def("CreateNewElement", [](ModelPart& rModelPart, const std::string ElementName, ModelPart::IndexType Id, std::vector<ModelPart::IndexType>& NodeIdList, ModelPart::PropertiesType::Pointer pProperties) {
            if (!KratosComponents<Element>::Has(ElementName)) {
                std::stringstream msg;
                KratosComponents<Element> instance;
                instance.PrintData(msg);
                KRATOS_ERROR << "The Element \"" << ElementName << "\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following Elements are registered:\n" << msg.str() << std::endl;
            }
            Geometry<Node>::PointsArrayType pElementNodeList;
            for(unsigned int i = 0; i < NodeIdList.size(); i++) {
                pElementNodeList.push_back(rModelPart.pGetNode(NodeIdList[i]));
            }
            return rModelPart.CreateNewElement(ElementName, Id, pElementNodeList, pProperties);
        }, py::return_value_policy::reference_internal)
        .def("CreateNewElement", [](ModelPart &rModelPart, const std::string &ElementName, ModelPart::IndexType Id, ModelPart::GeometryType::Pointer pGeometry, ModelPart::PropertiesType::Pointer pProperties)
             { return rModelPart.CreateNewElement(ElementName, Id, pGeometry, pProperties); }) // This one was already a lambda
        .def("CreateNewCondition", [](ModelPart& rModelPart, const std::string ConditionName, ModelPart::IndexType Id, std::vector<ModelPart::IndexType>& NodeIdList, ModelPart::PropertiesType::Pointer pProperties) {
            if (!KratosComponents<Condition>::Has(ConditionName)) {
                std::stringstream msg;
                KratosComponents<Condition> instance;
                instance.PrintData(msg);
                KRATOS_ERROR << "The Condition \"" << ConditionName << "\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following Conditions are registered:\n" << msg.str() << std::endl;
            }
            Geometry<Node>::PointsArrayType pConditionNodeList;
            for(unsigned int i = 0; i <NodeIdList.size(); i++) {
                pConditionNodeList.push_back(rModelPart.pGetNode(NodeIdList[i]));
            }
            return rModelPart.CreateNewCondition(ConditionName, Id, pConditionNodeList, pProperties);
        }, py::return_value_policy::reference_internal)
        .def("CreateNewCondition", [](ModelPart &rModelPart, const std::string &ConditionName, ModelPart::IndexType Id, ModelPart::GeometryType::Pointer pGeometry, ModelPart::PropertiesType::Pointer pProperties)
              { return rModelPart.CreateNewCondition(ConditionName, Id, pGeometry, pProperties); }) // This one was already a lambda
        .def("GetCommunicator", [](ModelPart& rModelPart) -> Communicator& { return rModelPart.GetCommunicator(); }, py::return_value_policy::reference_internal)
        .def("Check", &ModelPart::Check)
        .def("IsSubModelPart", &ModelPart::IsSubModelPart)
        .def("IsDistributed", &ModelPart::IsDistributed)
        .def("AddNodes", [](ModelPart& rModelPart, std::vector<ModelPart::IndexType>& NodesIds){ rModelPart.AddNodes(NodesIds); })
        .def("AddConditions", [](ModelPart& rModelPart, std::vector<ModelPart::IndexType>& ConditionsIds){ rModelPart.AddConditions(ConditionsIds); })
        .def("AddElements", [](ModelPart& rModelPart, std::vector<ModelPart::IndexType>& ElementsIds){ rModelPart.AddElements(ElementsIds); })
        .def("AddGeometries", [](ModelPart &rModelPart, std::vector<ModelPart::IndexType>& rGeometriesIds) {rModelPart.AddGeometries(rGeometriesIds);})
        .def("GetParentModelPart", [](ModelPart &self) -> ModelPart &
             { return self.GetParentModelPart(); }, py::return_value_policy::reference_internal)
        .def("GetRootModelPart", [](ModelPart &self) -> ModelPart &
             { return self.GetRootModelPart(); }, py::return_value_policy::reference_internal)
        .def("GetModel", [](ModelPart &self) -> Model &
             { return self.GetModel(); }, py::return_value_policy::reference_internal)
        .def_property_readonly("SubModelParts", [](ModelPart &self) 
            { return Kratos::Python::SubModelPartView(self.SubModelParts()); }, py::keep_alive<0, 1>())  // Keep self alive as long as SubModelPartView exists
        .def_property_readonly("MasterSlaveConstraints", [](ModelPart& rModelPart) -> const ModelPart::MasterSlaveConstraintContainerType& { return rModelPart.MasterSlaveConstraints(); })
        .def("GetHistoricalVariablesNames", [](ModelPart &rModelPart) -> std::unordered_set<std::string>
             {
            std::unordered_set<std::string> variable_names;
            for(auto & variable: rModelPart.GetNodalSolutionStepVariablesList()) {
                variable_names.insert(variable.Name());
            }
            return variable_names; })
        .def("GetNonHistoricalVariablesNames", [](ModelPart &rModelPart, ModelPart::NodesContainerType &rContainer, bool doFullSearch = false) -> std::unordered_set<std::string> {
            // Using the restored templated function
            return GetNonHistoricalVariablesNames(rModelPart, rContainer, doFullSearch);
        }, py::arg("container"), py::arg("doFullSearch") = false)
        .def("GetNonHistoricalVariablesNames", [](ModelPart &rModelPart, ModelPart::ElementsContainerType &rContainer, bool doFullSearch = false) -> std::unordered_set<std::string> {
            // Using the restored templated function
            return GetNonHistoricalVariablesNames(rModelPart, rContainer, doFullSearch);
        }, py::arg("container"), py::arg("doFullSearch") = false)
        .def("GetNonHistoricalVariablesNames", [](ModelPart &rModelPart, ModelPart::ConditionsContainerType &rContainer, bool doFullSearch = false) -> std::unordered_set<std::string> {
            // Using the restored templated function
            return GetNonHistoricalVariablesNames(rModelPart, rContainer, doFullSearch);
        }, py::arg("container"), py::arg("doFullSearch") = false)
        .def("HasMasterSlaveConstraint", [](ModelPart &rModelPart, ModelPart::IndexType MasterSlaveConstraintId) -> bool
             { return rModelPart.HasMasterSlaveConstraint(MasterSlaveConstraintId); })
        .def("GetMasterSlaveConstraint", [](ModelPart& rModelPart, ModelPart::IndexType MasterSlaveConstraintId){ return rModelPart.pGetMasterSlaveConstraint(MasterSlaveConstraintId); }, py::return_value_policy::reference_internal)
        .def("GetMasterSlaveConstraints", [](ModelPart& rModelPart) -> const ModelPart::MasterSlaveConstraintContainerType& { return rModelPart.MasterSlaveConstraints(); }) // Matches property
        .def("RemoveMasterSlaveConstraint", [](ModelPart& rModelPart, ModelPart::IndexType MasterSlaveConstraintId){ rModelPart.RemoveMasterSlaveConstraint(MasterSlaveConstraintId); })
        .def("RemoveMasterSlaveConstraint", [](ModelPart& rModelPart, ModelPart::MasterSlaveConstraintType& rOtherMasterSlaveConstraint){ rModelPart.RemoveMasterSlaveConstraint(rOtherMasterSlaveConstraint); })
        .def("RemoveMasterSlaveConstraintFromAllLevels", [](ModelPart& rModelPart, ModelPart::IndexType MasterSlaveConstraintId){ rModelPart.RemoveMasterSlaveConstraintFromAllLevels(MasterSlaveConstraintId); })
        .def("RemoveMasterSlaveConstraintFromAllLevels", [](ModelPart& rModelPart, ModelPart::MasterSlaveConstraintType& rMasterSlaveConstraint){ rModelPart.RemoveMasterSlaveConstraintFromAllLevels(rMasterSlaveConstraint); })
        .def("RemoveMasterSlaveConstraints", &ModelPart::RemoveMasterSlaveConstraints) // Keep direct member call
        .def("RemoveMasterSlaveConstraintsFromAllLevels", &ModelPart::RemoveMasterSlaveConstraintsFromAllLevels) // Keep direct member call
        .def("AddMasterSlaveConstraint", [](ModelPart& rModelPart, ModelPart::MasterSlaveConstraintType::Pointer pConstraint){ rModelPart.AddMasterSlaveConstraint(pConstraint); })
        .def("AddMasterSlaveConstraints", [](ModelPart& rModelPart, std::vector<ModelPart::IndexType>& ConstraintIds){ rModelPart.AddMasterSlaveConstraints(ConstraintIds); })
        .def("CreateNewMasterSlaveConstraint", [](ModelPart& rModelPart, std::string ConstraintName, ModelPart::IndexType Id, ModelPart::DofsVectorType& rMasterDofsVector, ModelPart::DofsVectorType& rSlaveDofsVector, ModelPart::MatrixType RelationMatrix, ModelPart::VectorType ConstantVector) {
            if (!KratosComponents<MasterSlaveConstraint>::Has(ConstraintName)) {
                std::stringstream msg;
                KratosComponents<MasterSlaveConstraint> instance;
                instance.PrintData(msg);
                KRATOS_ERROR << "The Constraint \"" << ConstraintName << "\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following Constraints are registered:\n" << msg.str() << std::endl;
            }
            return rModelPart.CreateNewMasterSlaveConstraint(ConstraintName, Id, rMasterDofsVector, rSlaveDofsVector, RelationMatrix, ConstantVector);
        }, py::return_value_policy::reference_internal)
        .def("CreateNewMasterSlaveConstraint", [](ModelPart& rModelPart, std::string ConstraintName, ModelPart::IndexType Id, ModelPart::NodeType& rMasterNode, ModelPart::DoubleVariableType& rMasterVariable, ModelPart::NodeType& rSlaveNode, ModelPart::DoubleVariableType& rSlaveVariable, double Weight, double Constant) {
            if (!KratosComponents<MasterSlaveConstraint>::Has(ConstraintName)) {
                std::stringstream msg;
                KratosComponents<MasterSlaveConstraint> instance;
                instance.PrintData(msg);
                KRATOS_ERROR << "The Constraint \"" << ConstraintName << "\" is not registered!\nMaybe you need to import the application where it is defined?\nThe following Constraints are registered:\n" << msg.str() << std::endl;
            }
            return rModelPart.CreateNewMasterSlaveConstraint(ConstraintName, Id, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
        }, py::return_value_policy::reference_internal)
        .def("__str__", PrintObject<ModelPart>)
        ;
}

} // namespace Kratos::Python.
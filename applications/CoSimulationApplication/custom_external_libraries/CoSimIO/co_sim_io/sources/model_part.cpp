//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes
#include <algorithm>

// Project includes
#include "includes/model_part.hpp"
#include "includes/utilities.hpp"

namespace CoSimIO {

Node::Node(
    const IdType I_Id,
    const double I_X,
    const double I_Y,
    const double I_Z)
: mId(I_Id),
    mX(I_X),
    mY(I_Y),
    mZ(I_Z)
{
    CO_SIM_IO_ERROR_IF(I_Id < 1) << "Id must be >= 1!" << std::endl;
}

void Node::Print(std::ostream& rOStream) const
{
    rOStream << "CoSimIO-Node; Id: " << Id() << "\n";
    rOStream << "    Coordinates: [ " << X() << " | " << Y() << " | " << Z() << " ]" << std::endl;
}

void Node::save(CoSimIO::Internals::Serializer& rSerializer) const
{
    rSerializer.save("mId", mId);
    rSerializer.save("mX", mX);
    rSerializer.save("mY", mY);
    rSerializer.save("mZ", mZ);
}

void Node::load(CoSimIO::Internals::Serializer& rSerializer)
{
    rSerializer.load("mId", mId);
    rSerializer.load("mX", mX);
    rSerializer.load("mY", mY);
    rSerializer.load("mZ", mZ);
}

Element::Element(
    const IdType I_Id,
    const ElementType I_Type,
    const NodesContainerType& I_Nodes)
    : mId(I_Id),
        mType(I_Type),
        mNodes(I_Nodes)
{
    CO_SIM_IO_ERROR_IF(I_Id < 1) << "Id must be >= 1!" << std::endl;
    CO_SIM_IO_ERROR_IF(NumberOfNodes() < 1) << "No nodes were passed!" << std::endl;
    const int num_nodes_elem_type = Utilities::GetNumberOfNodesForElementType(I_Type);
    CO_SIM_IO_ERROR_IF_NOT(num_nodes_elem_type == static_cast<int>(NumberOfNodes())) << "Number of nodes (" << NumberOfNodes() << ") does not match expected number for element type (" << num_nodes_elem_type << ")!" << std::endl;
}


void Element::Print(std::ostream& rOStream) const
{
    rOStream << "CoSimIO-Element; Id: " << Id() << "\n";
    rOStream << "    Number of Nodes: " << NumberOfNodes() << "\n";
    rOStream << "    Node Ids: ";
    if (NumberOfNodes() > 0) {
        rOStream << mNodes[0]->Id();
    }
    for (std::size_t i=1; i<NumberOfNodes(); ++i) {
        rOStream << ", " << mNodes[i]->Id();
    }
    rOStream << std::endl;
}

void Element::save(CoSimIO::Internals::Serializer& rSerializer) const
{
    rSerializer.save("mId", mId);
    rSerializer.save("mType", static_cast<int>(mType));
    rSerializer.save("mNodes", mNodes);
}

void Element::load(CoSimIO::Internals::Serializer& rSerializer)
{
    rSerializer.load("mId", mId);
    int tmp;
    rSerializer.load("mType", tmp);
    mType = static_cast<ElementType>(tmp);
    rSerializer.load("mNodes", mNodes);
}


ModelPart::ModelPart(const std::string& I_Name) : ModelPart(I_Name, true) {}

ModelPart::ModelPart(const std::string& I_Name, const bool InitInternalModelParts) : mName(I_Name)
{
    CO_SIM_IO_ERROR_IF(I_Name.empty()) << "Please don't use empty names (\"\") when creating a ModelPart" << std::endl;
    CO_SIM_IO_ERROR_IF_NOT(I_Name.find(".") == std::string::npos) << "Please don't use names containing (\".\") when creating a ModelPart (used in \"" << I_Name << "\")" << std::endl;

    if (InitInternalModelParts) {
        InitializeInternalModelParts();
    }
}

std::size_t ModelPart::NumberOfNodes() const
{
    return mNodes.size();
}

std::size_t ModelPart::NumberOfLocalNodes() const
{
    return GetLocalModelPart().NumberOfNodes();
}

std::size_t ModelPart::NumberOfGhostNodes() const
{
    return GetGhostModelPart().NumberOfNodes();
}


Node& ModelPart::CreateNewNode(
    const IdType I_Id,
    const double I_X,
    const double I_Y,
    const double I_Z)
{
    CO_SIM_IO_ERROR_IF(HasNode(I_Id)) << "The Node with Id " << I_Id << " exists already!" << std::endl;

    CoSimIO::intrusive_ptr<Node> new_node(CoSimIO::make_intrusive<Node>(I_Id, I_X, I_Y, I_Z));

    mNodes.push_back(new_node);
    GetLocalModelPart().mNodes.push_back(new_node);

    return *new_node;
}

Node& ModelPart::CreateNewGhostNode(
    const IdType I_Id,
    const double I_X,
    const double I_Y,
    const double I_Z,
    const int PartitionIndex)
{
    CO_SIM_IO_ERROR_IF(HasNode(I_Id)) << "The Node with Id " << I_Id << " exists already!" << std::endl;

    CoSimIO::intrusive_ptr<Node> new_node(CoSimIO::make_intrusive<Node>(I_Id, I_X, I_Y, I_Z));

    mNodes.push_back(new_node);
    GetGhostModelPart().mNodes.push_back(new_node);
    GetPartitionModelPart(PartitionIndex).mNodes.push_back(new_node);

    return *new_node;
}

Element& ModelPart::CreateNewElement(
    const IdType I_Id,
    const ElementType I_Type,
    const ConnectivitiesType& I_Connectivities)
{
    CO_SIM_IO_ERROR_IF(HasElement(I_Id)) << "The Element with Id " << I_Id << " exists already!" << std::endl;

    Element::NodesContainerType nodes;
    nodes.reserve(I_Connectivities.size());
    for (const IdType node_id : I_Connectivities) {
        nodes.push_back(pGetNode(node_id));
    }

    CoSimIO::intrusive_ptr<Element> new_element(CoSimIO::make_intrusive<Element>(I_Id, I_Type, nodes));

    mElements.push_back(new_element);
    GetLocalModelPart().mElements.push_back(new_element);

    return *new_element;
}

Node& ModelPart::GetNode(const IdType I_Id)
{
    auto it_node = FindNode(I_Id);
    CO_SIM_IO_ERROR_IF(it_node == mNodes.end()) << "Node with Id " << I_Id << " does not exist!" << std::endl;
    return **it_node;
}

const Node& ModelPart::GetNode(const IdType I_Id) const
{
    auto it_node = FindNode(I_Id);
    CO_SIM_IO_ERROR_IF(it_node == mNodes.end()) << "Node with Id " << I_Id << " does not exist!" << std::endl;
    return **it_node;
}

ModelPart::NodePointerType ModelPart::pGetNode(const IdType I_Id)
{
    auto it_node = FindNode(I_Id);
    CO_SIM_IO_ERROR_IF(it_node == mNodes.end()) << "Node with Id " << I_Id << " does not exist!" << std::endl;
    return *it_node;
}

Element& ModelPart::GetElement(const IdType I_Id)
{
    auto it_elem = FindElement(I_Id);
    CO_SIM_IO_ERROR_IF(it_elem == mElements.end()) << "Element with Id " << I_Id << " does not exist!" << std::endl;
    return **it_elem;
}

const Element& ModelPart::GetElement(const IdType I_Id) const
{
    auto it_elem = FindElement(I_Id);
    CO_SIM_IO_ERROR_IF(it_elem == mElements.end()) << "Element with Id " << I_Id << " does not exist!" << std::endl;
    return **it_elem;
}

ModelPart::ElementPointerType ModelPart::pGetElement(const IdType I_Id)
{
    auto it_elem = FindElement(I_Id);
    CO_SIM_IO_ERROR_IF(it_elem == mElements.end()) << "Element with Id " << I_Id << " does not exist!" << std::endl;
    return *it_elem;
}

void ModelPart::Print(std::ostream& rOStream) const
{
    rOStream << "CoSimIO-ModelPart \"" << mName << "\"\n";
    rOStream << "    Number of Nodes: " << NumberOfNodes() << "\n";
    rOStream << "    Number of Elements: " << NumberOfElements() << std::endl;
}

void ModelPart::Clear()
{
    InitializeInternalModelParts();

    mPartitionModelParts.clear();

    mElements.clear();
    mElements.shrink_to_fit();

    mNodes.clear();
    mNodes.shrink_to_fit();
}

ModelPart::NodesContainerType::const_iterator ModelPart::FindNode(const IdType I_Id) const
{
    return std::find_if(
        mNodes.begin(), mNodes.end(),
        [I_Id](const NodePointerType& rp_node) { return rp_node->Id() == I_Id;});
}

ModelPart::NodesContainerType::iterator ModelPart::FindNode(const IdType I_Id)
{
    return std::find_if(
        mNodes.begin(), mNodes.end(),
        [I_Id](const NodePointerType& rp_node) { return rp_node->Id() == I_Id;});
}

ModelPart::ElementsContainerType::const_iterator ModelPart::FindElement(const IdType I_Id) const
{
    return std::find_if(
        mElements.begin(), mElements.end(),
        [I_Id](const ElementPointerType& rp_elem) { return rp_elem->Id() == I_Id;});
}

ModelPart::ElementsContainerType::iterator ModelPart::FindElement(const IdType I_Id)
{
    return std::find_if(
        mElements.begin(), mElements.end(),
        [I_Id](const ElementPointerType& rp_elem) { return rp_elem->Id() == I_Id;});
}

bool ModelPart::HasNode(const IdType I_Id) const
{
    return FindNode(I_Id) != mNodes.end();
}

bool ModelPart::HasElement(const IdType I_Id) const
{
    return FindElement(I_Id) != mElements.end();
}

ModelPart& ModelPart::GetLocalModelPart()
{
    CO_SIM_IO_ERROR_IF_NOT(mpLocalModelPart) << "Internal ModelPart, access is not allowed!" << std::endl;
    return *mpLocalModelPart;
}

const ModelPart& ModelPart::GetLocalModelPart() const
{
    CO_SIM_IO_ERROR_IF_NOT(mpLocalModelPart) << "Internal ModelPart, access is not allowed!" << std::endl;
    return *mpLocalModelPart;
}

ModelPart& ModelPart::GetGhostModelPart()
{
    CO_SIM_IO_ERROR_IF_NOT(mpGhostModelPart) << "Internal ModelPart, access is not allowed!" << std::endl;
    return *mpGhostModelPart;
}

const ModelPart& ModelPart::GetGhostModelPart() const
{
    CO_SIM_IO_ERROR_IF_NOT(mpGhostModelPart) << "Internal ModelPart, access is not allowed!" << std::endl;
    return *mpGhostModelPart;
}

ModelPart& ModelPart::GetPartitionModelPart(const int PartitionIndex)
{
    ModelPart* p_partition_model_part;
    auto mp_iter = mPartitionModelParts.find(PartitionIndex);
    if (mp_iter != mPartitionModelParts.end()) { // a ModelPart for this partition exists already
        p_partition_model_part = mp_iter->second.get();
    } else {
        auto partition_mp = CoSimIO::make_unique<ModelPart>(std::to_string(PartitionIndex));
        p_partition_model_part = partition_mp.get();
        mPartitionModelParts[PartitionIndex] = std::move(partition_mp);
    }

    return *p_partition_model_part;
}

const ModelPart& ModelPart::GetPartitionModelPart(const int PartitionIndex) const
{
    ModelPart* p_partition_model_part;
    auto mp_iter = mPartitionModelParts.find(PartitionIndex);
    if (mp_iter != mPartitionModelParts.end()) { // a ModelPart for this partition exists already
        p_partition_model_part = mp_iter->second.get();
    } else {
        CO_SIM_IO_ERROR << "No ModelPart exists for partition index " << PartitionIndex << " and cannot be created in a const function!" << std::endl;
    }

    return *p_partition_model_part;
}

void ModelPart::InitializeInternalModelParts()
{
    mpLocalModelPart = std::unique_ptr<ModelPart>(new ModelPart("local", false));
    mpGhostModelPart = std::unique_ptr<ModelPart>(new ModelPart("ghost", false));
}

void ModelPart::save(CoSimIO::Internals::Serializer& rSerializer) const
{
    rSerializer.save("mName", mName);
    rSerializer.save("mNodes", mNodes);
    rSerializer.save("mElements", mElements);
    rSerializer.save("mpLocalModelPart", mpLocalModelPart);
    rSerializer.save("mpGhostModelPart", mpGhostModelPart);
    rSerializer.save("mPartitionModelParts", mPartitionModelParts);
}

void ModelPart::load(CoSimIO::Internals::Serializer& rSerializer)
{
    rSerializer.load("mName", mName);
    rSerializer.load("mNodes", mNodes);
    rSerializer.load("mElements", mElements);
    rSerializer.load("mpLocalModelPart", mpLocalModelPart);
    rSerializer.load("mpGhostModelPart", mpGhostModelPart);
    rSerializer.load("mPartitionModelParts", mPartitionModelParts);
}

} //namespace CoSimIO

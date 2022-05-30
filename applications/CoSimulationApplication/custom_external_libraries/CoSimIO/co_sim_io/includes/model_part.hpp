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

#ifndef CO_SIM_IO_MODEL_PART_INCLUDED
#define CO_SIM_IO_MODEL_PART_INCLUDED

/* This file contains the implementation of th  CoSimIO::ModelPart
It serves as a data container when exchanging data
Also it is used in order to be consistent with Kratos to reduce compatibility problems
This is a simplified version of Kratos::ModelPart
see https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/includes/model_part.h
*/

// System includes
#include <vector>
#include <string>
#include <unordered_map>
#include <atomic>
#include <ostream>

// Project includes
#include "define.hpp"
#include "serializer.hpp"

namespace CoSimIO {

namespace Internals {

template<class TDataType>
class PointerVector
{
public:

    using ContainerType = std::vector<TDataType>;
    using BaseType = typename TDataType::element_type; // to be used with smart pointers
    using const_iterator_type = typename ContainerType::const_iterator;

    struct const_iterator_adaptor
    {
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = TDataType;
        using pointer           = TDataType*;
        using reference         = TDataType&;

        const_iterator_adaptor(const_iterator_type ptr) : m_ptr(ptr) {}

        const BaseType& operator*() const { return **m_ptr; }
        // pointer operator->() { return m_ptr; }
        const_iterator_adaptor& operator++() { m_ptr++; return *this; }
        const_iterator_adaptor operator++(int) { const_iterator_adaptor tmp = *this; ++(*this); return tmp; }
        friend bool operator== (const const_iterator_adaptor& a, const const_iterator_adaptor& b) { return a.m_ptr == b.m_ptr; };
        friend bool operator!= (const const_iterator_adaptor& a, const const_iterator_adaptor& b) { return a.m_ptr != b.m_ptr; };

    private:
        const_iterator_type m_ptr;
    };

    PointerVector(const ContainerType& rPointerVector) : mPointerVector(rPointerVector) {}

    const_iterator_adaptor begin() const {return const_iterator_adaptor(mPointerVector.begin());}
    const_iterator_adaptor end()   const {return const_iterator_adaptor(mPointerVector.end());}

private:
    const ContainerType& mPointerVector;
};

} //namespace Internals

class CO_SIM_IO_API Node
{
public:
    Node(
        const IdType I_Id,
        const double I_X,
        const double I_Y,
        const double I_Z);

    Node(
        const IdType I_Id,
        const CoordinatesType& I_Coordinates)
    : Node(I_Id, I_Coordinates[0], I_Coordinates[1], I_Coordinates[2])
    { }

    // delete copy and assignment CTor
    Node(const Node&) = delete;
    Node& operator=(Node const&) = delete;

    IdType Id() const { return mId; }
    double X() const { return mX; }
    double Y() const { return mY; }
    double Z() const { return mZ; }
    CoordinatesType Coordinates() const { return {mX, mY, mZ}; }

    void Print(std::ostream& rOStream) const;

private:
    IdType mId;
    double mX;
    double mY;
    double mZ;

    //*********************************************
    //this block is needed for refcounting
    mutable std::atomic<int> mReferenceCounter{0};

    friend void intrusive_ptr_add_ref(const Node* x)
    {
        x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
    }

    friend void intrusive_ptr_release(const Node* x)
    {
        if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
            std::atomic_thread_fence(std::memory_order_acquire);
            delete x;
        }
    }
    //*********************************************

    Node() = default; // needed for Serializer

    friend class CoSimIO::Internals::Serializer; // needs "CoSimIO::Internals::" because it is in different namespace

    void save(CoSimIO::Internals::Serializer& rSerializer) const;

    void load(CoSimIO::Internals::Serializer& rSerializer);
};

/// output stream function
inline std::ostream & operator <<(
    std::ostream& rOStream,
    const Node& rThis)
{
    rThis.Print(rOStream);
    return rOStream;
}


class CO_SIM_IO_API Element
{
public:
    using NodePointerType = CoSimIO::intrusive_ptr<Node>;
    using NodesContainerType = std::vector<NodePointerType>;

    Element(
        const IdType I_Id,
        const ElementType I_Type,
        const NodesContainerType& I_Nodes);

    // delete copy and assignment CTor
    Element(const Element&) = delete;
    Element& operator=(Element const&) = delete;

    IdType Id() const { return mId; }
    ElementType Type() const { return mType; }
    std::size_t NumberOfNodes() const { return mNodes.size(); }
    const Internals::PointerVector<NodePointerType> Nodes() const {return Internals::PointerVector<NodePointerType>(mNodes);}
    NodesContainerType::const_iterator NodesBegin() const { return mNodes.begin(); }
    NodesContainerType::const_iterator NodesEnd() const { return mNodes.end(); }

    void Print(std::ostream& rOStream) const;

private:
    IdType mId;
    ElementType mType;
    NodesContainerType mNodes;

    //*********************************************
    //this block is needed for refcounting
    mutable std::atomic<int> mReferenceCounter{0};

    friend void intrusive_ptr_add_ref(const Element* x)
    {
        x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
    }

    friend void intrusive_ptr_release(const Element* x)
    {
        if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
            std::atomic_thread_fence(std::memory_order_acquire);
            delete x;
        }
    }
    //*********************************************

    Element() = default; // needed for Serializer

    friend class CoSimIO::Internals::Serializer; // needs "CoSimIO::Internals::" because it is in different namespace

    void save(CoSimIO::Internals::Serializer& rSerializer) const;

    void load(CoSimIO::Internals::Serializer& rSerializer);
};

/// output stream function
inline std::ostream & operator <<(
    std::ostream& rOStream,
    const Element& rThis)
{
    rThis.Print(rOStream);
    return rOStream;
}


class CO_SIM_IO_API ModelPart
{
public:

    using NodePointerType = CoSimIO::intrusive_ptr<Node>;
    using ElementPointerType = CoSimIO::intrusive_ptr<Element>;
    using NodesContainerType = std::vector<NodePointerType>;
    using ElementsContainerType = std::vector<ElementPointerType>;

    using PartitionModelPartsContainerType = std::unordered_map<int, std::unique_ptr<ModelPart>>;

    explicit ModelPart(const std::string& I_Name);

    // delete copy and assignment CTor
    ModelPart(const ModelPart&) = delete;
    ModelPart& operator=(ModelPart const&) = delete;

    const std::string& Name() const { return mName; }

    std::size_t NumberOfNodes() const;
    std::size_t NumberOfLocalNodes() const;
    std::size_t NumberOfGhostNodes() const;

    std::size_t NumberOfElements() const { return mElements.size(); }

    Node& CreateNewNode(
        const IdType I_Id,
        const double I_X,
        const double I_Y,
        const double I_Z);

    Node& CreateNewGhostNode(
        const IdType I_Id,
        const double I_X,
        const double I_Y,
        const double I_Z,
        const int PartitionIndex);

    Element& CreateNewElement(
        const IdType I_Id,
        const ElementType I_Type,
        const ConnectivitiesType& I_Connectivities);

    const Internals::PointerVector<NodePointerType> Nodes() const {return Internals::PointerVector<NodePointerType>(mNodes);}
    const Internals::PointerVector<NodePointerType> LocalNodes() const {return Internals::PointerVector<NodePointerType>(GetLocalModelPart().Nodes());}
    const Internals::PointerVector<NodePointerType> GhostNodes() const {return Internals::PointerVector<NodePointerType>(GetGhostModelPart().Nodes());}

    const Internals::PointerVector<ElementPointerType> Elements() const {return Internals::PointerVector<ElementPointerType>(mElements);}

    const ModelPart& GetLocalModelPart() const;
    const ModelPart& GetGhostModelPart() const;
    const PartitionModelPartsContainerType& GetPartitionModelParts() const {return mPartitionModelParts;}

    NodesContainerType::const_iterator NodesBegin() const { return mNodes.begin(); }
    NodesContainerType::const_iterator LocalNodesBegin() const { return GetLocalModelPart().mNodes.begin(); }
    NodesContainerType::const_iterator GhostNodesBegin() const { return GetGhostModelPart().mNodes.begin(); }

    ElementsContainerType::const_iterator ElementsBegin() const { return mElements.begin(); }

    NodesContainerType::const_iterator NodesEnd() const { return mNodes.end(); }
    NodesContainerType::const_iterator LocalNodesEnd() const { return GetLocalModelPart().mNodes.end(); }
    NodesContainerType::const_iterator GhostNodesEnd() const { return GetGhostModelPart().mNodes.end(); }

    ElementsContainerType::const_iterator ElementsEnd() const { return mElements.end(); }

    Node& GetNode(const IdType I_Id);

    const Node& GetNode(const IdType I_Id) const;

    NodePointerType pGetNode(const IdType I_Id);

    Element& GetElement(const IdType I_Id);

    const Element& GetElement(const IdType I_Id) const;

    ElementPointerType pGetElement(const IdType I_Id);

    void Print(std::ostream& rOStream) const;

    void Clear();

protected:
    friend class std::unique_ptr<ModelPart>;
    ModelPart(const std::string& I_Name, const bool InitInternalModelParts);

private:
    std::string mName;
    NodesContainerType mNodes; // contains all nodes, local and ghost
    ElementsContainerType mElements; // contains all elements, local and ghost

    std::unique_ptr<ModelPart> mpLocalModelPart;
    std::unique_ptr<ModelPart> mpGhostModelPart;

    // this contains the the map of ModelPart with Ghost Nodes
    // for communication / synchronization the same is needed for the local nodes, i.e.
    // the ModelPart needs to know where its local nodes exist as ghost nodes
    // in Kratos this is computed with the ParallelFillCommunicator
    // Here is can be added in the future if required
    PartitionModelPartsContainerType mPartitionModelParts;

    NodesContainerType::const_iterator FindNode(const IdType I_Id) const;
    NodesContainerType::iterator FindNode(const IdType I_Id);

    ElementsContainerType::const_iterator FindElement(const IdType I_Id) const;
    ElementsContainerType::iterator FindElement(const IdType I_Id);

    bool HasNode(const IdType I_Id) const;

    bool HasElement(const IdType I_Id) const;

    ModelPart& GetLocalModelPart();

    ModelPart& GetGhostModelPart();

    ModelPart& GetPartitionModelPart(const int PartitionIndex);
    const ModelPart& GetPartitionModelPart(const int PartitionIndex) const;

    void InitializeInternalModelParts();

    ModelPart() = default; // needed for Serializer

    friend class CoSimIO::Internals::Serializer; // needs "CoSimIO::Internals::" because it is in different namespace

    void save(CoSimIO::Internals::Serializer& rSerializer) const;

    void load(CoSimIO::Internals::Serializer& rSerializer);
};

/// output stream function
inline std::ostream & operator <<(
    std::ostream& rOStream,
    const ModelPart& rThis)
{
    rThis.Print(rOStream);
    return rOStream;
}

} //namespace CoSimIO

#endif // CO_SIM_IO_MODEL_PART_INCLUDED

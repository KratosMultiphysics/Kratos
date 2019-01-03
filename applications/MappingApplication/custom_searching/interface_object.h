//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_INTERFACE_SEARCH_OBJECT_INCLUDED_H_INCLUDED )
#define  KRATOS_INTERFACE_SEARCH_OBJECT_INCLUDED_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"


namespace Kratos
{
///@addtogroup MappingApplication
///@{

///@name Kratos Classes
///@{

/// Object used by the bin-search
/** This object is used by the bin search. It is the baseclass for objects that hold information
 * about geometric entities, e.g. Nodes or Geometries
*/
class InterfaceObject : public Point
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceObject
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceObject);

    typedef Point BaseType;

    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    typedef Node<3> NodeType;
    typedef NodeType* NodePointerType;

    typedef Geometry<NodeType> GeometryType;
    typedef GeometryType* GeometryPointerType;

    ///@}
    ///@name  Enum's
    ///@{

    enum class ConstructionType
    {
        Node_Coords,
        Geometry_Center,
        Element_Center,
        Condition_Center
    };

    ///@}
    ///@name Life Cycle
    ///@{

    explicit InterfaceObject(const CoordinatesArrayType& rCoordinates)
        : Point(rCoordinates) { }

    /// Destructor.
    virtual ~InterfaceObject() = default;


    ///@}
    ///@name Operations
    ///@{

    virtual void UpdateCoordinates()
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    ///@}
    ///@name Access
    ///@{

    virtual NodePointerType pGetBaseNode() const
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }

    virtual GeometryPointerType pGetBaseGeometry() const
    {
        KRATOS_ERROR << "Base class function called!" << std::endl;
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "InterfaceObject" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override{rOStream << "InterfaceObject";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override{}


    ///@}

protected:
    ///@name Protected Operations
    ///@{

    // This constructor is called by its derived classes
    InterfaceObject() : Point(0.0, 0.0, 0.0)
    {
    }

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_ERROR << "This object is not supposed to be used with serialization!" << std::endl;
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_ERROR << "This object is not supposed to be used with serialization!" << std::endl;
    }

    ///@}

}; // Class InterfaceObject

///@}

///@} addtogroup block


class InterfaceNode : public InterfaceObject
{
public:
    typedef InterfaceObject::NodePointerType NodePointerType;

    InterfaceNode() {}

    explicit InterfaceNode(NodePointerType pNode)
        : mpNode(pNode)
    {
        UpdateCoordinates();
    }

    void UpdateCoordinates() override
    {
        noalias(Coordinates()) = mpNode->Coordinates();
    }

    NodePointerType pGetBaseNode() const override
    {
        return mpNode;
    }

private:
    NodePointerType mpNode;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_ERROR << "This object is not supposed to be used with serialization!" << std::endl;
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_ERROR << "This object is not supposed to be used with serialization!" << std::endl;
    }
};

class InterfaceGeometryObject : public InterfaceObject
{
public:
    typedef InterfaceObject::GeometryPointerType GeometryPointerType;

    InterfaceGeometryObject() {}

    explicit InterfaceGeometryObject(GeometryPointerType pGeometry)
        : mpGeometry(pGeometry)
    {
        UpdateCoordinates();
    }

    void UpdateCoordinates() override
    {
        noalias(Coordinates()) = mpGeometry->Center();
    }

    GeometryPointerType pGetBaseGeometry() const override
    {
        return mpGeometry;
    }

private:
    GeometryPointerType mpGeometry;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_ERROR << "This object is not supposed to be used with serialization!" << std::endl;
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_ERROR << "This object is not supposed to be used with serialization!" << std::endl;
    }
};

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_SEARCH_OBJECT_INCLUDED_H_INCLUDED  defined



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

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
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
    typedef Kratos::shared_ptr<NodeType> NodePointerType;

    typedef Geometry<NodeType> GeometryType;
    typedef Kratos::shared_ptr<GeometryType> GeometryPointerType;

    ///@}
    ///@name  Enum's
    ///@{

    enum ConstructionType
    {
        Node_Coords,
        Geometry_Center,
        Element_Center,
        Condition_Center
    };

    ///@}
    ///@name Life Cycle
    ///@{

    InterfaceObject(const CoordinatesArrayType& rCoordinates) : Point(rCoordinates)
    { }

    /// Destructor.
    virtual ~InterfaceObject(){}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void UpdateCoordinates(const CoordinatesArrayType& rCoordinates)
    {
        noalias(Coordinates()) = rCoordinates;
    }

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
    ///@name Inquiry
    ///@{


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
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    // This constructor is called by its derived classes
    InterfaceObject() : Point(0.0f, 0.0f, 0.0f)
    {
    }

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    // /// Assignment operator.
    // InterfaceObject& operator=(InterfaceObject const& rOther){}

    // /// Copy constructor.
    // InterfaceObject(InterfaceObject const& rOther){}

    ///@}
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

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                 InterfaceObject& rThis){}

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                 const InterfaceObject& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block


class InterfaceNode : public InterfaceObject
{
public:
    typedef InterfaceObject::NodePointerType NodePointerType;

    InterfaceNode() {}

    InterfaceNode(NodePointerType pNode) : mpNode(pNode)
    {
        UpdateCoordinates();
    }

    NodePointerType pGetBaseNode() const override
    {
        return mpNode;
    }

private:
    NodePointerType mpNode;

    void UpdateCoordinates() override
    {
        noalias(Coordinates()) = mpNode->Coordinates();
    }

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

    InterfaceGeometryObject(GeometryPointerType pGeometry) : mpGeometry(pGeometry)
    {
        UpdateCoordinates();
    }

    GeometryPointerType pGetBaseGeometry() const override
    {
        return mpGeometry;
    }

private:
    GeometryPointerType mpGeometry;

    void UpdateCoordinates() override
    {
        noalias(Coordinates()) = mpGeometry->Center();
    }

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



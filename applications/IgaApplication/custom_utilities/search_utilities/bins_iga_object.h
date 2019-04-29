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
#include "includes/element.h"
#include "iga_application.h"
#include "iga_application_variables.h"


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
class BinsIgaObject : public Point
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BinsIgaObject
    KRATOS_CLASS_POINTER_DEFINITION(BinsIgaObject);

    using BaseType = Point;

    using CoordinatesArrayType = typename BaseType::CoordinatesArrayType;

    using ElementType = Element;
    using ElementPointerType = ElementType::Pointer;

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    BinsIgaObject() : Point()
    { };

    BinsIgaObject(ElementPointerType pElement) : mpElement(pElement)
    {
        array_1d<double, 3> coords = ZeroVector(3);
        const Kratos::ProcessInfo process_info;
        pElement->Calculate(COORDINATES, coords, process_info);
        noalias(Coordinates()) = coords;
    };

    //BinsIgaObject(const CoordinatesArrayType& rCoordinates, ElementPointerType pElement) : Point(rCoordinates), mpElement(pElement)
    //{
    //    Vector coords = ZeroVector(3);
    //    const Kratos::ProcessInfo process_info;
    //    pElement->Calculate(COORDINATES, coords, process_info);
    //    noalias(Coordinates()) = coords;
    //};

    BinsIgaObject(const CoordinatesArrayType& rCoordinates) : Point(rCoordinates)
    {
    };

    /// Destructor.
    virtual ~BinsIgaObject(){}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void UpdateCoordinates(const CoordinatesArrayType& rCoordinates)
    {
        noalias(Coordinates()) = rCoordinates;
    };

    ///@}
    ///@name Access
    ///@{

    ElementPointerType pGetBaseElement() const
    {
        return mpElement;
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
        buffer << "BinsIgaObject" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override{rOStream << "BinsIgaObject";}

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

    ElementPointerType mpElement;


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
    // BinsIgaObject& operator=(BinsIgaObject const& rOther){}

    // /// Copy constructor.
    // BinsIgaObject(BinsIgaObject const& rOther){}

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

}; // Class BinsIgaObject

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


// /// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                 BinsIgaObject& rThis){}

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                 const BinsIgaObject& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_SEARCH_OBJECT_INCLUDED_H_INCLUDED  defined
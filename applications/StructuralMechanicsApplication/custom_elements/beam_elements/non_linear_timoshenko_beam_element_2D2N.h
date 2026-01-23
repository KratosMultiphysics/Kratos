// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "custom_elements/beam_elements/timoshenko_beam_element_2D2N.h"

namespace Kratos
{

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

/**
 * @class NonLinearTimoshenkoBeamElement2D2N
 * @ingroup StructuralMechanicsApplication
 * @brief This is a non-linear (geometric) Timoshenko beam element of 2 nodes. It extends the LinearTimoshenkoBeamElement2D2N
 * with geometrically exact capabilities.
 * Reference: "Nonlinear finite element methods", P. Wriggers, Springer 2008.
 * @author Alejandro Cornejo
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) NonLinearTimoshenkoBeamElement2D2N
    : public LinearTimoshenkoBeamElement2D2N
{

public:

    ///@name Type Definitions
    ///@{

    /// The base element type
    using BaseType = LinearTimoshenkoBeamElement2D2N;

    // Counted pointer of NonLinearTimoshenkoBeamElement2D2N
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(NonLinearTimoshenkoBeamElement2D2N);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    NonLinearTimoshenkoBeamElement2D2N()
    {
    }

    // Constructor using an array of nodes
    NonLinearTimoshenkoBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry) 
        : LinearTimoshenkoBeamElement2D2N(NewId, pGeometry)
    {
    }

    // Constructor using an array of nodes with properties
    NonLinearTimoshenkoBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : LinearTimoshenkoBeamElement2D2N(NewId, pGeometry, pProperties)
    {
    }

    // Copy constructor
    NonLinearTimoshenkoBeamElement2D2N(NonLinearTimoshenkoBeamElement2D2N const& rOther)
        : BaseType(rOther)
    {
    }

    // Create method
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<NonLinearTimoshenkoBeamElement2D2N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    // Create method
    Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const override
    {
        return Kratos::make_intrusive<NonLinearTimoshenkoBeamElement2D2N>(NewId, pGeom, pProperties);
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "NonLinearTimoshenkoBeamElement2D2N";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NonLinearTimoshenkoBeamElement2D2N";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

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
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer &rSerializer) const override;

    void load(Serializer &rSerializer) override;

}; // class NonLinearTimoshenkoBeamElement2D2N.

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

} // namespace Kratos.

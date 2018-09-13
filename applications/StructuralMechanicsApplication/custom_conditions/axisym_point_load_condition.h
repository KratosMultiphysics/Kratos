// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_AXISYM_POINT_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_AXISYM_POINT_LOAD_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/point_load_condition.h"

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

/// Axisymmetric point load condition

/**
 * Implements a point load condition for structural analysis.
 */

class AxisymPointLoadCondition
    : public PointLoadCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of AxisymPointLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION(AxisymPointLoadCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AxisymPointLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    AxisymPointLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    ~AxisymPointLoadCondition() override;

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    //std::string Info() const;

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
    //      virtual String Info() const;

    /// Print information about this object.
    //      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    //      virtual void PrintData(std::ostream& rOStream) const;
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

    AxisymPointLoadCondition() : PointLoadCondition()
    {
    }

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

     /**
     * It calcules the integration load for the point load
     */
    double GetPointLoadIntegrationWeight() override;

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    /// Assignment operator.
    //AxisymPointLoadCondition& operator=(const AxisymPointLoadCondition& rOther);
    /// Copy constructor.
    //AxisymPointLoadCondition(const AxisymPointLoadCondition& rOther);
    ///@}

}; // Class AxisymPointLoadCondition

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_AXISYM_POINT_LOAD_CONDITION_H_INCLUDED  defined

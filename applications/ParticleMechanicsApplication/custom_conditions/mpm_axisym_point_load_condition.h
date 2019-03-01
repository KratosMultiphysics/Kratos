//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#if !defined(KRATOS_MPM_AXISYM_POINT_LOAD_CONDITION_H_INCLUDED )
#define      KRATOS_MPM_AXISYM_POINT_LOAD_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/mpm_point_load_condition.h"

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

class MPMAxisymPointLoadCondition
    : public MPMPointLoadCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MPMAxisymPointLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION(MPMAxisymPointLoadCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MPMAxisymPointLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    MPMAxisymPointLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    ~MPMAxisymPointLoadCondition() override;

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

    MPMAxisymPointLoadCondition() : MPMPointLoadCondition()
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
    //MPMAxisymPointLoadCondition& operator=(const MPMAxisymPointLoadCondition& rOther);
    /// Copy constructor.
    //MPMAxisymPointLoadCondition(const MPMAxisymPointLoadCondition& rOther);
    ///@}

}; // Class MPMAxisymPointLoadCondition

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_AXISYM_POINT_LOAD_CONDITION_H_INCLUDED  defined

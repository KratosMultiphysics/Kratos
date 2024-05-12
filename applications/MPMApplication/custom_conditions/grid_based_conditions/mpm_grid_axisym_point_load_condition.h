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


#if !defined(KRATOS_MPM_GRID_AXISYM_POINT_LOAD_CONDITION_H_INCLUDED )
#define      KRATOS_MPM_GRID_AXISYM_POINT_LOAD_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/grid_based_conditions/mpm_grid_point_load_condition.h"

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

class MPMGridAxisymPointLoadCondition
    : public MPMGridPointLoadCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MPMGridAxisymPointLoadCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MPMGridAxisymPointLoadCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MPMGridAxisymPointLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    MPMGridAxisymPointLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    ~MPMGridAxisymPointLoadCondition() override;

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
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

    MPMGridAxisymPointLoadCondition() : MPMGridPointLoadCondition()
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
    //MPMGridAxisymPointLoadCondition& operator=(const MPMGridAxisymPointLoadCondition& rOther);
    /// Copy constructor.
    //MPMGridAxisymPointLoadCondition(const MPMGridAxisymPointLoadCondition& rOther);
    ///@}

}; // Class MPMGridAxisymPointLoadCondition

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_MPM_GRID_AXISYM_POINT_LOAD_CONDITION_H_INCLUDED  defined

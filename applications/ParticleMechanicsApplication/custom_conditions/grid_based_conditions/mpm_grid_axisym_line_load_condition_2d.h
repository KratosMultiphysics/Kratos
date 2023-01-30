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


#if !defined(KRATOS_MPM_GRID_AXISYM_LINE_LOAD_CONDITION_2D_H_INCLUDED )
#define      KRATOS_MPM_GRID_AXISYM_LINE_LOAD_CONDITION_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/grid_based_conditions/mpm_grid_line_load_condition_2d.h"

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

/// Axisymmetric line load condition

/**
 * Implements a line load condition for structural analysis.
 */

class MPMGridAxisymLineLoadCondition2D
    : public MPMGridLineLoadCondition2D
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MPMGridAxisymLineLoadCondition2D
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MPMGridAxisymLineLoadCondition2D);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MPMGridAxisymLineLoadCondition2D(IndexType NewId, GeometryType::Pointer pGeometry);
    MPMGridAxisymLineLoadCondition2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    ~MPMGridAxisymLineLoadCondition2D() override;

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

    MPMGridAxisymLineLoadCondition2D() : MPMGridLineLoadCondition2D()
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
     * This functions computes the integration weight to consider
     * @param IntegrationPoints: The array containing the integration points
     * @param PointNumber: The id of the integration point considered
     * @param detJ: The determinant of the jacobian of the element
     */
    double GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber,
        const double detJ
        ) override;

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
    //MPMGridAxisymLineLoadCondition2D& operator=(const MPMGridAxisymLineLoadCondition2D& rOther);
    /// Copy constructor.
    //MPMGridAxisymLineLoadCondition2D(const MPMGridAxisymLineLoadCondition2D& rOther);
    ///@}

}; // Class MPMGridAxisymLineLoadCondition2D

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_MPM_GRID_AXISYM_LINE_LOAD_CONDITION_2D_H_INCLUDED  defined

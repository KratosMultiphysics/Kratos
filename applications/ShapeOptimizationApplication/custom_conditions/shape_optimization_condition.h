// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#if !defined(KRATOS_SHAPE_OPTIMIZATION_CONDITION_H_INCLUDED )
#define  KRATOS_SHAPE_OPTIMIZATION_CONDITION_H_INCLUDED

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/condition.h"

// ==============================================================================

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

/// Shape optimization condition

/**
 * Implements a shape definition to allow for shape optimization.
 * This works currently only works for 3D geometries composed of 3-noded triangles
 */
class ShapeOptimizationCondition
    : public Condition
{
public:

    ///@name Type Definitions
    ///@{
    // Counted pointer of ShapeOptimizationCondition
    KRATOS_CLASS_POINTER_DEFINITION( ShapeOptimizationCondition );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ShapeOptimizationCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    ShapeOptimizationCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Copy constructor
    ShapeOptimizationCondition( ShapeOptimizationCondition const& rOther);

    /// Destructor
    virtual ~ShapeOptimizationCondition();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new condition pointer
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId,
			      NodesArrayType const& ThisNodes,
			      PropertiesType::Pointer pProperties ) const override;


    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId,
			     NodesArrayType const& ThisNodes) const override;



    //************* COMPUTING  METHODS

    /**
     * Calculate a condition variable
     */

    // void Calculate(const Variable< array_1d<double,3> > &rVariable, double &Output , ProcessInfo &rCurrentProcessInfo );

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) override;

    ///@}
    ///@name Access
    ///@{
    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
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
    ShapeOptimizationCondition() {};
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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // class ShapeOptimizationCondition.

} // namespace Kratos.

#endif // KRATOS_SHAPE_OPTIMIZATION_CONDITION_H_INCLUDED defined

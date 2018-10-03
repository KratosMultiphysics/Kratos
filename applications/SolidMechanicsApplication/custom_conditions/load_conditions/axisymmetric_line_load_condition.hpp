//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_AXISYMMETRIC_LINE_LOAD_CONDITION_H_INCLUDED)
#define  KRATOS_AXISYMMETRIC_LINE_LOAD_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/load_conditions/line_load_condition.hpp"

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

/// Load Condition for 2D axisymmetric geometries. (base class)

/**
 * Implements a Load definition for structural analysis.
 * This works for arbitrary geometries in 2D (base class)
 */
class KRATOS_API(SOLID_MECHANICS_APPLICATION) AxisymmetricLineLoadCondition
    : public LineLoadCondition
{
public:

    ///@name Type Definitions
    ///@{
    // Counted pointer of AxisymmetricLineLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION( AxisymmetricLineLoadCondition );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AxisymmetricLineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    AxisymmetricLineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Copy constructor
    AxisymmetricLineLoadCondition( AxisymmetricLineLoadCondition const& rOther);

    /// Destructor
    ~AxisymmetricLineLoadCondition() override;

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
    AxisymmetricLineLoadCondition() {};
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculate Condition Kinematics
     */
    void CalculateKinematics(ConditionVariables& rVariables,
				     const double& rPointNumber) override;

    /**
     * Calculation and addition of the matrices of the LHS
     */
    void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    ConditionVariables& rVariables,
                                    double& rIntegrationWeight) override;

    /**
     * Calculation and addition of the vectors of the RHS
     */
    void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                                    ConditionVariables& rVariables,
                                    double& rIntegrationWeight) override;

    /**
     * Calculation of the contidion radius (axisymmetry)
     */
    void CalculateRadius(double & rCurrentRadius,
			 double & rReferenceRadius,
			 const Vector& rN);

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


}; // class AxisymmetricLineLoadCondition.

} // namespace Kratos.

#endif // KRATOS_AXISYMMETRIC_LINE_LOAD_CONDITION_H_INCLUDED defined

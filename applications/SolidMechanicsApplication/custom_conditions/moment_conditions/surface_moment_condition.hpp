//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SURFACE_MOMENT_CONDITION_H_INCLUDED)
#define  KRATOS_SURFACE_MOMENT_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/moment_conditions/moment_condition.hpp"

namespace Kratos
{
///@name Kratos Globals///@{
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

// Surface moment condition for 3D geometries.

class KRATOS_API(SOLID_MECHANICS_APPLICATION) SurfaceMomentCondition
    : public MomentCondition
{
public:

    ///@name Type Definitions
    ///@{
    // Counted pointer of SurfaceMomentCondition
    KRATOS_CLASS_POINTER_DEFINITION( SurfaceMomentCondition );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SurfaceMomentCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    SurfaceMomentCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Copy constructor
    SurfaceMomentCondition( SurfaceMomentCondition const& rOther);

    /// Destructor
    ~SurfaceMomentCondition() override;

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
    SurfaceMomentCondition() {};
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Initialize System Matrices
     */
    void InitializeConditionVariables(ConditionVariables& rVariables,
					    const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Calculate Condition Kinematics
     */
    void CalculateKinematics(ConditionVariables& rVariables,
				     const double& rPointNumber) override;

    /**
     * Calculate the External Moment of the Condition
     */
    void CalculateExternalMoment(ConditionVariables& rVariables) override;

    /**
     * Calculation of the Moment Stiffness Matrix which usually is subtracted to the global stiffness matrix
     */
    void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
				     ConditionVariables& rVariables,
				     double& rIntegrationWeight) override;


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

}; // class SurfaceMomentCondition.

} // namespace Kratos.

#endif // KRATOS_SURFACE_MOMENT_CONDITION_H_INCLUDED defined

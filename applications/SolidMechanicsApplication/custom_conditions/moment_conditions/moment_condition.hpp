//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_MOMENT_CONDITION_H_INCLUDED)
#define  KRATOS_MOMENT_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/boundary_condition.hpp"

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

/// Load Condition for 3D and 2D geometries. (base class)

/**
 * Implements a General Load definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */
class KRATOS_API(SOLID_MECHANICS_APPLICATION) MomentCondition
    : public BoundaryCondition
{
public:

    ///@name Type Definitions
    ///@{

    ///Type for size
    typedef GeometryData::SizeType SizeType;

    // Counted pointer of MomentCondition
    KRATOS_CLASS_POINTER_DEFINITION( MomentCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    MomentCondition();

    /// Default constructor.
    MomentCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    MomentCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Copy constructor
    MomentCondition( MomentCondition const& rOther);

    /// Destructor
    ~MomentCondition() override;

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


    /**
     * Sets on rConditionDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rConditionDofList,
		    ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(EquationIdVectorType& rResult,
			  ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues,
			 int Step = 0) override;

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues,
				   int Step = 0) override;

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues,
				    int Step = 0) override;

    /**
     * this function is designed to make the element to assemble an rRHS vector
     * identified by a variable rRHSVariable by assembling it to the nodes on the variable
     * rDestinationVariable.
     * @param rRHSVector: input variable containing the RHS vector to be assembled
     * @param rRHSVariable: variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable: variable in the database to which the rRHSvector will be assembled
      * @param rCurrentProcessInfo: the current process info instance
     */
    void AddExplicitContribution(const VectorType& rRHS,
					 const Variable<VectorType>& rRHSVariable,
					 Variable<array_1d<double,3> >& rDestinationVariable,
					 const ProcessInfo& rCurrentProcessInfo) override;

    //************************************************************************************
    //************************************************************************************
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
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Initialize Explicit Contributions
     */
    void InitializeExplicitContributions();


    /**
     * Initialize System Matrices
     */

    unsigned int GetDofsSize() override;


    /**
     * Calculate the External Load of the Condition
     */
    virtual void CalculateExternalMoment(ConditionVariables& rVariables);


    /**
     * Calculation of the External Forces Vector for a force or pressure vector
     */
    void CalculateAndAddExternalForces(Vector& rRightHandSideVector,
					       ConditionVariables& rVariables,
					       double& rIntegrationWeight) override;


    /**
     * Calculation of the External Forces Vector for a force or pressure vector
     */
    double& CalculateAndAddExternalEnergy(double& rEnergy,
						  ConditionVariables& rVariables,
						  double& rIntegrationWeight,
						  const ProcessInfo& rCurrentProcessInfo) override;

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


}; // class MomentCondition.

} // namespace Kratos.

#endif // KRATOS_MOMENT_CONDITION_H_INCLUDED defined

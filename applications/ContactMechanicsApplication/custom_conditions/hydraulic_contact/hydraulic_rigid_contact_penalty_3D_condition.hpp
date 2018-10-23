//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:                LMonforte $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:               January 2018 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_HYDRAULIC_RIGID_CONTACT_PENALTY_CONDITION_H_INCLUDED )
#define  KRATOS_HYDRAULIC_RIGID_CONTACT_PENALTY_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/rigid_contact/point_rigid_contact_condition.hpp"

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

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(CONTACT_MECHANICS_APPLICATION) HydraulicRigidContactPenalty3DCondition
    : public PointRigidContactCondition
{
public:

   ///@name Type Definitions

    ///Tensor order 1 definition
    typedef BoundedVector<double, 3>     PointType;

    ///@{
    // Counted pointer of WaterPointRigidContactCondition
    KRATOS_CLASS_POINTER_DEFINITION( HydraulicRigidContactPenalty3DCondition );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Serialization constructor
    HydraulicRigidContactPenalty3DCondition(){};

    /// Default constructor.
    HydraulicRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    HydraulicRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    HydraulicRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall);


    /// Copy constructor
    HydraulicRigidContactPenalty3DCondition( HydraulicRigidContactPenalty3DCondition const& rOther);


    /// Destructor.
    virtual ~HydraulicRigidContactPenalty3DCondition();


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
    Condition::Pointer Create(IndexType NewId, NodesArrayType const&
                              ThisNodes,  PropertiesType::Pointer pProperties) const override;

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
     * Called at the beginning of each iteration
     */
    virtual void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    //************* GETTING METHODS

    /**
     * Sets on rConditionDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rConditionDofList,
		    ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(EquationIdVectorType& rResult,
			  ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues,
			 int Step = 0 ) override;

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues,
				   int Step = 0 ) override;

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues,
				    int Step = 0 ) override;



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
     * Calculate Condition Kinematics
     */
    virtual void CalculateKinematics(ConditionVariables& rVariables,
				     const ProcessInfo& rCurrentProcessInfo,
				     const double& rPointNumber) override;



    /**
     * Calculation of the Load Stiffness Matrix which usually is subtracted to the global stiffness matrix
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
				     ConditionVariables& rVariables,
				     double& rIntegrationWeight) override;

    /**
     * Calculation of the External Forces Vector for a force or pressure vector
     */
    virtual void CalculateAndAddContactForces(Vector& rRightHandSideVector,
					      ConditionVariables& rVariables,
					      double& rIntegrationWeight ) override;


    void CalculateAndAddNormalContactForce(Vector& rRightHandSideVector, ConditionVariables& rVariables, double& rIntegrationWeight);

    double& CalculateNormalForceModulus( double& rNormalForceModulus, ConditionVariables& rVariables );


    /**
     * Calculation of the Contact Force Factors
     */
    virtual void CalculateContactFactors(ConditionVariables &rContact);

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

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PointRigidContactCondition )
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PointRigidContactCondition )
    }


}; // Class HydraulicRigidContactPenalty3DCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    HydraulicRigidContactPenalty3DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const HydraulicRigidContactPenalty3DCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_HYDRAULIC_RIGID_CONTACT_PENALTY_CONDITION_H_INCLUDED  defined

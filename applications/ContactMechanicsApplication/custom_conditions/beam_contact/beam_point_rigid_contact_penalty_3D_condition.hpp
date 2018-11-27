//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2013 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_BEAM_POINT_RIGID_CONTACT_PENALTY_3D_CONDITION_H_INCLUDED )
#define  KRATOS_BEAM_POINT_RIGID_CONTACT_PENALTY_3D_CONDITION_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_conditions/thermal_contact/beam_point_rigid_contact_condition.hpp"

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
class KRATOS_API(CONTACT_MECHANICS_APPLICATION) BeamPointRigidContactPenalty3DCondition
    : public BeamPointRigidContactCondition
{
protected:

    typedef struct
     {
        bool   Slip;
        double Sign;

        double DeltaTime;
        double PreviousTangentForceModulus;

        double FrictionCoefficient;
        double DynamicFrictionCoefficient;
        double StaticFrictionCoefficient;

     } TangentialContactVariables;



public:

    ///@name Type Definitions

    ///Tensor order 1 definition
    //typedef BoundedVector<double, 3>     PointType;
    typedef array_1d<double, 3>             PointType;

    ///@{
    // Counted pointer of BeamPointRigidContactCondition
    KRATOS_CLASS_POINTER_DEFINITION( BeamPointRigidContactPenalty3DCondition );
    ///@}

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BeamPointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    BeamPointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    BeamPointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall);


    /// Copy constructor
    BeamPointRigidContactPenalty3DCondition( BeamPointRigidContactPenalty3DCondition const& rOther);


    /// Destructor.
    virtual ~BeamPointRigidContactPenalty3DCondition();


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
                              ThisNodes,  PropertiesType::Pointer pProperties) const;

    /**
     * Called at the beginning of each step
     */

    virtual void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);


    /**
     * Called at the beginning of each iteration
     */
    virtual void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

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

    // A protected default constructor necessary for serialization
    BeamPointRigidContactPenalty3DCondition() {};

    ///@}
    ///@name Protected member Variables
    ///@{

    TangentialContactVariables mTangentialVariables;

    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Initialize General Variables
     */
    void InitializeConditionVariables(ConditionVariables& rVariables,
				    const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate Condition Kinematics
     */
    virtual void CalculateKinematics(ConditionVariables& rVariables,
				     const ProcessInfo& rCurrentProcessInfo,
				     const double& rPointNumber);


    /**
     * Calculation of the Load Stiffness Matrix which usually is subtracted to the global stiffness matrix
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
				     ConditionVariables& rVariables,
				     double& rIntegrationWeight);

    virtual void CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix,
				     ConditionVariables& rVariables,
				     double& rIntegrationWeight);

    /**
     * Calculation of the External Forces Vector for a force or pressure vector
     */
    virtual void CalculateAndAddContactForces(Vector& rRightHandSideVector,
					      ConditionVariables& rVariables,
					      double& rIntegrationWeight );


    virtual void CalculateAndAddNormalContactForce(Vector& rRightHandSideVector, ConditionVariables& rVariables, double& rIntegrationWeight);


    virtual void CalculateAndAddTangentContactForce(Vector& rRightHandSideVector, ConditionVariables& rVariables, double& rIntegrationWeight);


    double& CalculateNormalForceModulus( double& rNormalForceModulus, ConditionVariables& rVariables );

    double& CalculateTangentRelativeMovement( double& rTangentRelativeMovement, ConditionVariables& rVariables );

    double CalculateCoulombsFrictionLaw( double& rTangentForceModulus, double& rNormalForceModulus, ConditionVariables& rVariables );

    double CalculateFrictionCoefficient(double & rTangentRelativeMovement);


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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BeamPointRigidContactCondition )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BeamPointRigidContactCondition )
    }


}; // Class BeamPointRigidContactPenalty3DCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    BeamPointRigidContactPenalty3DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const BeamPointRigidContactPenalty3DCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_POINT_RIGID_CONTACT_PENALTY_3D_CONDITION_H_INCLUDED  defined

//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_POINT_RIGID_CONTACT_PENALTY_3D_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_RIGID_CONTACT_PENALTY_3D_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/point_rigid_contact_condition.hpp"

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
class PointRigidContactPenalty3DCondition
    : public PointRigidContactCondition
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
        double IntegrationWeight;

        double Cohesion;
        double Neighb_distance; 

     } TangentialContactVariables;



public:

   ///@name Type Definitions

    ///Tensor order 1 definition
    typedef bounded_vector<double, 3>     PointType;

    ///@{
    // Counted pointer of PointRigidContactCondition
    KRATOS_CLASS_POINTER_DEFINITION( PointRigidContactPenalty3DCondition );
    ///@}
 
    ///@}
    ///@name Life Cycle
    ///@{

    /// Serialization constructor
    PointRigidContactPenalty3DCondition(){};

    /// Default constructor.
    PointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    PointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    PointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall);
  

    /// Copy constructor
    PointRigidContactPenalty3DCondition( PointRigidContactPenalty3DCondition const& rOther);


    /// Destructor.
    virtual ~PointRigidContactPenalty3DCondition();


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
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, 
			     NodesArrayType const& ThisNodes) const;

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
    void InitializeGeneralVariables(GeneralVariables& rVariables, 
				    const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate Condition Kinematics
     */
    virtual void CalculateKinematics(GeneralVariables& rVariables,
				     const ProcessInfo& rCurrentProcessInfo,
				     const double& rPointNumber);



    /**
     * Calculation of the Load Stiffness Matrix which usually is subtracted to the global stiffness matrix
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
				     GeneralVariables& rVariables,
				     double& rIntegrationWeight);

    virtual void CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix,
				     GeneralVariables& rVariables,
				     double& rIntegrationWeight);

    /**
     * Calculation of the External Forces Vector for a force or pressure vector 
     */
    virtual void CalculateAndAddContactForces(Vector& rRightHandSideVector,
					      GeneralVariables& rVariables,
					      double& rIntegrationWeight );


    virtual void CalculateAndAddNormalContactForce(Vector& rRightHandSideVector, GeneralVariables& rVariables, double& rIntegrationWeight);


    virtual void CalculateAndAddTangentContactForce(Vector& rRightHandSideVector, GeneralVariables& rVariables, double& rIntegrationWeight);


    double& CalculateNormalForceModulus( double& rNormalForceModulus, GeneralVariables& rVariables );
    

    double CalculateCoulombsFrictionLaw( double& rTangentForceModulus, double& rNormalForceModulus, GeneralVariables& rVariables );

    double CalculateEffectiveNormalForceModulus(const double&  rNormalForceModulus);

    double CalculateFrictionCoefficient(double & rTangentRelativeMovement);


    /**
     * Calculation of the Contact Force Factors
     */
    virtual void CalculateContactFactors(GeneralVariables &rContact);

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PointRigidContactCondition )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PointRigidContactCondition )
    }


}; // Class PointRigidContactPenalty3DCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    PointRigidContactPenalty3DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const PointRigidContactPenalty3DCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_POINT_RIGID_CONTACT_PENALTY_3D_CONDITION_H_INCLUDED  defined 




//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_EP_POINT_RIGID_CONTACT_PENALTY_3D_CONDITION_H_INCLUDED )
#define  KRATOS_EP_POINT_RIGID_CONTACT_PENALTY_3D_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/point_rigid_contact_penalty_3D_condition.hpp"

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
class KRATOS_API(CONTACT_MECHANICS_APPLICATION) EPPointRigidContactPenalty3DCondition
    : public PointRigidContactPenalty3DCondition
{

   protected:
      typedef struct
      {
         Vector PreviousStepForceVector;
         Vector t1;
         Vector t2;
         Vector n;
      } GeometricalInformation;

      typedef struct
      {
         // ConstitutiveInformation
         double TangentForceRatio;
         double  NormalTangentMatrix;
         double TangentTangentMatrix;
         Vector ForceDirection;
      }  ConstitutiveVariables;

   public:

   ///@name Type Definitions

    ///Tensor order 1 definition
    typedef BoundedVector<double, 3>     PointType;

    ///@{
    // Counted pointer of PointRigidContactCondition
    KRATOS_CLASS_POINTER_DEFINITION( EPPointRigidContactPenalty3DCondition );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Serialization constructor
    EPPointRigidContactPenalty3DCondition(){};

    /// Default constructor.
    EPPointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    EPPointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    EPPointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall);


    /// Copy constructor
    EPPointRigidContactPenalty3DCondition( EPPointRigidContactPenalty3DCondition const& rOther);


    /// Destructor.
    virtual ~EPPointRigidContactPenalty3DCondition();


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
     * Called at the end of each solution step
     */
    virtual void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Called at the end of each solution step
     */
    virtual void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Called at the beginning of each iteration
     */
    virtual void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

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

    GeometricalInformation mCurrentInfo;
    GeometricalInformation mSavedInfo;

    double mElasticYoungModulus; // using MCC + IMPLEX, at the finalizeSolutionStep the Contact Forces are correctly ingrated with the finalized of the Young Modulus at the continuum elements, that is quite different from the one used during the implex step.

    bool mImplex;

    ///@}
    ///@name Protected Operators
    ///@{

  ///@}
    ///@name Protected Operations
    ///@{

    void CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix,
                                    ConditionVariables& rVariables,
                                    double& rIntegrationWeight) override;




    virtual void CalculateAndAddTangentContactForce(Vector& rRightHandSideVector, ConditionVariables& rVariables, double& rIntegrationWeight) override;


    bool CalculateFrictionLaw( ConditionVariables & rVariables, ConstitutiveVariables & rConstitutiveVariables, Vector & rTangentForce);

    virtual double CalculateSomeSortOfArea();

    double CalculateEffectiveNormalForceModulus( const double& rNormalForceModulus, const double & rArea);

    Matrix ConvertToTheAppropriateSize( const Matrix & rForceMatrix);

    /**
     * Calculation of the Contact Force Factors
     */
    virtual void CalculateContactFactors(ConditionVariables &rContact) override;
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, PointRigidContactPenalty3DCondition )
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, PointRigidContactPenalty3DCondition )
    }


}; // Class EPPointRigidContactPenalty3DCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    EPPointRigidContactPenalty3DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const EPPointRigidContactPenalty3DCondition& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_POINT_RIGID_CONTACT_PENALTY_3D_CONDITION_H_INCLUDED  defined

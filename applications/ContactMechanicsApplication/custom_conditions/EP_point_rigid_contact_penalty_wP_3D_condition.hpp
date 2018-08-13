//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:                LMonforte $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2013 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_EP_POINT_RIGID_CONTACT_PENALTY_WP_3D_CONDITION_H_INCLUDED)
#define      KRATOS_EP_POINT_RIGID_CONTACT_PENALTY_WP_3D_CONDITION_H_INCLUDED

#include "custom_conditions/EP_point_rigid_contact_penalty_3D_condition.hpp"

namespace Kratos
{

class KRATOS_API(CONTACT_MECHANICS_APPLICATION) EPPointRigidContactPenaltywP3DCondition
    : public EPPointRigidContactPenalty3DCondition
{
 public:
  typedef Vector VectorType;

  KRATOS_CLASS_POINTER_DEFINITION( EPPointRigidContactPenaltywP3DCondition );

  /// Serialization constructor
  EPPointRigidContactPenaltywP3DCondition(){};

  /// Default constructor.
  EPPointRigidContactPenaltywP3DCondition(IndexType NewId, GeometryType::Pointer pGeometry);

  EPPointRigidContactPenaltywP3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

  EPPointRigidContactPenaltywP3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall);


  /// Copy constructor
  EPPointRigidContactPenaltywP3DCondition( EPPointRigidContactPenaltywP3DCondition const& rOther);


  /// Destructor.
  virtual ~EPPointRigidContactPenaltywP3DCondition();

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

 protected:
  ///@name Protected static Member Variables
  ///@{
  void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
                                VectorType& rRightHandSideVector,
                                Flags& rCalculationFlags) override;


  /**
   * Calculation of the Load Stiffness Matrix which usually is subtracted to the global stiffness matrix
   */
  void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
                           ConditionVariables& rVariables,
                           double& rIntegrationWeight) override;

  void CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix,
                                  ConditionVariables& rVariables,
                                  double& rIntegrationWeight) override;

 private:
  friend class Serializer;

  void save(Serializer& rSerializer) const override
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, EPPointRigidContactPenalty3DCondition )
        }

  void load(Serializer& rSerializer) override
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, EPPointRigidContactPenalty3DCondition )
        }



}; // end class


} // end Namespace Kratos



#endif

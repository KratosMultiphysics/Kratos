//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:             September 2018 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_RIGID_BODY_POINT_LINK_SEGREGATED_V_CONDITION_H_INCLUDED )
#define  KRATOS_RIGID_BODY_POINT_LINK_SEGREGATED_V_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/rigid_body_links/rigid_body_point_link_condition.hpp"

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

/// Rigid Body Point Rigid Contact Condition for 3D and 2D geometries. (base class)

/**
 * Implements a Contact Point Load definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */
class RigidBodyPointLinkSegregatedVCondition
    : public RigidBodyPointLinkCondition
{
 public:

  ///@name Type Definitions

  typedef Vector                             VectorType;
  typedef Element                           ElementType;
  typedef Node<3>::Pointer             PointPointerType;
  typedef Quaternion<double>             QuaternionType;
  typedef Node<3>::DofsContainerType  DofsContainerType;
  typedef GeometryData::SizeType               SizeType;

  ///@{

  // Counted pointer of RigidBodyPointLinkSegregatedVCondition
  KRATOS_CLASS_POINTER_DEFINITION( RigidBodyPointLinkSegregatedVCondition );

  enum StepType{VELOCITY_STEP = 0, PRESSURE_STEP = 1};

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  RigidBodyPointLinkSegregatedVCondition( IndexType NewId, GeometryType::Pointer pGeometry );

  RigidBodyPointLinkSegregatedVCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

  /// Copy constructor
  RigidBodyPointLinkSegregatedVCondition( RigidBodyPointLinkSegregatedVCondition const& rOther);

  /// Destructor
  virtual ~RigidBodyPointLinkSegregatedVCondition();

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

  //************* STARTING - ENDING  METHODS

  /**
   * is called to initialize the condition
   * if the condition needs to perform any operation before any calculation is done
   * the condition variables will be initialized and set using this method
   */
  void Initialize() override;

  /**
   * Called at the beginning of each solution step
   */
  void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;


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


  //************* COMPUTING  METHODS

  /**
   * this is called during the assembling process in order
   * to calculate all condition contributions to the global system
   * matrix and the right hand side
   * @param rLeftHandSideMatrix: the condition left hand side matrix
   * @param rRightHandSideVector: the condition right hand side
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                            VectorType& rRightHandSideVector,
                            ProcessInfo& rCurrentProcessInfo) override;

  /**
   * this is called during the assembling process in order
   * to calculate the condition right hand side vector only
   * @param rRightHandSideVector: the condition right hand side vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateRightHandSide(VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override;

  /**
   * this is called during the assembling process in order
   * to calculate the second derivatives contributions for the LHS and RHS
   * @param rLeftHandSideMatrix: the condition left hand side matrix
   * @param rRightHandSideVector: the condition right hand side
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                                               VectorType& rRightHandSideVector,
                                               ProcessInfo& rCurrentProcessInfo) override;

  /**
   * this is called during the assembling process in order
   * to calculate the condition left hand side matrix for the second derivatives constributions
   * @param rLeftHandSideMatrix: the condition left hand side matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                     ProcessInfo& rCurrentProcessInfo) override;


  /**
   * this is called during the assembling process in order
   * to calculate the condition right hand side vector for the second derivatives constributions
   * @param rRightHandSideVector: the condition right hand side vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
                                     ProcessInfo& rCurrentProcessInfo) override;


  /**
   * This function provides the place to perform checks on the completeness of the input.
   * It is designed to be called only once (or anyway, not often) typically at the beginning
   * of the calculations, so to verify that nothing is missing from the input
   * or that no common error is found.
   * @param rCurrentProcessInfo
   */
  virtual int Check(const ProcessInfo& rCurrentProcessInfo) override;


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

  StepType mStepVariable;

  ///@}
  ///@name Protected Operators
  ///@{

  ///@}
  ///@name Protected Operations
  ///@{

  /**
   * Sets process information to set member variables like mStepVariable
   */
  void SetProcessInformation(const ProcessInfo& rCurrentProcessInfo);


  /**
   * Get element size from the dofs
   */
  SizeType GetDofsSize() override;

  ///@}
  ///@name Protected  Access
  ///@{

  RigidBodyPointLinkSegregatedVCondition() {};

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

  void save( Serializer& rSerializer ) const override
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, RigidBodyPointLinkCondition )
  }

  void load( Serializer& rSerializer ) override
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, RigidBodyPointLinkCondition )
  }

}; // class RigidBodyPointLinkSegregatedVCondition.

} // namespace Kratos.

#endif // KRATOS_RIGID_BODY_POINT_LINK_SEGREGATED_V_CONDITION_H_INCLUDED defined

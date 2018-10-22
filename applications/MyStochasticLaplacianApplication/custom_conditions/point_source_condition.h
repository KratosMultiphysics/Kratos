//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    , KratosAppGenerator
//

#if !defined(KRATOS_POINT_SOURCE_CONDITION_H_INCLUDED )
#define KRATOS_POINT_SOURCE_CONDITION_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/condition.h"
#include "my_stochastic_laplacian_application_variables.h"


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


class PointSourceCondition : public Condition {
public:



  ///@name Type Definitions
  ///@{

  typedef Condition BaseType;

  ///@}
  ///@name Pointer Definitions
  /// Pointer definition of PointSourceCondition
  KRATOS_CLASS_POINTER_DEFINITION(PointSourceCondition);

  ///@}
  ///@name Life Cycle
  ///@{

  /**
   * Constructor.
   */
  PointSourceCondition(IndexType NewId = 0);

  /**
   * Constructor using an array of nodes
   */
  PointSourceCondition(IndexType NewId, const NodesArrayType& ThisNodes);

  /**
   * Constructor using Geometry
   */
  PointSourceCondition(IndexType NewId, GeometryType::Pointer pGeometry);

  /**
   * Constructor using Properties
   */
  PointSourceCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

  /**
   * Copy Constructor
   */
  PointSourceCondition(PointSourceCondition const& rOther);

  /**
   * Destructor
   */
  virtual ~PointSourceCondition();

  ///@}
  ///@name Operators
  ///@{

  /// Assignment operator.
  PointSourceCondition & operator=(PointSourceCondition const& rOther);

  ///@}
  ///@name Operations
  ///@{

  /**
   * CONDITIONS inherited from this class have to implement next
   * Create and Clone methods: MANDATORY
   */

  /**
   * creates a new condition pointer
   * @param NewId: the ID of the new condition
   * @param ThisNodes: the nodes of the new condition
   * @param pProperties: the properties assigned to the new condition
   * @return a Pointer to the new condition
   */
  Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

  /**
   * creates a new condition pointer
   * @param NewId: the ID of the new condition
   * @param pGeom: the geometry to be employed
   * @param pProperties: the properties assigned to the new condition
   * @return a Pointer to the new condition
   */
  Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const;

  /**
   * creates a new condition pointer and clones the previous condition data
   * @param NewId: the ID of the new condition
   * @param ThisNodes: the nodes of the new condition
   * @param pProperties: the properties assigned to the new condition
   * @return a Pointer to the new condition
   */
  Condition::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;

  /**
   * this determines the condition equation ID vector for all condition
   * DOFs
   * @param rResult: the condition equation ID vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo);

  /**
   * determines the condition list of DOFs
   * @param rConditionDofList: the list of DOFs
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& CurrentProcessInfo);

  /**
   * CONDITIONS inherited from this class have to implement next
   * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
   * they can be managed internally with a private method to do the same calculations
   * only once: MANDATORY
   */

  /**
   * this is called during the assembling process in order
   * to calculate all condition contributions to the global system
   * matrix and the right hand side
   * @param rLeftHandSideMatrix: the condition left hand side matrix
   * @param rRightHandSideVector: the condition right hand side
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateLocalSystem(
      MatrixType& rLeftHandSideMatrix,
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo);

  /**
   * this is called during the assembling process in order
   * to calculate the condition left hand side matrix only
   * @param rLeftHandSideMatrix: the condition left hand side matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

  /**
   * this is called during the assembling process in order
   * to calculate the condition right hand side vector only
   * @param rRightHandSideVector: the condition right hand side vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

  /**
   * this is called during the assembling process in order
   * to calculate the first derivatives contributions for the LHS and RHS
   * @param rLeftHandSideMatrix: the condition left hand side matrix
   * @param rRightHandSideVector: the condition right hand side
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateFirstDerivativesContributions(
      MatrixType& rLeftHandSideMatrix,
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo);

  /**
   * this is called during the assembling process in order
   * to calculate the condition left hand side matrix for the first derivatives constributions
   * @param rLeftHandSideMatrix: the condition left hand side matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

  /**
   * this is called during the assembling process in order
   * to calculate the condition right hand side vector for the first derivatives constributions
   * @param rRightHandSideVector: the condition right hand side vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

  /**
   * CONDITIONS inherited from this class must implement this methods
   * if they need to add dynamic condition contributions
   * note: second derivatives means the accelerations if the displacements are the dof of the analysis
   * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
   * CalculateSecondDerivativesContributions,
   * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
   */


 /**
   * this is called during the assembling process in order
   * to calculate the second derivative contributions for the LHS and RHS
   * @param rLeftHandSideMatrix: the condition left hand side matrix
   * @param rRightHandSideVector: the condition right hand side
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateSecondDerivativesContributions(
      MatrixType& rLeftHandSideMatrix,
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo);

  /**
   * this is called during the assembling process in order
   * to calculate the condition left hand side matrix for the second derivatives constributions
   * @param rLeftHandSideMatrix: the condition left hand side matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateSecondDerivativesLHS(
      MatrixType& rLeftHandSideMatrix,
      ProcessInfo& rCurrentProcessInfo);

  /**
   * this is called during the assembling process in order
   * to calculate the condition right hand side vector for the second derivatives constributions
   * @param rRightHandSideVector: the condition right hand side vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateSecondDerivativesRHS(
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo);

  /**
   * this is called during the assembling process in order
   * to calculate the condition mass matrix
   * @param rMassMatrix: the condition mass matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

  /**
   * this is called during the assembling process in order
   * to calculate the condition damping matrix
   * @param rDampingMatrix: the condition damping matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo);

  /**
   * This method provides the place to perform checks on the completeness of the input
   * and the compatibility with the problem options as well as the contitutive laws selected
   * It is designed to be called only once (or anyway, not often) typically at the beginning
   * of the calculations, so to verify that nothing is missing from the input
   * or that no common error is found.
   * @param rCurrentProcessInfo
   * this method is: MANDATORY
   */
  virtual int Check(const ProcessInfo& rCurrentProcessInfo);

  ///@}
  ///@name Access
  ///@{


  ///@}
  ///@name Inquiry
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  virtual std::string Info() const;

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const;

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const;

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
  ///@name Serialization
  ///@{

  friend class Serializer;

  virtual void save(Serializer& rSerializer) const;
  virtual void load(Serializer& rSerializer);

  ///@}
  ///@name Private  Access
  ///@{

  ///@}
  ///@name Private Inquiry
  ///@{

  ///@}
  ///@name Un accessible methods
  ///@{

  ///@}

}; // Class PointSourceCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_POINT_SOURCE_CONDITION_H_INCLUDED  defined

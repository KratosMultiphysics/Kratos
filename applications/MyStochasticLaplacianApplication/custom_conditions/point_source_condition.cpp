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

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Include Base h
#include "includes/checks.h"
#include "custom_conditions/point_source_condition.h"


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

/**
 * Constructor.
 */
PointSourceCondition::PointSourceCondition(IndexType NewId)
    : Condition(NewId)  {
}

/**
 * Constructor using an array of nodes
 */
PointSourceCondition::PointSourceCondition(IndexType NewId, const NodesArrayType& ThisNodes)
    : Condition(NewId, ThisNodes)  {
}

/**
 * Constructor using Geometry
 */
PointSourceCondition::PointSourceCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)  {
}

/**
 * Constructor using Properties
 */
PointSourceCondition::PointSourceCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)  {
}

/**
 * Copy Constructor
 */
PointSourceCondition::PointSourceCondition(PointSourceCondition const& rOther)
    : Condition(rOther)  {
}

/**
 * Destructor
 */
PointSourceCondition::~PointSourceCondition() {
}

///@}
///@name Operators
///@{

/// Assignment operator.
PointSourceCondition & PointSourceCondition::operator=(PointSourceCondition const& rOther) {
  BaseType::operator=(rOther);
  Flags::operator =(rOther);
  // mpProperties = rOther.mpProperties;
  return *this;
}

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
Condition::Pointer PointSourceCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const {

  KRATOS_TRY
  return Condition::Pointer(new PointSourceCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
  KRATOS_CATCH("");
}

/**
 * creates a new condition pointer
 * @param NewId: the ID of the new condition
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new condition
 * @return a Pointer to the new condition
 */
Condition::Pointer PointSourceCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const {

  KRATOS_TRY
  return Condition::Pointer(new PointSourceCondition(NewId, pGeom, pProperties));
  KRATOS_CATCH("");
}

/**
 * creates a new condition pointer and clones the previous condition data
 * @param NewId: the ID of the new condition
 * @param ThisNodes: the nodes of the new condition
 * @param pProperties: the properties assigned to the new condition
 * @return a Pointer to the new condition
 */
Condition::Pointer PointSourceCondition::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const {
  KRATOS_TRY
  return Condition::Pointer(new PointSourceCondition(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
  KRATOS_CATCH("");
}

/**
 * this determines the condition equation ID vector for all conditional
 * DOFs
 * @param rResult: the condition equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) {
  unsigned int number_of_nodes = GetGeometry().PointsNumber();
  if (rResult.size() != number_of_nodes)
    rResult.resize(number_of_nodes, false);

  for (unsigned int i = 0; i < number_of_nodes; i++)
    rResult[i] = GetGeometry()[i].GetDof(SOLUTION).EquationId();


}

/**
 * determines the condition equation list of DOFs
 * @param ConditionDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& CurrentProcessInfo) {
  unsigned int number_of_nodes = GetGeometry().PointsNumber();
  if (rConditionDofList.size() != number_of_nodes)
    rConditionDofList.resize(number_of_nodes);

  for (unsigned int i = 0; i < number_of_nodes; i++)
    rConditionDofList[i] = GetGeometry()[i].pGetDof(SOLUTION);


}

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
void PointSourceCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

    if(rLeftHandSideMatrix.size1() != 1)
        rLeftHandSideMatrix.resize(1,1,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(1,1);

    if(rRightHandSideVector.size() != 1)
        rRightHandSideVector.resize(1,false);
    double load = 0.0; // GetGeometry()[0].GetSolutionStepValue(FORCING);
    //double load = GetGeometry()[0].GetSolutionStepValue(FORCING);
    //std::cout<<load;
    rRightHandSideVector[0] = load;
    KRATOS_CATCH("")
}

/**
 * this is called during the assembling process in order
 * to calculate the condition left hand side matrix only
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) {
}

/**
 * this is called during the assembling process in order
 * to calculate the condition right hand side vector only
 * @param rRightHandSideVector: the condition right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
}

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rRightHandSideVector: the condition right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::CalculateFirstDerivativesContributions(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo) {

  if (rLeftHandSideMatrix.size1() != 0)
    rLeftHandSideMatrix.resize(0, 0, false);
  if (rRightHandSideVector.size() != 0)
    rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the condition left hand side matrix for the first derivatives constributions
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) {
  if (rLeftHandSideMatrix.size1() != 0)
    rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the condition right hand side vector for the first derivatives constributions
 * @param rRightHandSideVector: the condition right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
  if (rRightHandSideVector.size() != 0)
    rRightHandSideVector.resize(0, false);
}

/**
 * CONDITION inherited from this class must implement this methods
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
void PointSourceCondition::CalculateSecondDerivativesContributions(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo) {

  if (rLeftHandSideMatrix.size1() != 0)
    rLeftHandSideMatrix.resize(0, 0, false);
  if (rRightHandSideVector.size() != 0)
    rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the condition left hand side matrix for the second derivatives constributions
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo) {

  if (rLeftHandSideMatrix.size1() != 0)
    rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the condition right hand side vector for the second derivatives constributions
 * @param rRightHandSideVector: the condition right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::CalculateSecondDerivativesRHS(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo) {

  if (rRightHandSideVector.size() != 0)
    rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the condition mass matrix
 * @param rMassMatrix: the condition mass matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {
  if (rMassMatrix.size1() != 0)
    rMassMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the condition damping matrix
 * @param rDampingMatrix: the condition damping matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void PointSourceCondition::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {
  if (rDampingMatrix.size1() != 0)
    rDampingMatrix.resize(0, 0, false);
}

/**
 * This method provides the place to perform checks on the completeness of the input
 * and the compatibility with the problem options as well as the contitutive laws selected
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 * this method is: MANDATORY
 */
int PointSourceCondition::Check(const ProcessInfo& rCurrentProcessInfo) {

  KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Condition::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(SOLUTION)
    KRATOS_CHECK_VARIABLE_KEY(FORCING)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    const unsigned int number_of_points = GetGeometry().size();
    for ( unsigned int i = 0; i < number_of_points; i++ )
    {
        Node<3> &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(SOLUTION,rnode)
        KRATOS_CHECK_DOF_IN_NODE(SOLUTION,rnode)
    }

    return ierr;

    KRATOS_CATCH("");
}

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

std::string PointSourceCondition::Info() const {
  std::stringstream buffer;
  buffer << "PointSourceCondition #" << Id();
  return buffer.str();
}

/// Print information about this object.

void PointSourceCondition::PrintInfo(std::ostream& rOStream) const {
  rOStream << "PointSourceCondition #" << Id();
}

/// Print object's data.

void PointSourceCondition::PrintData(std::ostream& rOStream) const {
  pGetGeometry()->PrintData(rOStream);
}

///@}
///@name Friends
///@{

///@}

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

void PointSourceCondition::save(Serializer& rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );

  // List
  // To be completed with the class member list
}

void PointSourceCondition::load(Serializer& rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );

  // List
  // To be completed with the class member list
}

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
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >> (std::istream& rIStream, PointSourceCondition& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const PointSourceCondition& rThis) {
  rThis.PrintInfo(rOStream);
  rOStream << " : " << std::endl;
  rThis.PrintData(rOStream);
  return rOStream;
}

} // namespace Kratos.

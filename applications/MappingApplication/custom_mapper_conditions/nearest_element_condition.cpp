//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Include Base h
#include "nearest_element_condition.h"


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
 * Constructor using an array of nodes
 */
NearestElementCondition::NearestElementCondition(IndexType NewId, const NodesArrayType& ThisNodes)
    : Condition(NewId, ThisNodes) {
}

/**
 * Constructor using Geometry
 */
NearestElementCondition::NearestElementCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry) {
}

/**
 * Constructor using Properties
 */
NearestElementCondition::NearestElementCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties) {
}

/**
 * Copy Constructor
 */
NearestElementCondition::NearestElementCondition(NearestElementCondition const& rOther)
    : Condition(rOther) {
}

/**
 * Destructor
 */
NearestElementCondition::~NearestElementCondition() {
}

///@}
///@name Operators
///@{

/// Assignment operator.
NearestElementCondition & NearestElementCondition::operator=(NearestElementCondition const& rOther) {
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
Condition::Pointer NearestElementCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const {

  KRATOS_TRY
  return Condition::Pointer(new NearestElementCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
  KRATOS_CATCH("");
}

/**
 * creates a new condition pointer
 * @param NewId: the ID of the new condition
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new condition
 * @return a Pointer to the new condition
 */
Condition::Pointer NearestElementCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const {

  KRATOS_TRY
  return Condition::Pointer(new NearestElementCondition(NewId, pGeom, pProperties));
  KRATOS_CATCH("");
}

/**
 * creates a new condition pointer and clones the previous condition data
 * @param NewId: the ID of the new condition
 * @param ThisNodes: the nodes of the new condition
 * @param pProperties: the properties assigned to the new condition
 * @return a Pointer to the new condition
 */
Condition::Pointer NearestElementCondition::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const {
  KRATOS_TRY
  return Condition::Pointer(new NearestElementCondition(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
  KRATOS_CATCH("");
}

/**
 * this determines the condition equation ID vector for all conditional
 * DOFs
 * @param rResult: the condition equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) {
  unsigned int number_of_nodes = GetGeometry().PointsNumber();
  if (rResult.size() != number_of_nodes)
    rResult.resize(number_of_nodes, false);

}

/**
 * determines the condition equation list of DOFs
 * @param ConditionDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementCondition::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& CurrentProcessInfo) {
  unsigned int number_of_nodes = GetGeometry().PointsNumber();
  if (rConditionDofList.size() != number_of_nodes)
    rConditionDofList.resize(number_of_nodes);

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
void NearestElementCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo) {
}

/**
 * this is called during the assembling process in order
 * to calculate the condition left hand side matrix only
 * @param rLeftHandSideMatrix: the condition left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) {
}

/**
 * this is called during the assembling process in order
 * to calculate the condition right hand side vector only
 * @param rRightHandSideVector: the condition right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void NearestElementCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
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
int NearestElementCondition::Check(const ProcessInfo& rCurrentProcessInfo) {

  KRATOS_TRY

  if (this->Id() < 1) {
    KRATOS_THROW_ERROR(std::logic_error, "NearestElementCondition found with Id 0 or negative","")
  }

  return 0;

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

std::string NearestElementCondition::Info() const {
  std::stringstream buffer;
  buffer << "NearestElementCondition #" << Id();
  return buffer.str();
}

/// Print information about this object.

void NearestElementCondition::PrintInfo(std::ostream& rOStream) const {
  rOStream << "NearestElementCondition #" << Id();
}

/// Print object's data.

void NearestElementCondition::PrintData(std::ostream& rOStream) const {
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

void NearestElementCondition::save(Serializer& rSerializer) const {
    KRATOS_ERROR << "This Object cannot be serialized!" << std::endl;
}

void NearestElementCondition::load(Serializer& rSerializer) {
    KRATOS_ERROR << "This Object cannot be serialized!" << std::endl;
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
inline std::istream & operator >> (std::istream& rIStream, NearestElementCondition& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const NearestElementCondition& rThis) {
  rThis.PrintInfo(rOStream);
  rOStream << " : " << std::endl;
  rThis.PrintData(rOStream);
  return rOStream;
}

} // namespace Kratos.

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand@{KRATOS_APP_AUTHOR}
//

#if !defined(KRATOS_@{KRATOS_NAME_UPPER}_H_INCLUDED )
#define KRATOS_@{KRATOS_NAME_UPPER}_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/element.h"


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

@{KRATOS_CLASS_TEMPLATE}
class @{KRATOS_NAME_CAMEL} @{KRATOS_CLASS_BASE_HEADER} {
public:

@{KRATOS_CLASS_LOCAL_FLAGS}
  // KRATOS_DEFINE_LOCAL_FLAG(PERFORM_STEP1);
  // KRATOS_DEFINE_LOCAL_FLAG(DO_EXPENSIVE_CHECKS);

  ///@name Type Definitions
  ///@{

  typedef BaseType @{KRATOS_CLASS_BASE}

  ///@}
  ///@name Pointer Definitions
  /// Pointer definition of @{KRATOS_NAME_CAMEL}
  KRATOS_CLASS_POINTER_DEFINITION(@{KRATOS_NAME_CAMEL});

  ///@}
  ///@name Life Cycle
  ///@{

  /**
   * Constructor.
   */
  @{KRATOS_NAME_CAMEL}(IndexType NewId = 0)
      : BaseType(NewId)
      , Flags()
      , mpProperties(new PropertiesType) @{KRATOS_INIT_MEMBER_LIST} {
  }

  /**
   * Constructor using an array of nodes
   */
  @{KRATOS_NAME_CAMEL}(IndexType NewId, const NodesArrayType& ThisNodes)
      : BaseType(NewId,GeometryType::Pointer(new GeometryType(ThisNodes)))
      , Flags()
      , mpProperties(new PropertiesType) @{KRATOS_INIT_MEMBER_LIST} {
  }

  /**
   * Constructor using Geometry
   */
  @{KRATOS_NAME_CAMEL}(IndexType NewId, GeometryType::Pointer pGeometry)
      : BaseType(NewId,pGeometry)
      , Flags()
      , mpProperties(new PropertiesType) @{KRATOS_INIT_MEMBER_LIST} {
  }

  /**
   * Constructor using Properties
   */
  @{KRATOS_NAME_CAMEL}(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
      : BaseType(NewId,pGeometry)
      , Flags()
      , mpProperties(pProperties) @{KRATOS_INIT_MEMBER_LIST} {
  }

  /**
   * Copy Constructor
   */
  @{KRATOS_NAME_CAMEL}(@{KRATOS_NAME_CAMEL} const& rOther)
      : BaseType(rOther)
      , Flags(rOther)
      , mpProperties(rOther.mpProperties) @{KRATOS_CC_INIT_MEMBER_LIST} {
  }

  /**
   * Destructor
   */
  virtual ~@{KRATOS_NAME_CAMEL}() {
  }

  ///@}
  ///@name Operators
  ///@{

  /// Assignment operator.
  @{KRATOS_NAME_CAMEL} & operator=(@{KRATOS_NAME_CAMEL} const& rOther) {
    BaseType::operator=(rOther);
    Flags::operator =(rOther);
    mpProperties = rOther.mpProperties;
    return *this;
  }

  ///@}
  ///@name Operations
  ///@{

  /**
   * ELEMENTS inherited from this class have to implement next
   * Create and Clone methods: MANDATORY
   */

  /**
   * creates a new element pointer
   * @param NewId: the ID of the new element
   * @param ThisNodes: the nodes of the new element
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  virtual Pointer Create(
      IndexType NewId,
      NodesArrayType const& ThisNodes,
      PropertiesType::Pointer pProperties) const {

    KRATOS_TRY
    return Element::Pointer(new @{KRATOS_NAME_CAMEL}(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
  }

  /**
   * creates a new element pointer
   * @param NewId: the ID of the new element
   * @param pGeom: the geometry to be employed
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  virtual Pointer Create(
      IndexType NewId,
      GeometryType::Pointer pGeom,
      PropertiesType::Pointer pProperties) const {

    KRATOS_TRY
    return Element::Pointer(new @{KRATOS_NAME_CAMEL}(NewId, pGeom, pProperties));
    KRATOS_CATCH("");
  }

  /**
   * creates a new element pointer and clones the previous element data
   * @param NewId: the ID of the new element
   * @param ThisNodes: the nodes of the new element
   * @param pProperties: the properties assigned to the new element
   * @return a Pointer to the new element
   */
  virtual Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const {
    KRATOS_TRY
    return Element::Pointer(new @{KRATOS_NAME_CAMEL}(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
    KRATOS_CATCH("");
  }

  /**
   * this determines the elemental equation ID vector for all elemental
   * DOFs
   * @param rResult: the elemental equation ID vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) {
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if (rResult.size() != number_of_nodes)
      rResult.resize(number_of_nodes, false);

@{KRATOS_ELEMENT_ECUATION_ID_DOFS}
  }

  /**
   * determines the elemental list of DOFs
   * @param ElementalDofList: the list of DOFs
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if (ElementalDofList.size() != number_of_nodes)
      ElementalDofList.resize(number_of_nodes);

@{KRATOS_ELEMENT_LIST_DOFS}
  }

  /**
   * ELEMENTS inherited from this class have to implement next
   * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
   * they can be managed internally with a private method to do the same calculations
   * only once: MANDATORY
   */

  /**
   * this is called during the assembling process in order
   * to calculate all elemental contributions to the global system
   * matrix and the right hand side
   * @param rLeftHandSideMatrix: the elemental left hand side matrix
   * @param rRightHandSideVector: the elemental right hand side
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateLocalSystem(
      MatrixType& rLeftHandSideMatrix,
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo) {
  }

  /**
   * this function provides a more general interface to the element.
   * it is designed so that rLHSvariables and rRHSvariables are passed TO the element
   * thus telling what is the desired output
   * @param rLeftHandSideMatrices: container with the output left hand side matrices
   * @param rLHSVariables: paramter describing the expected LHSs
   * @param rRightHandSideVectors: container for the desired RHS output
   * @param rRHSVariables: parameter describing the expected RHSs
   */
  virtual void CalculateLocalSystem(
      std::vector< MatrixType >& rLeftHandSideMatrices, const std::vector< Variable< MatrixType > >& rLHSVariables,
      std::vector< VectorType >& rRightHandSideVectors, const std::vector< Variable< VectorType > >& rRHSVariables,
      ProcessInfo& rCurrentProcessInfo) {
  }

  /**
   * this is called during the assembling process in order
   * to calculate the elemental left hand side matrix only
   * @param rLeftHandSideMatrix: the elemental left hand side matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) {
  }

  /**
   * this function provides a more general interface to the element.
   * it is designed so that rLHSvariables are passed TO the element
   * thus telling what is the desired output
   * @param rLeftHandSideMatrices: container for the desired LHS output
   * @param rLHSVariables: parameter describing the expected LHSs
   */
  virtual void CalculateLeftHandSide(
      std::vector< MatrixType >& rLeftHandSideMatrices,
      const std::vector< Variable< MatrixType > >& rLHSVariables,
      ProcessInfo& rCurrentProcessInfo) {
  }

  /**
   * this is called during the assembling process in order
   * to calculate the elemental right hand side vector only
   * @param rRightHandSideVector: the elemental right hand side vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
  }

  /**
   * this function provides a more general interface to the element.
   * it is designed so that rRHSvariables are passed TO the element
   * thus telling what is the desired output
   * @param rRightHandSideVectors: container for the desired RHS output
   * @param rRHSVariables: parameter describing the expected RHSs
   */
  virtual void CalculateRightHandSide(
      std::vector< VectorType >& rRightHandSideVectors,
      const std::vector< Variable< VectorType > >& rRHSVariables,
      ProcessInfo& rCurrentProcessInfo) {
  }

  /**
   * this is called during the assembling process in order
   * to calculate the first derivatives contributions for the LHS and RHS
   * @param rLeftHandSideMatrix: the elemental left hand side matrix
   * @param rRightHandSideVector: the elemental right hand side
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateFirstDerivativesContributions(
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
   * to calculate the elemental left hand side matrix for the first derivatives constributions
   * @param rLeftHandSideMatrix: the elemental left hand side matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) {
    if (rLeftHandSideMatrix.size1() != 0)
      rLeftHandSideMatrix.resize(0, 0, false);
  }

  /**
   * this is called during the assembling process in order
   * to calculate the elemental right hand side vector for the first derivatives constributions
   * @param rRightHandSideVector: the elemental right hand side vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
    if (rRightHandSideVector.size() != 0)
      rRightHandSideVector.resize(0, false);
  }

  /**
   * ELEMENTS inherited from this class must implement this methods
   * if they need to add dynamic element contributions
   * note: second derivatives means the accelerations if the displacements are the dof of the analysis
   * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
   * CalculateSecondDerivativesContributions,
   * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
   */


 /**
   * this is called during the assembling process in order
   * to calculate the second derivative contributions for the LHS and RHS
   * @param rLeftHandSideMatrix: the elemental left hand side matrix
   * @param rRightHandSideVector: the elemental right hand side
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateSecondDerivativesContributions(
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
   * to calculate the elemental left hand side matrix for the second derivatives constributions
   * @param rLeftHandSideMatrix: the elemental left hand side matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateSecondDerivativesLHS(
      MatrixType& rLeftHandSideMatrix,
      ProcessInfo& rCurrentProcessInfo) {

    if (rLeftHandSideMatrix.size1() != 0)
      rLeftHandSideMatrix.resize(0, 0, false);
  }

  /**
   * this is called during the assembling process in order
   * to calculate the elemental right hand side vector for the second derivatives constributions
   * @param rRightHandSideVector: the elemental right hand side vector
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateSecondDerivativesRHS(
      VectorType& rRightHandSideVector,
      ProcessInfo& rCurrentProcessInfo) {

    if (rRightHandSideVector.size() != 0)
      rRightHandSideVector.resize(0, false);
  }

  /**
   * this is called during the assembling process in order
   * to calculate the elemental mass matrix
   * @param rMassMatrix: the elemental mass matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {
    if (rMassMatrix.size1() != 0)
      rMassMatrix.resize(0, 0, false);
  }

  /**
   * this is called during the assembling process in order
   * to calculate the elemental damping matrix
   * @param rDampingMatrix: the elemental damping matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  virtual void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {
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
  virtual int Check(const ProcessInfo& rCurrentProcessInfo) {

    KRATOS_TRY

    if (this->Id() < 1) {
      KRATOS_THROW_ERROR(std::logic_error, "@{KRATOS_NAME_CAMEL} found with Id 0 or negative","")
    }

    if (this->GetGeometry().Area() <= 0) {
      std::cout << "error on @{KRATOS_NAME_CAMEL} -> " << this->Id() << std::endl;
      KRATOS_THROW_ERROR(std::logic_error, "Area cannot be less than or equal to 0","")
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

  virtual std::string Info() const {
    std::stringstream buffer;
    buffer << "@{KRATOS_NAME_CAMEL} #" << Id();
    return buffer.str();
  }

  /// Print information about this object.

  virtual void PrintInfo(std::ostream& rOStream) const {
    rOStream << "@{KRATOS_NAME_CAMEL} #" << Id();
  }

  /// Print object's data.

  virtual void PrintData(std::ostream& rOStream) const {
    pGetGeometry()->PrintData(rOStream);
  }

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

@{KRATOS_STATIC_MEMBERS_LIST}

  ///@}
  ///@name Member Variables
  ///@{

@{KRATOS_MEMBERS_LIST}

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

  virtual void save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, @{KRATOS_CLASS_BASE} );

    // List
    rSerializer.save("Data", mData);
  }

  virtual void load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, @{KRATOS_CLASS_BASE} );

    // List
    rSerializer.load("Data", mData);
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

  /// Copy constructor.
  //@{KRATOS_NAME_CAMEL}(@{KRATOS_NAME_CAMEL} const& rOther);

  ///@}

}; // Class @{KRATOS_NAME_CAMEL}

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >> (std::istream& rIStream, @{KRATOS_NAME_CAMEL}& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const @{KRATOS_NAME_CAMEL}& rThis) {
  rThis.PrintInfo(rOStream);
  rOStream << " : " << std::endl;
  rThis.PrintData(rOStream);
  return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_@{KRATOS_NAME_UPPER}_H_INCLUDED  defined

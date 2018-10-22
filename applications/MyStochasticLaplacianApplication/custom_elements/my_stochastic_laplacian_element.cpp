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
#include "utilities/geometry_utilities.h"
#include "custom_elements/my_stochastic_laplacian_element.h"


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
MyStochasticLaplacianElement::MyStochasticLaplacianElement(IndexType NewId)
    : Element(NewId)  {
}

/**
 * Constructor using an array of nodes
 */
MyStochasticLaplacianElement::MyStochasticLaplacianElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)  {
}

/**
 * Constructor using Geometry
 */
MyStochasticLaplacianElement::MyStochasticLaplacianElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)  {
}

/**
 * Constructor using Properties
 */
MyStochasticLaplacianElement::MyStochasticLaplacianElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)  {
}

/**
 * Copy Constructor
 */
MyStochasticLaplacianElement::MyStochasticLaplacianElement(MyStochasticLaplacianElement const& rOther)
    : Element(rOther)  {
}

/**
 * Destructor
 */
MyStochasticLaplacianElement::~MyStochasticLaplacianElement() {
}

///@}
///@name Operators
///@{

/// Assignment operator.
MyStochasticLaplacianElement & MyStochasticLaplacianElement::operator=(MyStochasticLaplacianElement const& rOther) {
  BaseType::operator=(rOther);
  Flags::operator =(rOther);
  // mpProperties = rOther.mpProperties;
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
Element::Pointer MyStochasticLaplacianElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const {

  KRATOS_TRY
  return Element::Pointer(new MyStochasticLaplacianElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
  KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer MyStochasticLaplacianElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const {

  KRATOS_TRY
  return Element::Pointer(new MyStochasticLaplacianElement(NewId, pGeom, pProperties));
  KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer MyStochasticLaplacianElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const {
  KRATOS_TRY
  return Element::Pointer(new MyStochasticLaplacianElement(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
  KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void MyStochasticLaplacianElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) {
  unsigned int number_of_nodes = GetGeometry().PointsNumber();
  if (rResult.size() != number_of_nodes)
    rResult.resize(number_of_nodes, false);

  for (unsigned int i = 0; i < number_of_nodes; i++)
    rResult[i] = GetGeometry()[i].GetDof(SOLUTION).EquationId();


}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void MyStochasticLaplacianElement::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) {
  unsigned int number_of_nodes = GetGeometry().PointsNumber();
  if (rElementalDofList.size() != number_of_nodes)
    rElementalDofList.resize(number_of_nodes);

  for (unsigned int i = 0; i < number_of_nodes; i++)
    rElementalDofList[i] = GetGeometry()[i].pGetDof(SOLUTION);


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
void MyStochasticLaplacianElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

    BoundedMatrix<double,3,2> DN_DX;  // Gradients matrix 
    BoundedMatrix<double,2,2> D;      // Conductivity matrix 
    D = ZeroMatrix(2,2); //initializing the matrix as zero
    array_1d<double,3> N; //dimension = number of nodes . Position of the gauss point 
    array_1d<double,3> temp; //dimension = number of nodes . . since we are using a residualbased approach 

    const unsigned int number_of_points = GetGeometry().size();
    // we have: number_of_points = 3 // we have a triangular mesh

    if(rLeftHandSideMatrix.size1() != number_of_points)
        rLeftHandSideMatrix.resize(number_of_points,number_of_points,false); //resizing the system in case it does not have the right size 

    if(rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points,false);

    // Getting data for the given geometry
    double area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area); //asking for gradients and other info
    // DN_DX : a stack matrix of size 3*2 to store the shape function's derivatives
    // N : an array_1d to store the shape functions at baricenter
    // area : the volume of the element

    // Reading properties and conditions
    const double integrated_permittivity = area * GetProperties()[CONDUCTIVITY];
    D(0,0)=integrated_permittivity;
    D(1,1)=integrated_permittivity;

    // Main loop (one Gauss point)
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());
    
    noalias(rLeftHandSideMatrix) = prod(DN_DX, Matrix(prod(D, trans(DN_DX))));  // Bt D B
    
    // Compute the forcing solving integral using ONE Gauss point
    // f_i = int_{Triang} forcing * shapefunction ~= sum_{l=1}**{ngauss} gauss_weight_l * forcing(x_l) * shapefunction(x_l)
    // forcing(x_l) = N0(x_l) * f0 + N1(x_l) * f1 + N2(x_l) * f2
    // 
    // RHS = forcing * N * det(J) * Gauss_weight = forcing * N * area
    // recall in case number of Gauss integration points is 1 we have det(J) = 2 * area_element and Gauss_weight = 0.5
    double forcing_gauss_point;
    forcing_gauss_point = N[0] * GetGeometry()[0].GetSolutionStepValue(FORCING) + N[1] * GetGeometry()[1].GetSolutionStepValue(FORCING) + N[2] * GetGeometry()[2].GetSolutionStepValue(FORCING);
    rRightHandSideVector[0] = area * forcing_gauss_point * N[0];
    rRightHandSideVector[1] = area * forcing_gauss_point * N[1];
    rRightHandSideVector[2] = area * forcing_gauss_point * N[2];

    //KRATOS_WATCH(gauss_weight)
    //KRATOS_WATCH(rRightHandSideVector)
    

    // Subtracting the dirichlet term
    // RHS -= LHS*DUMMY_UNKNOWNs
    for(unsigned int iii = 0; iii<number_of_points; iii++)
        temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(SOLUTION);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

    KRATOS_CATCH("");
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix only
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void MyStochasticLaplacianElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) {
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void MyStochasticLaplacianElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
}

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void MyStochasticLaplacianElement::CalculateFirstDerivativesContributions(
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
void MyStochasticLaplacianElement::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) {
  if (rLeftHandSideMatrix.size1() != 0)
    rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the first derivatives constributions
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void MyStochasticLaplacianElement::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
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
void MyStochasticLaplacianElement::CalculateSecondDerivativesContributions(
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
void MyStochasticLaplacianElement::CalculateSecondDerivativesLHS(
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
void MyStochasticLaplacianElement::CalculateSecondDerivativesRHS(
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
void MyStochasticLaplacianElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {
  if (rMassMatrix.size1() != 0)
    rMassMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental damping matrix
 * @param rDampingMatrix: the elemental damping matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void MyStochasticLaplacianElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {
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
int MyStochasticLaplacianElement::Check(const ProcessInfo& rCurrentProcessInfo) {

  KRATOS_TRY
  
      // Base class checks for positive Jacobian and Id > 0
      int ierr = Element::Check(rCurrentProcessInfo);
      if(ierr != 0) return ierr;
  
      // Check that all required variables have been registered
      KRATOS_CHECK_VARIABLE_KEY(SOLUTION)
  
      unsigned const int number_of_points = GetGeometry().size();  //added cornejo
      // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
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

std::string MyStochasticLaplacianElement::Info() const {
  std::stringstream buffer;
  buffer << "MyStochasticLaplacianElement #" << Id();
  return buffer.str();
}

/// Print information about this object.

void MyStochasticLaplacianElement::PrintInfo(std::ostream& rOStream) const {
  rOStream << "MyStochasticLaplacianElement #" << Id();
}

/// Print object's data.

void MyStochasticLaplacianElement::PrintData(std::ostream& rOStream) const {
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

void MyStochasticLaplacianElement::save(Serializer& rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

  // List
  // To be completed with the class member list
}

void MyStochasticLaplacianElement::load(Serializer& rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );

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
inline std::istream & operator >> (std::istream& rIStream, MyStochasticLaplacianElement& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const MyStochasticLaplacianElement& rThis) {
  rThis.PrintInfo(rOStream);
  rOStream << " : " << std::endl;
  rThis.PrintData(rOStream);
  return rOStream;
}

} // namespace Kratos.

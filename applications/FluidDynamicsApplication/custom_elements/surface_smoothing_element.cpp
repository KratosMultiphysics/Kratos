//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    KratosAppGenerator
//

// System includes
// see the header file

// External includes
// see the header file

// Include Base headers
// see the header file

// Project includes
#include "custom_elements/surface_smoothing_element.h"

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
SurfaceSmoothingElement::SurfaceSmoothingElement(IndexType NewId)
    : Element(NewId) 
{
}

/**
 * Constructor using an array of nodes
 */
SurfaceSmoothingElement::SurfaceSmoothingElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes) 
{
}

/**
 * Constructor using Geometry
 */
SurfaceSmoothingElement::SurfaceSmoothingElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) 
{
}

/**
 * Constructor using Properties
 */
SurfaceSmoothingElement::SurfaceSmoothingElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) 
{
}

/**
 * Copy Constructor
 */
SurfaceSmoothingElement::SurfaceSmoothingElement(SurfaceSmoothingElement const& rOther)
    : Element(rOther) 
{
}

/**
 * Destructor
 */
SurfaceSmoothingElement::~SurfaceSmoothingElement()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
SurfaceSmoothingElement & SurfaceSmoothingElement::operator=(SurfaceSmoothingElement const& rOther)
{
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
Element::Pointer SurfaceSmoothingElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<SurfaceSmoothingElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer SurfaceSmoothingElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<SurfaceSmoothingElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer SurfaceSmoothingElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<SurfaceSmoothingElement>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;

    // num_dof = num_nodes
    if (rResult.size() != num_nodes){
        rResult.resize(num_nodes, false);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rResult[i]  =  this->GetGeometry()[i].GetDof(DISTANCE_AUX).EquationId();
    }
}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;

    // num_dof = num_nodes
    if (rElementalDofList.size() != num_nodes){
        rElementalDofList.resize(num_nodes);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rElementalDofList[i] = this->GetGeometry()[i].pGetDof(DISTANCE_AUX);
    }
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
void SurfaceSmoothingElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;

    const double epsilon = 1.0e-8;

    BoundedMatrix<double,num_nodes,num_dim> DN_DX;  // Gradients matrix 
    array_1d<double,num_nodes> N; //dimension = number of nodes . Position of the gauss point 
    array_1d<double,num_nodes> tempVdof; //dimension = number of DOFs . . since we are using a residualbased approach
    array_1d<double,num_nodes> tempVold; //dimension = number of DOFs . . since we are using a residualbased approach
    //array_1d<double,num_nodes> tempRHS = ZeroVector(num_nodes);  
    BoundedMatrix<double,num_nodes,num_nodes> tempM;
    tempM = ZeroMatrix(num_nodes,num_nodes);
    BoundedMatrix<double,num_nodes,num_nodes> tempA;
    tempA = ZeroMatrix(num_nodes,num_nodes);

    // num_dof = num_nodes
    if(rLeftHandSideMatrix.size1() != num_nodes)
        rLeftHandSideMatrix.resize(num_nodes,num_nodes,false); //resizing the system in case it does not have the right size 

    if(rRightHandSideVector.size() != num_nodes)
        rRightHandSideVector.resize(num_nodes,false);

    // Getting data for the given geometry
    double area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area); //asking for gradients and other info 

    // Main loop (one Gauss point)
    //const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
    //noalias(rLeftHandSideMatrix) = prod(DN_DX, Matrix(prod(D, trans(DN_DX))));  // Bt D B

    // Subtracting the dirichlet term
    // RHS -= LHS*DUMMY_UNKNOWNs

    for(unsigned int i = 0; i<num_nodes; i++){
        tempVold[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        tempVdof[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE_AUX);

        for(unsigned int j = 0; j<num_nodes; j++){
            tempM(i,j) = area*N[i]*N[j];

            for (unsigned int k = 0; k<num_dim; k++){
                tempA(i,j) += area*epsilon*DN_DX(i,k)*DN_DX(j,k);
            }

            //tempRHS(i) += epsilon*DN_DX(j,0)*GetGeometry()[j].FastGetSolutionStepValue(DISTANCE_GRADIENT_X)*area*N[i];
            //tempRHS(i) += epsilon*DN_DX(j,1)*GetGeometry()[j].FastGetSolutionStepValue(DISTANCE_GRADIENT_Y)*area*N[i];
            //tempRHS(i) += epsilon*DN_DX(j,2)*GetGeometry()[j].FastGetSolutionStepValue(DISTANCE_GRADIENT_Z)*area*N[i];
        }
    }
    noalias(rLeftHandSideMatrix) = tempM + tempA;
    noalias(rRightHandSideVector) = prod(tempM,tempVold) /* + prod(tempA,tempVold) */ /* + tempRHS */ - prod(rLeftHandSideMatrix,tempVdof); //Phi_smooth + epsilon*laplacian(Phi_smooth) = Phi_old

    KRATOS_CATCH("");
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix only
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::CalculateFirstDerivativesContributions(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
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
void SurfaceSmoothingElement::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the first derivatives constributions
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
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
void SurfaceSmoothingElement::CalculateSecondDerivativesContributions(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
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
void SurfaceSmoothingElement::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the second derivatives constributions
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::CalculateSecondDerivativesRHS(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental mass matrix
 * @param rMassMatrix: the elemental mass matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    if (rMassMatrix.size1() != 0)
        rMassMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental damping matrix
 * @param rDampingMatrix: the elemental damping matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void SurfaceSmoothingElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
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
int SurfaceSmoothingElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1) <<"SurfaceSmoothingElement found with Id 0 or negative" << std::endl;

    KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0) << "On SurfaceSmoothingElement -> "
        << this->Id() <<  "; Area cannot be less than or equal to 0" << std::endl;
  
      // Base class checks for positive Jacobian and Id > 0
      int ierr = Element::Check(rCurrentProcessInfo);
      if(ierr != 0) return ierr;
  
      // Check that all required variables have been registered
      KRATOS_CHECK_VARIABLE_KEY(DISTANCE)
  
      unsigned const int number_of_points = GetGeometry().size(); 
      // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
      for ( unsigned int i = 0; i < number_of_points; i++ )
      {
          Node<3> &rnode = this->GetGeometry()[i];
          KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE,rnode)
          KRATOS_CHECK_DOF_IN_NODE(DISTANCE,rnode)
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

std::string SurfaceSmoothingElement::Info() const {
    std::stringstream buffer;
    buffer << "SurfaceSmoothingElement #" << Id();
    return buffer.str();
}

/// Print information about this object.

void SurfaceSmoothingElement::PrintInfo(std::ostream& rOStream) const {
    rOStream << "SurfaceSmoothingElement #" << Id();
}

/// Print object's data.

void SurfaceSmoothingElement::PrintData(std::ostream& rOStream) const {
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

void SurfaceSmoothingElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

    // List
    // To be completed with the class member list
}

void SurfaceSmoothingElement::load(Serializer& rSerializer)
{
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
inline std::istream & operator >> (std::istream& rIStream, SurfaceSmoothingElement& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const SurfaceSmoothingElement& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.
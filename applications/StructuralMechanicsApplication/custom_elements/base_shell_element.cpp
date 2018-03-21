//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (based on the work of Massimo Petracca and Peter Wilson)
//

// System includes

// External includes


// Project includes
#include "custom_elements/base_shell_element.h"


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
 * Constructor using Geometry
 */
BaseShellElement::BaseShellElement(IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    SetBaseMembers();
}

/**
 * Constructor using Properties
 */
BaseShellElement::BaseShellElement(IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    SetBaseMembers();
}

/**
 * Destructor
 */
BaseShellElement::~BaseShellElement() {
}

///@}
///@name Operators
///@{


///@}
///@name Operations
///@{

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void BaseShellElement::EquationIdVector(EquationIdVectorType& rResult,
                                        ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != mNumDofs)
        rResult.resize(mNumDofs, false);

    GeometryType& geom = GetGeometry();

    for (SizeType i = 0; i < geom.size(); ++i)
    {
        const SizeType index = i * 6;
        NodeType& i_node = geom[i];

        rResult[index]     = i_node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = i_node.GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = i_node.GetDof(DISPLACEMENT_Z).EquationId();

        rResult[index + 3] = i_node.GetDof(ROTATION_X).EquationId();
        rResult[index + 4] = i_node.GetDof(ROTATION_Y).EquationId();
        rResult[index + 5] = i_node.GetDof(ROTATION_Z).EquationId();
    }
}

/**
 * determines the elemental list of DOFs
 * @param rElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void BaseShellElement::GetDofList(DofsVectorType& rElementalDofList,
                                  ProcessInfo& rCurrentProcessInfo)
{
    rElementalDofList.resize(0);
    rElementalDofList.reserve(mNumDofs);

    GeometryType& geom = GetGeometry();

    for (SizeType i = 0; i < geom.size(); ++i)
    {
        NodeType& i_node = geom[i];

        rElementalDofList.push_back(i_node.pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(i_node.pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(i_node.pGetDof(DISPLACEMENT_Z));

        rElementalDofList.push_back(i_node.pGetDof(ROTATION_X));
        rElementalDofList.push_back(i_node.pGetDof(ROTATION_Y));
        rElementalDofList.push_back(i_node.pGetDof(ROTATION_Z));
    }
}

void BaseShellElement::GetValuesVector(Vector& rValues, int Step)
{
    if (rValues.size() != mNumDofs)
        rValues.resize(mNumDofs, false);

    const GeometryType& geom = GetGeometry();

    for (SizeType i = 0; i < geom.size(); ++i)
    {
        const NodeType& i_node = geom[i];
        const array_1d<double, 3>& disp = i_node.FastGetSolutionStepValue(DISPLACEMENT, Step);
        const array_1d<double, 3>& rot = i_node.FastGetSolutionStepValue(ROTATION, Step);

        const SizeType index = i * 6;
        rValues[index]     = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];

        rValues[index + 3] = rot[0];
        rValues[index + 4] = rot[1];
        rValues[index + 5] = rot[2];
    }
}

void BaseShellElement::GetFirstDerivativesVector(Vector& rValues, int Step)
{
    if (rValues.size() != mNumDofs)
        rValues.resize(mNumDofs, false);

    const GeometryType& geom = GetGeometry();

    for (SizeType i = 0; i < geom.size(); ++i)
    {
        const NodeType& i_node = geom[i];
        const array_1d<double, 3>& vel = i_node.FastGetSolutionStepValue(VELOCITY, Step);
        // TODO also include the angular velocity

        const SizeType index = i * 6;
        rValues[index]     = vel[0];
        rValues[index + 1] = vel[1];
        rValues[index + 2] = vel[2];
        rValues[index + 3] = 0.0;
        rValues[index + 4] = 0.0;
        rValues[index + 5] = 0.0;
    }
}

void BaseShellElement::GetSecondDerivativesVector(Vector& rValues, int Step)
{
    if (rValues.size() != mNumDofs)
        rValues.resize(mNumDofs, false);

    const GeometryType & geom = GetGeometry();

    for (SizeType i = 0; i < geom.size(); ++i)
    {
        const NodeType& i_node = geom[i];
        const array_1d<double, 3>& acc = i_node.FastGetSolutionStepValue(ACCELERATION, Step);
        // TODO also include the angular acceleration

        const SizeType index = i * 6;
        rValues[index]     = acc[0];
        rValues[index + 1] = acc[1];
        rValues[index + 2] = acc[2];
        rValues[index + 3] = 0.0;
        rValues[index + 4] = 0.0;
        rValues[index + 5] = 0.0;
    }
}

void BaseShellElement::ResetConstitutiveLaw()
{
    KRATOS_TRY

    const GeometryType& geom = GetGeometry();
    const Matrix& shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

    const Properties& props = GetProperties();
    for(SizeType i = 0; i < mSections.size(); i++)
        mSections[i]->ResetCrossSection(props, geom, row(shapeFunctionsValues, i));

    KRATOS_CATCH("")
}


void BaseShellElement::BaseInitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& geom = this->GetGeometry();
    const Matrix& shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());
    for (SizeType i = 0; i < mSections.size(); ++i)
        mSections[i]->InitializeNonLinearIteration(GetProperties(), geom, row(shapeFunctionsValues, i), rCurrentProcessInfo);
}

void BaseShellElement::BaseFinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& geom = this->GetGeometry();
    const Matrix& shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());
    for (SizeType i = 0; i < mSections.size(); ++i)
        mSections[i]->FinalizeNonLinearIteration(GetProperties(), geom, row(shapeFunctionsValues, i), rCurrentProcessInfo);
}

void BaseShellElement::BaseInitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
	const PropertiesType& props = GetProperties();
	const GeometryType & geom = GetGeometry();
	const Matrix& shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

	for (SizeType i = 0; i < mSections.size(); ++i)
		mSections[i]->InitializeSolutionStep(props, geom, row(shapeFunctionsValues, i), rCurrentProcessInfo);
}


void BaseShellElement::BaseFinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    const PropertiesType& props = GetProperties();
    const GeometryType& geom = GetGeometry();
    const Matrix& shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

    for (SizeType i = 0; i < mSections.size(); i++)
        mSections[i]->FinalizeSolutionStep(props, geom, row(shapeFunctionsValues, i), rCurrentProcessInfo);
}

void BaseShellElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
	VectorType& rRightHandSideVector,
	ProcessInfo& rCurrentProcessInfo)
{
    // Calculation flags
    const bool calculate_stiffness_matrix_flag = true;
    const bool calculate_residual_vector_flag = true;

	CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                 calculate_stiffness_matrix_flag, calculate_residual_vector_flag);
}

void BaseShellElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
	ProcessInfo& rCurrentProcessInfo)
{
    // Calculation flags
    const bool CalculateStiffnessMatrixFlag = true; // TODO check is this can be false => see solids
    const bool CalculateResidualVectorFlag = true;

	Matrix dummy;
	CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo,
                 calculate_stiffness_matrix_flag, calculate_residual_vector_flag);
}

// /**
//  * ELEMENTS inherited from this class have to implement next
//  * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
//  * they can be managed internally with a private method to do the same calculations
//  * only once: MANDATORY
//  */

// /**
//  * this is called during the assembling process in order
//  * to calculate all elemental contributions to the global system
//  * matrix and the right hand side
//  * @param rLeftHandSideMatrix: the elemental left hand side matrix
//  * @param rRightHandSideVector: the elemental right hand side
//  * @param rCurrentProcessInfo: the current process info instance
//  */
// void BaseShellElement::CalculateLocalSystem(
//     MatrixType& rLeftHandSideMatrix,
//     VectorType& rRightHandSideVector,
//     ProcessInfo& rCurrentProcessInfo) {
// }

// /**
//  * this is called during the assembling process in order
//  * to calculate the elemental left hand side matrix only
//  * @param rLeftHandSideMatrix: the elemental left hand side matrix
//  * @param rCurrentProcessInfo: the current process info instance
//  */
// void BaseShellElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) {
// }

// /**
//  * this is called during the assembling process in order
//  * to calculate the elemental right hand side vector only
//  * @param rRightHandSideVector: the elemental right hand side vector
//  * @param rCurrentProcessInfo: the current process info instance
//  */
// void BaseShellElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
// }

// /**
//  * this is called during the assembling process in order
//  * to calculate the first derivatives contributions for the LHS and RHS
//  * @param rLeftHandSideMatrix: the elemental left hand side matrix
//  * @param rRightHandSideVector: the elemental right hand side
//  * @param rCurrentProcessInfo: the current process info instance
//  */
// void BaseShellElement::CalculateFirstDerivativesContributions(
//     MatrixType& rLeftHandSideMatrix,
//     VectorType& rRightHandSideVector,
//     ProcessInfo& rCurrentProcessInfo) {

//   if (rLeftHandSideMatrix.size1() != 0)
//     rLeftHandSideMatrix.resize(0, 0, false);
//   if (rRightHandSideVector.size() != 0)
//     rRightHandSideVector.resize(0, false);
// }

// /**
//  * this is called during the assembling process in order
//  * to calculate the elemental left hand side matrix for the first derivatives constributions
//  * @param rLeftHandSideMatrix: the elemental left hand side matrix
//  * @param rCurrentProcessInfo: the current process info instance
//  */
// void BaseShellElement::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) {
//   if (rLeftHandSideMatrix.size1() != 0)
//     rLeftHandSideMatrix.resize(0, 0, false);
// }

// /**
//  * this is called during the assembling process in order
//  * to calculate the elemental right hand side vector for the first derivatives constributions
//  * @param rRightHandSideVector: the elemental right hand side vector
//  * @param rCurrentProcessInfo: the current process info instance
//  */
// void BaseShellElement::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {
//   if (rRightHandSideVector.size() != 0)
//     rRightHandSideVector.resize(0, false);
// }

// /**
//  * ELEMENTS inherited from this class must implement this methods
//  * if they need to add dynamic element contributions
//  * note: second derivatives means the accelerations if the displacements are the dof of the analysis
//  * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
//  * CalculateSecondDerivativesContributions,
//  * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
//  */

// /**
//  * this is called during the assembling process in order
//  * to calculate the elemental right hand side vector for the second derivatives constributions
//  * @param rRightHandSideVector: the elemental right hand side vector
//  * @param rCurrentProcessInfo: the current process info instance
//  */
// void BaseShellElement::CalculateSecondDerivativesRHS(
//     VectorType& rRightHandSideVector,
//     ProcessInfo& rCurrentProcessInfo) {

//   if (rRightHandSideVector.size() != 0)
//     rRightHandSideVector.resize(0, false);
// }

// /**
//  * this is called during the assembling process in order
//  * to calculate the elemental mass matrix
//  * @param rMassMatrix: the elemental mass matrix
//  * @param rCurrentProcessInfo: the current process info instance
//  */
// void BaseShellElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {
//   if (rMassMatrix.size1() != 0)
//     rMassMatrix.resize(0, 0, false);
// }

void BaseShellElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    if ( rDampingMatrix.size1() != mNumDofs )
        rDampingMatrix.resize( mNumDofs, mNumDofs, false );

    noalias( rDampingMatrix ) = ZeroMatrix( mNumDofs, mNumDofs );

    // 1.-Calculate StiffnessMatrix:

    MatrixType StiffnessMatrix  = Matrix();
    VectorType ResidualVector  = Vector();

    CalculateAll(StiffnessMatrix, ResidualVector, rCurrentProcessInfo, true, false);

    // 2.-Calculate MassMatrix:

    MatrixType MassMatrix  = Matrix();

    CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

    // 3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0.0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) )
        alpha = GetProperties()[RAYLEIGH_ALPHA];
    else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) )
        alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];

    double beta  = 0.0;
    if( GetProperties().Has(RAYLEIGH_BETA) )
        beta = GetProperties()[RAYLEIGH_BETA];
    else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
        beta = rCurrentProcessInfo[RAYLEIGH_BETA];

    // 4.-Compose the Damping Matrix:

    // Rayleigh Damping Matrix: alpha*M + beta*K
    noalias( rDampingMatrix ) += alpha * MassMatrix;
    noalias( rDampingMatrix ) += beta  * StiffnessMatrix;

    KRATOS_CATCH( "" )
}

// /**
//  * This method provides the place to perform checks on the completeness of the input
//  * and the compatibility with the problem options as well as the contitutive laws selected
//  * It is designed to be called only once (or anyway, not often) typically at the beginning
//  * of the calculations, so to verify that nothing is missing from the input
//  * or that no common error is found.
//  * @param rCurrentProcessInfo
//  * this method is: MANDATORY
//  */
// int BaseShellElement::Check(const ProcessInfo& rCurrentProcessInfo) {

//   KRATOS_TRY

//   if (this->Id() < 1) {
//     KRATOS_THROW_ERROR(std::logic_error, "BaseShellElement found with Id 0 or negative","")
//   }

//   if (this->GetGeometry().Area() <= 0) {
//     std::cout << "error on BaseShellElement -> " << this->Id() << std::endl;
//     KRATOS_THROW_ERROR(std::logic_error, "Area cannot be less than or equal to 0","")
//   }

//   return 0;

//   KRATOS_CATCH("");
// }

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

std::string BaseShellElement::Info() const {
  std::stringstream buffer;
  buffer << "BaseShellElement #" << Id();
  return buffer.str();
}

/// Print information about this object.

void BaseShellElement::PrintInfo(std::ostream& rOStream) const {
  rOStream << "BaseShellElement #" << Id();
}

/// Print object's data.

void BaseShellElement::PrintData(std::ostream& rOStream) const {
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

void BaseShellElement::SetBaseMembers()
{
    mNumDofs = 6 * GetGeometry().PointsNumber(); // 6 dofs per node

    const GeometryType::IntegrationPointsArrayType& integrationPoints =
    GetGeometry().IntegrationPoints(mIntegrationMethod);

    mNumGPs = integrationPoints.size();
}

void BaseShellElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_ERROR << "You have called to the CalculateAll from the base class for shell elements" << std::endl;
}

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

void BaseShellElement::save(Serializer& rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
  rSerializer.save("Sections", mSections);
  rSerializer.save("NumDofs", mNumDofs);
  rSerializer.save("NumGPs", mNumGPs);
  rSerializer.save("IntM", (int)mIntegrationMethod);
}

void BaseShellElement::load(Serializer& rSerializer) {
  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
  rSerializer.load("Sections", mSections);
  rSerializer.load("NumDofs", mNumDofs);
  rSerializer.load("NumGPs", mNumGPs);
  int temp;
  rSerializer.load("IntM", temp);
  mIntegrationMethod = (IntegrationMethod)temp;
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
inline std::istream & operator >> (std::istream& rIStream, BaseShellElement& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const BaseShellElement& rThis) {
  rThis.PrintInfo(rOStream);
  rOStream << " : " << std::endl;
  rThis.PrintData(rOStream);
  return rOStream;
}

} // namespace Kratos.

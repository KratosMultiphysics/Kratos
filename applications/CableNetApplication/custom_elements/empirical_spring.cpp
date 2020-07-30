//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//  given an empirical obtained polynomial (u,f)

// System includes

// External includes

// Project includes
#include "custom_elements/empirical_spring.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include "cable_net_application_variables.h"
#include "includes/checks.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"


namespace Kratos {
EmpiricalSpringElement3D2N::EmpiricalSpringElement3D2N(IndexType NewId,
                                   GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

EmpiricalSpringElement3D2N::EmpiricalSpringElement3D2N(IndexType NewId,
                                   GeometryType::Pointer pGeometry,
                                   PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

Element::Pointer
EmpiricalSpringElement3D2N::Create(IndexType NewId, NodesArrayType const& rThisNodes,
                         PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<EmpiricalSpringElement3D2N>(NewId, rGeom.Create(rThisNodes),
            pProperties);
}

Element::Pointer
EmpiricalSpringElement3D2N::Create(IndexType NewId, GeometryType::Pointer pGeom,
                         PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<EmpiricalSpringElement3D2N>(NewId, pGeom,
            pProperties);
}

EmpiricalSpringElement3D2N::~EmpiricalSpringElement3D2N() {}

void EmpiricalSpringElement3D2N::EquationIdVector(EquationIdVectorType& rResult,
                                        ProcessInfo& rCurrentProcessInfo)
{

    if (rResult.size() != msLocalSize) {
        rResult.resize(msLocalSize);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * 3;
        rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] =
            GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] =
            GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }
}
void EmpiricalSpringElement3D2N::GetDofList(DofsVectorType& rElementalDofList,
                                  ProcessInfo& rCurrentProcessInfo)
{

    if (rElementalDofList.size() != msLocalSize) {
        rElementalDofList.resize(msLocalSize);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * 3;
        rElementalDofList[index] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] =
            GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index + 2] =
            GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
    }
}

BoundedMatrix<double, EmpiricalSpringElement3D2N::msLocalSize,
EmpiricalSpringElement3D2N::msLocalSize>
EmpiricalSpringElement3D2N::CreateElementStiffnessMatrix(
    ProcessInfo& rCurrentProcessInfo)
{
    BoundedMatrix<double, msLocalSize, msLocalSize> stiffness_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);

    const Vector polynomial = GetProperties()[SPRING_DEFORMATION_EMPIRICAL_POLYNOMIAL];
    const double current_spring_stiffness = EvaluatePolynomialFirstDerivative(polynomial);

    stiffness_matrix(0,0) = current_spring_stiffness;
    stiffness_matrix(0,msDimension) = -current_spring_stiffness;
    stiffness_matrix(msDimension,0) = -current_spring_stiffness;
    stiffness_matrix(msDimension,msDimension) = current_spring_stiffness;
    GlobalizeMatrix(stiffness_matrix);
    return stiffness_matrix;
}

double EmpiricalSpringElement3D2N::GetElementElongation() const {
    return StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this)
        -StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
}

double EmpiricalSpringElement3D2N::EvaluatePolynomial(const Vector& rPolynomial) const{
    const double current_disp = GetElementElongation();
    double current_int_force(0.00);
    for (SizeType i=0;i<rPolynomial.size();++i) current_int_force += rPolynomial[i]*std::pow(current_disp,rPolynomial.size()-1-i);
    return current_int_force;
}


double EmpiricalSpringElement3D2N::EvaluatePolynomialFirstDerivative(const Vector& rPolynomial) const{
    const double current_disp = GetElementElongation();
    double current_int_force_derivative(0.00);
    for (SizeType i=0;i<rPolynomial.size()-1;++i)
        current_int_force_derivative += rPolynomial[i]*std::pow(current_disp,rPolynomial.size()-2-i)*(rPolynomial.size()-1-i);
    return current_int_force_derivative;
}


void EmpiricalSpringElement3D2N::GetValuesVector(Vector& rValues, int Step)
{

    KRATOS_TRY
    if (rValues.size() != msLocalSize) {
        rValues.resize(msLocalSize, false);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msDimension;
        const auto& disp =
            GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);

        rValues[index] = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];
    }
    KRATOS_CATCH("")
}

void EmpiricalSpringElement3D2N::GetFirstDerivativesVector(Vector& rValues, int Step)
{

    KRATOS_TRY
    if (rValues.size() != msLocalSize) {
        rValues.resize(msLocalSize, false);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msDimension;
        const auto& vel =
            GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);

        rValues[index] = vel[0];
        rValues[index + 1] = vel[1];
        rValues[index + 2] = vel[2];
    }
    KRATOS_CATCH("")
}

void EmpiricalSpringElement3D2N::GetSecondDerivativesVector(Vector& rValues, int Step)
{

    KRATOS_TRY
    if (rValues.size() != msLocalSize) {
        rValues.resize(msLocalSize, false);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msDimension;
        const auto& acc =
            GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);

        rValues[index] = acc[0];
        rValues[index + 1] = acc[1];
        rValues[index + 2] = acc[2];
    }

    KRATOS_CATCH("")
}

void EmpiricalSpringElement3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    CalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    CalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void EmpiricalSpringElement3D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    rRightHandSideVector = ZeroVector(msLocalSize);
    BoundedVector<double, msLocalSize> internal_forces = ZeroVector(msLocalSize);

    const Vector polynomial = GetProperties()[SPRING_DEFORMATION_EMPIRICAL_POLYNOMIAL];
    const double scalar_local_force = EvaluatePolynomial(polynomial);

    internal_forces[0] = -scalar_local_force;
    internal_forces[msDimension] = scalar_local_force;

    GlobalizeVector(internal_forces);

    noalias(rRightHandSideVector) -= internal_forces;
    KRATOS_CATCH("")
}

void EmpiricalSpringElement3D2N::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY;
    // resizing the matrices + create memory for LHS
    rLeftHandSideMatrix = ZeroMatrix(msLocalSize, msLocalSize);
    // creating LHS
    noalias(rLeftHandSideMatrix) =
        CreateElementStiffnessMatrix(rCurrentProcessInfo);
    KRATOS_CATCH("")
}

int EmpiricalSpringElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if (dimension != msDimension ||number_of_nodes != msNumberOfNodes) {
        KRATOS_ERROR << "The element works only in 3D and with 2 nodes" << std::endl;
    }
    // verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        NodeType& rnode = GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, rnode);

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode);
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode);
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode);
    }

    if (GetProperties().Has(DENSITY) == false ||
            GetProperties()[DENSITY] <= numerical_limit) {
        KRATOS_ERROR << "DENSITY not provided for this element " << Id()
                     << std::endl;
    }

    if (GetProperties().Has(CROSS_AREA) == false ||
            GetProperties()[CROSS_AREA] <= numerical_limit) {
        KRATOS_ERROR << "CROSS_AREA not provided for this element " << Id()
                     << std::endl;
    }


    if (GetProperties().Has(SPRING_DEFORMATION_EMPIRICAL_POLYNOMIAL) == false) {
        KRATOS_ERROR << "SPRING_DEFORMATION_EMPIRICAL_POLYNOMIAL not provided for this element " << Id()
                     << std::endl;
    }
    else
    {
        KRATOS_ERROR_IF(GetProperties()[SPRING_DEFORMATION_EMPIRICAL_POLYNOMIAL].size()<2)
        << "polynomial must be of order higher than 0 for element " << Id() << std::endl;
    }




    return 0;

    KRATOS_CATCH("")
}


void EmpiricalSpringElement3D2N::LocalizeVector(
    BoundedVector<double, EmpiricalSpringElement3D2N::msLocalSize>& rInputVector)
{
    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);
    CreateTransformationMatrix(transformation_matrix);
    noalias(rInputVector) = prod(transformation_matrix, rInputVector);
}

void EmpiricalSpringElement3D2N::GlobalizeVector(
    BoundedVector<double, EmpiricalSpringElement3D2N::msLocalSize>& rInputVector)
{
    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);
    CreateTransformationMatrix(transformation_matrix);
    noalias(rInputVector) = prod(Matrix(trans(transformation_matrix)), rInputVector);
}

void EmpiricalSpringElement3D2N::GlobalizeMatrix(
    BoundedMatrix<double, EmpiricalSpringElement3D2N::msLocalSize,
    EmpiricalSpringElement3D2N::msLocalSize>& rInputMatrix)
{
    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);
    CreateTransformationMatrix(transformation_matrix);
    rInputMatrix = prod(rInputMatrix,Matrix(trans(transformation_matrix)));
    rInputMatrix = prod(transformation_matrix,rInputMatrix);
}


void EmpiricalSpringElement3D2N::CreateTransformationMatrix(
    BoundedMatrix<double, EmpiricalSpringElement3D2N::msLocalSize,
    EmpiricalSpringElement3D2N::msLocalSize>& rRotationMatrix)
{

    KRATOS_TRY
    const double numeric_limit = std::numeric_limits<double>::epsilon();
    // 1st calculate transformation matrix
    typedef BoundedVector<double, msDimension> arraydim;
    typedef BoundedVector<double, msLocalSize> arraylocal;
    arraydim direction_vector_x = ZeroVector(msDimension);
    arraydim direction_vector_y = ZeroVector(msDimension);
    arraydim direction_vector_z = ZeroVector(msDimension);
    arraylocal reference_coordinates = ZeroVector(msLocalSize);
    arraydim global_z_vector = ZeroVector(msDimension);
    global_z_vector[2] = 1.0;

    WriteTransformationCoordinates(reference_coordinates);

    for (unsigned int i = 0; i < msDimension; ++i) {
        direction_vector_x[i] =
            (reference_coordinates[i + msDimension] - reference_coordinates[i]);
    }
    // local x-axis (e1_local) is the beam axis  (in GID is e3_local)
    double VectorNorm;
    VectorNorm = MathUtils<double>::Norm(direction_vector_x);
    if (VectorNorm > numeric_limit) {
        direction_vector_x /= VectorNorm;
    }

    else {
        KRATOS_ERROR << "length of element" << Id() << "~ zero" << std::endl;
    }

    if (std::abs(direction_vector_x[2]-1.00) <= numeric_limit) {
        direction_vector_y[1] = 1.0;
        direction_vector_z[0] = -1.0;
    }

    else if (std::abs(direction_vector_x[2]+1.00) <= numeric_limit) {
        direction_vector_y[1] = 1.0;
        direction_vector_z[0] = 1.0;
    }

    else {
        MathUtils<double>::UnitCrossProduct(direction_vector_y, direction_vector_x,
                                            global_z_vector);
        MathUtils<double>::UnitCrossProduct(direction_vector_z, direction_vector_y,
                                            direction_vector_x);
    }

    // 2nd fill big rotation matrix
    BoundedMatrix<double, msDimension, msDimension> current_coordinate_system =
        ZeroMatrix(msDimension, msDimension);
    for (unsigned int i = 0; i < msDimension; ++i) {
        current_coordinate_system(i, 0) = direction_vector_x[i];
        current_coordinate_system(i, 1) = direction_vector_y[i];
        current_coordinate_system(i, 2) = direction_vector_z[i];
    }

    rRotationMatrix = ZeroMatrix(msLocalSize, msLocalSize);
    // Building the rotation matrix for the local element matrix
    for (unsigned int kk = 0; kk < msLocalSize; kk += msDimension) {
        for (unsigned int i = 0; i < msDimension; ++i) {
            for (unsigned int j = 0; j < msDimension; ++j) {
                rRotationMatrix(i + kk, j + kk) = current_coordinate_system(i, j);
            }
        }
    }
    KRATOS_CATCH("")
}

void EmpiricalSpringElement3D2N::WriteTransformationCoordinates(
    BoundedVector<double, EmpiricalSpringElement3D2N::msLocalSize>
    & rReferenceCoordinates)
{
    KRATOS_TRY;
    rReferenceCoordinates = ZeroVector(msLocalSize);
    Vector current_displacement = ZeroVector(msLocalSize);
    GetValuesVector(current_displacement, 0);
    rReferenceCoordinates[0] =
        GetGeometry()[0].X0() + current_displacement[0];
    rReferenceCoordinates[1] =
        GetGeometry()[0].Y0() + current_displacement[1];
    rReferenceCoordinates[2] =
        GetGeometry()[0].Z0() + current_displacement[2];
    rReferenceCoordinates[3] =
        GetGeometry()[1].X0() + current_displacement[3];
    rReferenceCoordinates[4] =
        GetGeometry()[1].Y0() + current_displacement[4];
    rReferenceCoordinates[5] =
        GetGeometry()[1].Z0() + current_displacement[5];
    KRATOS_CATCH("");
}


void EmpiricalSpringElement3D2N::CalculateLumpedMassVector(VectorType& rMassVector)
{
    KRATOS_TRY

    // Clear matrix
    if (rMassVector.size() != msLocalSize) {
        rMassVector.resize(msLocalSize, false);
    }

    const double A = GetProperties()[CROSS_AREA];
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double rho = GetProperties()[DENSITY];

    const double total_mass = A * L * rho;

    for (int i = 0; i < msNumberOfNodes; ++i) {
        for (int j = 0; j < msDimension; ++j) {
            int index = i * msDimension + j;

            rMassVector[index] = total_mass * 0.50;
        }
    }

    KRATOS_CATCH("")
}


void EmpiricalSpringElement3D2N::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    // Compute lumped mass matrix
    VectorType temp_vector(msLocalSize);
    CalculateLumpedMassVector(temp_vector);

    // Clear matrix
    if (rMassMatrix.size1() != msLocalSize || rMassMatrix.size2() != msLocalSize) {
        rMassMatrix.resize(msLocalSize, msLocalSize, false);
    }
    rMassMatrix = ZeroMatrix(msLocalSize, msLocalSize);

    // Fill the matrix
    for (IndexType i = 0; i < msLocalSize; ++i) {
        rMassMatrix(i, i) = temp_vector[i];
    }

    KRATOS_CATCH("")
}

void EmpiricalSpringElement3D2N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        msLocalSize);
}


void EmpiricalSpringElement3D2N::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    Variable<double >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    auto& r_geom = GetGeometry();

    if (rDestinationVariable == NODAL_MASS) {
        VectorType element_mass_vector(msLocalSize);
        CalculateLumpedMassVector(element_mass_vector);

        for (SizeType i = 0; i < msNumberOfNodes; ++i) {
            double& r_nodal_mass = r_geom[i].GetValue(NODAL_MASS);
            int index = i * msDimension;

            #pragma omp atomic
            r_nodal_mass += element_mass_vector(index);
        }
    }

    KRATOS_CATCH("")
}

void EmpiricalSpringElement3D2N::AddExplicitContribution(
    const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
    Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {

        BoundedVector<double, msLocalSize> damping_residual_contribution = ZeroVector(msLocalSize);
        Vector current_nodal_velocities = ZeroVector(msLocalSize);
        GetFirstDerivativesVector(current_nodal_velocities);
        Matrix damping_matrix;
        ProcessInfo temp_process_information; // cant pass const ProcessInfo
        CalculateDampingMatrix(damping_matrix, temp_process_information);
        // current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);

        for (size_t i = 0; i < msNumberOfNodes; ++i) {
            size_t index = msDimension * i;
            array_1d<double, 3>& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            for (size_t j = 0; j < msDimension; ++j) {
                #pragma omp atomic
                r_force_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
            }
        }
    } else if (rDestinationVariable == NODAL_INERTIA) {

        // Getting the vector mass
        VectorType mass_vector(msLocalSize);
        CalculateLumpedMassVector(mass_vector);

        for (int i = 0; i < msNumberOfNodes; ++i) {
            double& r_nodal_mass = GetGeometry()[i].GetValue(NODAL_MASS);
            array_1d<double, msDimension>& r_nodal_inertia = GetGeometry()[i].GetValue(NODAL_INERTIA);
            int index = i * msDimension;

            #pragma omp atomic
            r_nodal_mass += mass_vector[index];

            for (int k = 0; k < msDimension; ++k) {
                #pragma omp atomic
                r_nodal_inertia[k] += 0.0;
            }
        }
    }

    KRATOS_CATCH("")
}

void EmpiricalSpringElement3D2N::save(Serializer &rSerializer) const {
  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

void EmpiricalSpringElement3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

} // namespace Kratos.

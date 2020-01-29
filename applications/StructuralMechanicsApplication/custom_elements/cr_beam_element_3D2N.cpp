// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//
// System includes

// External includes

// Project includes
#include "custom_elements/cr_beam_element_3D2N.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{
CrBeamElement3D2N::CrBeamElement3D2N(IndexType NewId,
                                     GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

CrBeamElement3D2N::CrBeamElement3D2N(IndexType NewId,
                                     GeometryType::Pointer pGeometry,
                                     PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

Element::Pointer
CrBeamElement3D2N::Create(IndexType NewId, NodesArrayType const& rThisNodes,
                          PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<CrBeamElement3D2N>(NewId, rGeom.Create(rThisNodes),
            pProperties);
}

Element::Pointer
CrBeamElement3D2N::Create(IndexType NewId, GeometryType::Pointer pGeom,
                          PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<CrBeamElement3D2N>(NewId, pGeom,
            pProperties);
}

CrBeamElement3D2N::~CrBeamElement3D2N() {}

void CrBeamElement3D2N::EquationIdVector(EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != msElementSize) {
        rResult.resize(msElementSize);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msNumberOfNodes * msDimension;
        rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] =
            GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] =
            GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();

        rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
        rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
        rResult[index + 5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
    }
}

void CrBeamElement3D2N::GetDofList(DofsVectorType& rElementalDofList,
                                   ProcessInfo& rCurrentProcessInfo)
{

    if (rElementalDofList.size() != msElementSize) {
        rElementalDofList.resize(msElementSize);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msNumberOfNodes * msDimension;
        rElementalDofList[index] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] =
            GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index + 2] =
            GetGeometry()[i].pGetDof(DISPLACEMENT_Z);

        rElementalDofList[index + 3] = GetGeometry()[i].pGetDof(ROTATION_X);
        rElementalDofList[index + 4] = GetGeometry()[i].pGetDof(ROTATION_Y);
        rElementalDofList[index + 5] = GetGeometry()[i].pGetDof(ROTATION_Z);
    }
}

void CrBeamElement3D2N::Initialize()
{

    KRATOS_TRY;
    KRATOS_CATCH("")
}

void CrBeamElement3D2N::GetSecondDerivativesVector(Vector& rValues, int Step)
{

    KRATOS_TRY
    if (rValues.size() != msElementSize) {
        rValues.resize(msElementSize, false);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msDimension * 2;
        const auto& acc =
            GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const auto& ang_acc = GetGeometry()[i].FastGetSolutionStepValue(
                                  ANGULAR_ACCELERATION, Step);

        rValues[index] = acc[0];
        rValues[index + 1] = acc[1];
        rValues[index + 2] = acc[2];

        rValues[index + 3] = ang_acc[0];
        rValues[index + 4] = ang_acc[1];
        rValues[index + 5] = ang_acc[2];
    }
    KRATOS_CATCH("")
}

void CrBeamElement3D2N::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    mDeformationPreviousIteration = mDeformationCurrentIteration;
    GetValuesVector(mDeformationCurrentIteration, 0);
    KRATOS_CATCH("")
}

void CrBeamElement3D2N::GetFirstDerivativesVector(Vector& rValues, int Step)
{

    KRATOS_TRY
    if (rValues.size() != msElementSize) {
        rValues.resize(msElementSize, false);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msDimension * 2;
        const auto& vel =
            GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        const auto& ang_vel =
            GetGeometry()[i].FastGetSolutionStepValue(ANGULAR_VELOCITY, Step);

        rValues[index] = vel[0];
        rValues[index + 1] = vel[1];
        rValues[index + 2] = vel[2];

        rValues[index + 3] = ang_vel[0];
        rValues[index + 4] = ang_vel[1];
        rValues[index + 5] = ang_vel[2];
    }
    KRATOS_CATCH("")
}

void CrBeamElement3D2N::GetValuesVector(Vector& rValues, int Step)
{
    KRATOS_TRY
    if (rValues.size() != msElementSize) {
        rValues.resize(msElementSize, false);
    }

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msDimension * 2;
        const auto& disp =
            GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const auto& rot =
            GetGeometry()[i].FastGetSolutionStepValue(ROTATION, Step);

        rValues[index] = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];

        rValues[index + 3] = rot[0];
        rValues[index + 4] = rot[1];
        rValues[index + 5] = rot[2];
    }
    KRATOS_CATCH("")
}

BoundedVector<double, CrBeamElement3D2N::msElementSize>
CrBeamElement3D2N::CalculateBodyForces() const
{
    KRATOS_TRY
    // getting shapefunctionvalues for linear SF
    const Matrix& Ncontainer =
        GetGeometry().ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

    BoundedVector<double, msDimension> equivalent_line_load =
        ZeroVector(msDimension);
    BoundedVector<double, msElementSize> body_forces_global =
        ZeroVector(msElementSize);

    const double A = GetProperties()[CROSS_AREA];
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    const double rho = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);

    // calculating equivalent line load
    for (int i = 0; i < msNumberOfNodes; ++i) {
        noalias(equivalent_line_load) +=
            (A * rho * Ncontainer(0, i)) *
            GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    // adding the nodal forces
    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msLocalSize;
        for (int j = 0; j < msDimension; ++j) {
            body_forces_global[j + index] =
                equivalent_line_load[j] * Ncontainer(0, i) * l;
        }
    }

    // adding the nodal moments
    CalculateAndAddWorkEquivalentNodalForcesLineLoad(equivalent_line_load,
            body_forces_global, l);

    // return the total ForceVector
    return body_forces_global;
    KRATOS_CATCH("")
}

void CrBeamElement3D2N::CalculateAndAddWorkEquivalentNodalForcesLineLoad(
    const BoundedVector<double, CrBeamElement3D2N::msDimension> ForceInput,
    BoundedVector<double, CrBeamElement3D2N::msElementSize>
    & rRightHandSideVector,
    const double GeometryLength) const
{
    KRATOS_TRY;
    // calculate orthogonal load vector
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    Vector geometric_orientation = ZeroVector(msDimension);
    geometric_orientation[0] =
        GetGeometry()[1].X() - GetGeometry()[0].X();
    geometric_orientation[1] =
        GetGeometry()[1].Y() - GetGeometry()[0].Y();
    if (msDimension == 3) {
        geometric_orientation[2] =
            GetGeometry()[1].Z() - GetGeometry()[0].Z();
    }

    const double vector_norm_a = MathUtils<double>::Norm(geometric_orientation);
    if (vector_norm_a > numerical_limit) {
        geometric_orientation /= vector_norm_a;
    }

    Vector line_load_direction = ZeroVector(msDimension);
    for (int i = 0; i < msDimension; ++i) {
        line_load_direction[i] = ForceInput[i];
    }

    const double vector_norm_b = MathUtils<double>::Norm(line_load_direction);
    if (vector_norm_b > numerical_limit) {
        line_load_direction /= vector_norm_b;
    }

    double cos_angle = 0.00;
    for (int i = 0; i < msDimension; ++i) {
        cos_angle += line_load_direction[i] * geometric_orientation[i];
    }

    const double sin_angle = std::sqrt(1.00 - (cos_angle * cos_angle));
    const double norm_force_vector_orthogonal = sin_angle * vector_norm_b;

    Vector node_a = ZeroVector(msDimension);
    node_a[0] = GetGeometry()[0].X();
    node_a[1] = GetGeometry()[0].Y();
    if (msDimension == 3) {
        node_a[2] = GetGeometry()[0].Z();
    }

    Vector node_b = ZeroVector(msDimension);
    node_b = node_a + line_load_direction;

    Vector node_c = ZeroVector(msDimension);
    node_c = node_a + (geometric_orientation * cos_angle);

    Vector load_orthogonal_direction = ZeroVector(msDimension);
    load_orthogonal_direction = node_b - node_c;
    const double vector_norm_c =
        MathUtils<double>::Norm(load_orthogonal_direction);
    if (vector_norm_c > numerical_limit) {
        load_orthogonal_direction /= vector_norm_c;
    }

    // now caluclate respective work equivilent nodal moments

    const double custom_moment =
        norm_force_vector_orthogonal * GeometryLength * GeometryLength / 12.00;

    Vector moment_a = ZeroVector(msDimension);
    moment_a = MathUtils<double>::CrossProduct(geometric_orientation,
               load_orthogonal_direction);
    moment_a *= custom_moment;

    for (int i = 0; i < msDimension; ++i) {
        rRightHandSideVector[(1 * msDimension) + i] += moment_a[i];
        rRightHandSideVector[(3 * msDimension) + i] -= moment_a[i];
    }

    KRATOS_CATCH("")
}

void CrBeamElement3D2N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        msElementSize);
}

Vector CrBeamElement3D2N::CalculateLocalNodalForces() const
{
    // deformation modes
    Vector element_forces_t = CalculateElementForces();
    // Nodal element forces local
    Matrix transformation_matrix_s = CalculateTransformationS();
    Vector nodal_forces_local_qe =
        prod(transformation_matrix_s, element_forces_t);
    return nodal_forces_local_qe;
}

BoundedMatrix<double, CrBeamElement3D2N::msElementSize,
CrBeamElement3D2N::msElementSize>
CrBeamElement3D2N::CreateElementStiffnessMatrix_Material() const
{

    KRATOS_TRY;
    const double E = GetProperties()[YOUNG_MODULUS];
    const double G = CalculateShearModulus();
    const double A = GetProperties()[CROSS_AREA];
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);

    const double It = GetProperties()[TORSIONAL_INERTIA];
    const double Iy = GetProperties()[I22];
    const double Iz = GetProperties()[I33];

    double Ay = 0.00;
    if (GetProperties().Has(AREA_EFFECTIVE_Y)) {
        Ay = GetProperties()[AREA_EFFECTIVE_Y];
    }

    double Az = 0.00;
    if (GetProperties().Has(AREA_EFFECTIVE_Z)) {
        Az = GetProperties()[AREA_EFFECTIVE_Z];
    }
    const double Psi_y = CalculatePsi(Iy, Az);
    const double Psi_z = CalculatePsi(Iz, Ay);

    BoundedMatrix<double, msElementSize, msElementSize> local_stiffness_matrix =
        ZeroMatrix(msElementSize, msElementSize);
    const double L3 = L * L * L;
    const double L2 = L * L;

    local_stiffness_matrix(0, 0) = E * A / L;
    local_stiffness_matrix(6, 0) = -1.0 * local_stiffness_matrix(0, 0);
    local_stiffness_matrix(0, 6) = local_stiffness_matrix(6, 0);
    local_stiffness_matrix(6, 6) = local_stiffness_matrix(0, 0);

    local_stiffness_matrix(1, 1) = 12.0 * E * Iz * Psi_z / L3;
    local_stiffness_matrix(1, 7) = -1.0 * local_stiffness_matrix(1, 1);
    local_stiffness_matrix(1, 5) = 6.0 * E * Iz * Psi_z / L2;
    local_stiffness_matrix(1, 11) = local_stiffness_matrix(1, 5);

    local_stiffness_matrix(2, 2) = 12.0 * E * Iy * Psi_y / L3;
    local_stiffness_matrix(2, 8) = -1.0 * local_stiffness_matrix(2, 2);
    local_stiffness_matrix(2, 4) = -6.0 * E * Iy * Psi_y / L2;
    local_stiffness_matrix(2, 10) = local_stiffness_matrix(2, 4);

    local_stiffness_matrix(4, 2) = local_stiffness_matrix(2, 4);
    local_stiffness_matrix(5, 1) = local_stiffness_matrix(1, 5);
    local_stiffness_matrix(3, 3) = G * It / L;
    local_stiffness_matrix(4, 4) = E * Iy * (3.0 * Psi_y + 1.0) / L;
    local_stiffness_matrix(5, 5) = E * Iz * (3.0 * Psi_z + 1.0) / L;
    local_stiffness_matrix(4, 8) = -1.0 * local_stiffness_matrix(4, 2);
    local_stiffness_matrix(5, 7) = -1.0 * local_stiffness_matrix(5, 1);
    local_stiffness_matrix(3, 9) = -1.0 * local_stiffness_matrix(3, 3);
    local_stiffness_matrix(4, 10) = E * Iy * (3.0 * Psi_y - 1.0) / L;
    local_stiffness_matrix(5, 11) = E * Iz * (3.0 * Psi_z - 1.0) / L;

    local_stiffness_matrix(7, 1) = local_stiffness_matrix(1, 7);
    local_stiffness_matrix(7, 5) = local_stiffness_matrix(5, 7);
    local_stiffness_matrix(7, 7) = local_stiffness_matrix(1, 1);
    local_stiffness_matrix(7, 11) = local_stiffness_matrix(7, 5);

    local_stiffness_matrix(8, 2) = local_stiffness_matrix(2, 8);
    local_stiffness_matrix(8, 4) = local_stiffness_matrix(4, 8);
    local_stiffness_matrix(8, 8) = local_stiffness_matrix(2, 2);
    local_stiffness_matrix(8, 10) = local_stiffness_matrix(8, 4);

    local_stiffness_matrix(9, 3) = local_stiffness_matrix(3, 9);
    local_stiffness_matrix(9, 9) = local_stiffness_matrix(3, 3);

    local_stiffness_matrix(10, 2) = local_stiffness_matrix(2, 10);
    local_stiffness_matrix(10, 4) = local_stiffness_matrix(4, 10);
    local_stiffness_matrix(10, 8) = local_stiffness_matrix(8, 10);
    local_stiffness_matrix(10, 10) = local_stiffness_matrix(4, 4);

    local_stiffness_matrix(11, 1) = local_stiffness_matrix(1, 11);
    local_stiffness_matrix(11, 5) = local_stiffness_matrix(5, 11);
    local_stiffness_matrix(11, 7) = local_stiffness_matrix(7, 11);
    local_stiffness_matrix(11, 11) = local_stiffness_matrix(5, 5);

    return local_stiffness_matrix;
    KRATOS_CATCH("")
}

BoundedMatrix<double, CrBeamElement3D2N::msElementSize,
CrBeamElement3D2N::msElementSize>
CrBeamElement3D2N::CreateElementStiffnessMatrix_Geometry() const
{

    KRATOS_TRY;

    Vector nodal_forces_local_qe = CalculateLocalNodalForces();

    const double N = nodal_forces_local_qe[6];
    const double Mt   = nodal_forces_local_qe[9];
    const double my_A = nodal_forces_local_qe[4];
    const double mz_A = nodal_forces_local_qe[5];
    const double my_B = nodal_forces_local_qe[10];
    const double mz_B = nodal_forces_local_qe[11];

    const double L = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    const double Qy = -1.00 * (mz_A + mz_B) / L;
    const double Qz = (my_A + my_B) / L;


    BoundedMatrix<double, msElementSize, msElementSize> local_stiffness_matrix =
        ZeroMatrix(msElementSize, msElementSize);

    local_stiffness_matrix(0, 1) = -Qy / L;
    local_stiffness_matrix(0, 2) = -Qz / L;
    local_stiffness_matrix(0, 7) = -1.0 * local_stiffness_matrix(0, 1);
    local_stiffness_matrix(0, 8) = -1.0 * local_stiffness_matrix(0, 2);

    local_stiffness_matrix(1, 0) = local_stiffness_matrix(0, 1);

    local_stiffness_matrix(1, 1) = 1.2 * N / L;

    local_stiffness_matrix(1, 3) = my_A / L;
    local_stiffness_matrix(1, 4) = Mt / L;

    local_stiffness_matrix(1, 5) = N / 10.0;

    local_stiffness_matrix(1, 6) = local_stiffness_matrix(0, 7);
    local_stiffness_matrix(1, 7) = -1.00 * local_stiffness_matrix(1, 1);
    local_stiffness_matrix(1, 9) = my_B / L;
    local_stiffness_matrix(1, 10) = -1.00 * local_stiffness_matrix(1, 4);
    local_stiffness_matrix(1, 11) = local_stiffness_matrix(1, 5);

    local_stiffness_matrix(2, 0) = local_stiffness_matrix(0, 2);
    local_stiffness_matrix(2, 2) = local_stiffness_matrix(1, 1);
    local_stiffness_matrix(2, 3) = mz_A / L;
    local_stiffness_matrix(2, 4) = -1.00 * local_stiffness_matrix(1, 5);
    local_stiffness_matrix(2, 5) = local_stiffness_matrix(1, 4);
    local_stiffness_matrix(2, 6) = local_stiffness_matrix(0, 8);
    local_stiffness_matrix(2, 8) = local_stiffness_matrix(1, 7);
    local_stiffness_matrix(2, 9) = mz_B / L;
    local_stiffness_matrix(2, 10) = local_stiffness_matrix(2, 4);
    local_stiffness_matrix(2, 11) = local_stiffness_matrix(1, 10);

    for (int i = 0; i < 3; ++i) {
        local_stiffness_matrix(3, i) = local_stiffness_matrix(i, 3);
    }
    local_stiffness_matrix(3, 4) = (-mz_A / 3.00) + (mz_B / 6.00);
    local_stiffness_matrix(3, 5) = (my_A / 3.00) - (my_B / 6.00);
    local_stiffness_matrix(3, 7) = -my_A / L;
    local_stiffness_matrix(3, 8) = -mz_A / L;
    local_stiffness_matrix(3, 10) = L * Qy / 6.00;
    local_stiffness_matrix(3, 11) = L * Qz / 6.00;

    for (int i = 0; i < 4; ++i) {
        local_stiffness_matrix(4, i) = local_stiffness_matrix(i, 4);
    }
    local_stiffness_matrix(4, 4) = 2.00 * L * N / 15.00;
    local_stiffness_matrix(4, 7) = -Mt / L;
    local_stiffness_matrix(4, 8) = N / 10.00;
    local_stiffness_matrix(4, 9) = local_stiffness_matrix(3, 10);
    local_stiffness_matrix(4, 10) = -L * N / 30.00;
    local_stiffness_matrix(4, 11) = Mt / 2.00;

    for (int i = 0; i < 5; ++i) {
        local_stiffness_matrix(5, i) = local_stiffness_matrix(i, 5);
    }
    local_stiffness_matrix(5, 5) = local_stiffness_matrix(4, 4);
    local_stiffness_matrix(5, 7) = -N / 10.0;
    local_stiffness_matrix(5, 8) = -Mt / L;
    local_stiffness_matrix(5, 9) = local_stiffness_matrix(3, 11);
    local_stiffness_matrix(5, 10) = -1.00 * local_stiffness_matrix(4, 11);
    local_stiffness_matrix(5, 11) = local_stiffness_matrix(4, 10);

    for (int i = 0; i < 6; ++i) {
        local_stiffness_matrix(6, i) = local_stiffness_matrix(i, 6);
    }
    local_stiffness_matrix(6, 7) = local_stiffness_matrix(0, 1);
    local_stiffness_matrix(6, 8) = local_stiffness_matrix(0, 2);

    for (int i = 0; i < 7; ++i) {
        local_stiffness_matrix(7, i) = local_stiffness_matrix(i, 7);
    }
    local_stiffness_matrix(7, 7) = local_stiffness_matrix(1, 1);
    local_stiffness_matrix(7, 9) = -1.00 * local_stiffness_matrix(1, 9);
    local_stiffness_matrix(7, 10) = local_stiffness_matrix(4, 1);
    local_stiffness_matrix(7, 11) = local_stiffness_matrix(2, 4);

    for (int i = 0; i < 8; ++i) {
        local_stiffness_matrix(8, i) = local_stiffness_matrix(i, 8);
    }
    local_stiffness_matrix(8, 8) = local_stiffness_matrix(1, 1);
    local_stiffness_matrix(8, 9) = -1.00 * local_stiffness_matrix(2, 9);
    local_stiffness_matrix(8, 10) = local_stiffness_matrix(1, 5);
    local_stiffness_matrix(8, 11) = local_stiffness_matrix(1, 4);

    for (int i = 0; i < 9; ++i) {
        local_stiffness_matrix(9, i) = local_stiffness_matrix(i, 9);
    }
    local_stiffness_matrix(9, 10) = (mz_A / 6.00) - (mz_B / 3.00);
    local_stiffness_matrix(9, 11) = (-my_A / 6.00) + (my_B / 3.00);

    for (int i = 0; i < 10; ++i) {
        local_stiffness_matrix(10, i) = local_stiffness_matrix(i, 10);
    }
    local_stiffness_matrix(10, 10) = local_stiffness_matrix(4, 4);

    for (int i = 0; i < 11; ++i) {
        local_stiffness_matrix(11, i) = local_stiffness_matrix(i, 11);
    }
    local_stiffness_matrix(11, 11) = local_stiffness_matrix(4, 4);

    return local_stiffness_matrix;
    KRATOS_CATCH("")
}

BoundedMatrix<double, CrBeamElement3D2N::msLocalSize,
CrBeamElement3D2N::msLocalSize>
CrBeamElement3D2N::CalculateDeformationStiffness() const
{

    KRATOS_TRY
    BoundedMatrix<double, msLocalSize, msLocalSize> Kd =
        ZeroMatrix(msLocalSize, msLocalSize);
    const double E = GetProperties()[YOUNG_MODULUS];
    const double G = CalculateShearModulus();
    const double A = GetProperties()[CROSS_AREA];
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);

    const double It = GetProperties()[TORSIONAL_INERTIA];
    const double Iy = GetProperties()[I22];
    const double Iz = GetProperties()[I33];

    double Ay = 0.00;
    if (GetProperties().Has(AREA_EFFECTIVE_Y)) {
        Ay = GetProperties()[AREA_EFFECTIVE_Y];
    }

    double Az = 0.00;
    if (GetProperties().Has(AREA_EFFECTIVE_Z)) {
        Az = GetProperties()[AREA_EFFECTIVE_Z];
    }
    const double Psi_y = CalculatePsi(Iy, Az);
    const double Psi_z = CalculatePsi(Iz, Ay);

    Kd(0, 0) = G * It / L;
    Kd(1, 1) = E * Iy / L;
    Kd(2, 2) = E * Iz / L;
    Kd(3, 3) = E * A / L;
    Kd(4, 4) = 3.0 * E * Iy * Psi_y / L;
    Kd(5, 5) = 3.0 * E * Iz * Psi_z / L;


    const double N = (E*A/L) * (l-L);
    const double N1 = l * N / 12.00;
    const double N2 = l * N / 20.00;

    Kd(1, 1) += N1;
    Kd(2, 2) += N1;
    Kd(4, 4) += N2;
    Kd(5, 5) += N2;


    return Kd;
    KRATOS_CATCH("")
}

BoundedMatrix<double, CrBeamElement3D2N::msElementSize,
CrBeamElement3D2N::msElementSize>
CrBeamElement3D2N::CalculateInitialLocalCS() const
{

    KRATOS_TRY
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    array_1d<double, msDimension> direction_vector_x = ZeroVector(msDimension);
    array_1d<double, msDimension> direction_vector_y = ZeroVector(msDimension);
    array_1d<double, msDimension> direction_vector_z = ZeroVector(msDimension);
    array_1d<double, msLocalSize> reference_coordinates = ZeroVector(msLocalSize);

    reference_coordinates[0] = GetGeometry()[0].X0();
    reference_coordinates[1] = GetGeometry()[0].Y0();
    reference_coordinates[2] = GetGeometry()[0].Z0();
    reference_coordinates[3] = GetGeometry()[1].X0();
    reference_coordinates[4] = GetGeometry()[1].Y0();
    reference_coordinates[5] = GetGeometry()[1].Z0();

    for (unsigned int i = 0; i < msDimension; ++i) {
        direction_vector_x[i] =
            (reference_coordinates[i + msDimension] - reference_coordinates[i]);
    }
    Matrix temp_matrix = ZeroMatrix(msDimension);

    // take user defined local axis 2 from GID input
    if (Has(LOCAL_AXIS_2)) {
        double vector_norm = MathUtils<double>::Norm(direction_vector_x);
        if (vector_norm > numerical_limit) {
            direction_vector_x /= vector_norm;
        }

        direction_vector_y = GetValue(LOCAL_AXIS_2);

        direction_vector_z[0] = direction_vector_x[1] * direction_vector_y[2] -
                                direction_vector_x[2] * direction_vector_y[1];
        direction_vector_z[1] = direction_vector_x[2] * direction_vector_y[0] -
                                direction_vector_x[0] * direction_vector_y[2];
        direction_vector_z[2] = direction_vector_x[0] * direction_vector_y[1] -
                                direction_vector_x[1] * direction_vector_y[0];

        vector_norm = MathUtils<double>::Norm(direction_vector_z);
        if (vector_norm > numerical_limit) {
            direction_vector_z /= vector_norm;
        } else
            KRATOS_ERROR << "LOCAL_AXIS_3 has length 0 for element " << Id()
                         << std::endl;

        for (int i = 0; i < msDimension; ++i) {
            temp_matrix(i, 0) = direction_vector_x[i];
            temp_matrix(i, 1) = direction_vector_y[i];
            temp_matrix(i, 2) = direction_vector_z[i];
        }
    }

    // if no user defined local axis 2 input available use this
    else {
        // use orientation class 1st constructor
        double theta_custom = 0.00;
        if (GetProperties().Has(ANG_ROT)) {
            theta_custom = GetProperties()[ANG_ROT];
        }

        typedef array_1d<double, msDimension> arraydim;
        arraydim global_z = ZeroVector(msDimension);
        global_z[2] = 1.0;

        arraydim v2 = ZeroVector(msDimension);
        arraydim v3 = ZeroVector(msDimension);

        double vector_norm;
        vector_norm = MathUtils<double>::Norm(direction_vector_x);
        if (vector_norm > numerical_limit) {
            direction_vector_x /= vector_norm;
        }

        if (std::abs(direction_vector_x[2] - 1.00) < numerical_limit) {
            v2[1] = 1.0;
            v3[0] = -1.0;
        }

        else if (std::abs(direction_vector_x[2] + 1.00) < numerical_limit) {
            v2[1] = 1.0;
            v3[0] = 1.0;
        }

        else {
            MathUtils<double>::UnitCrossProduct(v2, global_z, direction_vector_x);
            MathUtils<double>::UnitCrossProduct(v3, direction_vector_x, v2);
        }

        // manual rotation around the beam axis
        if (std::abs(theta_custom) > numerical_limit) {
            const Vector nz_temp = v3;
            const Vector ny_temp = v2;
            const double cos_theta = std::cos(theta_custom);
            const double sin_theta = std::sin(theta_custom);

            v2 = ny_temp * cos_theta + nz_temp * sin_theta;
            vector_norm = MathUtils<double>::Norm(v2);
            if (vector_norm > numerical_limit) {
                v2 /= vector_norm;
            }

            v3 = nz_temp * cos_theta - ny_temp * sin_theta;
            vector_norm = MathUtils<double>::Norm(v3);
            if (vector_norm > numerical_limit) {
                v3 /= vector_norm;
            }
        }

        for (int i = 0; i < msDimension; ++i) {
            temp_matrix(i, 0) = direction_vector_x[i];
            temp_matrix(i, 1) = v2[i];
            temp_matrix(i, 2) = v3[i];
        }
    }

    BoundedMatrix<double, msElementSize, msElementSize> reference_transformation;

    // Create big rotation Matrix
    AssembleSmallInBigMatrix(temp_matrix, reference_transformation);

    return reference_transformation;
    KRATOS_CATCH("")
}

Vector CrBeamElement3D2N::GetIncrementDeformation() const
{
    KRATOS_TRY;
    return mDeformationCurrentIteration - mDeformationPreviousIteration;
    KRATOS_CATCH("")
}


void CrBeamElement3D2N::UpdateQuaternionParameters(
    double& rScalNodeA,double& rScalNodeB,
    Vector& rVecNodeA,Vector& rVecNodeB) const
{
    KRATOS_TRY;
    BoundedVector<double, msDimension> d_phi_a = ZeroVector(msDimension);
    BoundedVector<double, msDimension> d_phi_b = ZeroVector(msDimension);
    Vector increment_deformation = GetIncrementDeformation();



    for (unsigned int i = 0; i < msDimension; ++i) {
        d_phi_a[i] = increment_deformation[i + 3];
        d_phi_b[i] = increment_deformation[i + 9];
    }

    // calculating quaternions
    Vector drA_vec = ZeroVector(msDimension);
    Vector drB_vec = ZeroVector(msDimension);
    double drA_sca, drB_sca;

    drA_vec = 0.50 * d_phi_a;
    drB_vec = 0.50 * d_phi_b;

    drA_sca = 0.00;
    drB_sca = 0.00;
    for (unsigned int i = 0; i < msDimension; ++i) {
        drA_sca += drA_vec[i] * drA_vec[i];
        drB_sca += drB_vec[i] * drB_vec[i];
    }
    drA_sca = 1.00 - drA_sca;
    drB_sca = 1.00 - drB_sca;

    drA_sca = std::sqrt(drA_sca);
    drB_sca = std::sqrt(drB_sca);

    Vector temp_vector = ZeroVector(msDimension);
    double temp_scalar = 0.00;

    // Node A
    temp_vector = mQuaternionVEC_A;
    temp_scalar = mQuaternionSCA_A;

    rScalNodeA = drA_sca * temp_scalar;
    for (unsigned int i = 0; i < msDimension; ++i) {
        rScalNodeA -= drA_vec[i] * temp_vector[i];
    }

    rVecNodeA = drA_sca * temp_vector;
    rVecNodeA += temp_scalar * drA_vec;
    rVecNodeA += MathUtils<double>::CrossProduct(drA_vec, temp_vector);

    // Node B
    temp_vector = mQuaternionVEC_B;
    temp_scalar = mQuaternionSCA_B;

   rScalNodeB = drB_sca * temp_scalar;
    for (unsigned int i = 0; i < msDimension; ++i) {
        rScalNodeB -= drB_vec[i] * temp_vector[i];
    }

    rVecNodeB = drB_sca * temp_vector;
    rVecNodeB += temp_scalar * drB_vec;
    rVecNodeB += MathUtils<double>::CrossProduct(drB_vec, temp_vector);
    KRATOS_CATCH("");
}

void CrBeamElement3D2N::SaveQuaternionParameters()
{
    KRATOS_TRY;
    double quaternion_sca_a = 0.0;
    double quaternion_sca_b = 0.0;
    Vector quaternion_vec_a = ZeroVector(msDimension);
    Vector quaternion_vec_b = ZeroVector(msDimension);

    UpdateQuaternionParameters(quaternion_sca_a,
    quaternion_sca_b,quaternion_vec_a,quaternion_vec_b);

    mQuaternionVEC_A = quaternion_vec_a;
    mQuaternionVEC_B = quaternion_vec_b;
    mQuaternionSCA_A = quaternion_sca_a;
    mQuaternionSCA_B = quaternion_sca_b;
    KRATOS_CATCH("");
}

BoundedMatrix<double, CrBeamElement3D2N::msDimension,
CrBeamElement3D2N::msDimension>
CrBeamElement3D2N::UpdateRotationMatrixLocal(Vector& Bisectrix,
        Vector& VectorDifference) const
{

    KRATOS_TRY
    const double numerical_limit = std::numeric_limits<double>::epsilon();

    double quaternion_sca_a = 0.0;
    double quaternion_sca_b = 0.0;
    Vector quaternion_vec_a = ZeroVector(msDimension);
    Vector quaternion_vec_b = ZeroVector(msDimension);

    UpdateQuaternionParameters(quaternion_sca_a,
    quaternion_sca_b,quaternion_vec_a,quaternion_vec_b);


    Vector temp_vector = ZeroVector(msDimension);
    // scalar part of difference quaternion
    double scalar_diff;
    scalar_diff = (quaternion_sca_a + quaternion_sca_b) *
                  (quaternion_sca_a + quaternion_sca_b);

    temp_vector = quaternion_vec_a + quaternion_vec_b;
    scalar_diff += MathUtils<double>::Norm(temp_vector) *
                   MathUtils<double>::Norm(temp_vector);

    scalar_diff = 0.50 * std::sqrt(scalar_diff);

    // mean rotation quaternion
    double mean_rotation_scalar;
    mean_rotation_scalar =
        (quaternion_sca_a + quaternion_sca_b) * 0.50;
    mean_rotation_scalar = mean_rotation_scalar / scalar_diff;

    BoundedVector<double, msDimension> mean_rotation_vector =
        ZeroVector(msDimension);
    mean_rotation_vector =
        (quaternion_vec_a + quaternion_vec_b) * 0.50;
    mean_rotation_vector = mean_rotation_vector / scalar_diff;

    // vector part of difference quaternion
    VectorDifference = ZeroVector(msDimension);
    VectorDifference = quaternion_sca_a * quaternion_vec_b;
    VectorDifference -= quaternion_sca_b * quaternion_vec_a;
    VectorDifference += MathUtils<double>::CrossProduct(quaternion_vec_a,
                        quaternion_vec_b);

    VectorDifference = 0.50 * VectorDifference / scalar_diff;

    // rotate inital element basis
    const double r0 = mean_rotation_scalar;
    const double r1 = mean_rotation_vector[0];
    const double r2 = mean_rotation_vector[1];
    const double r3 = mean_rotation_vector[2];

    BoundedMatrix<double, msElementSize, msElementSize>
    reference_transformation = CalculateInitialLocalCS();
    Vector rotated_nx0 = ZeroVector(msDimension);
    Vector rotated_ny0 = ZeroVector(msDimension);
    Vector rotated_nz0 = ZeroVector(msDimension);
    for (IndexType i = 0; i < msDimension; ++i) {
        rotated_nx0[i] = reference_transformation(i, 0);
        rotated_ny0[i] = reference_transformation(i, 1);
        rotated_nz0[i] = reference_transformation(i, 2);
    }

    Quaternion<double> q(r0, r1, r2, r3);
    q.RotateVector3(rotated_nx0);
    q.RotateVector3(rotated_ny0);
    q.RotateVector3(rotated_nz0);

    BoundedMatrix<double, msDimension, msDimension> rotated_coordinate_system =
        ZeroMatrix(msDimension, msDimension);
    for (unsigned int i = 0; i < msDimension; ++i) {
        rotated_coordinate_system(i, 0) = rotated_nx0[i];
        rotated_coordinate_system(i, 1) = rotated_ny0[i];
        rotated_coordinate_system(i, 2) = rotated_nz0[i];
    }

    // rotate basis to element axis + redefine R
    Bisectrix = ZeroVector(msDimension);
    Vector delta_x = ZeroVector(msDimension);
    double vector_norm;

    BoundedVector<double, msLocalSize> current_nodal_position =
        GetCurrentNodalPosition();
    for (unsigned int i = 0; i < msDimension; ++i)
        delta_x[i] =
            current_nodal_position[msDimension + i] - current_nodal_position[i];

    vector_norm = MathUtils<double>::Norm(delta_x);
    if (vector_norm > numerical_limit) {
        delta_x /= vector_norm;
    }

    Bisectrix = rotated_nx0 + delta_x;
    vector_norm = MathUtils<double>::Norm(Bisectrix);
    if (vector_norm > numerical_limit) {
        Bisectrix /= vector_norm;
    }

    BoundedMatrix<double, msDimension, msDimension> n_xyz =
        ZeroMatrix(msDimension);
    for (unsigned int i = 0; i < msDimension; ++i) {
        n_xyz(i, 0) = -rotated_coordinate_system(i, 0);
        n_xyz(i, 1) = rotated_coordinate_system(i, 1);
        n_xyz(i, 2) = rotated_coordinate_system(i, 2);
    }

    BoundedMatrix<double, msDimension, msDimension> Identity =
        ZeroMatrix(msDimension);
    for (unsigned int i = 0; i < msDimension; ++i) {
        Identity(i, i) = 1.0;
    }
    Identity -= 2.0 * outer_prod(Bisectrix, Bisectrix);
    n_xyz = prod(Identity, n_xyz);

    return n_xyz;
    KRATOS_CATCH("")
}

Vector CrBeamElement3D2N::CalculateSymmetricDeformationMode() const
{

    Vector phi_s = ZeroVector(msDimension);
    Vector vector_difference = ZeroVector(msDimension);
    Vector bisectrix = ZeroVector(msDimension);
    Matrix rotation_matrix = UpdateRotationMatrixLocal(bisectrix, vector_difference);
    phi_s = prod(Matrix(trans(rotation_matrix)), vector_difference);
    phi_s *= 4.00;
    return phi_s;
}

Vector CrBeamElement3D2N::CalculateAntiSymmetricDeformationMode() const
{
    Vector phi_a = ZeroVector(msDimension);
    Vector vector_difference = ZeroVector(msDimension);
    Vector bisectrix = ZeroVector(msDimension);
    Matrix rotation_matrix = UpdateRotationMatrixLocal(bisectrix, vector_difference);


    Vector rotated_nx0 = ZeroVector(msDimension);
    for (unsigned int i = 0; i < msDimension; ++i) {
        rotated_nx0[i] = rotation_matrix(i, 0);
    }
    Vector temp_vector = ZeroVector(msDimension);
    MathUtils<double>::CrossProduct(temp_vector, rotated_nx0, bisectrix);
    phi_a = prod(Matrix(trans(rotation_matrix)), temp_vector);
    phi_a *= 4.00;

    return phi_a;
}

BoundedMatrix<double, CrBeamElement3D2N::msElementSize,
CrBeamElement3D2N::msLocalSize>
CrBeamElement3D2N::CalculateTransformationS() const
{

    KRATOS_TRY
    const double L = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    BoundedMatrix<double, msElementSize, msLocalSize> S =
        ZeroMatrix(msElementSize, msLocalSize);
    S(0, 3) = -1.00;
    S(1, 5) = 2.00 / L;
    S(2, 4) = -2.00 / L;
    S(3, 0) = -1.00;
    S(4, 1) = -1.00;
    S(4, 4) = 1.00;
    S(5, 2) = -1.00;
    S(5, 5) = 1.00;
    S(6, 3) = 1.00;
    S(7, 5) = -2.00 / L;
    S(8, 4) = 2.00 / L;
    S(9, 0) = 1.00;
    S(10, 1) = 1.00;
    S(10, 4) = 1.00;
    S(11, 2) = 1.00;
    S(11, 5) = 1.00;

    return S;
    KRATOS_CATCH("")
}

void CrBeamElement3D2N::CalculateMassMatrix(MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    if (rMassMatrix.size1() != msElementSize) {
        rMassMatrix.resize(msElementSize, msElementSize, false);
    }
    rMassMatrix = ZeroMatrix(msElementSize, msElementSize);

    bool use_consistent_mass_matrix = false;

    if (GetProperties().Has(USE_CONSISTENT_MASS_MATRIX)) {
        use_consistent_mass_matrix = GetProperties()[USE_CONSISTENT_MASS_MATRIX];
    }

    if (!use_consistent_mass_matrix) {
        CalculateLumpedMassMatrix(rMassMatrix, rCurrentProcessInfo);
    } else {
        CalculateConsistentMassMatrix(rMassMatrix, rCurrentProcessInfo);
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix = GetTransformationMatrixGlobal();

        BoundedMatrix<double, msElementSize, msElementSize> aux_matrix =
            prod(rotation_matrix, rMassMatrix);
        rMassMatrix = prod(aux_matrix, Matrix(trans(rotation_matrix)));
    }
    KRATOS_CATCH("")
}

void CrBeamElement3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY;
    ConstCalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void CrBeamElement3D2N::ConstCalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    ConstCalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    ConstCalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    KRATOS_CATCH("");
}

void CrBeamElement3D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    ConstCalculateRightHandSide(rRightHandSideVector,rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void CrBeamElement3D2N::ConstCalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    Vector internal_forces = CalculateGlobalNodalForces();
    rRightHandSideVector = ZeroVector(msElementSize);
    noalias(rRightHandSideVector) -= internal_forces;
    // add bodyforces
    noalias(rRightHandSideVector) += CalculateBodyForces();
    KRATOS_CATCH("")
}



void CrBeamElement3D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY;
    ConstCalculateLeftHandSide(rLeftHandSideMatrix,rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void CrBeamElement3D2N::ConstCalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) const
{

    KRATOS_TRY;
    // resizing the matrices + create memory for LHS
    rLeftHandSideMatrix = ZeroMatrix(msElementSize, msElementSize);
    // creating LHS
    noalias(rLeftHandSideMatrix) += CreateElementStiffnessMatrix_Material();
    noalias(rLeftHandSideMatrix) += CreateElementStiffnessMatrix_Geometry();

    BoundedMatrix<double, msElementSize, msElementSize> total_rotation_matrix = GetTransformationMatrixGlobal();


    BoundedMatrix<double, msElementSize, msElementSize> aux_matrix =
        prod(total_rotation_matrix, rLeftHandSideMatrix);
    noalias(rLeftHandSideMatrix) = prod(aux_matrix, trans(total_rotation_matrix));

    KRATOS_CATCH("")
}

Vector CrBeamElement3D2N::CalculateGlobalNodalForces() const
{
    Vector nodal_forces_local_qe = CalculateLocalNodalForces();

    BoundedMatrix<double, msElementSize, msElementSize> total_rotation_matrix = GetTransformationMatrixGlobal();

    // Nodal element forces global
    BoundedVector<double, msElementSize> nodal_forces_global_q =
        ZeroVector(msElementSize);
    nodal_forces_global_q = prod(total_rotation_matrix, nodal_forces_local_qe);
    return nodal_forces_global_q;
}

BoundedVector<double, CrBeamElement3D2N::msLocalSize>
CrBeamElement3D2N::CalculateElementForces() const
{
    KRATOS_TRY;
    BoundedVector<double, msLocalSize> deformation_modes_total_v =
        ZeroVector(msLocalSize);
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);

    Vector phi_s = CalculateSymmetricDeformationMode();
    Vector phi_a = CalculateAntiSymmetricDeformationMode();

    deformation_modes_total_v[3] = l - L;
    for (int i = 0; i < 3; ++i) {
        deformation_modes_total_v[i] = phi_s[i];
    }
    for (int i = 0; i < 2; ++i) {
        deformation_modes_total_v[i + 4] = phi_a[i + 1];
    }

    // calculate element forces
    BoundedVector<double, msLocalSize> element_forces_t =
        ZeroVector(msLocalSize);
    BoundedMatrix<double, msLocalSize, msLocalSize> deformation_stiffness_Kd =
        ZeroMatrix(msLocalSize);

    deformation_stiffness_Kd = CalculateDeformationStiffness();
    element_forces_t = prod(deformation_stiffness_Kd, deformation_modes_total_v);

    return element_forces_t;
    KRATOS_CATCH("")
}

double CrBeamElement3D2N::CalculatePsi(const double I, const double A_eff) const
{

    KRATOS_TRY;
    const double E = GetProperties()[YOUNG_MODULUS];
    const double L = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    const double G = CalculateShearModulus();

    const double phi = (12.0 * E * I) / (L * L * G * A_eff);
    double psi;
    // interpret input A_eff == 0 as shearstiff -> psi = 1.0
    if (A_eff == 0.00) {
        psi = 1.00;
    } else {
        psi = 1.0 / (1.0 + phi);
    }

    return psi;
    KRATOS_CATCH("")
}


BoundedVector<double, CrBeamElement3D2N::msLocalSize>
CrBeamElement3D2N::GetCurrentNodalPosition() const
{
    BoundedVector<double, msLocalSize> current_nodal_position =
        ZeroVector(msLocalSize);
    for (unsigned int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msDimension;
        current_nodal_position[index] =
            GetGeometry()[i].X0() +
            GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X, 0);
        current_nodal_position[index + 1] =
            GetGeometry()[i].Y0() +
            GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y, 0);
        current_nodal_position[index + 2] =
            GetGeometry()[i].Z0() +
            GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Z, 0);
    }

    return current_nodal_position;
}

void CrBeamElement3D2N::Calculate(const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == LOCAL_ELEMENT_ORIENTATION) {
        if(rOutput.size1() != msDimension || rOutput.size2() != msDimension) {
            rOutput.resize(msDimension, msDimension, false);
        }

        Matrix  transformation_matrix = GetTransformationMatrixGlobal();

        Vector base_1 = ZeroVector(3);
        Vector base_2 = ZeroVector(3);
        Vector base_3 = ZeroVector(3);

        for (SizeType i=0;i<msDimension;++i){
            base_1[i] = transformation_matrix(i,0);
            base_2[i] = transformation_matrix(i,1);
            base_3[i] = transformation_matrix(i,2);
        }

        column(rOutput,0) = base_1;
        column(rOutput,1) = base_2;
        column(rOutput,2) = base_3;
    }
}

void CrBeamElement3D2N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    // element with two nodes can only represent results at one node
    const unsigned int& write_points_number =
        GetGeometry().IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);
    if (rOutput.size() != write_points_number) {
        rOutput.resize(write_points_number);
    }


    // rOutput[GP 1,2,3][x,y,z]

    if (rVariable == MOMENT) {
        Vector nodal_forces_local_qe = CalculateLocalNodalForces();
        rOutput[0][0] = -1.0 * nodal_forces_local_qe[3] * 0.75 + nodal_forces_local_qe[9] * 0.25;
        rOutput[1][0] = -1.0 * nodal_forces_local_qe[3] * 0.50 + nodal_forces_local_qe[9] * 0.50;
        rOutput[2][0] = -1.0 * nodal_forces_local_qe[3] * 0.25 + nodal_forces_local_qe[9] * 0.75;

        rOutput[0][1] = -1.0 * nodal_forces_local_qe[4] * 0.75 + nodal_forces_local_qe[10] * 0.25;
        rOutput[1][1] = -1.0 * nodal_forces_local_qe[4] * 0.50 + nodal_forces_local_qe[10] * 0.50;
        rOutput[2][1] = -1.0 * nodal_forces_local_qe[4] * 0.25 + nodal_forces_local_qe[10] * 0.75;

        rOutput[0][2] = 1.0 * nodal_forces_local_qe[5] * 0.75 - nodal_forces_local_qe[11] * 0.25;
        rOutput[1][2] = 1.0 * nodal_forces_local_qe[5] * 0.50 - nodal_forces_local_qe[11] * 0.50;
        rOutput[2][2] = 1.0 * nodal_forces_local_qe[5] * 0.25 - nodal_forces_local_qe[11] * 0.75;
    } else if (rVariable == FORCE) {
        Vector nodal_forces_local_qe = CalculateLocalNodalForces();
        rOutput[0][0] = -1.0 * nodal_forces_local_qe[0] * 0.75 + nodal_forces_local_qe[6] * 0.25;
        rOutput[1][0] = -1.0 * nodal_forces_local_qe[0] * 0.50 + nodal_forces_local_qe[6] * 0.50;
        rOutput[2][0] = -1.0 * nodal_forces_local_qe[0] * 0.25 + nodal_forces_local_qe[6] * 0.75;

        rOutput[0][1] = -1.0 * nodal_forces_local_qe[1] * 0.75 + nodal_forces_local_qe[7] * 0.25;
        rOutput[1][1] = -1.0 * nodal_forces_local_qe[1] * 0.50 + nodal_forces_local_qe[7] * 0.50;
        rOutput[2][1] = -1.0 * nodal_forces_local_qe[1] * 0.25 + nodal_forces_local_qe[7] * 0.75;

        rOutput[0][2] = -1.0 * nodal_forces_local_qe[2] * 0.75 + nodal_forces_local_qe[8] * 0.25;
        rOutput[1][2] = -1.0 * nodal_forces_local_qe[2] * 0.50 + nodal_forces_local_qe[8] * 0.50;
        rOutput[2][2] = -1.0 * nodal_forces_local_qe[2] * 0.25 + nodal_forces_local_qe[8] * 0.75;
    }

    else if (rVariable == LOCAL_AXIS_1) {
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix = GetTransformationMatrixGlobal();
        for (SizeType i =0; i<msDimension; ++i) {
            rOutput[1][i] = column(rotation_matrix, 0)[i];
        }
    } else if (rVariable == LOCAL_AXIS_2) {
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix = GetTransformationMatrixGlobal();
        for (SizeType i =0; i<msDimension; ++i) {
            rOutput[1][i] = column(rotation_matrix, 1)[i];
        }
    } else if (rVariable == LOCAL_AXIS_3) {
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix = GetTransformationMatrixGlobal();
        for (SizeType i =0; i<msDimension; ++i) {
            rOutput[1][i] = column(rotation_matrix, 2)[i];
        }
    }


    KRATOS_CATCH("")
}

void CrBeamElement3D2N::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    KRATOS_CATCH("")
}


void CrBeamElement3D2N::AssembleSmallInBigMatrix(
    Matrix SmallMatrix,
    BoundedMatrix<double, CrBeamElement3D2N::msElementSize,
    CrBeamElement3D2N::msElementSize>& BigMatrix) const
{
    KRATOS_TRY
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    BigMatrix.clear();
    for (unsigned int kk = 0; kk < msElementSize; kk += msDimension) {
        for (int i = 0; i < msDimension; ++i) {
            for (int j = 0; j < msDimension; ++j) {
                if (std::abs(SmallMatrix(i, j)) <= numerical_limit) {
                    BigMatrix(i + kk, j + kk) = 0.00;
                } else {
                    BigMatrix(i + kk, j + kk) = SmallMatrix(i, j);
                }
            }
        }
    }
    KRATOS_CATCH("")
}

BoundedMatrix<double, CrBeamElement3D2N::msElementSize, CrBeamElement3D2N::msElementSize>
  CrBeamElement3D2N::GetTransformationMatrixGlobal() const
{
    KRATOS_TRY;
    Vector vector_difference = ZeroVector(msDimension);
    Vector bisectrix = ZeroVector(msDimension);
    Matrix rotation_matrix = UpdateRotationMatrixLocal(bisectrix, vector_difference);
    BoundedMatrix<double, msElementSize, msElementSize> total_rotation_matrix;
    AssembleSmallInBigMatrix(rotation_matrix, total_rotation_matrix);
    return total_rotation_matrix;
    KRATOS_CATCH("");
}


void CrBeamElement3D2N::BuildSingleMassMatrix(MatrixType& rMassMatrix,
        double Phi, double CT, double CR,
        double L, double dir) const
{
    KRATOS_TRY;
    const SizeType MatSize = msNumberOfNodes * 2;

    if (rMassMatrix.size1() != MatSize) {
        rMassMatrix.resize(MatSize, MatSize, false);
    }
    rMassMatrix = ZeroMatrix(MatSize, MatSize);
    BoundedMatrix<double, MatSize, MatSize> temp_mass_matrix =
        ZeroMatrix(MatSize, MatSize);
    const double Phi2 = Phi * Phi;
    const double L2 = L * L;

    temp_mass_matrix(0, 0) =
        (13.00 / 35.00) + (7.00 / 10.00) * Phi + (1.00 / 3.00) * Phi2;
    temp_mass_matrix(0, 1) =
        dir *
        ((11.00 / 210.00) + (11.00 / 210.00) * Phi + (1.00 / 24.00) * Phi2) * L;
    temp_mass_matrix(0, 2) =
        (9.00 / 70.00) + (3.00 / 10.00) * Phi + (1.00 / 6.00) * Phi2;
    temp_mass_matrix(0, 3) =
        -((13.00 / 420.00) + (3.00 / 40.00) * Phi + (1.00 / 24.00) * Phi2) * L *
        dir;
    temp_mass_matrix(1, 0) = temp_mass_matrix(0, 1);
    temp_mass_matrix(1, 1) =
        ((1.00 / 105.00) + (1.00 / 60.00) * Phi + (1.00 / 120.00) * Phi2) * L2;
    temp_mass_matrix(1, 2) =
        dir * ((13.00 / 420.00) + (3.00 / 40.00) * Phi + (1.00 / 24.00) * Phi2) *
        L;
    temp_mass_matrix(1, 3) =
        -((1.00 / 140.00) + (1.00 / 60.00) * Phi + (1.00 / 120.00) * Phi2) * L2;
    temp_mass_matrix(2, 0) = temp_mass_matrix(0, 2);
    temp_mass_matrix(2, 1) = temp_mass_matrix(1, 2);
    temp_mass_matrix(2, 2) =
        (13.00 / 35.00) + (7.00 / 10.00) * Phi + (1.00 / 3.00) * Phi2;
    temp_mass_matrix(2, 3) =
        -((11.00 / 210.00) + (11.00 / 210.00) * Phi + (1.00 / 24.00) * Phi2) * L *
        dir;
    temp_mass_matrix(3, 0) = temp_mass_matrix(0, 3);
    temp_mass_matrix(3, 1) = temp_mass_matrix(1, 3);
    temp_mass_matrix(3, 2) = temp_mass_matrix(2, 3);
    temp_mass_matrix(3, 3) =
        ((1.00 / 105.00) + (1.00 / 60.00) * Phi + (1.00 / 120.00) * Phi2) * L2;

    temp_mass_matrix *= CT;
    rMassMatrix += temp_mass_matrix;

    temp_mass_matrix = ZeroMatrix(MatSize, MatSize);

    temp_mass_matrix(0, 0) = 6.00 / 5.00;
    temp_mass_matrix(0, 1) = dir * ((1.00 / 10.00) - (1.00 / 2.00) * Phi) * L;
    temp_mass_matrix(0, 2) = -6.00 / 5.00;
    temp_mass_matrix(0, 3) = dir * ((1.00 / 10.00) - (1.00 / 2.00) * Phi) * L;
    temp_mass_matrix(1, 0) = temp_mass_matrix(0, 1);
    temp_mass_matrix(1, 1) =
        ((2.00 / 15.00) + (1.00 / 6.00) * Phi + (1.00 / 3.00) * Phi2) * L2;
    temp_mass_matrix(1, 2) = dir * ((-1.00 / 10.00) + (1.00 / 2.00) * Phi) * L;
    temp_mass_matrix(1, 3) =
        -((1.00 / 30.00) + (1.00 / 6.00) * Phi - (1.00 / 6.00) * Phi2) * L2;
    temp_mass_matrix(2, 0) = temp_mass_matrix(0, 2);
    temp_mass_matrix(2, 1) = temp_mass_matrix(1, 2);
    temp_mass_matrix(2, 2) = 6.00 / 5.00;
    temp_mass_matrix(2, 3) = dir * ((-1.00 / 10.00) + (1.00 / 2.00) * Phi) * L;
    temp_mass_matrix(3, 0) = temp_mass_matrix(0, 3);
    temp_mass_matrix(3, 1) = temp_mass_matrix(1, 3);
    temp_mass_matrix(3, 2) = temp_mass_matrix(2, 3);
    temp_mass_matrix(3, 3) =
        ((2.00 / 15.00) + (1.00 / 6.00) * Phi + (1.00 / 3.00) * Phi2) * L2;

    temp_mass_matrix *= CR;
    rMassMatrix += temp_mass_matrix;
    KRATOS_CATCH("")
}

void CrBeamElement3D2N::CalculateConsistentMassMatrix(
    MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    const int smallMatSize = msNumberOfNodes * 2;

    if (rMassMatrix.size1() != msElementSize) {
        rMassMatrix.resize(msElementSize, msElementSize, false);
    }
    rMassMatrix = ZeroMatrix(msElementSize, msElementSize);

    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double L2 = L * L;
    const double rho = GetProperties()[DENSITY];
    const double A = GetProperties()[CROSS_AREA];
    const double E = GetProperties()[YOUNG_MODULUS];

    const double Iy = GetProperties()[I22];
    const double Iz = GetProperties()[I33];

    double Ip = 0.00;
    if (GetProperties().Has(MASS_MOMENT_OF_INERTIA)){
        Ip = GetProperties()[MASS_MOMENT_OF_INERTIA];
    }
    else {
        //This is an approximation for a circular cross section
        Ip = Iy + Iz;
    }
    const double G = CalculateShearModulus();

    double Ay = 0.00;
    if (GetProperties().Has(AREA_EFFECTIVE_Y)) {
        Ay = GetProperties()[AREA_EFFECTIVE_Y];
    }

    double Az = 0.00;
    if (GetProperties().Has(AREA_EFFECTIVE_Z)) {
        Az = GetProperties()[AREA_EFFECTIVE_Z];
    }

    double IRy = Iy;
    if (GetProperties().Has(INERTIA_ROT_Y)) {
        IRy = GetProperties()[INERTIA_ROT_Y];
    }

    double IRz = Iz;
    if (GetProperties().Has(INERTIA_ROT_Z)) {
        IRz = GetProperties()[INERTIA_ROT_Z];
    }

    double Phiy = 0.00;
    double Phiz = 0.00;

    if (Ay != 0.00) {
        Phiz = (12.00 * E * Iz) / (L2 * G * Ay);
    }
    if (Az != 0.00) {
        Phiy = (12.00 * E * Iy) / (L2 * G * Az);
    }

    const double CTy = (rho * A * L) / ((1 + Phiy) * (1 + Phiy));
    const double CTz = (rho * A * L) / ((1 + Phiz) * (1 + Phiz));

    const double CRy = (rho * IRy) / ((1 + Phiy) * (1 + Phiy) * L);
    const double CRz = (rho * IRz) / ((1 + Phiz) * (1 + Phiz) * L);

    // longitudinal forces + torsional moment
    const double M00 = (1.00 / 3.00) * A * rho * L;
    const double M06 = M00 / 2.00;
    const double M33 = (Ip * L * rho) / 3.00;
    const double M39 = M33 / 2.00;

    rMassMatrix(0, 0) = M00;
    rMassMatrix(0, 6) = M06;
    rMassMatrix(6, 6) = M00;
    rMassMatrix(3, 3) = M33;
    rMassMatrix(3, 9) = M39;
    rMassMatrix(9, 9) = M33;

    Matrix temp_bending_mass_matrix = ZeroMatrix(smallMatSize, smallMatSize);
    BuildSingleMassMatrix(temp_bending_mass_matrix, Phiz, CTz, CRz, L, +1);

    rMassMatrix(1, 1) = temp_bending_mass_matrix(0, 0);
    rMassMatrix(1, 5) = temp_bending_mass_matrix(0, 1);
    rMassMatrix(1, 7) = temp_bending_mass_matrix(0, 2);
    rMassMatrix(1, 11) = temp_bending_mass_matrix(0, 3);
    rMassMatrix(5, 5) = temp_bending_mass_matrix(1, 1);
    rMassMatrix(5, 7) = temp_bending_mass_matrix(1, 2);
    rMassMatrix(5, 11) = temp_bending_mass_matrix(1, 3);
    rMassMatrix(7, 7) = temp_bending_mass_matrix(2, 2);
    rMassMatrix(7, 11) = temp_bending_mass_matrix(2, 3);
    rMassMatrix(11, 11) = temp_bending_mass_matrix(3, 3);

    temp_bending_mass_matrix = ZeroMatrix(smallMatSize, smallMatSize);
    BuildSingleMassMatrix(temp_bending_mass_matrix, Phiy, CTy, CRy, L, -1);

    rMassMatrix(2, 2) = temp_bending_mass_matrix(0, 0);
    rMassMatrix(2, 4) = temp_bending_mass_matrix(0, 1);
    rMassMatrix(2, 8) = temp_bending_mass_matrix(0, 2);
    rMassMatrix(2, 10) = temp_bending_mass_matrix(0, 3);
    rMassMatrix(4, 4) = temp_bending_mass_matrix(1, 1);
    rMassMatrix(4, 8) = temp_bending_mass_matrix(1, 2);
    rMassMatrix(4, 10) = temp_bending_mass_matrix(1, 3);
    rMassMatrix(8, 8) = temp_bending_mass_matrix(2, 2);
    rMassMatrix(8, 10) = temp_bending_mass_matrix(2, 3);
    rMassMatrix(10, 10) = temp_bending_mass_matrix(3, 3);

    for (int j = 1; j < 12; ++j) {
        for (int i = 0; i < j; ++i) {
            rMassMatrix(j, i) = rMassMatrix(i, j);
        }
    }

    KRATOS_CATCH("")
}

void CrBeamElement3D2N::CalculateLumpedMassMatrix(
    MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    if (rMassMatrix.size1() != msElementSize) {
        rMassMatrix.resize(msElementSize, msElementSize, false);
    }
    rMassMatrix = ZeroMatrix(msElementSize, msElementSize);
    const double A = GetProperties()[CROSS_AREA];
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double rho = GetProperties()[DENSITY];

    const double total_mass = A * L * rho;
    const double temp = 0.50 * total_mass;

    // w.r.t. Felippa - Chapter 31: LUMPED AND CONSISTENT MASS MATRICES - p.31–10
    const double rotational_inertia_lumped = total_mass * L * L * GetProperties()[LUMPED_MASS_ROTATION_COEFFICIENT];

    // translatonal mass
    for (int i = 0; i < msNumberOfNodes; ++i) {
        for (int j = 0; j < msDimension; ++j) {
            int index = i * (msDimension * 2) + j;
            rMassMatrix(index, index) = temp;

            // add rotational inertia
            rMassMatrix(index+msDimension, index+msDimension) = rotational_inertia_lumped;
        }
    }

    KRATOS_CATCH("")
}

CrBeamElement3D2N::IntegrationMethod
CrBeamElement3D2N::GetIntegrationMethod() const
{
    // do this to have 3GP as an output in GID
    return Kratos::GeometryData::GI_GAUSS_3;
}

void CrBeamElement3D2N::AddExplicitContribution(
    const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
    Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    BoundedVector<double, msElementSize> damping_residual_contribution = ZeroVector(msElementSize);
    // Calculate damping contribution to residual -->
    if ((GetProperties().Has(RAYLEIGH_ALPHA) ||
            GetProperties().Has(RAYLEIGH_BETA)) &&
            (rDestinationVariable != NODAL_INERTIA)) {
        Vector current_nodal_velocities = ZeroVector(msElementSize);
        GetFirstDerivativesVector(current_nodal_velocities);
        Matrix damping_matrix = ZeroMatrix(msElementSize, msElementSize);
        ProcessInfo temp_process_information; // cant pass const ProcessInfo
        CalculateDampingMatrix(damping_matrix, temp_process_information);
        // current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);
    }

    if (rRHSVariable == RESIDUAL_VECTOR &&
            rDestinationVariable == FORCE_RESIDUAL) {

        for (IndexType i = 0; i < msNumberOfNodes; ++i) {
            const SizeType index = msLocalSize * i;

            array_1d<double, 3>& r_force_residual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

            for (IndexType j = 0; j < msDimension; ++j) {
                #pragma omp atomic
                r_force_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
            }
        }
    }

    if (rRHSVariable == RESIDUAL_VECTOR &&
            rDestinationVariable == MOMENT_RESIDUAL) {

        for (IndexType i = 0; i < msNumberOfNodes; ++i) {
            const SizeType index = (msLocalSize * i) + msDimension;

            array_1d<double, 3>& r_moment_residual = GetGeometry()[i].FastGetSolutionStepValue(MOMENT_RESIDUAL);

            for (IndexType j = 0; j < msDimension; ++j) {
                #pragma omp atomic
                r_moment_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
            }
        }
    }

    if (rDestinationVariable == NODAL_INERTIA) {
        Matrix element_mass_matrix = ZeroMatrix(msElementSize, msElementSize);
        ProcessInfo temp_info; // Dummy
        CalculateMassMatrix(element_mass_matrix, temp_info);

        for (IndexType i = 0; i < msNumberOfNodes; ++i) {
            double aux_nodal_mass = 0.0;
            array_1d<double, 3> aux_nodal_inertia(3, 0.0);

            const SizeType index = i * msLocalSize;

            for (IndexType j = 0; j < msElementSize; ++j) {
                aux_nodal_mass += element_mass_matrix(index, j);
                for (IndexType k = 0; k < msDimension; ++k) {
                    aux_nodal_inertia[k] += element_mass_matrix(index + msDimension + k, j);
                }
            }

            #pragma omp atomic
            GetGeometry()[i].GetValue(NODAL_MASS) += aux_nodal_mass;

            array_1d<double, 3>& r_nodal_inertia = GetGeometry()[i].GetValue(NODAL_INERTIA);
            for (IndexType k = 0; k < msDimension; ++k) {
                #pragma omp atomic
                r_nodal_inertia[k] += std::abs(aux_nodal_inertia[k]);
            }
        }
    }

    KRATOS_CATCH("")
}

double CrBeamElement3D2N::CalculateShearModulus() const
{
    KRATOS_TRY;
    const double nu = GetProperties()[POISSON_RATIO];
    const double E = GetProperties()[YOUNG_MODULUS];
    const double G = E / (2.0 * (1.0 + nu));
    return G;
    KRATOS_CATCH("")
}

int CrBeamElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const double numerical_limit = std::numeric_limits<double>::epsilon();

    KRATOS_ERROR_IF(GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size() != 2)
            << "The beam element works only in 3D and with 2 noded elements" << std::endl;

    // verify that the variables are correctly initialized
    if (VELOCITY.Key() == 0) {
        KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is "
                     "correctly registered"
                     << "" << std::endl;
    }
    if (DISPLACEMENT.Key() == 0) {
        KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is "
                     "correctly registered"
                     << "" << std::endl;
    }
    if (ACCELERATION.Key() == 0) {
        KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is "
                     "correctly registered"
                     << "" << std::endl;
    }
    if (DENSITY.Key() == 0) {
        KRATOS_ERROR << "DENSITY has Key zero! (check if the application is "
                     "correctly registered"
                     << "" << std::endl;
    }
    if (CROSS_AREA.Key() == 0) {
        KRATOS_ERROR << "CROSS_AREA has Key zero! (check if the application is "
                     "correctly registered"
                     << "" << std::endl;
    }
    // verify that the dofs exist
    for (unsigned int i = 0; i < GetGeometry().size(); ++i) {
        if (GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false) {
            KRATOS_ERROR << "missing variable DISPLACEMENT on node "
                         << GetGeometry()[i].Id() << std::endl;
        }
        if (GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false ||
                GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false ||
                GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false) {
            KRATOS_ERROR
                    << "missing one of the dofs for the variable DISPLACEMENT on node "
                    << GetGeometry()[i].Id() << std::endl;
        }
    }

    KRATOS_ERROR_IF(!GetProperties().Has(CROSS_AREA) ||
                    GetProperties()[CROSS_AREA] <= numerical_limit)
            << "Please provide a reasonable value for \"CROSS_AREA\" for element #"
            << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(YOUNG_MODULUS) ||
                    GetProperties()[YOUNG_MODULUS] <= numerical_limit)
            << "Please provide a reasonable value for \"YOUNG_MODULUS\" for element #"
            << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(DENSITY) ||
                    GetProperties()[DENSITY] <= numerical_limit)
            << "Please provide a reasonable value for \"DENSITY\" for element #"
            << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(I22) ||
                    GetProperties()[I22] <= numerical_limit)
            << "Please provide a reasonable value for \"I22\" for element #"
            << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(I33) ||
                    GetProperties()[I33] <= numerical_limit)
            << "Please provide a reasonable value for \"I33\" for element #"
            << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(TORSIONAL_INERTIA) ||
                    GetProperties()[TORSIONAL_INERTIA] <= numerical_limit)
            << "Please provide a reasonable value for \"TORSIONAL_INERTIA\" for element #"
            << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(POISSON_RATIO))
            << "\"POISSON_RATIO\" not provided for element #" << Id() << std::endl;

    if (Has(LOCAL_AXIS_2)) {
        array_1d<double, msDimension> direction_vector_x = ZeroVector(msDimension);
        array_1d<double, msDimension> direction_vector_y = ZeroVector(msDimension);
        array_1d<double, msLocalSize> reference_coordinates = ZeroVector(msLocalSize);

        reference_coordinates[0] = GetGeometry()[0].X0();
        reference_coordinates[1] = GetGeometry()[0].Y0();
        reference_coordinates[2] = GetGeometry()[0].Z0();
        reference_coordinates[3] = GetGeometry()[1].X0();
        reference_coordinates[4] = GetGeometry()[1].Y0();
        reference_coordinates[5] = GetGeometry()[1].Z0();

        for (unsigned int i = 0; i < msDimension; ++i) {
            direction_vector_x[i] = (reference_coordinates[i + msDimension] - reference_coordinates[i]);
        }

        const double vector_norm = MathUtils<double>::Norm(direction_vector_x);
        if (vector_norm > numerical_limit) {
            direction_vector_x /= vector_norm;
        }

        direction_vector_y = GetValue(LOCAL_AXIS_2);

        KRATOS_ERROR_IF(MathUtils<double>::Norm(direction_vector_y)<numerical_limit)
                << "Given LOCAL_AXIS_2 has length 0 for element " << Id() << std::endl;

        // a tollerance of 1e-3 allows for a rough deviation of 0.06 degrees from 90.0 degrees
        KRATOS_ERROR_IF(std::abs(MathUtils<double>::Dot(direction_vector_x,direction_vector_y))>1e-3)
                << "LOCAL_AXIS_1 is not perpendicular to LOCAL_AXIS_2 for element " << Id() << std::endl;
    }

    KRATOS_ERROR_IF(StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this)
                    < std::numeric_limits<double>::epsilon())
            << "Element #" << Id() << " has a length of zero!" << std::endl;

    return 0;

    KRATOS_CATCH("")
}


void CrBeamElement3D2N::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    SaveQuaternionParameters();
}

void CrBeamElement3D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    rSerializer.save("NodalDeformationCurrent", mDeformationCurrentIteration);
    rSerializer.save("NodalDeformationPrevious", mDeformationPreviousIteration);
    rSerializer.save("QuaternionVecA", mQuaternionVEC_A);
    rSerializer.save("QuaternionVecB", mQuaternionVEC_B);
    rSerializer.save("QuaternionScaA", mQuaternionSCA_A);
    rSerializer.save("QuaternionScaB", mQuaternionSCA_B);
}

void CrBeamElement3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    rSerializer.load("NodalDeformationCurrent", mDeformationCurrentIteration);
    rSerializer.load("NodalDeformationPrevious", mDeformationPreviousIteration);
    rSerializer.load("QuaternionVecA", mQuaternionVEC_A);
    rSerializer.load("QuaternionVecB", mQuaternionVEC_B);
    rSerializer.load("QuaternionScaA", mQuaternionSCA_A);
    rSerializer.load("QuaternionScaB", mQuaternionSCA_B);
}

} // namespace Kratos.

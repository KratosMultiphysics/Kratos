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
#include "custom_elements/new_beam.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{
NewBeamElement3D2N::NewBeamElement3D2N(IndexType NewId,
                                     GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

NewBeamElement3D2N::NewBeamElement3D2N(IndexType NewId,
                                     GeometryType::Pointer pGeometry,
                                     PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

Element::Pointer
NewBeamElement3D2N::Create(IndexType NewId, NodesArrayType const& rThisNodes,
                          PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_shared<NewBeamElement3D2N>(NewId, rGeom.Create(rThisNodes),
            pProperties);
}

Element::Pointer
NewBeamElement3D2N::Create(IndexType NewId, GeometryType::Pointer pGeom,
                          PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<NewBeamElement3D2N>(NewId, pGeom,
            pProperties);
}

NewBeamElement3D2N::~NewBeamElement3D2N() {}

void NewBeamElement3D2N::EquationIdVector(EquationIdVectorType& rResult,
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

void NewBeamElement3D2N::GetDofList(DofsVectorType& rElementalDofList,
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

void NewBeamElement3D2N::Initialize()
{
    KRATOS_TRY;
    KRATOS_CATCH("")
}

void NewBeamElement3D2N::GetSecondDerivativesVector(Vector& rValues, int Step)
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

void NewBeamElement3D2N::GetFirstDerivativesVector(Vector& rValues, int Step)
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

void NewBeamElement3D2N::GetValuesVector(Vector& rValues, int Step)
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

void NewBeamElement3D2N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void NewBeamElement3D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    Vector internal_forces = CalculateGlobalNodalForces();
    rRightHandSideVector = ZeroVector(msElementSize);
    noalias(rRightHandSideVector) -= internal_forces;
    KRATOS_CATCH("")
}

void NewBeamElement3D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY;
    // resizing the matrices + create memory for LHS
    rLeftHandSideMatrix = ZeroMatrix(msElementSize, msElementSize);
    // creating LHS
    noalias(rLeftHandSideMatrix) += CreateElementStiffnessMatrix_Material();

    Vector vector_difference = ZeroVector(msDimension);
    Vector bisectrix = ZeroVector(msDimension);
    Matrix rotation_matrix = UpdateRotationMatrixLocal(bisectrix, vector_difference, false);

    BoundedMatrix<double, msElementSize, msElementSize> total_rotation_matrix;
    // Create big rotation Matrix
    AssembleSmallInBigMatrix(rotation_matrix, total_rotation_matrix);


    BoundedMatrix<double, msElementSize, msElementSize> aux_matrix =
        prod(total_rotation_matrix, rLeftHandSideMatrix);
    noalias(rLeftHandSideMatrix) = prod(aux_matrix, trans(total_rotation_matrix));

    KRATOS_CATCH("")
}

BoundedMatrix<double, NewBeamElement3D2N::msLocalSize,
NewBeamElement3D2N::msLocalSize>
NewBeamElement3D2N::CalculateDeformationStiffness() const
{

    KRATOS_TRY
    BoundedMatrix<double, msLocalSize, msLocalSize> Kd =
        ZeroMatrix(msLocalSize, msLocalSize);
    const double E = GetProperties()[YOUNG_MODULUS];
    const double G = CalculateShearModulus();
    const double A = GetProperties()[CROSS_AREA];
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    //const double L = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);

    const double J = GetProperties()[TORSIONAL_INERTIA];
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

    Kd(0, 0) = G * J / L;
    Kd(1, 1) = E * Iy / L;
    Kd(2, 2) = E * Iz / L;
    Kd(3, 3) = E * A / L;
    Kd(4, 4) = 3.0 * E * Iy * Psi_y / L;
    Kd(5, 5) = 3.0 * E * Iz * Psi_z / L;


    return Kd;
    KRATOS_CATCH("")
}

double NewBeamElement3D2N::CalculatePsi(const double I, const double A_eff) const
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

double NewBeamElement3D2N::CalculateShearModulus() const
{
    KRATOS_TRY;
    const double nu = GetProperties()[POISSON_RATIO];
    const double E = GetProperties()[YOUNG_MODULUS];
    const double G = E / (2.0 * (1.0 + nu));
    return G;
    KRATOS_CATCH("")
}

BoundedMatrix<double, NewBeamElement3D2N::msElementSize,
NewBeamElement3D2N::msElementSize>
NewBeamElement3D2N::CreateElementStiffnessMatrix_Material() const
{

    KRATOS_TRY;
    const double E = GetProperties()[YOUNG_MODULUS];
    const double G = CalculateShearModulus();
    const double A = GetProperties()[CROSS_AREA];
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);

    const double J = GetProperties()[TORSIONAL_INERTIA];
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
    local_stiffness_matrix(3, 3) = G * J / L;
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

BoundedMatrix<double, NewBeamElement3D2N::msElementSize,
NewBeamElement3D2N::msElementSize>
NewBeamElement3D2N::CalculateInitialLocalCS() const
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


void NewBeamElement3D2N::AssembleSmallInBigMatrix(
    Matrix SmallMatrix,
    BoundedMatrix<double, NewBeamElement3D2N::msElementSize,
    NewBeamElement3D2N::msElementSize>& BigMatrix) const
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

BoundedMatrix<double, NewBeamElement3D2N::msElementSize,
NewBeamElement3D2N::msLocalSize>
NewBeamElement3D2N::CalculateTransformationS() const
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

Vector NewBeamElement3D2N::UpdateIncrementDeformation()
{
    KRATOS_TRY
    Vector actual_deformation = ZeroVector(msElementSize);
    GetValuesVector(actual_deformation, 0);
    Vector increment_deformation =
        actual_deformation - mTotalNodalDeformation;

    //mTotalNodalDeformation = actual_deformation;
    return increment_deformation;
    KRATOS_CATCH("")
}

BoundedVector<double, NewBeamElement3D2N::msLocalSize>
    NewBeamElement3D2N::GetCurrentNodalPosition() const
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

BoundedMatrix<double, NewBeamElement3D2N::msDimension,
NewBeamElement3D2N::msDimension>
NewBeamElement3D2N::UpdateRotationMatrixLocal(Vector& Bisectrix,
        Vector& VectorDifference,const bool& rFinalize)
{

    KRATOS_TRY
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    BoundedVector<double, msDimension> d_phi_a = ZeroVector(msDimension);
    BoundedVector<double, msDimension> d_phi_b = ZeroVector(msDimension);
    Vector increment_deformation = UpdateIncrementDeformation();

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

    double quaternion_sca_a = drA_sca * temp_scalar;
    for (unsigned int i = 0; i < msDimension; ++i) {
        quaternion_sca_a -= drA_vec[i] * temp_vector[i];
    }

    Vector quaternion_vec_a = drA_sca * temp_vector;
    quaternion_vec_a += temp_scalar * drA_vec;
    quaternion_vec_a +=
        MathUtils<double>::CrossProduct(drA_vec, temp_vector);

    // Node B
    temp_vector = mQuaternionVEC_B;
    temp_scalar = mQuaternionSCA_B;

    double quaternion_sca_b = drB_sca * temp_scalar;
    for (unsigned int i = 0; i < msDimension; ++i) {
        quaternion_sca_b -= drB_vec[i] * temp_vector[i];
    }

    Vector quaternion_vec_b = drB_sca * temp_vector;
    quaternion_vec_b += temp_scalar * drB_vec;
    quaternion_vec_b +=
        MathUtils<double>::CrossProduct(drB_vec, temp_vector);

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


    if (rFinalize)
    {
        mQuaternionVEC_A = quaternion_vec_a;
        mQuaternionVEC_B = quaternion_vec_b;
        mQuaternionSCA_A = quaternion_sca_a;
        mQuaternionSCA_B = quaternion_sca_b;
    }

    return n_xyz;
    KRATOS_CATCH("")
}

Vector NewBeamElement3D2N::CalculateSymmetricDeformationMode() {

    Vector phi_s = ZeroVector(msDimension);
    Vector vector_difference = ZeroVector(msDimension);
    Vector bisectrix = ZeroVector(msDimension);
    Matrix rotation_matrix = UpdateRotationMatrixLocal(bisectrix, vector_difference, false);
    phi_s = prod(Matrix(trans(rotation_matrix)), vector_difference);
    phi_s *= 4.00;
    return phi_s;
}

Vector NewBeamElement3D2N::CalculateAntiSymmetricDeformationMode()
{
    Vector phi_a = ZeroVector(msDimension);
    Vector vector_difference = ZeroVector(msDimension);
    Vector bisectrix = ZeroVector(msDimension);
    Matrix rotation_matrix = UpdateRotationMatrixLocal(bisectrix, vector_difference, false);


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


BoundedVector<double, NewBeamElement3D2N::msLocalSize>
NewBeamElement3D2N::CalculateElementForces()
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

Vector NewBeamElement3D2N::CalculateLocalNodalForces()
{
    // deformation modes
    Vector element_forces_t = CalculateElementForces();
    // Nodal element forces local
    Matrix transformation_matrix_s = CalculateTransformationS();
    Vector nodal_forces_local_qe =
        prod(transformation_matrix_s, element_forces_t);
    return nodal_forces_local_qe;
}

Vector NewBeamElement3D2N::CalculateGlobalNodalForces()
{
    Vector nodal_forces_local_qe = CalculateLocalNodalForces();
    Vector vector_difference = ZeroVector(msDimension);
    Vector bisectrix = ZeroVector(msDimension);
    Matrix rotation_matrix = UpdateRotationMatrixLocal(bisectrix, vector_difference, false);

    BoundedMatrix<double, msElementSize, msElementSize> total_rotation_matrix;
    // Create big rotation Matrix
    AssembleSmallInBigMatrix(rotation_matrix, total_rotation_matrix);

    // Nodal element forces global
    BoundedVector<double, msElementSize> nodal_forces_global_q =
        ZeroVector(msElementSize);
    nodal_forces_global_q = prod(total_rotation_matrix, nodal_forces_local_qe);
    return nodal_forces_global_q;
}


void NewBeamElement3D2N::FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    GetValuesVector(mTotalNodalDeformation,1);
    Vector vector_difference = ZeroVector(msDimension);
    Vector bisectrix = ZeroVector(msDimension);
    UpdateRotationMatrixLocal(bisectrix,vector_difference,true);
}

int NewBeamElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo)
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

void NewBeamElement3D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    rSerializer.save("NodalDeformation", mTotalNodalDeformation);
    rSerializer.save("QuaternionVecA", mQuaternionVEC_A);
    rSerializer.save("QuaternionVecB", mQuaternionVEC_B);
    rSerializer.save("QuaternionScaA", mQuaternionSCA_A);
    rSerializer.save("QuaternionScaB", mQuaternionSCA_B);
}

void NewBeamElement3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    rSerializer.load("NodalDeformation", mTotalNodalDeformation);
    rSerializer.load("QuaternionVecA", mQuaternionVEC_A);
    rSerializer.load("QuaternionVecB", mQuaternionVEC_B);
    rSerializer.load("QuaternionScaA", mQuaternionSCA_A);
    rSerializer.load("QuaternionScaB", mQuaternionSCA_B);
}

} // namespace Kratos.

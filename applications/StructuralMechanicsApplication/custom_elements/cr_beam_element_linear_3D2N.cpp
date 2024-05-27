// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Klaus B. Sautter
//

// System includes

// External includes

// Project includes
#include "custom_elements/cr_beam_element_linear_3D2N.hpp"
#include "custom_utilities/static_condensation_utility.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{
CrBeamElementLinear3D2N::CrBeamElementLinear3D2N(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : CrBeamElement3D2N(NewId, pGeometry) {}

CrBeamElementLinear3D2N::CrBeamElementLinear3D2N(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : CrBeamElement3D2N(NewId, pGeometry, pProperties) {}

Element::Pointer
CrBeamElementLinear3D2N::Create(IndexType NewId,
                                NodesArrayType const& rThisNodes,
                                PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<CrBeamElementLinear3D2N>(
               NewId, rGeom.Create(rThisNodes), pProperties);
}

Element::Pointer
CrBeamElementLinear3D2N::Create(IndexType NewId,
                                GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<CrBeamElementLinear3D2N>(
               NewId, pGeom, pProperties);
}

CrBeamElementLinear3D2N::~CrBeamElementLinear3D2N() {}


int CrBeamElementLinear3D2N::Check(const ProcessInfo& rCurrentProcessInfo) const {

    KRATOS_TRY

    const GeometryType& r_geom = GetGeometry();
    // verify that the rotational stiffness around axis 2 and 3 are greater or equal than 0.0 if they are given.
    for (unsigned int i = 0; i < r_geom.size(); ++i) {

        if (r_geom[i].Has(ROTATIONAL_STIFFNESS_AXIS_2)) {
            KRATOS_ERROR_IF(r_geom[i].GetValue(ROTATIONAL_STIFFNESS_AXIS_2) < 0.0) 
                << " ROTATIONAL_STIFFNESS_AXIS_2 should be greater or equal than 0.0" << std::endl;
        }
        if (r_geom[i].Has(ROTATIONAL_STIFFNESS_AXIS_3)) {
            KRATOS_ERROR_IF(r_geom[i].GetValue(ROTATIONAL_STIFFNESS_AXIS_3) < 0.0) 
                << " ROTATIONAL_STIFFNESS_AXIS_3 should be greater or equal than 0.0" << std::endl;
        }
    }

    return CrBeamElement3D2N::Check(rCurrentProcessInfo);

    KRATOS_CATCH("")
}

void CrBeamElementLinear3D2N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    Vector nodal_deformation = ZeroVector(msElementSize);
    GetValuesVector(nodal_deformation);
    rRightHandSideVector = ZeroVector(msElementSize);
    rRightHandSideVector -= prod(rLeftHandSideMatrix, nodal_deformation);

    // add bodyforces
    rRightHandSideVector += CalculateBodyForces();
    KRATOS_CATCH("")
}

void CrBeamElementLinear3D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    rRightHandSideVector = ZeroVector(msElementSize);

    Matrix left_hand_side_matrix = ZeroMatrix(msElementSize, msElementSize);
    CalculateLeftHandSide(left_hand_side_matrix, rCurrentProcessInfo);
    Vector nodal_deformation = ZeroVector(msElementSize);
    GetValuesVector(nodal_deformation);
    rRightHandSideVector = ZeroVector(msElementSize);
    noalias(rRightHandSideVector) -=
        prod(left_hand_side_matrix, nodal_deformation);

    // add bodyforces
    noalias(rRightHandSideVector) += CalculateBodyForces();
    KRATOS_CATCH("")
}

void CrBeamElementLinear3D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY;
    BoundedMatrix<double, msElementSize, msElementSize> transformation_matrix =
        CalculateInitialLocalCS();
    rLeftHandSideMatrix = ZeroMatrix(msElementSize, msElementSize);
    noalias(rLeftHandSideMatrix) += CreateElementStiffnessMatrix_Material();


    const GeometryType& r_geometry = GetGeometry();

    // reduce rigidity if required
    for (IndexType i = 0; i < r_geometry.size(); ++i) {
        if (r_geometry[i].Has(ROTATIONAL_STIFFNESS_AXIS_2) || r_geometry[i].Has(ROTATIONAL_STIFFNESS_AXIS_3)) {
            BoundedMatrix<double, msElementSize, msElementSize> rigidity_reduction_matrix = ZeroMatrix(msElementSize, msElementSize);
            CalculateRigidityReductionMatrix(rigidity_reduction_matrix);
            rLeftHandSideMatrix = element_prod(rLeftHandSideMatrix, rigidity_reduction_matrix);
            break;
        }
    }
    
    //// start static condensation
    if (Has(CONDENSED_DOF_LIST)) {
        Vector dof_list_input = GetValue(CONDENSED_DOF_LIST);
        std::vector<int> dofList(dof_list_input.size());
        for (SizeType i = 0; i < dof_list_input.size(); ++i) {
            dofList[i] = dof_list_input[i];
        }
        StaticCondensationUtility::CondenseLeftHandSide(*this, rLeftHandSideMatrix,
                dofList);
    }
    //// end static condensation

    BoundedMatrix<double, msElementSize, msElementSize> aux_matrix =
        ZeroMatrix(msElementSize);
    aux_matrix = prod(transformation_matrix, rLeftHandSideMatrix);
    noalias(rLeftHandSideMatrix) =
        prod(aux_matrix, Matrix(trans(transformation_matrix)));

    KRATOS_CATCH("")
}

void CrBeamElementLinear3D2N::CalculateMassMatrix(MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    if (rMassMatrix.size1() != msElementSize) {
        rMassMatrix.resize(msElementSize, msElementSize, false);
    }
    rMassMatrix = ZeroMatrix(msElementSize, msElementSize);

    const bool compute_lumped_mass_matrix = StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(GetProperties(), rCurrentProcessInfo);

    if (compute_lumped_mass_matrix) {
        CalculateLumpedMassMatrix(rMassMatrix, rCurrentProcessInfo);
    } else {
        CalculateConsistentMassMatrix(rMassMatrix, rCurrentProcessInfo);
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix =
            CalculateInitialLocalCS();

        BoundedMatrix<double, msElementSize, msElementSize> aux_matrix =
            prod(rotation_matrix, rMassMatrix);
        rMassMatrix = prod(aux_matrix, Matrix(trans(rotation_matrix)));
    }

    KRATOS_CATCH("")
}


BoundedMatrix<double, CrBeamElement3D2N::msLocalSize,
CrBeamElement3D2N::msLocalSize>
CrBeamElementLinear3D2N::CalculateDeformationStiffness() const
{

    KRATOS_TRY
    BoundedMatrix<double, msLocalSize, msLocalSize> Kd =
        ZeroMatrix(msLocalSize, msLocalSize);
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

    Kd(0, 0) = G * J / L;
    Kd(1, 1) = E * Iy / L;
    Kd(2, 2) = E * Iz / L;
    Kd(3, 3) = E * A / L;
    Kd(4, 4) = 3.0 * E * Iy * Psi_y / L;
    Kd(5, 5) = 3.0 * E * Iz * Psi_z / L;

    return Kd;
    KRATOS_CATCH("")
}

void CrBeamElementLinear3D2N::Calculate(const Variable<Matrix>& rVariable,
     Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == LOCAL_ELEMENT_ORIENTATION) {
        if(rOutput.size1() != msElementSize || rOutput.size2() != msElementSize) {
            rOutput.resize(msElementSize, msElementSize, false);
        }
        noalias(rOutput) = CalculateInitialLocalCS();
    } else {
        CrBeamElement3D2N::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }

}

void CrBeamElementLinear3D2N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    // element with two nodes can only represent results at one node
    const unsigned int& write_points_number =
        GetGeometry().IntegrationPointsNumber(Kratos::GeometryData::IntegrationMethod::GI_GAUSS_3);
    if (rOutput.size() != write_points_number) {
        rOutput.resize(write_points_number);
    }

    Matrix left_hand_side_matrix = CreateElementStiffnessMatrix_Material();

    Vector nodal_deformation = ZeroVector(msElementSize);
    GetValuesVector(nodal_deformation);

    BoundedMatrix<double, msElementSize, msElementSize> transformation_matrix =
        CalculateInitialLocalCS();
    nodal_deformation =
        prod(Matrix(trans(transformation_matrix)), nodal_deformation);

    //// start static back condensation
    if (Has(CONDENSED_DOF_LIST)) {
        Vector dof_list_input = GetValue(CONDENSED_DOF_LIST);
        std::vector<int> dofList(dof_list_input.size());
        for (SizeType i = 0; i < dof_list_input.size(); ++i) {
            dofList[i] = dof_list_input[i];
        }
        Vector nodal_deformation_temp = nodal_deformation;
        StaticCondensationUtility::ConvertingCondensation(
            *this, nodal_deformation_temp, nodal_deformation, dofList,
            left_hand_side_matrix);
    }
    //// end static back condensation

    Vector stress = prod(left_hand_side_matrix, nodal_deformation);

    // rOutput[GP 1,2,3][x,y,z]

    if (rVariable == MOMENT) {
        rOutput[0][0] = -1.0 * stress[3] * 0.75 + stress[9] * 0.25;
        rOutput[1][0] = -1.0 * stress[3] * 0.50 + stress[9] * 0.50;
        rOutput[2][0] = -1.0 * stress[3] * 0.25 + stress[9] * 0.75;

        rOutput[0][1] = -1.0 * stress[4] * 0.75 + stress[10] * 0.25;
        rOutput[1][1] = -1.0 * stress[4] * 0.50 + stress[10] * 0.50;
        rOutput[2][1] = -1.0 * stress[4] * 0.25 + stress[10] * 0.75;

        rOutput[0][2] = -1.0 * stress[5] * 0.75 + stress[11] * 0.25;
        rOutput[1][2] = -1.0 * stress[5] * 0.50 + stress[11] * 0.50;
        rOutput[2][2] = -1.0 * stress[5] * 0.25 + stress[11] * 0.75;
    }
    if (rVariable == FORCE) {
        rOutput[0][0] = -1.0 * stress[0] * 0.75 + stress[6] * 0.25;
        rOutput[1][0] = -1.0 * stress[0] * 0.50 + stress[6] * 0.50;
        rOutput[2][0] = -1.0 * stress[0] * 0.25 + stress[6] * 0.75;

        rOutput[0][1] = -1.0 * stress[1] * 0.75 + stress[7] * 0.25;
        rOutput[1][1] = -1.0 * stress[1] * 0.50 + stress[7] * 0.50;
        rOutput[2][1] = -1.0 * stress[1] * 0.25 + stress[7] * 0.75;

        rOutput[0][2] = -1.0 * stress[2] * 0.75 + stress[8] * 0.25;
        rOutput[1][2] = -1.0 * stress[2] * 0.50 + stress[8] * 0.50;
        rOutput[2][2] = -1.0 * stress[2] * 0.25 + stress[8] * 0.75;
    }

    KRATOS_CATCH("")
}

void CrBeamElementLinear3D2N::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable, std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rVariable == LOCAL_AXES_VECTOR) {
        BoundedMatrix<double, msElementSize, msElementSize> transformation_matrix;
        transformation_matrix = CalculateInitialLocalCS();

        rOutput.resize(3);
        for (int i = 0; i < 3; ++i) {
            rOutput[i] = ZeroVector(3);
        }

        for (SizeType i = 0; i < 3; ++i) {
            for (SizeType j = 0; j < 3; ++j) {
                rOutput[i][j] = transformation_matrix(j, i);
            }
        }
    }

    KRATOS_CATCH("");
}


void CrBeamElementLinear3D2N::CalculateRigidityReductionMatrix(
    BoundedMatrix<double, msElementSize, msElementSize>& rRigidityReductionMatrix) const
{
    KRATOS_TRY


    rRigidityReductionMatrix = ZeroMatrix(msElementSize, msElementSize);
    const double E = GetProperties()[YOUNG_MODULUS];
    const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);


    const double Iy = GetProperties()[I22];
    const double Iz = GetProperties()[I33];

    const GeometryType& r_geometry = GetGeometry();
    
    // inititalise fixity factors as completely fixed
    double alpha1_yy = 1;
    double alpha2_yy = 1;

    double alpha1_zz = 1;
    double alpha2_zz = 1;

    // set fixity factors, 0 is no fixity, i.e. a hinge; 1 is completely fixed, the resulting stiffness matrix will remain unchanged
    if (r_geometry[0].Has(ROTATIONAL_STIFFNESS_AXIS_2)) {
        const double rotational_stiffness_y_1 = r_geometry[0].GetValue(ROTATIONAL_STIFFNESS_AXIS_2);
        alpha1_zz = (rotational_stiffness_y_1 > 0) ? 1 / (1 + 3 * E * Iy / (rotational_stiffness_y_1 * L)) : 0.0;
    }
    if (r_geometry[1].Has(ROTATIONAL_STIFFNESS_AXIS_2)) {
        const double rotational_stiffness_y_2 = r_geometry[1].GetValue(ROTATIONAL_STIFFNESS_AXIS_2);
        alpha2_zz = (rotational_stiffness_y_2 > 0) ? 1 / (1 + 3 * E * Iy / (rotational_stiffness_y_2 * L)) : 0.0;
    }
    if (r_geometry[0].Has(ROTATIONAL_STIFFNESS_AXIS_3)) {
        const double rotational_stiffness_z_1 = r_geometry[0].GetValue(ROTATIONAL_STIFFNESS_AXIS_3);
        alpha1_yy = (rotational_stiffness_z_1 > 0) ? 1 / (1 + 3 * E * Iz / (rotational_stiffness_z_1 * L)) : 0.0;

    }
    if (r_geometry[1].Has(ROTATIONAL_STIFFNESS_AXIS_3)) {
        const double rotational_stiffness_z_2 = r_geometry[1].GetValue(ROTATIONAL_STIFFNESS_AXIS_3);
        alpha2_yy = (rotational_stiffness_z_2 > 0) ? 1 / (1 + 3 * E * Iz / (rotational_stiffness_z_2 * L)) : 0.0;
    }

    const double yy_denominator = (4 - alpha1_yy * alpha2_yy);
    const double zz_denominator = (4 - alpha1_zz * alpha2_zz);

    // normal
    rRigidityReductionMatrix(0, 0) = 1;
    rRigidityReductionMatrix(0, 6) = 1;
    rRigidityReductionMatrix(6, 6) = 1;
    rRigidityReductionMatrix(6, 0) = 1;

    // yy
    rRigidityReductionMatrix(1, 1) = (alpha1_yy + alpha2_yy + alpha1_yy * alpha2_yy) / yy_denominator;
    rRigidityReductionMatrix(7, 7) = rRigidityReductionMatrix(1, 1);
    rRigidityReductionMatrix(1, 7) = rRigidityReductionMatrix(1, 1);
    rRigidityReductionMatrix(7, 1) = rRigidityReductionMatrix(1, 1);

    // zz
    rRigidityReductionMatrix(2, 2) = (alpha1_zz + alpha2_zz + alpha1_zz * alpha2_zz) / zz_denominator;
    rRigidityReductionMatrix(8, 8) = rRigidityReductionMatrix(2, 2);
    rRigidityReductionMatrix(2, 8) = rRigidityReductionMatrix(2, 2);
    rRigidityReductionMatrix(8, 2) = rRigidityReductionMatrix(2, 2);

    // tortion
    rRigidityReductionMatrix(3, 3) = 1;
    rRigidityReductionMatrix(9, 9) = 1;

    // yy 1,2 - rot z 1
    rRigidityReductionMatrix(1, 5) = (2 * alpha1_yy + alpha1_yy * alpha2_yy) / yy_denominator;
    rRigidityReductionMatrix(5, 1) = rRigidityReductionMatrix(1, 5);
    rRigidityReductionMatrix(5, 7) = rRigidityReductionMatrix(1, 5);
    rRigidityReductionMatrix(7, 5) = rRigidityReductionMatrix(1, 5);

    // zz 1,2 - rot y 1
    rRigidityReductionMatrix(2, 4) = (2 * alpha1_zz + alpha1_zz * alpha2_zz) / zz_denominator;
    rRigidityReductionMatrix(4, 2) = rRigidityReductionMatrix(2, 4);
    rRigidityReductionMatrix(4, 8) = rRigidityReductionMatrix(2, 4);
    rRigidityReductionMatrix(8, 4) = rRigidityReductionMatrix(2, 4);

    // yy 1,2 - rot z 2
    rRigidityReductionMatrix(1, 11) = (2 * alpha2_yy + alpha1_yy * alpha2_yy) / yy_denominator;
    rRigidityReductionMatrix(11, 1) = rRigidityReductionMatrix(1, 11);
    rRigidityReductionMatrix(11, 7) = rRigidityReductionMatrix(1, 11);
    rRigidityReductionMatrix(7, 11) = rRigidityReductionMatrix(1, 11);

    // zz 1,2 - rot y 2
    rRigidityReductionMatrix(2, 10) = (2 * alpha2_zz + alpha1_zz * alpha2_zz) / zz_denominator;
    rRigidityReductionMatrix(10, 2) = rRigidityReductionMatrix(2, 10);
    rRigidityReductionMatrix(10, 8) = rRigidityReductionMatrix(2, 10);
    rRigidityReductionMatrix(8, 10) = rRigidityReductionMatrix(2, 10);

    // rot z 1,2 
    rRigidityReductionMatrix(5, 5) = 3 * alpha1_yy / yy_denominator;
    rRigidityReductionMatrix(11, 11) = 3 * alpha2_yy / yy_denominator;
    rRigidityReductionMatrix(5, 11) = 3 * alpha1_yy * alpha2_yy / yy_denominator;
    rRigidityReductionMatrix(11, 5) = rRigidityReductionMatrix(5, 11);

    // rot y 1, 2
    rRigidityReductionMatrix(4, 4) = 3 * alpha1_zz / zz_denominator;
    rRigidityReductionMatrix(10, 10) = 3 * alpha2_zz / zz_denominator;
    rRigidityReductionMatrix(4, 10) = 3 * alpha1_zz * alpha2_zz / zz_denominator;
    rRigidityReductionMatrix(10, 4) = rRigidityReductionMatrix(4, 10);

    KRATOS_CATCH("");

}

void CrBeamElementLinear3D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, CrBeamElement3D2N);
}

void CrBeamElementLinear3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, CrBeamElement3D2N);
}

} // namespace Kratos.

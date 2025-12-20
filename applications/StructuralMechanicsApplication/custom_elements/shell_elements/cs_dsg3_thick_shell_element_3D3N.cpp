// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "cs_dsg3_thick_shell_element_3D3N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/EICR.hpp"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::Initialize(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        const auto& r_geometry = GetGeometry();

        const auto& r_integration_points = CustomTriangleAreaCoordinatesQuadrature(r_geometry.Area());

        // Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_integration_points.size())
            mConstitutiveLawVector.resize(r_integration_points.size());
        InitializeMaterial();

        if constexpr (is_corotational) {
            mpCoordinateTransformation = Kratos::make_unique<ShellT3_CorotationalCoordinateTransformation>(pGetGeometry());
            mpCoordinateTransformation->Initialize();
        }
    }
    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::Initialize")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::InitializeMaterial()
{
    KRATOS_TRY

    // TODO: ensure retro-compatibility with older SHELLS using Section class

    const auto& r_properties = GetProperties();
    const auto& r_geometry = GetGeometry();
    if (r_properties[CONSTITUTIVE_LAW] != nullptr) {
        auto N_values = Vector();
        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number] = r_properties[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, N_values);
        }
    } else {
        KRATOS_ERROR << "A constitutive law needs to be specified for the CS-DSG3 shell element with ID " << this->Id() << std::endl;
    }

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::InitializeMaterial")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
Element::Pointer CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    CSDSG3ThickShellElement3D3N::Pointer p_new_elem = Kratos::make_intrusive<CSDSG3ThickShellElement3D3N>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::Clone")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dofs_per_node = GetDoFsPerNode();

    IndexType local_index = 0;

    if (rResult.size() != dofs_per_node * number_of_nodes)
        rResult.resize(dofs_per_node * number_of_nodes, false);

    const IndexType xpos    = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType rot_pos = r_geometry[0].GetDofPosition(ROTATION_X);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_X, xpos    ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Y, xpos + 1).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Z, xpos + 2).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_X, rot_pos     ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Y, rot_pos + 1 ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Z, rot_pos + 2 ).EquationId();
    }
    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::EquationIdVector")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dofs_per_node = GetDoFsPerNode(); // u, v, w, theta_x, theta_y, theta_z
    rElementalDofList.resize(dofs_per_node * number_of_nodes);
    IndexType index = 0;

    const IndexType xpos    = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType rot_pos = r_geometry[0].GetDofPosition(ROTATION_X);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rElementalDofList[index++] = r_geometry[i].pGetDof(DISPLACEMENT_X, xpos    );
        rElementalDofList[index++] = r_geometry[i].pGetDof(DISPLACEMENT_Y, xpos + 1);
        rElementalDofList[index++] = r_geometry[i].pGetDof(DISPLACEMENT_Z, xpos + 2);
        rElementalDofList[index++] = r_geometry[i].pGetDof(ROTATION_X, rot_pos     );
        rElementalDofList[index++] = r_geometry[i].pGetDof(ROTATION_Y, rot_pos + 1 );
        rElementalDofList[index++] = r_geometry[i].pGetDof(ROTATION_Z, rot_pos + 2 );
    }
    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::GetDofList")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
double CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::CalculateArea(
    const array_3& r_coord_1,
    const array_3& r_coord_2,
    const array_3& r_coord_3
) const
{
    KRATOS_TRY
    const double x21 = r_coord_2[0] - r_coord_1[0];
    const double y21 = r_coord_2[1] - r_coord_1[1];
    const double x31 = r_coord_3[0] - r_coord_1[0];
    const double y31 = r_coord_3[1] - r_coord_1[1];
    return 0.5 * (x21 * y31 - y21 * x31);
    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::CalculateArea")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::CalculateRotationMatrixGlobalToLocal(
    bounded_3_matrix& rRotationMatrix,
    const bool UseInitialConfiguration
) const
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    array_3 v1, v2, v3; // basis vectors
    array_3 aux_0, aux_1, aux_2;

    if (UseInitialConfiguration) {
        noalias(aux_0) = r_geometry[0].GetInitialPosition();
        noalias(aux_1) = r_geometry[1].GetInitialPosition();
        noalias(aux_2) = r_geometry[2].GetInitialPosition();
    } else {
        noalias(aux_0) = r_geometry[0].Coordinates();
        noalias(aux_1) = r_geometry[1].Coordinates();
        noalias(aux_2) = r_geometry[2].Coordinates();
    }

    if (this->Has(LOCAL_AXIS_1)) {
        noalias(v1) = this->GetValue(LOCAL_AXIS_1); // We assume that the user has set a unit vector
        noalias(v2) = aux_2 - aux_0;
        v2 -= inner_prod(v1, v2) * v1; // v2 orthogonal to v1
        const double norm_v2 = norm_2(v2);
        if (norm_v2 <= 1.0e-8) { // colineal
            noalias(v2) = aux_1 - aux_0;
            v2 -= inner_prod(v1, v2) * v1; // v2 orthogonal to v1
        }
        v2 /= norm_2(v2);
    } else {
        noalias(v1) = aux_1 - aux_0;
        const double norm_v1 = norm_2(v1);
        KRATOS_DEBUG_ERROR_IF_NOT(norm_v1 > 0.0) << "Zero length local axis 1 for CSDSG3ThickShellElement3D3N " << this->Id() << std::endl;
        v1 /= norm_v1;
        noalias(v2) = aux_2 - aux_0;
        v2 -= inner_prod(v1, v2) * v1; // v2 orthogonal to v1
        const double norm_v2 = norm_2(v2);
        KRATOS_DEBUG_ERROR_IF_NOT(norm_v2 > 0.0) << "Zero length local axis 2 for CSDSG3ThickShellElement3D3N " << this->Id() << std::endl;
        v2 /= norm_v2;
    }
    noalias(v3) = MathUtils<double>::CrossProduct(v1, v2);

    row(rRotationMatrix, 0) = v1;
    row(rRotationMatrix, 1) = v2;
    row(rRotationMatrix, 2) = v3;

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::CalculateRotationMatrixGlobalToLocal")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::RotateLHSToGlobal(
    MatrixType& rLHS,
    const bounded_3_matrix& rRotationMatrix // provided
) const
{
    KRATOS_TRY

    // We will perform the rotation by blocks of 3x3 for efficiency
    bounded_3_matrix local_LHS_block;
    bounded_3_matrix temp_matrix;
    bounded_3_matrix temp_matrix2;

    // Loop over the six blocks in each direction
    for (IndexType i = 0; i < 6; ++i) {
        for (IndexType j = 0; j < 6; ++j) {
            // We extract the 3x3 block
            for (IndexType k = 0; k < 3; ++k) {
                for (IndexType l = 0; l < 3; ++l) {
                    local_LHS_block(k, l) = rLHS(i * 3 + k, j * 3 + l);
                }
            }

            // We perform the rotation to the block
            noalias(temp_matrix) = prod(trans(rRotationMatrix), local_LHS_block);
            noalias(temp_matrix2) = prod(temp_matrix, rRotationMatrix);

            // We store the values back in the LHS
            for (IndexType k = 0; k < 3; ++k) {
                for (IndexType l = 0; l < 3; ++l) {
                    rLHS(i * 3 + k, j * 3 + l) = temp_matrix2(k, l);
                }
            }
        }
    }

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::RotateLHSToGlobal")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::RotateRHSToGlobal(
    VectorType& rRHS,
    const bounded_3_matrix& rRotationMatrix // provided
) const
{
    KRATOS_TRY

    // We will perform the rotation by blocks of 3 for efficiency
    array_3 local_RHS_block;
    array_3 temp_vector;

    for (IndexType i = 0; i < 6; ++i) {
        // We extract the 3 components of the block
        for (IndexType k = 0; k < 3; ++k) {
            local_RHS_block[k] = rRHS[i * 3 + k];
        }
        // We perform the rotation to the block
        noalias(temp_vector) = prod(trans(rRotationMatrix), local_RHS_block);
        // We store the values back in the RHS
        for (IndexType k = 0; k < 3; ++k) {
            rRHS[i * 3 + k] = temp_vector[k];
        }
    }
    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::RotateRHSToGlobal")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::RotateRHSToLocal(
    VectorType& rRHS,
    const bounded_3_matrix& rRotationMatrix // provided
) const
{
    KRATOS_TRY

    // We will perform the rotation by blocks of 3 for efficiency
    array_3 local_RHS_block;
    array_3 temp_vector;

    for (IndexType i = 0; i < 6; ++i) {
        // We extract the 3 components of the block
        for (IndexType k = 0; k < 3; ++k) {
            local_RHS_block[k] = rRHS[i * 3 + k];
        }
        // We perform the rotation to the block
        noalias(temp_vector) = prod(rRotationMatrix, local_RHS_block);
        // We store the values back in the RHS
        for (IndexType k = 0; k < 3; ++k) {
            rRHS[i * 3 + k] = temp_vector[k];
        }
    }
    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::RotateRHSToLocal")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::CalculateSmoothedBendingShearB(
    MatrixType& rB,
    const double Area,
    const array_3& r_coord_1,
    const array_3& r_coord_2,
    const array_3& r_coord_3
)
{
    const IndexType strain_size = GetStrainSize();
    const auto& r_geometry = GetGeometry();
    const auto& r_props = GetProperties();
    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    if (rB.size1() != strain_size || rB.size2() != system_size)
        rB.resize(strain_size, system_size, false);
    rB.clear();

    MatrixType B1(strain_size, system_size), B2(strain_size, system_size), B3(strain_size, system_size);

    const double one_third = 1.0 / 3.0;

    const array_3 initial_center = one_third * (r_coord_1 + r_coord_2 + r_coord_3);
    const double area_1 = CalculateArea(initial_center, r_coord_1, r_coord_2);
    const double area_2 = CalculateArea(initial_center, r_coord_2, r_coord_3);
    const double area_3 = CalculateArea(initial_center, r_coord_3, r_coord_1);

    CalculateBbendingShearTriangle(B1, area_1, initial_center, r_coord_1, r_coord_2);
    CalculateBbendingShearTriangle(B2, area_2, initial_center, r_coord_2, r_coord_3);
    CalculateBbendingShearTriangle(B3, area_3, initial_center, r_coord_3, r_coord_1);

    const double w1 = area_1 / Area;
    const double w2 = area_2 / Area;
    const double w3 = area_3 / Area;

    for (IndexType row = 3; row < 8; ++row) { // bending + shear rows
        for (IndexType col = 0; col < 6; ++col) {
            const IndexType c1 = col;
            const IndexType c2 = col + 6;
            const IndexType c3 = col + 12;
            // --- Block 1 ---
            rB(row, c1) += w1 * (B1(row, c1) * one_third + B1(row, c2)) + w2 * (B2(row, c1) * one_third) + w3 * (B3(row, c1) * one_third + B3(row, c3));
            // --- Block 2 ---
            rB(row, c2) += w1 * (B1(row, c1) * one_third + B1(row, c3)) + w2 * (B2(row, c1) * one_third + B2(row, c2)) + w3 * (B3(row, c1) * one_third);
            // --- Block 3 ---
            rB(row, c3) += w1 * (B1(row, c1) * one_third) + w2 * (B2(row, c1) * one_third + B2(row, c3)) + w3 * (B3(row, c1) * one_third + B3(row, c2));
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::CalculateBbendingShearTriangle(
    MatrixType& rB,
    const double Area,
    const array_3& r_local_coord_1,
    const array_3& r_local_coord_2,
    const array_3& r_local_coord_3
)
{
    KRATOS_TRY
    const IndexType strain_size = GetStrainSize();
    const IndexType number_of_nodes = GetGeometry().PointsNumber();
    const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    if (rB.size1() != strain_size || rB.size2() != system_size)
        rB.resize(strain_size, system_size, false);
    rB.clear();

    const double x1 = r_local_coord_1[0];
    const double y1 = r_local_coord_1[1];
    const double x2 = r_local_coord_2[0];
    const double y2 = r_local_coord_2[1];
    const double x3 = r_local_coord_3[0];
    const double y3 = r_local_coord_3[1];

    const double x12 = x1 - x2;
    const double x23 = x2 - x3;
    const double x31 = x3 - x1;
    const double x32 = -x23;
    const double x13 = -x31;
    const double x21 = -x12;

    const double y12 = y1 - y2;
    const double y23 = y2 - y3;
    const double y31 = y3 - y1;
    const double y32 = -y23;
    const double y13 = -y31;
    const double y21 = -y12;

    const double aux_prod = 0.5 / Area;

    // Bending components
    rB(3, 4) = y23;
    rB(4, 3) = x23;
    rB(5, 3) = y32;
    rB(5, 4) = x32;

    rB(3, 10) = y31;
    rB(4, 9) = x31;
    rB(5, 9) = y13;
    rB(5, 10) = x13;

    rB(3, 16) = y12;
    rB(4, 15) = x12;
    rB(5, 15) = y21;
    rB(5, 16) = x21;

    // Shear components
    // Bs1
    rB(6, 2) = y21 + y13;
    rB(6, 4) = Area;
    rB(7, 2) = x31 + x12;
    rB(7, 3) = -Area;
    // Bs2
    rB(6, 8) = y31;
    rB(6, 9) = y12 * y31 * 0.5;
    rB(6, 10) = x21 * y31 * 0.5;
    rB(7, 8) = x13;
    rB(7, 9) = y21 * x31 * 0.5;
    rB(7, 10) = x12 * x31 * 0.5;
    // Bs3
    rB(6, 14) = y12;
    rB(6, 15) = y21 * y31 * 0.5;
    rB(6, 16) = y12 * x31 * 0.5;
    rB(7, 14) = x21;
    rB(7, 15) = x12 * y31 * 0.5;
    rB(7, 16) = x21 * x31 * 0.5;

    rB *= aux_prod;

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::CalculateBbendingShearTriangle")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::CalculateBmTriangle(
    MatrixType& rB,
    const double Area,
    const array_3& r_local_coord_1,
    const array_3& r_local_coord_2,
    const array_3& r_local_coord_3,
    const double area_coords_1,
    const double area_coords_2,
    const double area_coords_3
)
{
    KRATOS_TRY
    const IndexType strain_size = GetStrainSize();
    const IndexType number_of_nodes = GetGeometry().PointsNumber();
    const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    if (rB.size1() != strain_size || rB.size2() != system_size)
        rB.resize(strain_size, system_size, false);
    rB.clear();

    const double x1 = r_local_coord_1[0];
    const double y1 = r_local_coord_1[1];
    const double x2 = r_local_coord_2[0];
    const double y2 = r_local_coord_2[1];
    const double x3 = r_local_coord_3[0];
    const double y3 = r_local_coord_3[1];

    const double x12 = x1 - x2;
    const double x23 = x2 - x3;
    const double x31 = x3 - x1;
    const double x32 = -x23;
    const double x13 = -x31;
    const double x21 = -x12;

    const double y12 = y1 - y2;
    const double y23 = y2 - y3;
    const double y31 = y3 - y1;
    const double y32 = -y23;
    const double y13 = -y31;
    const double y21 = -y12;

    const double aux_prod = 0.5 / Area;

    // Membrane components with drilling rotations (Zhang et al 2011)
    // CST membrane part
    rB(0, 0) = y23;
    rB(1, 1) = x32;
    rB(2, 0) = x32;
    rB(2, 1) = y23;

    rB(0, 6) = y31;
    rB(1, 7) = x13;
    rB(2, 6) = x13;
    rB(2, 7) = y31;

    rB(0, 12) = y12;
    rB(1, 13) = x21;
    rB(2, 12) = x21;
    rB(2, 13) = y12;

    const double alpha = 1.5;
    const double temp = alpha / 6.0;
    const double temp2 = 2.0 * temp;

    // Drilling rotation membrane components
    rB(0, 5) = temp * y23 * (y13 + y12);
    rB(1, 5) = temp * x32 * (x31 + x21);
    rB(2, 5) = temp2 * (x31 * y13 + x21 * y21);

    rB(0, 11) = temp * y31 * (y21 + y23);
    rB(1, 11) = temp * x13 * (x12 + x32);
    rB(2, 11) = temp2 * (x12 * y21 + x32 * y32);

    rB(0, 17) = temp * y12 * (y32 + y31);
    rB(1, 17) = temp * x21 * (x23 + x13);
    rB(2, 17) = temp2 * (x23 * y32 + x13 * y13);

    // Beta parameters for the higher order membrane part
    const double b1 = 1.0;
    const double b2 = 2.0;
    const double b3 = 1.0;
    const double b4 = 0.0;
    const double b5 = 1.0;
    const double b6 = -1.0;
    const double b7 = -1.0;
    const double b8 = -1.0;
    const double b9 = -2.0;

    bounded_3_matrix Q, Q1, Q2, Q3, Te, aux3_by_3;
    BoundedMatrix<double, 9, 3> TTu;
    BoundedMatrix<double, 3, 9> B_m_high_order;
    TTu.clear();

    const double LL21 = x21*x21 + y21*y21;
    const double LL32 = x32*x32 + y32*y32;
    const double LL13 = x13*x13 + y13*y13;
    const double A2 = 2.0 * Area;
    const double A4 = 4.0 * Area;
    Te(0, 0) = y23 * y13 * LL21;
    Te(0, 1) = y31 * y21 * LL32;
    Te(0, 2) = y12 * y32 * LL13;
    Te(1, 0) = x23 * x13 * LL21;
    Te(1, 1) = x31 * x21 * LL32;
    Te(1, 2) = x12 * x32 * LL13;
    Te(2, 0) = (y23 * x31 + x32 * y13) * LL21;
    Te(2, 1) = (y31 * x12 + x13 * y21) * LL32;
    Te(2, 2) = (y12 * x23 + x21 * y32) * LL13;
    Te /= (Area * A4);

    const double c21 = A2 / (LL21 * 3.0);
    const double c32 = A2 / (LL32 * 3.0);
    const double c13 = A2 / (LL13 * 3.0);

    Q1(0, 0) = b1 * c21;
    Q1(0, 1) = b2 * c21;
    Q1(0, 2) = b3 * c21;
    Q1(1, 0) = b4 * c32;
    Q1(1, 1) = b5 * c32;
    Q1(1, 2) = b6 * c32;
    Q1(2, 0) = b7 * c13;
    Q1(2, 1) = b8 * c13;
    Q1(2, 2) = b9 * c13;

    Q2(0, 0) = b9 * c21;
    Q2(0, 1) = b7 * c21;
    Q2(0, 2) = b8 * c21;
    Q2(1, 0) = b3 * c32;
    Q2(1, 1) = b1 * c32;
    Q2(1, 2) = b2 * c32;
    Q2(2, 0) = b6 * c13;
    Q2(2, 1) = b4 * c13;
    Q2(2, 2) = b5 * c13;

    Q3(0, 0) = b5 * c21;
    Q3(0, 1) = b6 * c21;
    Q3(0, 2) = b4 * c21;
    Q3(1, 0) = b8 * c32;
    Q3(1, 1) = b9 * c32;
    Q3(1, 2) = b7 * c32;
    Q3(2, 0) = b2 * c13;
    Q3(2, 1) = b3 * c13;
    Q3(2, 2) = b1 * c13;

    noalias(Q) = area_coords_1 * Q1;
    noalias(Q) += area_coords_2 * Q2;
    noalias(Q) += area_coords_3 * Q3;

    for (IndexType i = 0; i < 3; ++i) {
        TTu(0, i) = x32;
        TTu(1, i) = y32;

        TTu(3, i) = x13;
        TTu(4, i) = y13;

        TTu(6, i) = x21;
        TTu(7, i) = y21;
    }
    TTu(2, 0) = A4;
    TTu(5, 1) = A4;
    TTu(8, 2) = A4;
    TTu *= (1.0 / A4);

    const double beta0 = 0.5 * (1.0 - 4.0 * std::pow(GetMaterialProperty<double>(POISSON_RATIO, GetProperties()), 2));
    noalias(aux3_by_3) = (1.5 * std::sqrt(beta0)) * prod(trans(Q), trans(Te));
    noalias(B_m_high_order) = trans(prod(TTu, aux3_by_3)) / aux_prod;

    // We assemble into the global one
    BoundedVector<IndexType, 9> local_indices; // u,v,theta_z to global size
    local_indices[0] = 0;
    local_indices[1] = 1;
    local_indices[2] = 5;
    local_indices[3] = 6;
    local_indices[4] = 7;
    local_indices[5] = 11;
    local_indices[6] = 12;
    local_indices[7] = 13;
    local_indices[8] = 17;

    for (IndexType row = 0; row < 3; ++row) {
        for (IndexType col = 0; col < 9; ++col) {
            rB(row, local_indices[col]) += B_m_high_order(row, col);
        }
    }

    rB *= aux_prod;

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::CalculateBmTriangle")
}
/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::GetNodalValuesVector(
    VectorType& rNodalValues,
    const bounded_3_matrix& rT) const
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();

    if constexpr (is_corotational) {
        noalias(rNodalValues) = mpCoordinateTransformation->CalculateLocalDisplacements(
            mpCoordinateTransformation->CreateLocalCoordinateSystem(), Vector());
    } else { // Linear
        IndexType index = 0;
        array_3 displacement;
        array_3 rotation;
        for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
            noalias(displacement) = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
            noalias(rotation) = r_geometry[i].FastGetSolutionStepValue(ROTATION);

            rNodalValues[index++] = displacement[0];
            rNodalValues[index++] = displacement[1];
            rNodalValues[index++] = displacement[2];
            rNodalValues[index++] = rotation[0];
            rNodalValues[index++] = rotation[1];
            rNodalValues[index++] = rotation[2];
        }
        RotateRHSToLocal(rNodalValues, rT);
    }

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::GetNodalValuesVector")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY

    const IndexType strain_size = GetStrainSize();
    const auto& r_geometry = GetGeometry();
    const auto& r_props = GetProperties();
    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    if (rLHS.size1() != system_size || rLHS.size2() != system_size)
        rLHS.resize(system_size, system_size, false);
    rLHS.clear();

    if (rRHS.size() != system_size)
        rRHS.resize(system_size, false);
    rRHS.clear();

    bounded_3_matrix rotation_matrix;
    CalculateRotationMatrixGlobalToLocal(rotation_matrix, true);

    array_3 local_coords_1, local_coords_2, local_coords_3;
    noalias(local_coords_1) = ZeroVector(3);
    noalias(local_coords_2) = prod(rotation_matrix, r_geometry[1].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    noalias(local_coords_3) = prod(rotation_matrix, r_geometry[2].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    const double area = CalculateArea(local_coords_1, local_coords_2, local_coords_3);

    VectorType nodal_values(system_size);
    GetNodalValuesVector(nodal_values, rotation_matrix);

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    // Let's initialize the constitutive law's values
    VectorType gen_strain_vector(strain_size), gen_stress_vector(strain_size); // Generalized
    MatrixType gen_constitutive_matrix(strain_size, strain_size);
    cl_values.SetStrainVector(gen_strain_vector);
    cl_values.SetStressVector(gen_stress_vector);
    cl_values.SetConstitutiveMatrix(gen_constitutive_matrix);

    const auto& r_integration_points = CustomTriangleAreaCoordinatesQuadrature(area);
    double zeta1, zeta2, zeta3, weight;
    MatrixType B(strain_size, system_size);
    MatrixType temporal(strain_size, system_size);
    MatrixType Bm(strain_size, system_size);
    MatrixType B_bs_smoothed(strain_size, system_size);
    CalculateSmoothedBendingShearB(B_bs_smoothed, area, local_coords_1, local_coords_2, local_coords_3); // constant for all points

    for (SizeType i_point = 0; i_point < r_integration_points.size(); ++i_point) {
        zeta1 = r_integration_points[i_point].X();
        zeta2 = r_integration_points[i_point].Y();
        zeta3 = r_integration_points[i_point].Z();
        weight = r_integration_points[i_point].Weight();

        CalculateBmTriangle(Bm, area, local_coords_1, local_coords_2, local_coords_3, zeta1, zeta2, zeta3);

        noalias(B) = Bm + B_bs_smoothed;

        // We compute the strain at the integration point
        noalias(gen_strain_vector) = prod(B, nodal_values);

        // We call the constitutive law to compute the stress
        cl_values.SetStrainVector(gen_strain_vector);
        mConstitutiveLawVector[i_point]->CalculateMaterialResponseCauchy(cl_values);
        noalias(gen_stress_vector) = cl_values.GetStressVector();
        noalias(gen_constitutive_matrix) = cl_values.GetConstitutiveMatrix();

        // We integrate the LHS and RHS
        noalias(temporal) = prod(gen_constitutive_matrix, B);
        noalias(rLHS) += weight * prod(trans(B), temporal);
        noalias(rRHS) -= weight * prod(trans(B), gen_stress_vector);

    }
    if constexpr (is_corotational) {
        this->mpCoordinateTransformation->FinalizeCalculations(mpCoordinateTransformation->CreateLocalCoordinateSystem(),
                                                               Vector(),
                                                               nodal_values,
                                                               rLHS,
                                                               rRHS,
                                                               true,
                                                               true);
    } else {
        RotateLHSToGlobal(rLHS, rotation_matrix);
        RotateRHSToGlobal(rRHS, rotation_matrix);
    }
    AddBodyForces(area, rRHS);

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::CalculateLocalSystem")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY

    const IndexType strain_size = GetStrainSize();
    const auto& r_geometry = GetGeometry();
    const auto& r_props = GetProperties();
    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    if (rLHS.size1() != system_size || rLHS.size2() != system_size)
        rLHS.resize(system_size, system_size, false);
    rLHS.clear();

    bounded_3_matrix rotation_matrix;
    CalculateRotationMatrixGlobalToLocal(rotation_matrix, true);

    array_3 local_coords_1, local_coords_2, local_coords_3;
    noalias(local_coords_1) = ZeroVector(3);
    noalias(local_coords_2) = prod(rotation_matrix, r_geometry[1].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    noalias(local_coords_3) = prod(rotation_matrix, r_geometry[2].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    const double area = CalculateArea(local_coords_1, local_coords_2, local_coords_3);

    VectorType nodal_values(system_size);
    GetNodalValuesVector(nodal_values, rotation_matrix);

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    // Let's initialize the constitutive law's values
    VectorType gen_strain_vector(strain_size), gen_stress_vector(strain_size); // Generalized
    MatrixType gen_constitutive_matrix(strain_size, strain_size);
    cl_values.SetStrainVector(gen_strain_vector);
    cl_values.SetStressVector(gen_stress_vector);
    cl_values.SetConstitutiveMatrix(gen_constitutive_matrix);

    const auto& r_integration_points = CustomTriangleAreaCoordinatesQuadrature(area);
    double zeta1, zeta2, zeta3, weight;
    MatrixType B(strain_size, system_size);
    MatrixType temporal(strain_size, system_size);
    MatrixType Bm(strain_size, system_size);
    MatrixType B_bs_smoothed(strain_size, system_size);
    CalculateSmoothedBendingShearB(B_bs_smoothed, area, local_coords_1, local_coords_2, local_coords_3); // constant for all points

    for (SizeType i_point = 0; i_point < r_integration_points.size(); ++i_point) {
        zeta1 = r_integration_points[i_point].X();
        zeta2 = r_integration_points[i_point].Y();
        zeta3 = r_integration_points[i_point].Z();
        weight = r_integration_points[i_point].Weight();

        CalculateBmTriangle(Bm, area, local_coords_1, local_coords_2, local_coords_3, zeta1, zeta2, zeta3);

        noalias(B) = Bm + B_bs_smoothed;

        // We compute the strain at the integration point
        noalias(gen_strain_vector) = prod(B, nodal_values);

        // We call the constitutive law to compute the stress
        cl_values.SetStrainVector(gen_strain_vector);
        mConstitutiveLawVector[i_point]->CalculateMaterialResponseCauchy(cl_values);
        noalias(gen_constitutive_matrix) = cl_values.GetConstitutiveMatrix();

        // We integrate the LHS and RHS
        noalias(temporal) = prod(gen_constitutive_matrix, B);
        noalias(rLHS) += weight * prod(trans(B), temporal);
    }
    if constexpr (is_corotational) {
        Vector empty_vector;
        this->mpCoordinateTransformation->FinalizeCalculations(mpCoordinateTransformation->CreateLocalCoordinateSystem(),
                                                               Vector(),
                                                               nodal_values,
                                                               rLHS,
                                                               empty_vector,
                                                               true,
                                                               true);
    } else {
        RotateLHSToGlobal(rLHS, rotation_matrix);
    }

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::CalculateLeftHandSide")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY

    const IndexType strain_size = GetStrainSize();
    const auto& r_geometry = GetGeometry();
    const auto& r_props = GetProperties();
    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    if (rRHS.size() != system_size)
        rRHS.resize(system_size, false);
    rRHS.clear();

    bounded_3_matrix rotation_matrix;
    CalculateRotationMatrixGlobalToLocal(rotation_matrix, true);

    array_3 local_coords_1, local_coords_2, local_coords_3;
    noalias(local_coords_1) = ZeroVector(3);
    noalias(local_coords_2) = prod(rotation_matrix, r_geometry[1].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    noalias(local_coords_3) = prod(rotation_matrix, r_geometry[2].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    const double area = CalculateArea(local_coords_1, local_coords_2, local_coords_3);

    VectorType nodal_values(system_size);
    GetNodalValuesVector(nodal_values, rotation_matrix);

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    // Let's initialize the constitutive law's values
    VectorType gen_strain_vector(strain_size), gen_stress_vector(strain_size); // Generalized
    MatrixType gen_constitutive_matrix(strain_size, strain_size);
    cl_values.SetStrainVector(gen_strain_vector);
    cl_values.SetStressVector(gen_stress_vector);
    cl_values.SetConstitutiveMatrix(gen_constitutive_matrix);

    const auto& r_integration_points = CustomTriangleAreaCoordinatesQuadrature(area);
    double zeta1, zeta2, zeta3, weight;
    MatrixType B(strain_size, system_size);
    MatrixType Bm(strain_size, system_size);
    MatrixType B_bs_smoothed(strain_size, system_size);
    CalculateSmoothedBendingShearB(B_bs_smoothed, area, local_coords_1, local_coords_2, local_coords_3); // constant for all points

    for (SizeType i_point = 0; i_point < r_integration_points.size(); ++i_point) {
        zeta1 = r_integration_points[i_point].X();
        zeta2 = r_integration_points[i_point].Y();
        zeta3 = r_integration_points[i_point].Z();
        weight = r_integration_points[i_point].Weight();

        CalculateBmTriangle(Bm, area, local_coords_1, local_coords_2, local_coords_3, zeta1, zeta2, zeta3);

        noalias(B) = Bm + B_bs_smoothed;

        // We compute the strain at the integration point
        noalias(gen_strain_vector) = prod(B, nodal_values);

        // We call the constitutive law to compute the stress
        cl_values.SetStrainVector(gen_strain_vector);
        mConstitutiveLawVector[i_point]->CalculateMaterialResponseCauchy(cl_values);
        noalias(gen_stress_vector) = cl_values.GetStressVector();
        noalias(gen_constitutive_matrix) = cl_values.GetConstitutiveMatrix();

        // We integrate the LHS and RHS
        noalias(rRHS) -= weight * prod(trans(B), gen_stress_vector);

    }
    if constexpr (is_corotational) {
        Matrix empty_matrix;
        this->mpCoordinateTransformation->FinalizeCalculations(mpCoordinateTransformation->CreateLocalCoordinateSystem(),
                                                               Vector(),
                                                               nodal_values,
                                                               empty_matrix,
                                                               rRHS,
                                                               true,
                                                               false);
    } else {
        RotateRHSToGlobal(rRHS, rotation_matrix);
    }
    AddBodyForces(area, rRHS);

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::CalculateRightHandSide")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::AddBodyForces(
    const double Area,
    VectorType &rRightHandSideVector)
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const auto& r_props = GetProperties();

    array_3 body_forces = ZeroVector(3);
    constexpr double one_third = 1.0 / 3.0;

    double density = 0.0;
    if (r_props.Has(DENSITY)) {
        density = r_props[DENSITY];
    } else {
        // In composite shells the density is to be retrieved from the CL
        mConstitutiveLawVector[0]->GetValue(DENSITY, density);
        KRATOS_ERROR_IF(density <= 0.0) << "DENSITY is null as far as the shell CL is concerned... Please implement the GetValue(DENSITY) " <<  std::endl;
    }

    // We interpolate the volume acceleration at the center of the element
    for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
        noalias(body_forces) += r_geometry[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }
    body_forces *= density * one_third * Area; // total weight of the element
    body_forces * one_third; // we distribute it equally to the 3 nodes

    for (IndexType inode = 0; inode < r_geometry.size(); ++inode) {
        IndexType index = inode * 6;
        rRightHandSideVector[index + 0] += body_forces[0];
        rRightHandSideVector[index + 1] += body_forces[1];
        rRightHandSideVector[index + 2] += body_forces[2];
    }

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::AddBodyForces")
}


/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::FinalizeNonLinearIteration(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if constexpr (is_corotational) {
        this->mpCoordinateTransformation->FinalizeNonLinearIteration();
    }

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::FinalizeNonLinearIteration")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::FinalizeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    bool required = false;
    for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
        if (mConstitutiveLawVector[point_number]->RequiresFinalizeMaterialResponse()) {
            required = true;
            break;
        }
    }
    if (required) {
        const IndexType strain_size = GetStrainSize();
        const auto& r_geometry = GetGeometry();
        const auto& r_props = GetProperties();
        const IndexType number_of_nodes = r_geometry.PointsNumber();
        const IndexType system_size = number_of_nodes * GetDoFsPerNode();

        bounded_3_matrix rotation_matrix;
        CalculateRotationMatrixGlobalToLocal(rotation_matrix, true);

        array_3 local_coords_1, local_coords_2, local_coords_3, center;
        noalias(local_coords_1) = ZeroVector(3);
        noalias(local_coords_2) = prod(rotation_matrix, r_geometry[1].GetInitialPosition() - r_geometry[0].GetInitialPosition());
        noalias(local_coords_3) = prod(rotation_matrix, r_geometry[2].GetInitialPosition() - r_geometry[0].GetInitialPosition());
        const double area = CalculateArea(local_coords_1, local_coords_2, local_coords_3);

        VectorType nodal_values(system_size);
        GetNodalValuesVector(nodal_values, rotation_matrix);

        ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rCurrentProcessInfo);
        auto &r_cl_options = cl_values.GetOptions();
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

        // Let's initialize the constitutive law's values
        VectorType gen_strain_vector(strain_size), gen_stress_vector(strain_size);
        MatrixType gen_constitutive_matrix(strain_size, strain_size);
        cl_values.SetStrainVector(gen_strain_vector);
        cl_values.SetStressVector(gen_stress_vector);
        cl_values.SetConstitutiveMatrix(gen_constitutive_matrix);

        const auto& r_integration_points = CustomTriangleAreaCoordinatesQuadrature(area);
        double zeta1, zeta2, zeta3;
        MatrixType B(strain_size, system_size);
        MatrixType Bm(strain_size, system_size);
        MatrixType B_bs_smoothed(strain_size, system_size);
        CalculateSmoothedBendingShearB(B_bs_smoothed, area, local_coords_1, local_coords_2, local_coords_3); // constant for all points

        for (SizeType i_point = 0; i_point < r_integration_points.size(); ++i_point) {
            zeta1 = r_integration_points[i_point].X();
            zeta2 = r_integration_points[i_point].Y();
            zeta3 = r_integration_points[i_point].Z();

            CalculateBmTriangle(Bm, area, local_coords_1, local_coords_2, local_coords_3, zeta1, zeta2, zeta3);
            noalias(B) = Bm + B_bs_smoothed;

            // We compute the strain at the integration point
            noalias(gen_strain_vector) = prod(B, nodal_values);

            // We call the constitutive law to compute the stress
            cl_values.SetStrainVector(gen_strain_vector);
            mConstitutiveLawVector[i_point]->FinalizeMaterialResponseCauchy(cl_values);
        }
    }
    if constexpr (is_corotational) {
        this->mpCoordinateTransformation->FinalizeSolutionStep();
    }
    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::FinalizeSolutionStep")
}


/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::InitializeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    bool required = false;
    for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
        if (mConstitutiveLawVector[point_number]->RequiresInitializeMaterialResponse()) {
            required = true;
            break;
        }
    }
    if (required) {
        const IndexType strain_size = GetStrainSize();
        const auto& r_geometry = GetGeometry();
        const auto& r_props = GetProperties();
        const IndexType number_of_nodes = r_geometry.PointsNumber();
        const IndexType system_size = number_of_nodes * GetDoFsPerNode();

        bounded_3_matrix rotation_matrix;
        CalculateRotationMatrixGlobalToLocal(rotation_matrix, true);

        array_3 local_coords_1, local_coords_2, local_coords_3, center;
        noalias(local_coords_1) = ZeroVector(3);
        noalias(local_coords_2) = prod(rotation_matrix, r_geometry[1].GetInitialPosition() - r_geometry[0].GetInitialPosition());
        noalias(local_coords_3) = prod(rotation_matrix, r_geometry[2].GetInitialPosition() - r_geometry[0].GetInitialPosition());
        const double area = CalculateArea(local_coords_1, local_coords_2, local_coords_3);

        VectorType nodal_values(system_size);
        GetNodalValuesVector(nodal_values, rotation_matrix);

        ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rCurrentProcessInfo);
        auto &r_cl_options = cl_values.GetOptions();
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        // Let's initialize the constitutive law's values
        VectorType gen_strain_vector(strain_size), gen_stress_vector(strain_size);
        MatrixType gen_constitutive_matrix(strain_size, strain_size);
        cl_values.SetStrainVector(gen_strain_vector);
        cl_values.SetStressVector(gen_stress_vector);
        cl_values.SetConstitutiveMatrix(gen_constitutive_matrix);

        const auto& r_integration_points = CustomTriangleAreaCoordinatesQuadrature(area);
        double zeta1, zeta2, zeta3;

        MatrixType B(strain_size, system_size);
        MatrixType Bm(strain_size, system_size);
        MatrixType B_bs_smoothed(strain_size, system_size);
        CalculateSmoothedBendingShearB(B_bs_smoothed, area, local_coords_1, local_coords_2, local_coords_3); // constant for all points

        for (SizeType i_point = 0; i_point < r_integration_points.size(); ++i_point) {
            zeta1 = r_integration_points[i_point].X();
            zeta2 = r_integration_points[i_point].Y();
            zeta3 = r_integration_points[i_point].Z();

            CalculateBmTriangle(Bm, area, local_coords_1, local_coords_2, local_coords_3, zeta1, zeta2, zeta3);
            noalias(B) = Bm + B_bs_smoothed;

            // We compute the strain at the integration point
            noalias(gen_strain_vector) = prod(B, nodal_values);

            // We call the constitutive law to compute the stress
            cl_values.SetStrainVector(gen_strain_vector);
            mConstitutiveLawVector[i_point]->InitializeMaterialResponse(cl_values, ConstitutiveLaw::StressMeasure_Cauchy);
        }
    }

    if constexpr (is_corotational) {
        this->mpCoordinateTransformation->InitializeSolutionStep();
    }

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::InitializeSolutionStep")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
int CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::Check(
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto& r_properties = GetProperties();
    KRATOS_ERROR_IF_NOT(mConstitutiveLawVector[0]->GetStrainSize() == 8) << "The constitutive law used is not suitable for shell calculations, the StrainSize is NOT 8..." << std::endl;
    KRATOS_ERROR_IF_NOT(r_properties.Has(THICKNESS)) << "THICKNESS not provided for CSDSG3ThickShellElement3D3N " << this->Id() << std::endl;
    KRATOS_ERROR_IF_NOT(r_properties.GetValue(THICKNESS) > 0.0) << "Wrong value for THICKNESS in the  CSDSG3ThickShellElement3D3N " << this->Id() << std::endl;
    return mConstitutiveLawVector[0]->Check(r_properties, GetGeometry(), rCurrentProcessInfo);

    KRATOS_CATCH("CSDSG3ThickShellElement3D3N::Check")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    if (rMassMatrix.size1() != system_size || rMassMatrix.size2() != system_size)
        rMassMatrix.resize(system_size, system_size, false);
    rMassMatrix.clear();
    const auto& r_props = GetProperties();

    const bool compute_lumped_mass_matrix = StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(r_props, rCurrentProcessInfo);

    double density = 0.0;
    if (r_props.Has(DENSITY)) {
        density = r_props[DENSITY];
    } else {
        // In composite shells the density is to be retrieved from the CL
        mConstitutiveLawVector[0]->GetValue(DENSITY, density);
        KRATOS_ERROR_IF(density <= 0.0) << "DENSITY is null as far as the shell CL is concerned... Please implement the GetValue(DENSITY) " <<  std::endl;
    }
    const double thickness = r_props[THICKNESS];

    bounded_3_matrix rotation_matrix;
    CalculateRotationMatrixGlobalToLocal(rotation_matrix, true);

    array_3 local_coords_1, local_coords_2, local_coords_3;
    noalias(local_coords_1) = ZeroVector(3);
    noalias(local_coords_2) = prod(rotation_matrix, r_geometry[1].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    noalias(local_coords_3) = prod(rotation_matrix, r_geometry[2].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    const double ref_area = CalculateArea(local_coords_1, local_coords_2, local_coords_3);

    if (compute_lumped_mass_matrix) { // Lumped
        const double nodal_mass = density * thickness * ref_area / 3.0;

        for (SizeType i=0; i < number_of_nodes; i++) {
            SizeType index = i * GetDoFsPerNode();
            rMassMatrix(index, index) = nodal_mass;
            rMassMatrix(index + 1, index + 1) = nodal_mass;
            rMassMatrix(index + 2, index + 2) = nodal_mass;
        }
    } else { // Consistent mass matrix
        // General matrix form as per Felippa plane stress CST
        // The shell assumes that the density is provided by the CL if composite shell is used
        const double factor = thickness * thickness / 12.0;

        for (SizeType row = 0; row < system_size; row++) {
            if (row % 6 < 3) { // Translational entry
                for (SizeType col = 0; col < 3; col++) {
                    rMassMatrix(row, 6 * col + row % 6) = 1.0;
                }
            } else { // Rotational entry
                for (SizeType col = 0; col < 3; col++) {
                    rMassMatrix(row, 6 * col + row % 6) = factor;
                }
            }
            // Diagonal entry
            rMassMatrix(row, row) *= 2.0;
        }
        rMassMatrix *= thickness * density * ref_area / 12.0;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    const auto& r_geometry = GetGeometry();
    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(*this, rDampingMatrix, rCurrentProcessInfo, system_size);
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::save(
    Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    int IntMethod = int(GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
    rSerializer.save("pCoordinateTransformation", mpCoordinateTransformation);
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void CSDSG3ThickShellElement3D3N<IS_COROTATIONAL>::load(
    Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
    rSerializer.load("pCoordinateTransformation", mpCoordinateTransformation);
}

/***********************************************************************************/
/***********************************************************************************/

// Either corotational or linear
template class CSDSG3ThickShellElement3D3N<true>;
template class CSDSG3ThickShellElement3D3N<false>;

} // Namespace Kratos

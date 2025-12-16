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
#include "mitc4_andes_shell_thick_element_3D4N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/EICR.hpp"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::Initialize(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        const auto& r_integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

        // Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_integration_points.size())
            mConstitutiveLawVector.resize(r_integration_points.size());
        InitializeMaterial();

        if constexpr (is_corotational) {
            mpCoordinateTransformation = Kratos::make_unique<ShellQ4_CorotationalCoordinateTransformation>(pGetGeometry());
            mpCoordinateTransformation->Initialize();
        }
    }
    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::Initialize")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::InitializeMaterial()
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

    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::InitializeMaterial")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
Element::Pointer MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    MITC4AndesShellThickElement3D4N::Pointer p_new_elem = Kratos::make_intrusive<MITC4AndesShellThickElement3D4N>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::Clone")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::EquationIdVector(
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
    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::EquationIdVector")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::GetDofList(
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
        rElementalDofList[index++] = r_geometry[i].pGetDof(ROTATION_Y, rot_pos +  1);
        rElementalDofList[index++] = r_geometry[i].pGetDof(ROTATION_Z, rot_pos +  2);
    }
    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::GetDofList")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
double MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::CalculateArea(
    const array_3& r_coord_1, 
    const array_3& r_coord_2, 
    const array_3& r_coord_3,
    const array_3& r_coord_4
) const
{
    KRATOS_TRY
    // const double x21 = r_coord_2[0] - r_coord_1[0];
    // const double y21 = r_coord_2[1] - r_coord_1[1];
    // const double x31 = r_coord_3[0] - r_coord_1[0];
    // const double y31 = r_coord_3[1] - r_coord_1[1];
    // return 0.5 * (x21 * y31 - y21 * x31);
    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::CalculateArea")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::CalculateRotationMatrixGlobalToLocal(
    bounded_3_matrix& rRotationMatrix,
    const bool UseInitialConfiguration
) const
{
    KRATOS_TRY
    // const auto& r_geometry = GetGeometry();
    // array_3 v1, v2, v3; // basis vectors
    // array_3 aux_0, aux_1, aux_2;

    // if (UseInitialConfiguration) {
    //     noalias(aux_0) = r_geometry[0].GetInitialPosition();
    //     noalias(aux_1) = r_geometry[1].GetInitialPosition();
    //     noalias(aux_2) = r_geometry[2].GetInitialPosition();
    // } else {
    //     noalias(aux_0) = r_geometry[0].Coordinates();
    //     noalias(aux_1) = r_geometry[1].Coordinates();
    //     noalias(aux_2) = r_geometry[2].Coordinates();
    // }

    // if (this->Has(LOCAL_AXIS_1)) {
    //     noalias(v1) = this->GetValue(LOCAL_AXIS_1); // We assume that the user has set a unit vector
    //     noalias(v2) = aux_2 - aux_0;
    //     v2 -= inner_prod(v1, v2) * v1; // v2 orthogonal to v1
    //     const double norm_v2 = norm_2(v2);
    //     if (norm_v2 <= 1.0e-8) { // colineal
    //         noalias(v2) = aux_1 - aux_0;
    //         v2 -= inner_prod(v1, v2) * v1; // v2 orthogonal to v1
    //     }
    //     v2 /= norm_2(v2);
    // } else {
    //     noalias(v1) = aux_1 - aux_0;
    //     const double norm_v1 = norm_2(v1);
    //     KRATOS_DEBUG_ERROR_IF_NOT(norm_v1 > 0.0) << "Zero length local axis 1 for MITC4AndesShellThickElement3D4N " << this->Id() << std::endl;
    //     v1 /= norm_v1;
    //     noalias(v2) = aux_2 - aux_0;
    //     v2 -= inner_prod(v1, v2) * v1; // v2 orthogonal to v1
    //     const double norm_v2 = norm_2(v2);
    //     KRATOS_DEBUG_ERROR_IF_NOT(norm_v2 > 0.0) << "Zero length local axis 2 for MITC4AndesShellThickElement3D4N " << this->Id() << std::endl;
    //     v2 /= norm_v2; 
    // }
    // noalias(v3) = MathUtils<double>::CrossProduct(v1, v2);

    // row(rRotationMatrix, 0) = v1;
    // row(rRotationMatrix, 1) = v2;
    // row(rRotationMatrix, 2) = v3;

    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::CalculateRotationMatrixGlobalToLocal")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::RotateLHSToGlobal(
    MatrixType& rLHS,
    const bounded_3_matrix& rRotationMatrix // provided
) const
{
    KRATOS_TRY

    // We will perform the rotation by blocks of 3x3 for efficiency
    // bounded_3_matrix local_LHS_block;
    // bounded_3_matrix temp_matrix;
    // bounded_3_matrix temp_matrix2;

    // // Loop over the six blocks in each direction
    // for (IndexType i = 0; i < 6; ++i) {
    //     for (IndexType j = 0; j < 6; ++j) {
    //         // We extract the 3x3 block
    //         for (IndexType k = 0; k < 3; ++k) {
    //             for (IndexType l = 0; l < 3; ++l) {
    //                 local_LHS_block(k, l) = rLHS(i * 3 + k, j * 3 + l);
    //             }
    //         }

    //         // We perform the rotation to the block
    //         noalias(temp_matrix) = prod(trans(rRotationMatrix), local_LHS_block);
    //         noalias(temp_matrix2) = prod(temp_matrix, rRotationMatrix);

    //         // We store the values back in the LHS
    //         for (IndexType k = 0; k < 3; ++k) {
    //             for (IndexType l = 0; l < 3; ++l) {
    //                 rLHS(i * 3 + k, j * 3 + l) = temp_matrix2(k, l);
    //             }
    //         }
    //     }
    // }

    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::RotateLHSToGlobal")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::RotateRHSToGlobal(
    VectorType& rRHS,
    const bounded_3_matrix& rRotationMatrix // provided
) const
{
    KRATOS_TRY

    // We will perform the rotation by blocks of 3 for efficiency
    // array_3 local_RHS_block;
    // array_3 temp_vector;

    // for (IndexType i = 0; i < 6; ++i) {
    //     // We extract the 3 components of the block
    //     for (IndexType k = 0; k < 3; ++k) {
    //         local_RHS_block[k] = rRHS[i * 3 + k];
    //     }
    //     // We perform the rotation to the block
    //     noalias(temp_vector) = prod(trans(rRotationMatrix), local_RHS_block);
    //     // We store the values back in the RHS
    //     for (IndexType k = 0; k < 3; ++k) {
    //         rRHS[i * 3 + k] = temp_vector[k];
    //     }
    // }
    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::RotateRHSToGlobal")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::RotateRHSToLocal(
    VectorType& rRHS,
    const bounded_3_matrix& rRotationMatrix // provided
) const
{
    KRATOS_TRY

    // We will perform the rotation by blocks of 3 for efficiency
    // array_3 local_RHS_block;
    // array_3 temp_vector;

    // for (IndexType i = 0; i < 6; ++i) {
    //     // We extract the 3 components of the block
    //     for (IndexType k = 0; k < 3; ++k) {
    //         local_RHS_block[k] = rRHS[i * 3 + k];
    //     }
    //     // We perform the rotation to the block
    //     noalias(temp_vector) = prod(rRotationMatrix, local_RHS_block);
    //     // We store the values back in the RHS
    //     for (IndexType k = 0; k < 3; ++k) {
    //         rRHS[i * 3 + k] = temp_vector[k];
    //     }
    // }
    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::RotateRHSToLocal")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::CalculateShearBendingB(
    MatrixType& rB,
    const double Area,
    const array_3& r_local_coord_1,
    const array_3& r_local_coord_2,
    const array_3& r_local_coord_3,
    const array_3& r_coord_4
)
{
    KRATOS_TRY



    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::CalculateBbendingShearTriangle")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::CalculateMembraneB(
    MatrixType& rB,
    const double Area,
    const array_3& r_local_coord_1,
    const array_3& r_local_coord_2,
    const array_3& r_local_coord_3,
    const array_3& r_local_coord_4
)
{
    KRATOS_TRY


    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::CalculateBmTriangle")
}
/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::GetNodalValuesVector(
    VectorType& rNodalValues,
    const bounded_3_matrix& rT) const
{
    KRATOS_TRY
    // const auto& r_geometry = GetGeometry();

    // if constexpr (is_corotational) {
    //     noalias(rNodalValues) = mpCoordinateTransformation->CalculateLocalDisplacements(
    //         mpCoordinateTransformation->CreateLocalCoordinateSystem(), Vector());
    // } else { // Linear
    //     IndexType index = 0;
    //     array_3 displacement;
    //     array_3 rotation;
    //     for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
    //         noalias(displacement) = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
    //         noalias(rotation) = r_geometry[i].FastGetSolutionStepValue(ROTATION);

    //         rNodalValues[index++] = displacement[0];
    //         rNodalValues[index++] = displacement[1];
    //         rNodalValues[index++] = displacement[2];
    //         rNodalValues[index++] = rotation[0];
    //         rNodalValues[index++] = rotation[1];
    //         rNodalValues[index++] = rotation[2];
    //     }
    //     RotateRHSToLocal(rNodalValues, rT);
    // }

    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::GetNodalValuesVector")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY

    // const IndexType strain_size = GetStrainSize();
    // const auto& r_geometry = GetGeometry();
    // const auto& r_props = GetProperties();
    // const IndexType number_of_nodes = r_geometry.PointsNumber();
    // const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    // if (rLHS.size1() != system_size || rLHS.size2() != system_size)
    //     rLHS.resize(system_size, system_size, false);
    // rLHS.clear();

    // if (rRHS.size() != system_size)
    //     rRHS.resize(system_size, false);
    // rRHS.clear();

    // bounded_3_matrix rotation_matrix;
    // CalculateRotationMatrixGlobalToLocal(rotation_matrix, true);

    // array_3 local_coords_1, local_coords_2, local_coords_3;
    // noalias(local_coords_1) = ZeroVector(3);
    // noalias(local_coords_2) = prod(rotation_matrix, r_geometry[1].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    // noalias(local_coords_3) = prod(rotation_matrix, r_geometry[2].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    // const double area = CalculateArea(local_coords_1, local_coords_2, local_coords_3);

    // VectorType nodal_values(system_size);
    // GetNodalValuesVector(nodal_values, rotation_matrix);

    // ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    // auto &r_cl_options = cl_values.GetOptions();
    // r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    // r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    // r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    
    // // Let's initialize the constitutive law's values
    // VectorType gen_strain_vector(strain_size), gen_stress_vector(strain_size); // Generalized
    // MatrixType gen_constitutive_matrix(strain_size, strain_size);
    // cl_values.SetStrainVector(gen_strain_vector);
    // cl_values.SetStressVector(gen_stress_vector);
    // cl_values.SetConstitutiveMatrix(gen_constitutive_matrix);

    // const auto& r_integration_points = CustomTriangleAreaCoordinatesQuadrature(area);
    // double zeta1, zeta2, zeta3, weight;
    // MatrixType B(strain_size, system_size);
    // MatrixType temporal(strain_size, system_size);
    // MatrixType Bm(strain_size, system_size);
    // MatrixType B_bs_smoothed(strain_size, system_size);
    // CalculateSmoothedBendingShearB(B_bs_smoothed, area, local_coords_1, local_coords_2, local_coords_3); // constant for all points

    // for (SizeType i_point = 0; i_point < r_integration_points.size(); ++i_point) {
    //     zeta1 = r_integration_points[i_point].X();
    //     zeta2 = r_integration_points[i_point].Y();
    //     zeta3 = r_integration_points[i_point].Z();
    //     weight = r_integration_points[i_point].Weight();

    //     CalculateBmTriangle(Bm, area, local_coords_1, local_coords_2, local_coords_3, zeta1, zeta2, zeta3);

    //     noalias(B) = Bm + B_bs_smoothed;

    //     // We compute the strain at the integration point
    //     noalias(gen_strain_vector) = prod(B, nodal_values);

    //     // We call the constitutive law to compute the stress
    //     cl_values.SetStrainVector(gen_strain_vector);
    //     mConstitutiveLawVector[i_point]->CalculateMaterialResponseCauchy(cl_values);
    //     noalias(gen_stress_vector) = cl_values.GetStressVector();
    //     noalias(gen_constitutive_matrix) = cl_values.GetConstitutiveMatrix();

    //     // We integrate the LHS and RHS
    //     noalias(temporal) = prod(gen_constitutive_matrix, B);
    //     noalias(rLHS) += weight * prod(trans(B), temporal);
    //     noalias(rRHS) -= weight * prod(trans(B), gen_stress_vector);

    // }
    // if constexpr (is_corotational) {
    //     this->mpCoordinateTransformation->FinalizeCalculations(mpCoordinateTransformation->CreateLocalCoordinateSystem(),
    //                                                            Vector(),
    //                                                            nodal_values,
    //                                                            rLHS,
    //                                                            rRHS,
    //                                                            true,
    //                                                            true);
    // } else {
    //     RotateLHSToGlobal(rLHS, rotation_matrix);
    //     RotateRHSToGlobal(rRHS, rotation_matrix);
    // }
    // AddBodyForces(area, rRHS);

    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::CalculateLocalSystem")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY

    // const IndexType strain_size = GetStrainSize();
    // const auto& r_geometry = GetGeometry();
    // const auto& r_props = GetProperties();
    // const IndexType number_of_nodes = r_geometry.PointsNumber();
    // const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    // if (rLHS.size1() != system_size || rLHS.size2() != system_size)
    //     rLHS.resize(system_size, system_size, false);
    // rLHS.clear();

    // bounded_3_matrix rotation_matrix;
    // CalculateRotationMatrixGlobalToLocal(rotation_matrix, true);

    // array_3 local_coords_1, local_coords_2, local_coords_3;
    // noalias(local_coords_1) = ZeroVector(3);
    // noalias(local_coords_2) = prod(rotation_matrix, r_geometry[1].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    // noalias(local_coords_3) = prod(rotation_matrix, r_geometry[2].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    // const double area = CalculateArea(local_coords_1, local_coords_2, local_coords_3);

    // VectorType nodal_values(system_size);
    // GetNodalValuesVector(nodal_values, rotation_matrix);

    // ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    // auto &r_cl_options = cl_values.GetOptions();
    // r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    // r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    // r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    
    // // Let's initialize the constitutive law's values
    // VectorType gen_strain_vector(strain_size), gen_stress_vector(strain_size); // Generalized
    // MatrixType gen_constitutive_matrix(strain_size, strain_size);
    // cl_values.SetStrainVector(gen_strain_vector);
    // cl_values.SetStressVector(gen_stress_vector);
    // cl_values.SetConstitutiveMatrix(gen_constitutive_matrix);

    // const auto& r_integration_points = CustomTriangleAreaCoordinatesQuadrature(area);
    // double zeta1, zeta2, zeta3, weight;
    // MatrixType B(strain_size, system_size);
    // MatrixType temporal(strain_size, system_size);
    // MatrixType Bm(strain_size, system_size);
    // MatrixType B_bs_smoothed(strain_size, system_size);
    // CalculateSmoothedBendingShearB(B_bs_smoothed, area, local_coords_1, local_coords_2, local_coords_3); // constant for all points

    // for (SizeType i_point = 0; i_point < r_integration_points.size(); ++i_point) {
    //     zeta1 = r_integration_points[i_point].X();
    //     zeta2 = r_integration_points[i_point].Y();
    //     zeta3 = r_integration_points[i_point].Z();
    //     weight = r_integration_points[i_point].Weight();

    //     CalculateBmTriangle(Bm, area, local_coords_1, local_coords_2, local_coords_3, zeta1, zeta2, zeta3);

    //     noalias(B) = Bm + B_bs_smoothed;

    //     // We compute the strain at the integration point
    //     noalias(gen_strain_vector) = prod(B, nodal_values);

    //     // We call the constitutive law to compute the stress
    //     cl_values.SetStrainVector(gen_strain_vector);
    //     mConstitutiveLawVector[i_point]->CalculateMaterialResponseCauchy(cl_values);
    //     noalias(gen_constitutive_matrix) = cl_values.GetConstitutiveMatrix();

    //     // We integrate the LHS and RHS
    //     noalias(temporal) = prod(gen_constitutive_matrix, B);
    //     noalias(rLHS) += weight * prod(trans(B), temporal);
    // }
    // if constexpr (is_corotational) {
    //     this->mpCoordinateTransformation->FinalizeCalculations(mpCoordinateTransformation->CreateLocalCoordinateSystem(),
    //                                                            Vector(),
    //                                                            nodal_values,
    //                                                            rLHS,
    //                                                            Vector(),
    //                                                            true,
    //                                                            true);
    // } else {
    //     RotateLHSToGlobal(rLHS, rotation_matrix);
    // }

    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::CalculateLeftHandSide")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY

    // const IndexType strain_size = GetStrainSize();
    // const auto& r_geometry = GetGeometry();
    // const auto& r_props = GetProperties();
    // const IndexType number_of_nodes = r_geometry.PointsNumber();
    // const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    // if (rRHS.size() != system_size)
    //     rRHS.resize(system_size, false);
    // rRHS.clear();

    // bounded_3_matrix rotation_matrix;
    // CalculateRotationMatrixGlobalToLocal(rotation_matrix, true);

    // array_3 local_coords_1, local_coords_2, local_coords_3;
    // noalias(local_coords_1) = ZeroVector(3);
    // noalias(local_coords_2) = prod(rotation_matrix, r_geometry[1].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    // noalias(local_coords_3) = prod(rotation_matrix, r_geometry[2].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    // const double area = CalculateArea(local_coords_1, local_coords_2, local_coords_3);

    // VectorType nodal_values(system_size);
    // GetNodalValuesVector(nodal_values, rotation_matrix);

    // ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    // auto &r_cl_options = cl_values.GetOptions();
    // r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    // r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    // r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    
    // // Let's initialize the constitutive law's values
    // VectorType gen_strain_vector(strain_size), gen_stress_vector(strain_size); // Generalized
    // MatrixType gen_constitutive_matrix(strain_size, strain_size);
    // cl_values.SetStrainVector(gen_strain_vector);
    // cl_values.SetStressVector(gen_stress_vector);
    // cl_values.SetConstitutiveMatrix(gen_constitutive_matrix);

    // const auto& r_integration_points = CustomTriangleAreaCoordinatesQuadrature(area);
    // double zeta1, zeta2, zeta3, weight;
    // MatrixType B(strain_size, system_size);
    // MatrixType Bm(strain_size, system_size);
    // MatrixType B_bs_smoothed(strain_size, system_size);
    // CalculateSmoothedBendingShearB(B_bs_smoothed, area, local_coords_1, local_coords_2, local_coords_3); // constant for all points

    // for (SizeType i_point = 0; i_point < r_integration_points.size(); ++i_point) {
    //     zeta1 = r_integration_points[i_point].X();
    //     zeta2 = r_integration_points[i_point].Y();
    //     zeta3 = r_integration_points[i_point].Z();
    //     weight = r_integration_points[i_point].Weight();

    //     CalculateBmTriangle(Bm, area, local_coords_1, local_coords_2, local_coords_3, zeta1, zeta2, zeta3);

    //     noalias(B) = Bm + B_bs_smoothed;

    //     // We compute the strain at the integration point
    //     noalias(gen_strain_vector) = prod(B, nodal_values);

    //     // We call the constitutive law to compute the stress
    //     cl_values.SetStrainVector(gen_strain_vector);
    //     mConstitutiveLawVector[i_point]->CalculateMaterialResponseCauchy(cl_values);
    //     noalias(gen_stress_vector) = cl_values.GetStressVector();
    //     noalias(gen_constitutive_matrix) = cl_values.GetConstitutiveMatrix();

    //     // We integrate the LHS and RHS
    //     noalias(rRHS) -= weight * prod(trans(B), gen_stress_vector);

    // }
    // if constexpr (is_corotational) {
    //     this->mpCoordinateTransformation->FinalizeCalculations(mpCoordinateTransformation->CreateLocalCoordinateSystem(),
    //                                                            Vector(),
    //                                                            nodal_values,
    //                                                            Matrix(),
    //                                                            rRHS,
    //                                                            true,
    //                                                            false);
    // } else {
    //     RotateRHSToGlobal(rRHS, rotation_matrix);
    // }
    // AddBodyForces(area, rRHS);

    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::CalculateRightHandSide")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::AddBodyForces(
    const double Area,
    VectorType &rRightHandSideVector)
{
    KRATOS_TRY
    // const auto& r_geometry = GetGeometry();
    // const auto& r_props = GetProperties();

    // array_3 body_forces = ZeroVector(3);
    // constexpr double one_third = 1.0 / 3.0;

    // double density = 0.0;
    // if (r_props.Has(DENSITY)) {
    //     density = r_props[DENSITY];
    // } else {
    //     // In composite shells the density is to be retrieved from the CL
    //     mConstitutiveLawVector[0]->GetValue(DENSITY, density);
    //     KRATOS_ERROR_IF(density <= 0.0) << "DENSITY is null as far as the shell CL is concerned... Please implement the GetValue(DENSITY) " <<  std::endl;
    // }

    // // We interpolate the volume acceleration at the center of the element
    // for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
    //     noalias(body_forces) += r_geometry[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    // }
    // body_forces *= density * one_third * Area; // total weight of the element
    // body_forces * one_third; // we distribute it equally to the 3 nodes

    // for (IndexType inode = 0; inode < r_geometry.size(); ++inode) {
    //     IndexType index = inode * 6;
    //     rRightHandSideVector[index + 0] += body_forces[0];
    //     rRightHandSideVector[index + 1] += body_forces[1];
    //     rRightHandSideVector[index + 2] += body_forces[2];
    // }

    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::AddBodyForces")
}


/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::FinalizeNonLinearIteration(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if constexpr (is_corotational) {
        this->mpCoordinateTransformation->FinalizeNonLinearIteration();
    }

    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::FinalizeNonLinearIteration")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::FinalizeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // bool required = false;
    // for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
    //     if (mConstitutiveLawVector[point_number]->RequiresFinalizeMaterialResponse()) {
    //         required = true;
    //         break;
    //     }
    // }
    // if (required) {
    //     const IndexType strain_size = GetStrainSize();
    //     const auto& r_geometry = GetGeometry();
    //     const auto& r_props = GetProperties();
    //     const IndexType number_of_nodes = r_geometry.PointsNumber();
    //     const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    //     bounded_3_matrix rotation_matrix;
    //     CalculateRotationMatrixGlobalToLocal(rotation_matrix, true);

    //     array_3 local_coords_1, local_coords_2, local_coords_3, center;
    //     noalias(local_coords_1) = ZeroVector(3);
    //     noalias(local_coords_2) = prod(rotation_matrix, r_geometry[1].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    //     noalias(local_coords_3) = prod(rotation_matrix, r_geometry[2].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    //     const double area = CalculateArea(local_coords_1, local_coords_2, local_coords_3);

    //     VectorType nodal_values(system_size);
    //     GetNodalValuesVector(nodal_values, rotation_matrix);

    //     ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rCurrentProcessInfo);
    //     auto &r_cl_options = cl_values.GetOptions();
    //     r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    //     r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    //     r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        
    //     // Let's initialize the constitutive law's values
    //     VectorType gen_strain_vector(strain_size), gen_stress_vector(strain_size);
    //     MatrixType gen_constitutive_matrix(strain_size, strain_size);
    //     cl_values.SetStrainVector(gen_strain_vector);
    //     cl_values.SetStressVector(gen_stress_vector);
    //     cl_values.SetConstitutiveMatrix(gen_constitutive_matrix);

    //     const auto& r_integration_points = CustomTriangleAreaCoordinatesQuadrature(area);
    //     double zeta1, zeta2, zeta3, weight;
    //     MatrixType B(strain_size, system_size);
    //     MatrixType Bm(strain_size, system_size);
    //     MatrixType B_bs_smoothed(strain_size, system_size);
    //     CalculateSmoothedBendingShearB(B_bs_smoothed, area, local_coords_1, local_coords_2, local_coords_3); // constant for all points

    //     for (SizeType i_point = 0; i_point < r_integration_points.size(); ++i_point) {
    //         zeta1 = r_integration_points[i_point].X();
    //         zeta2 = r_integration_points[i_point].Y();
    //         zeta3 = r_integration_points[i_point].Z();

    //         CalculateBmTriangle(Bm, area, local_coords_1, local_coords_2, local_coords_3, zeta1, zeta2, zeta3);
    //         noalias(B) = Bm + B_bs_smoothed;

    //         // We compute the strain at the integration point
    //         noalias(gen_strain_vector) = prod(B, nodal_values);

    //         // We call the constitutive law to compute the stress
    //         cl_values.SetStrainVector(gen_strain_vector);
    //         mConstitutiveLawVector[i_point]->FinalizeMaterialResponseCauchy(cl_values);
    //     }
    // }
    // if constexpr (is_corotational) {
    //     this->mpCoordinateTransformation->FinalizeSolutionStep();
    // }
    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::FinalizeSolutionStep")
}


/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::InitializeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // bool required = false;
    // for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
    //     if (mConstitutiveLawVector[point_number]->RequiresInitializeMaterialResponse()) {
    //         required = true;
    //         break;
    //     }
    // }
    // if (required) {
    //     const IndexType strain_size = GetStrainSize();
    //     const auto& r_geometry = GetGeometry();
    //     const auto& r_props = GetProperties();
    //     const IndexType number_of_nodes = r_geometry.PointsNumber();
    //     const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    //     bounded_3_matrix rotation_matrix;
    //     CalculateRotationMatrixGlobalToLocal(rotation_matrix, true);

    //     array_3 local_coords_1, local_coords_2, local_coords_3, center;
    //     noalias(local_coords_1) = ZeroVector(3);
    //     noalias(local_coords_2) = prod(rotation_matrix, r_geometry[1].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    //     noalias(local_coords_3) = prod(rotation_matrix, r_geometry[2].GetInitialPosition() - r_geometry[0].GetInitialPosition());
    //     const double area = CalculateArea(local_coords_1, local_coords_2, local_coords_3);

    //     VectorType nodal_values(system_size);
    //     GetNodalValuesVector(nodal_values, rotation_matrix);

    //     ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rCurrentProcessInfo);
    //     auto &r_cl_options = cl_values.GetOptions();
    //     r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    //     r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        
    //     // Let's initialize the constitutive law's values
    //     VectorType gen_strain_vector(strain_size), gen_stress_vector(strain_size);
    //     MatrixType gen_constitutive_matrix(strain_size, strain_size);
    //     cl_values.SetStrainVector(gen_strain_vector);
    //     cl_values.SetStressVector(gen_stress_vector);
    //     cl_values.SetConstitutiveMatrix(gen_constitutive_matrix);

    //     const auto& r_integration_points = CustomTriangleAreaCoordinatesQuadrature(area);
    //     double zeta1, zeta2, zeta3, weight;

    //     MatrixType B(strain_size, system_size);
    //     MatrixType Bm(strain_size, system_size);
    //     MatrixType B_bs_smoothed(strain_size, system_size);
    //     CalculateSmoothedBendingShearB(B_bs_smoothed, area, local_coords_1, local_coords_2, local_coords_3); // constant for all points

    //     for (SizeType i_point = 0; i_point < r_integration_points.size(); ++i_point) {
    //         zeta1 = r_integration_points[i_point].X();
    //         zeta2 = r_integration_points[i_point].Y();
    //         zeta3 = r_integration_points[i_point].Z();

    //         CalculateBmTriangle(Bm, area, local_coords_1, local_coords_2, local_coords_3, zeta1, zeta2, zeta3);
    //         noalias(B) = Bm + B_bs_smoothed;

    //         // We compute the strain at the integration point
    //         noalias(gen_strain_vector) = prod(B, nodal_values);

    //         // We call the constitutive law to compute the stress
    //         cl_values.SetStrainVector(gen_strain_vector);
    //         mConstitutiveLawVector[i_point]->InitializeMaterialResponse(cl_values, ConstitutiveLaw::StressMeasure_Cauchy);
    //     }
    // }

    // if constexpr (is_corotational) {
    //     this->mpCoordinateTransformation->InitializeSolutionStep();
    // }

    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::InitializeSolutionStep")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
int MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::Check(
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto& r_properties = GetProperties();
    KRATOS_ERROR_IF_NOT(mConstitutiveLawVector[0]->GetStrainSize() == 8) << "The constitutive law used is not suitable for shell calculations, the StrainSize is NOT 8..." << std::endl;
    KRATOS_ERROR_IF_NOT(r_properties.Has(THICKNESS)) << "THICKNESS not provided for MITC4AndesShellThickElement3D4N " << this->Id() << std::endl;
    KRATOS_ERROR_IF_NOT(r_properties.GetValue(THICKNESS) > 0.0) << "Wrong value for THICKNESS in the  MITC4AndesShellThickElement3D4N " << this->Id() << std::endl;
    return mConstitutiveLawVector[0]->Check(r_properties, GetGeometry(), rCurrentProcessInfo);

    KRATOS_CATCH("MITC4AndesShellThickElement3D4N::Check")
}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{

}

/***********************************************************************************/
/***********************************************************************************/

template <bool IS_COROTATIONAL>
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::CalculateDampingMatrix(
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
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::save(
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
void MITC4AndesShellThickElement3D4N<IS_COROTATIONAL>::load(
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

template class MITC4AndesShellThickElement3D4N<true>;
template class MITC4AndesShellThickElement3D4N<false>;

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos

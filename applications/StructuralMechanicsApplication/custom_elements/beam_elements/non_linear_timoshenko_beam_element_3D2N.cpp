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

#include "non_linear_timoshenko_beam_element_3D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement3D2N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dofs_per_node = GetDoFsPerNode(); // u, v, theta
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    IndexType local_index = 0;

    if (rResult.size() != dofs_per_node * number_of_nodes)
        rResult.resize(dofs_per_node * number_of_nodes, false);

    const IndexType xpos    = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType rot_pos = r_geometry[0].GetDofPosition(ROTATION_X);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_X, xpos).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Y, xpos + 1).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Z, xpos + 2).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_X, rot_pos).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Y, rot_pos + 1).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Z, rot_pos + 2).EquationId();
    }

    KRATOS_CATCH("EquationIdVector")
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement3D2N::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dofs_per_node = GetDoFsPerNode(); // u, v, theta
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    rElementalDofList.resize(dofs_per_node * number_of_nodes);
    SizeType index = 0;

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rElementalDofList[index++]   = r_geom[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index++] = r_geom[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index++] = r_geom[i].pGetDof(DISPLACEMENT_Z);
        rElementalDofList[index++] = r_geom[i].pGetDof(ROTATION_X);
        rElementalDofList[index++] = r_geom[i].pGetDof(ROTATION_Y);
        rElementalDofList[index++] = r_geom[i].pGetDof(ROTATION_Z);
    }

    KRATOS_CATCH("GetDofList")
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer NonLinearTimoshenkoBeamElement3D2N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    NonLinearTimoshenkoBeamElement3D2N::Pointer p_new_elem = Kratos::make_intrusive<NonLinearTimoshenkoBeamElement3D2N>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    // The rotation operators
    p_new_elem->SetRotationOperators(mRotationOperators);
    return p_new_elem;

    KRATOS_CATCH("Clone");
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement3D2N::Initialize(
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        const auto &r_geom = GetGeometry();
        const auto& r_integration_points = r_geom.IntegrationPoints(mThisIntegrationMethod); // Lobatto by default

        // Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_integration_points.size())
            mConstitutiveLawVector.resize(r_integration_points.size());
        InitializeMaterial();

        BoundedMatrix<double, 3, 3> T;
        noalias(T) = CalculateInitialRotationOperator();
        noalias(mRotationOperators[0]) = T;
        noalias(mRotationOperators[1]) = T;
    }

    KRATOS_CATCH("Initialize")
}

/***********************************************************************************/
/***********************************************************************************/

BoundedMatrix<double, 3, 3> NonLinearTimoshenkoBeamElement3D2N::CalculateInitialRotationOperator(
    const bool UseCurrentConfiguration
)
{
    BoundedMatrix<double, 3, 3> T, rot_operator;
    noalias(T) = StructuralMechanicsElementUtilities::GetFrenetSerretMatrix3D(GetGeometry(), UseCurrentConfiguration); // ref conf
    // We initialize the rotation operator with the reference local system
    for (IndexType i_node = 0; i_node < 2; ++i_node) {
        // In Romero and Armero, the directors convention is different from Kratos
        noalias(column(rot_operator, 0)) = row(T, 1);
        noalias(column(rot_operator, 1)) = row(T, 2);
        noalias(column(rot_operator, 2)) = row(T, 0);
    }
    return rot_operator;
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement3D2N::InitializeMaterial()
{
    KRATOS_TRY

    if (GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
        const auto& r_geometry   = GetGeometry();
        const auto& r_properties = GetProperties();
        auto N_values            = Vector();
        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number] = r_properties[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, N_values);
        }
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH("InitializeMaterial")
}

/***********************************************************************************/
/***********************************************************************************/

int NonLinearTimoshenkoBeamElement3D2N::Check(
    const ProcessInfo& rCurrentProcessInfo
) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mConstitutiveLawVector[0]->GetStrainSize() == 6) << "The strain size of the CL is not 6, hence is not compatible with this 3D beam element" << std::endl;

    KRATOS_ERROR_IF_NOT(CalculateReferenceLength() > 0.0) << "The element has null length..." << std::endl;

    return mConstitutiveLawVector[0]->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);

    KRATOS_CATCH("Check")
}

/***********************************************************************************/
/***********************************************************************************/

Vector NonLinearTimoshenkoBeamElement3D2N::CalculateStrainVector(        
    const double N1,
    const double N2,
    const double dN1,
    const double dN2
)
{
    const auto &r_geom = GetGeometry();
    Vector generalized_strain(6);

    // The current tangent vector to the beam axis r'
    array3 dr = r_geom[1].Coordinates() - r_geom[0].Coordinates();
    const double current_L = norm_2(dr);
    dr /= current_L;

    BoundedMatrix<double, 3, 3> current_rot, d_current_rot; // current rotation and its derivative w.r.t "s"
    // We interpolate the rotation operators and its derivative
    noalias(current_rot) = N1 * mRotationOperators[0] + N2 * mRotationOperators[1];
    noalias(d_current_rot) = dN1 * mRotationOperators[0] + dN2 * mRotationOperators[1];

    array3 d1, d2, d3, d1_s, d2_s, d3_s;
    noalias(d1) = column(current_rot, 0);
    noalias(d2) = column(current_rot, 1);
    noalias(d3) = column(current_rot, 2);

    noalias(d1_s) = column(d_current_rot, 0);
    noalias(d2_s) = column(d_current_rot, 1);
    noalias(d3_s) = column(d_current_rot, 2);

    generalized_strain[0] = inner_prod(d3, dr) - 1.0; // axial
    generalized_strain[1] = inner_prod(d1, dr); // shear y
    generalized_strain[2] = inner_prod(d2, dr); // shear z

    generalized_strain[3] = inner_prod(d2, d1_s) - inner_prod(d1, d2_s); // torsion
    generalized_strain[4] = inner_prod(d3, d2_s) - inner_prod(d2, d3_s); // bending y
    generalized_strain[5] = inner_prod(d1, d3_s) - inner_prod(d3, d1_s); // bending z

    return generalized_strain;
}
/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement3D2N::CalculateGeneralizedResponse(
    const IndexType IntegrationPoint,
    ConstitutiveLaw::Parameters rValues
)
{
    // Here the stress and constitutive matrix are filled according to Kratos convention
    mConstitutiveLawVector[IntegrationPoint]->CalculateMaterialResponseCauchy(rValues);

    // We swap order to be compatible with Romero and Armero
    if (rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        auto &r_stress = rValues.GetStressVector();
        Vector temp_stress = r_stress; // we make a copy
        r_stress[0] = temp_stress[1];
        r_stress[1] = temp_stress[2];
        r_stress[2] = temp_stress[0];
        r_stress[3] = temp_stress[4];
        r_stress[4] = temp_stress[5];
        r_stress[5] = temp_stress[3];
    }

    if (rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        auto &r_D = rValues.GetConstitutiveMatrix();
        Matrix temp_D = r_D; // we make a copy
        r_D(0, 0) = temp_D(1, 1);
        r_D(1, 1) = temp_D(2, 2);
        r_D(2, 2) = temp_D(0, 0);
        r_D(3, 3) = temp_D(4, 4);
        r_D(4, 4) = temp_D(5, 5);
        r_D(5, 5) = temp_D(3, 3);
    }
}

/***********************************************************************************/
/***********************************************************************************/

BoundedMatrix<double, 12, 6> NonLinearTimoshenkoBeamElement3D2N::CalculateDoFMappingMatrix(
        const Vector &rD1,
        const Vector &rD2,
        const Vector &rD3
)
{
    BoundedMatrix<double, 12, 6> matrix;
    matrix.clear();

    noalias(project(matrix, range(0, 3), range(0, 3)))  = IdentityMatrix(3);
    noalias(project(matrix, range(3, 6), range(3, 6)))  = ConstitutiveLawUtilities<6>::CalculateSpinMatrix(rD1);
    noalias(project(matrix, range(6, 9), range(3, 6)))  = ConstitutiveLawUtilities<6>::CalculateSpinMatrix(rD2);
    noalias(project(matrix, range(9, 12), range(3, 6))) = ConstitutiveLawUtilities<6>::CalculateSpinMatrix(rD3);

    return matrix;
}

/***********************************************************************************/
/***********************************************************************************/

double NonLinearTimoshenkoBeamElement3D2N::CalculateReferenceLength() const
{
    return StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
}

/***********************************************************************************/
/***********************************************************************************/

BoundedMatrix<double, 6, 12> NonLinearTimoshenkoBeamElement3D2N::CalculateB(
    const double N1, // 0.5 * (1.0 - xi)
    const double N2, // 0.5 * (1.0 + xi)
    const double dN1, // -1.0 / L0
    const double dN2 // 1 / L0
)
{
    const auto &r_geom = GetGeometry();

    BoundedMatrix<double, 6, 12> b_matrix;

    BoundedMatrix<double, 6, 12> B1, B2; // B = [B1, B2]
    B1.clear();
    B2.clear();

    // The current tangent vector to the beam axis r'
    array3 dr = r_geom[1].Coordinates() - r_geom[0].Coordinates();
    const double current_L = norm_2(dr);
    dr /= current_L;

    // We interpolate the directors
    const array3 d1 = N1 * column(mRotationOperators[0], 0) + N2 * column(mRotationOperators[1], 0);
    const array3 d2 = N1 * column(mRotationOperators[0], 1) + N2 * column(mRotationOperators[1], 1);
    const array3 d3 = N1 * column(mRotationOperators[0], 2) + N2 * column(mRotationOperators[1], 2);
    // We interpolate the derivative of the directors (d/ds)
    const array3 d1_s = dN1 * column(mRotationOperators[0], 0) + dN2 * column(mRotationOperators[1], 0);
    const array3 d2_s = dN1 * column(mRotationOperators[0], 1) + dN2 * column(mRotationOperators[1], 1);
    const array3 d3_s = dN1 * column(mRotationOperators[0], 2) + dN2 * column(mRotationOperators[1], 2);

    BoundedMatrix<double, 12, 6> dof_mapper_1, dof_mapper_2; // map from (u, d1, d2, d3) -> (u, rotation)
    noalias(dof_mapper_1) = CalculateDoFMappingMatrix(column(mRotationOperators[0], 0), column(mRotationOperators[0], 1), column(mRotationOperators[0], 2));
    noalias(dof_mapper_2) = CalculateDoFMappingMatrix(column(mRotationOperators[1], 0), column(mRotationOperators[1], 1), column(mRotationOperators[1], 2));

    // Let's start assigning components...
    for (IndexType i = 0; i < 3; ++i) {
        // Gamma components
        // node 1
        B1(0, i)     = dN1 * d1[i];
        B1(0, i + 3) = N1 * dr[i];

        B1(1, i)     = dN1 * d2[i];
        B1(1, i + 6) = N1 * dr[i];

        B1(2, i)     = dN1 * d3[i];
        B1(2, i + 9) = N1 * dr[i];
        // node 2
        B2(0, i)     = dN2 * d1[i];
        B2(0, i + 3) = N2 * dr[i];

        B2(1, i)     = dN2 * d2[i];
        B2(1, i + 6) = N2 * dr[i];

        B2(2, i)     = dN2 * d3[i];
        B2(2, i + 9) = N2 * dr[i];

        // Omega components
        // node 1
        B1(3, i + 6) = 0.5 * (dN1 * d3[i] - N1 * d3_s[i]);
        B1(3, i + 9) = 0.5 * (dN1 * d2[i] - N1 * d2_s[i]);

        B1(4, i + 3) = 0.5 * (dN1 * d3[i] - N1 * d3_s[i]);
        B1(4, i + 9) = 0.5 * (dN1 * d1[i] - N1 * d1_s[i]);

        B1(5, i + 3) = 0.5 * (dN1 * d2[i] - N1 * d2_s[i]);
        B1(5, i + 6) = 0.5 * (dN1 * d1[i] - N1 * d1_s[i]);

        // node 2
        B2(3, i + 6) = 0.5 * (dN2 * d3[i] - N2 * d3_s[i]);
        B2(3, i + 9) = 0.5 * (dN2 * d2[i] - N2 * d2_s[i]);

        B2(4, i + 3) = 0.5 * (dN2 * d3[i] - N2 * d3_s[i]);
        B2(4, i + 9) = 0.5 * (dN2 * d1[i] - N2 * d1_s[i]);

        B2(5, i + 3) = 0.5 * (dN2 * d2[i] - N2 * d2_s[i]);
        B2(5, i + 6) = 0.5 * (dN2 * d1[i] - N2 * d1_s[i]);
    }

    noalias(project(b_matrix, range(0, 6), range(0, 6)))  = prod(B1, dof_mapper_1);
    noalias(project(b_matrix, range(0, 6), range(6, 12))) = prod(B2, dof_mapper_2);

    return b_matrix; // 6x12 matrix
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement3D2N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement3D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    VectorType dummy_rhs;
    CalculateAll(rLeftHandSideMatrix, dummy_rhs, rCurrentProcessInfo, true, false);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement3D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
    MatrixType dummy_lhs;
    CalculateAll(dummy_lhs, rRightHandSideVector, rCurrentProcessInfo, false, true);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement3D2N::CalculateAll(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo,
    const bool ComputeLHS,
    const bool ComputeRHS
)
{
    const auto &r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType mat_size = GetDoFsPerNode() * number_of_nodes;

    if (ComputeLHS) {
        if (rLHS.size1() != mat_size || rLHS.size2() != mat_size) {
            rLHS.resize(mat_size, mat_size, false);
        }
        rLHS.clear();
    }

    if (ComputeRHS) {
        if (rRHS.size() != mat_size) {
            rRHS.resize(mat_size, false);
        }
        rRHS.clear();
    }

    
}

/***********************************************************************************/
/***********************************************************************************/

// void NonLinearTimoshenkoBeamElement3D2N::CalculateOnIntegrationPoints(
//     const Variable<double>& rVariable,
//     std::vector<double>& rOutput,
//     const ProcessInfo& rProcessInfo
//     )
// {
//     const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());
//     rOutput.resize(r_integration_points.size());
//     for (std::size_t i = 0; i < rOutput.size(); ++i) {
//         rOutput[i] = 0.0;
//     }
// }

/***********************************************************************************/
/***********************************************************************************/

} // namespace Kratos

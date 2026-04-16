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

void NonLinearTimoshenkoBeamElement3D2N::Initialize(const ProcessInfo& rCurrentProcessInfo)
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
        noalias(T) = StructuralMechanicsElementUtilities::GetFrenetSerretMatrix3D(r_geom, false); // ref conf
        // We initialize the IP rotation operators with the reference local system
        noalias(mRotationOperators[0]) = trans(T);
        noalias(mRotationOperators[1]) = trans(T);
    }


    KRATOS_CATCH("Initialize")
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

int NonLinearTimoshenkoBeamElement3D2N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mConstitutiveLawVector[0]->GetStrainSize() == 6) << "The strain size of the CL is not 6, hence is not compatible with this 3D beam element" << std::endl;

    return mConstitutiveLawVector[0]->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);

    KRATOS_CATCH("Check")
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

double NonLinearTimoshenkoBeamElement3D2N::CalculateReferenceLength()
{
    return StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
}

/***********************************************************************************/
/***********************************************************************************/

BoundedMatrix<double, 6, 12> NonLinearTimoshenkoBeamElement3D2N::CalculateB(
    const double xi, // [-1, 1] with Lobatto
    const double N1, // 0.5 * (1.0 - xi)
    const double N2, // 0.5 * (1.0 + xi)
    const double dN1, // -1.0 / L0
    const double dN2 // 1 / L0
    )
{
    const auto &r_geom = GetGeometry();

    BoundedMatrix<double, 6, 12> b_matrix;
    b_matrix.clear();

    BoundedMatrix<double, 12, 6> dof_mapper_1, dof_mapper_2;
    // noalias()

    BoundedMatrix<double, 6, 24> B; // without the DoF mapping operation
    B.clear();

    BoundedMatrix<double, 6, 12> B1, B2; // B = [B1, B2]
    B1.clear();
    B2.clear();

    // The tangent vector to the beam axis r'
    array3 dr = r_geom[1].Coordinates() - r_geom[0].Coordinates();
    const double current_L = norm_2(dr);
    KRATOS_ERROR_IF_NOT(current_L > 0.0) << "The length of the 3D beam element is null..." << std::endl;
    dr /= current_L;

    // We interpolate the directors
    const array3 d1 = N1 * column(mRotationOperators[0], 0) + N2 * column(mRotationOperators[1], 0);
    const array3 d2 = N1 * column(mRotationOperators[0], 1) + N2 * column(mRotationOperators[1], 1);
    const array3 d3 = N1 * column(mRotationOperators[0], 2) + N2 * column(mRotationOperators[1], 2);
    // We interpolate the derivative of the directors (d/ds)
    const array3 d1_s = dN1 * column(mRotationOperators[0], 0) + dN2 * column(mRotationOperators[1], 0);
    const array3 d2_s = dN1 * column(mRotationOperators[0], 1) + dN2 * column(mRotationOperators[1], 1);
    const array3 d3_s = dN1 * column(mRotationOperators[0], 2) + dN2 * column(mRotationOperators[1], 2);

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
    }

    return b_matrix;
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

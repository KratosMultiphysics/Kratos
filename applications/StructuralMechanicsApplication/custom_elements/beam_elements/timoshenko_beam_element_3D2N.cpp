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
#include "timoshenko_beam_element_3D2N.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer LinearTimoshenkoBeamElement3D2N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    LinearTimoshenkoBeamElement3D2N::Pointer p_new_elem = Kratos::make_intrusive<LinearTimoshenkoBeamElement3D2N>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}


/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::EquationIdVector(
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
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_X    , rot_pos + 1).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Y    , rot_pos + 2).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Z    , rot_pos + 3).EquationId();
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dofs_per_node = GetDoFsPerNode();
    rElementalDofList.resize(dofs_per_node * number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const SizeType index = i * dofs_per_node;
        rElementalDofList[index]     = r_geom[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] = r_geom[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index + 2] = r_geom[i].pGetDof(DISPLACEMENT_Z);
        rElementalDofList[index + 3] = r_geom[i].pGetDof(ROTATION_X    );
        rElementalDofList[index + 4] = r_geom[i].pGetDof(ROTATION_Y    );
        rElementalDofList[index + 5] = r_geom[i].pGetDof(ROTATION_Z    );
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> LinearTimoshenkoBeamElement3D2N::GetLocalAxesBodyForce(
    const Element &rElement,
    const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const IndexType PointNumber
    ) const
{
    return array_1d<double, 3>();
    // const double angle = GetAngle();
    // const auto body_force = StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);

    // const double c = std::cos(angle);
    // const double s = std::sin(angle);
    // array_1d<double, 3> local_body_force = ZeroVector(3);
    // local_body_force[0] = c * body_force[0] + s * body_force[1];
    // local_body_force[1] = -s * body_force[0] + c * body_force[1];
    // return local_body_force;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::GetNodalValuesVector(
    VectorType& rNodalValues
    ) const
{
    KRATOS_TRY

    const auto& r_geom = GetGeometry();
    const SizeType num_nodes = r_geom.size();
    const SizeType global_size = GetDoFsPerNode() * num_nodes;

    if (rNodalValues.size() != global_size)
        rNodalValues.resize(global_size, false);
    
    BoundedVector<double, 12> global_values;
    BoundedMatrix<double, 3, 3> T;
    noalias(T) = StructuralMechanicsElementUtilities::GetFrenetSerretMatrix3D(r_geom);

    const auto& r_displ_0    = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT);
    const auto& r_rotation_0 = r_geom[0].FastGetSolutionStepValue(ROTATION);

    // Here we rotate the vectors to local axes
    const VectorType& r_local_displ_0 = prod(T, r_displ_0);
    const VectorType& r_local_rot_0   = prod(T, r_rotation_0);

    global_values[0] = r_local_displ_0[0];
    global_values[1] = r_local_displ_0[1];
    global_values[2] = r_local_displ_0[2];

    global_values[3] = r_local_rot_0[0];
    global_values[4] = r_local_rot_0[1];
    global_values[5] = r_local_rot_0[2];

    const auto& r_displ_1    = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT);
    const auto& r_rotation_1 = r_geom[1].FastGetSolutionStepValue(ROTATION);

    const VectorType& r_local_displ_1 = prod(T, r_displ_1);
    const VectorType& r_local_rot_1   = prod(T, r_rotation_1);

    global_values[6] = r_rotation_1[0];
    global_values[7] = r_rotation_1[1];
    global_values[8] = r_rotation_1[2];

    global_values[9]  = r_local_rot_1[0];
    global_values[10] = r_local_rot_1[1];
    global_values[11] = r_local_rot_1[2];

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateAxialStrain(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    VectorType N_u_derivatives(2);
    GetFirstDerivativesNu0ShapeFunctionsValues(N_u_derivatives, Length, Phi, xi);
    return N_u_derivatives[0] * rNodalValues[0] + N_u_derivatives[1] * rNodalValues[6];
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateShearStrainXY(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    VectorType N_derivatives(4), N_theta(4);
    GetFirstDerivativesShapeFunctionsValues(N_derivatives, Length, Phi, xi);
    GetNThetaShapeFunctionsValues(N_theta, Length, Phi, xi);
    const VectorType N_s = N_derivatives - N_theta;
    return N_s[0] * rNodalValues[1] + N_s[1] * rNodalValues[5] + N_s[2] * rNodalValues[7] + N_s[3] * rNodalValues[11];
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateShearStrainXZ(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    VectorType N_derivatives(4), N_theta(4);
    GetFirstDerivativesShapeFunctionsValues(N_derivatives, Length, Phi, xi);
    GetNThetaShapeFunctionsValues(N_theta, Length, Phi, xi);
    const VectorType N_s = N_derivatives + N_theta;
    return N_s[0] * rNodalValues[2] + N_s[1] * rNodalValues[4] + N_s[2] * rNodalValues[8] + N_s[3] * rNodalValues[10];
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateBendingCurvatureX(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    VectorType N_theta_x_derivatives(2);
    GetFirstDerivativesNu0ShapeFunctionsValues(N_theta_x_derivatives, Length, Phi, xi);
    return N_theta_x_derivatives[0] * rNodalValues[3] + N_theta_x_derivatives[1] * rNodalValues[9];
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateBendingCurvatureY(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    // VectorType N_theta_derivatives(6);
    // GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, Length, Phi, xi);
    // return N_theta_derivatives[0] * rNodalValues[1] + N_theta_derivatives[1] * rNodalValues[2] +
    //        N_theta_derivatives[4] * rNodalValues[4] + N_theta_derivatives[5] * rNodalValues[5] +
    //        N_theta_derivatives[2] * rNodalValues[7] + N_theta_derivatives[3] * rNodalValues[8];
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateBendingCurvatureZ(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    // VectorType N_theta_derivatives(6);
    // GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, Length, Phi, xi);
    // return N_theta_derivatives[0] * rNodalValues[1] + N_theta_derivatives[1] * rNodalValues[2] +
    //        N_theta_derivatives[4] * rNodalValues[4] + N_theta_derivatives[5] * rNodalValues[5] +
    //        N_theta_derivatives[2] * rNodalValues[7] + N_theta_derivatives[3] * rNodalValues[8];
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::RotateLHS(
    MatrixType& rLHS,
    const GeometryType& rGeometry
    )
{
    // const double angle = GetAngle();

    // if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    //     BoundedMatrix<double, 3, 3> T, Tt;
    //     BoundedMatrix<double, 9, 9> global_size_T, aux_product;
    //     StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T, angle);
    //     StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);
    //     noalias(aux_product) = prod(rLHS, trans(global_size_T));
    //     noalias(rLHS) = prod(global_size_T, aux_product);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::RotateRHS(
    VectorType& rRHS,
    const GeometryType& rGeometry
    )
{
    // const double angle = GetAngle();
    // if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    //     BoundedMatrix<double, 3, 3> T;
    //     BoundedMatrix<double, 9, 9> global_size_T;
    //     BoundedVector<double, 9> local_rhs;
    //     noalias(local_rhs) = rRHS;
    //     StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T, angle);
    //     StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);

    //     noalias(rRHS) = prod(global_size_T, local_rhs);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::RotateAll(
    MatrixType& rLHS,
    VectorType& rRHS,
    const GeometryType& rGeometry
    )
{
    // const double angle = GetAngle();
    // if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    //     BoundedMatrix<double, 3, 3> T;
    //     BoundedMatrix<double, 9, 9> global_size_T, aux_product;
    //     BoundedVector<double, 9> local_rhs;
    //     StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T, angle);
    //     StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);

    //     noalias(local_rhs) = rRHS;
    //     noalias(rRHS) = prod(global_size_T, local_rhs);

    //     noalias(aux_product) = prod(rLHS, trans(global_size_T));
    //     noalias(rLHS) = prod(global_size_T, aux_product);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::CalculateGeneralizedStrainsVector(
    VectorType& rStrain,
    const double Length,
    const double Phi,
    const double xi,
    const VectorType &rNodalValues
    ) const
{



}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::save(
    Serializer& rSerializer
    ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::load(
    Serializer& rSerializer
    )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos

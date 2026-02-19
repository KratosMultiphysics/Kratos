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

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer NonLinearTimoshenkoBeamElement3D2N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    NonLinearTimoshenkoBeamElement3D2N::Pointer p_new_elem = Kratos::make_intrusive<NonLinearTimoshenkoBeamElement3D2N>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));
    return p_new_elem;
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement3D2N::BuildShapeFunctionsValuesVectors(
    const double xi, // Local coordinate in the beam axis direction
    const double L, // Reference length of the beam element
    const Vector& rNodalValues, // Vector containing the nodal values in local axes
    BeamElementData &rData)
{
    Vector local_N(2); // Shape functions in local axes
    Vector local_dN(2); // Derivatives of shape functions in local axes
    GetNu0ShapeFunctionsValues(local_N, L, 0.0, xi);
    GetFirstDerivativesNu0ShapeFunctionsValues(local_dN, L, 0.0, xi);

    rData.du = local_dN[0] * rNodalValues[0] + local_dN[1] * rNodalValues[6];
    rData.dv = local_dN[0] * rNodalValues[1] + local_dN[1] * rNodalValues[7];
    rData.dw = local_dN[0] * rNodalValues[2] + local_dN[1] * rNodalValues[8];

    rData.theta_x = local_N[0] * rNodalValues[3] + local_N[1] * rNodalValues[9];
    rData.theta_y = local_N[0] * rNodalValues[4] + local_N[1] * rNodalValues[10];
    rData.theta_z = local_N[0] * rNodalValues[5] + local_N[1] * rNodalValues[11];

    rData.dtheta_x = local_dN[0] * rNodalValues[3] + local_dN[1] * rNodalValues[9];
    rData.dtheta_y = local_dN[0] * rNodalValues[4] + local_dN[1] * rNodalValues[10];
    rData.dtheta_z = local_dN[0] * rNodalValues[5] + local_dN[1] * rNodalValues[11];

    rData.L = L;

    rData.cx = std::cos(rData.theta_x);
    rData.cy = std::cos(rData.theta_y);
    rData.cz = std::cos(rData.theta_z);

    rData.sx = std::sin(rData.theta_x);
    rData.sy = std::sin(rData.theta_y);
    rData.sz = std::sin(rData.theta_z);

    // Assemble the shape functions and their derivatives
    for (IndexType i = 0; i < 2; ++i) {
        const IndexType index = i * 6; // 6 DoFs per node

        rData.Nu[index    ] = local_N[i];
        rData.Nv[index + 1] = local_N[i];
        rData.Nw[index + 2] = local_N[i];

        rData.dNu[index    ] = local_dN[i];
        rData.dNv[index + 1] = local_dN[i];
        rData.dNw[index + 2] = local_dN[i];

        rData.Ntheta_x[index + 3] = local_N[i];
        rData.Ntheta_y[index + 4] = local_N[i];
        rData.Ntheta_z[index + 5] = local_N[i];

        rData.dNtheta_x[index + 3] = local_dN[i];
        rData.dNtheta_y[index + 4] = local_dN[i];
        rData.dNtheta_z[index + 5] = local_dN[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement3D2N::CalculateStrains(
    const double xi, // Local coordinate in the beam axis direction
    BeamElementData &rData)
{
    // rData.longitudinal_strain
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement3D2N::CalculateB(
    const double xi, // Local coordinate in the beam axis direction
    BeamElementData &rData)
{

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
        noalias(rLHS) = ZeroMatrix(mat_size, mat_size);
    }

    if (ComputeRHS) {
        if (rRHS.size() != mat_size) {
            rRHS.resize(mat_size, false);
        }
        noalias(rRHS) = ZeroVector(mat_size);
    }

    // Dummy implementation: no internal computations yet.
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

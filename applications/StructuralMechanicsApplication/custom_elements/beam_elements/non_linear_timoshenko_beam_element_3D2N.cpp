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

// void NonLinearTimoshenkoBeamElement3D2N::Initialize(const ProcessInfo& rCurrentProcessInfo)
// {
//     // Dummy/empty initialize
// }

// void NonLinearTimoshenkoBeamElement3D2N::GetShapeFunctionsValues(
//     VectorType& rN,
//     const double Length,
//     const double Phi,
//     const double xi) const
// {
//     rN.clear();
// }

// void NonLinearTimoshenkoBeamElement3D2N::GetFirstDerivativesShapeFunctionsValues(
//     VectorType& rN,
//     const double Length,
//     const double Phi,
//     const double xi) const
// {
//     rN.clear();
// }

// void NonLinearTimoshenkoBeamElement3D2N::GetNThetaShapeFunctionsValues(
//     VectorType& rN,
//     const double Length,
//     const double Phi,
//     const double xi) const
// {
//     rN.clear();
// }

// void NonLinearTimoshenkoBeamElement3D2N::GetFirstDerivativesNThetaShapeFunctionsValues(
//     VectorType& rN,
//     const double Length,
//     const double Phi,
//     const double xi) const
// {
//     rN.clear();
// }

// void NonLinearTimoshenkoBeamElement3D2N::GetNu0ShapeFunctionsValues(
//     VectorType& rN,
//     const double Length,
//     const double Phi,
//     const double xi) const
// {
//     rN.clear();
// }

// void NonLinearTimoshenkoBeamElement3D2N::GetFirstDerivativesNu0ShapeFunctionsValues(
//     VectorType& rN,
//     const double Length,
//     const double Phi,
//     const double xi) const
// {
//     rN.clear();
// }

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

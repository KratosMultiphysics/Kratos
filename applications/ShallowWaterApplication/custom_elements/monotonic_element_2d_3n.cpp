//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "monotonic_element_2d_3n.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

void MonotonicElement2D3N::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rMassMatrix.size1() != 9) {
        rMassMatrix.resize(9,9,false);
    }

    rMassMatrix = ZeroMatrix(9,9);

    const double area = this->GetGeometry().Area();
    const double one_twelve = 1.0 / 12.0;

    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        for (size_t j = 0; j < 3; ++j)
        {
            const size_t j_block = 3 * j;

            // Algebraic mass matrix
            const double n = (i == j)? 2*one_twelve : one_twelve;
            rMassMatrix(i_block,     j_block    ) = area * n;
            rMassMatrix(i_block + 1, j_block + 1) = area * n;
            rMassMatrix(i_block + 2, j_block + 2) = area * n;
        }
    }
}

void MonotonicElement2D3N::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rDampingMatrix.size1() != 9) {
        rDampingMatrix.resize(9,9,false);
    }

    rDampingMatrix = ZeroMatrix(9,9);

    const double area = this->GetGeometry().Area();
    const double one_third = 1.0 / 3.0;
    const double one_twelve = 1.0 / 12.0;
    const double g = rCurrentProcessInfo[GRAVITY_Z];

    // Computing the diffusion factor (inverse of characteristic time)
    array_1d<double,3> v = ZeroVector(3);
    double h = 0.0;
    for (auto& r_node : this->GetGeometry())
    {
        v += r_node.FastGetSolutionStepValue(VELOCITY);
        h += r_node.FastGetSolutionStepValue(HEIGHT);
    }
    h = std::max(h, 0.0);
    const double lambda = norm_2(v) + std::sqrt(g*h);
    const double c = lambda / this->GetGeometry().Length();

    // Construction of the scalar-based diffusion matrix
    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        for (size_t j = 0; j < 3; ++j)
        {
            const size_t j_block = 3 * j;

            const double d = (i == j)? one_third - 2*one_twelve : -one_twelve;
            rDampingMatrix(i_block,     j_block    ) = c * area * d;
            rDampingMatrix(i_block + 1, j_block + 1) = c * area * d;
            rDampingMatrix(i_block + 2, j_block + 2) = c * area * d;
        }
    }
}

void MonotonicElement2D3N::AddGradientTerms(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    const double c2 = rData.gravity * rData.effective_height; // c=sqrt(gh)
    const double u_1 = rData.velocity[0];
    const double u_2 = rData.velocity[1];

    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        for (size_t j = 0; j < 3; ++j)
        {
            const size_t j_block = 3 * j;

            /* First component
             * A_1 = {{2 * u_1   0   -u_1^2 + gh},
             *        {  u_2    u_1   -u_1 * u_2},
             *        {   1      0        0     }}
             */
            const double g1_ij = rN[i] * rDN_DX(j,0);
            rLHS(i_block,     j_block    ) += g1_ij * 2*u_1;
            rLHS(i_block,     j_block + 2) += g1_ij * (-u_1*u_1 + c2);
            rLHS(i_block + 1, j_block    ) += g1_ij * u_2;
            rLHS(i_block + 1, j_block + 1) += g1_ij * u_1;
            rLHS(i_block + 1, j_block + 2) -= g1_ij * u_1*u_2;
            rLHS(i_block + 2, j_block    ) += g1_ij;

            /* Second component
             * A_2 = {{u_2    u_1      -u_1 * u_2},
             *        { 0   2 * u_2   -u_2^2 + gh},
             *        { 0      1            0    }}
             */
            const double g2_ij = rN[i] * rDN_DX(j,1);
            rLHS(i_block,     j_block    ) += g2_ij * u_2;
            rLHS(i_block,     j_block + 1) += g2_ij * u_1;
            rLHS(i_block,     j_block + 2) -= g2_ij * u_1*u_2;
            rLHS(i_block + 1, j_block + 1) += g2_ij * 2*u_2;
            rLHS(i_block + 1, j_block + 2) += g2_ij * (-u_2*u_2 + c2);
            rLHS(i_block + 2, j_block + 1) += g2_ij;
        }
    }
}

void MonotonicElement2D3N::AddSourceTerms(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ElementData& rData,
    const array_1d<double,3>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    const double c2 = rData.gravity * rData.effective_height; // c=sqrt(gh)
    const double abs_vel = norm_2(rData.velocity);
    const double height4_3 = std::pow(rData.effective_height, 1.33333333333) + 1e-6;
    const double friction = rData.gravity * rData.manning2 * abs_vel / height4_3;
    const double one_third = 1.0 / 3.0;

    for (size_t i = 0; i < 3; ++i)
    {
        const size_t i_block = 3 * i;
        for (size_t j = 0; j < 3; ++j)
        {
            const size_t j_block = 3 * j;

            /* gravity term */
            rRHS[i_block] -= c2 * rN[i] * rDN_DX(j,0) * rData.depth[j_block + 2];
            rRHS[i_block + 1] -= c2 * rN[i] * rDN_DX(j,1) * rData.depth[j_block + 2];

        }
        /* friction term */
        rLHS(i_block, i_block) += friction * one_third;
        rLHS(i_block + 1, i_block + 1) += friction * one_third;

        /* rain source */
        rRHS[i_block + 2] += rData.rain[i_block + 2] * one_third; 
    }
}

} // namespace kratos

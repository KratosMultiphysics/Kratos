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

} // namespace kratos

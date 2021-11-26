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
#include "includes/checks.h"
#include "boussinesq_element.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/phase_function.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

template<std::size_t TNumNodes>
const Variable<double>& BoussinesqElement<TNumNodes>::GetUnknownComponent(int Index) const
{
    switch (Index) {
        case 0: return VELOCITY_X;
        case 1: return VELOCITY_Y;
        case 2: return FREE_SURFACE_ELEVATION;
        default: KRATOS_ERROR << "BoussinesqElement::GetUnknownComponent index out of bounds." << std::endl;
    }
}

template<std::size_t TNumNodes>
typename BoussinesqElement<TNumNodes>::LocalVectorType BoussinesqElement<TNumNodes>::GetUnknownVector(const ElementData& rData) const
{
    std::size_t index = 0;
    array_1d<double,mLocalSize> unknown;
    for (std::size_t i = 0; i < TNumNodes; ++i) {
        unknown[index++] = rData.nodal_v[i][0];
        unknown[index++] = rData.nodal_v[i][1];
        unknown[index++] = rData.nodal_f[i];
    }
    return unknown;
}

template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::CalculateGaussPointData(ElementData& rData, const array_1d<double,TNumNodes>& rN)
{
    const double eta = inner_prod(rData.nodal_f, rN);
    const double H = -inner_prod(rData.nodal_z, rN);
    const double g = rData.gravity;
    const double e = rData.amplitude / H;  // the non linearity ratio
    const array_1d<double,3> v = WaveElementType::VectorProduct(rData.nodal_v, rN);

    rData.depth = H;
    rData.height = H + e * eta;
    rData.velocity = v;

    /**
     * A_1 = {{ e*u_1      0     g  },
     *        {   0      e*u_1   0  },
     *        {H + e*eta   0   e*u_1}}
     */
    rData.A1(0,0) = e*v[0];
    rData.A1(0,1) = 0;
    rData.A1(0,2) = g;
    rData.A1(1,0) = 0;
    rData.A1(1,1) = e*v[0];
    rData.A1(1,2) = 0;
    rData.A1(2,0) = H + e*eta;
    rData.A1(2,1) = 0;
    rData.A1(2,2) = e*v[0];

    /*
     * A_2 = {{ e*u_2      0      0  },
     *        {   0      e*u_2    g  },
     *        {   0   H + e*eta e*u_2}}
     */
    rData.A2(0,0) = e*v[1];
    rData.A2(0,1) = 0;
    rData.A2(0,2) = 0;
    rData.A2(1,0) = 0;
    rData.A2(1,1) = e*v[1];
    rData.A2(1,2) = g;
    rData.A2(2,0) = 0;
    rData.A2(2,1) = H + e*eta;
    rData.A2(2,2) = e*v[1];

    /// b_1
    rData.b1[0] = 0;
    rData.b1[1] = 0;
    rData.b1[2] = v[0];

    /// b_2
    rData.b2[0] = 0;
    rData.b2[1] = 0;
    rData.b2[2] = v[1];
}

template<std::size_t TNumNodes>
double BoussinesqElement<TNumNodes>::StabilizationParameter(const ElementData& rData) const
{
    const double lambda = std::sqrt(rData.gravity * std::abs(rData.height)) + norm_2(rData.velocity);
    const double epsilon = 1e-6;
    const double threshold = rData.relative_dry_height * rData.length;
    const double w = PhaseFunction::WetFraction(rData.height, threshold);
    return w * rData.length * rData.stab_factor / (lambda + epsilon);
}

template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::AddDispersiveTerms(
    LocalMatrixType& rMatrix,
    LocalVectorType& rVector,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    array_1d<double,3> gradients_vector_i;
    array_1d<double,3> gradients_vector_j;
    BoundedMatrix<double,3,3> gradients_matrix_k = ZeroMatrix(3,3);
    BoundedMatrix<double,3,3> laplacian;
    BoundedMatrix<double,3,3> height_aux;
    BoundedMatrix<double,3,3> velocity_aux;
    BoundedMatrix<double,3,3> height_stab_aux;
    BoundedMatrix<double,3,3> velocity_stab_aux;
    LocalMatrixType height_dispersion = ZeroMatrix(mLocalSize, mLocalSize);
    LocalMatrixType velocity_dispersion = ZeroMatrix(mLocalSize, mLocalSize);

    const double dispersive_ratio = rData.depth / rData.wavelength;
    const double m2 = std::pow(dispersive_ratio, 2);
    const double beta = -0.531;
    const double C1 = 0.5 * std::pow(beta, 2) - 0.166666666666;
    const double C2 = 0.5 * std::pow(beta, 2);
    const double C3 = beta + 0.5;
    const double C4 = beta;
    const double H = rData.depth;
    const double inv_lumped_mass = double(TNumNodes);
    const double l = StabilizationParameter(rData);

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        gradients_vector_i[0] = rDN_DX(i,0);
        gradients_vector_i[1] = rDN_DX(i,1);
        gradients_vector_i[2] = 0.0;

        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            gradients_vector_j[0] = rDN_DX(j,0);
            gradients_vector_j[1] = rDN_DX(j,1);
            gradients_vector_j[2] = 0.0;

            laplacian = -outer_prod(gradients_vector_i, gradients_vector_j);
            gradients_matrix_k(2,0) = rDN_DX(i,0);
            gradients_matrix_k(2,1) = rDN_DX(i,1);

            // Dispersive fields
            height_aux = (C1 + C3) * std::pow(H,3) * inv_lumped_mass * prod(gradients_matrix_k, laplacian);
            velocity_aux = (C2 + C4) * std::pow(H,2) * laplacian;

            // Stabilization terms
            height_stab_aux  = l * rDN_DX(j,0) * prod(rData.A1, height_aux);
            height_stab_aux += l * rDN_DX(j,1) * prod(rData.A2, height_aux);
            velocity_stab_aux = l * rDN_DX(i,0) * inv_lumped_mass * prod(rData.A1, height_aux);
            velocity_stab_aux = l * rDN_DX(i,1) * inv_lumped_mass * prod(rData.A2, height_aux);

            // Contribution to the local matrix
            MathUtils<double>::AddMatrix(height_dispersion, Weight * height_aux, 3*i, 3*j);
            MathUtils<double>::AddMatrix(velocity_dispersion, Weight * velocity_aux, 3*i, 3*j);
            MathUtils<double>::AddMatrix(height_dispersion, Weight * height_stab_aux, 3*i, 3*j);
            MathUtils<double>::AddMatrix(velocity_dispersion, Weight * velocity_stab_aux, 3*i, 3*j);
        }
    }
    LocalVectorType acceleration = this->GetAccelerationsVector(rData);
    LocalVectorType velocity = this->GetUnknownVector(rData);

    noalias(rVector) -= m2 * prod(height_dispersion, velocity);
    // noalias(rVector) -= m2 * prod(velocity_dispersion, acceleration);
}

template class BoussinesqElement<3>;

} // namespace Kratos

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
const Parameters BoussinesqElement<TNumNodes>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "required_variables"         : ["VELOCITY","FREE_SURFACE_ELEVATION","TOPOGRAPHY","ACCELERATION","VERTICAL_VELOCITY","VELOCITY_LAPLACIAN","VELOCITY_H_LAPLACIAN"],
        "required_dofs"              : ["VELOCITY_X","VELOCITY_Y","FREE_SURFACE_ELEVATION"],
        "compatible_geometries"      : ["Triangle2D3"],
        "element_integrates_in_time" : false
    })");
    return specifications;
}

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
void BoussinesqElement<TNumNodes>::GetNodalData(ElementData& rData, const GeometryType& rGeometry, int Step)
{
    for (std::size_t i = 0; i < TNumNodes; i++)
    {
        rData.nodal_f[i] = rGeometry[i].FastGetSolutionStepValue(FREE_SURFACE_ELEVATION, Step);
        rData.nodal_w[i] = rGeometry[i].FastGetSolutionStepValue(VERTICAL_VELOCITY, Step);
        rData.nodal_z[i] = rGeometry[i].FastGetSolutionStepValue(TOPOGRAPHY, Step);
        rData.nodal_v[i] = rGeometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        rData.nodal_v_lap[i] = rGeometry[i].FastGetSolutionStepValue(VELOCITY_LAPLACIAN, Step);
        rData.nodal_q_lap[i] = rGeometry[i].FastGetSolutionStepValue(VELOCITY_H_LAPLACIAN, Step);
    }
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::UpdateGaussPointData(ElementData& rData, const array_1d<double,TNumNodes>& rN)
{
    const double eta = inner_prod(rData.nodal_f, rN);
    const double H = -inner_prod(rData.nodal_z, rN);
    const double g = rData.gravity;
    const array_1d<double,3> v = WaveElementType::VectorProduct(rData.nodal_v, rN);

    rData.depth = std::max(0.0, H);
    rData.height = H + eta;
    rData.velocity = v;

    /**
     * A_1 = {{ u_1     0    g  },
     *        {  0     u_1   0  },
     *        {H + eta  0   u_1}}
     */
    rData.A1(0,0) = v[0];
    rData.A1(0,1) = 0;
    rData.A1(0,2) = g;
    rData.A1(1,0) = 0;
    rData.A1(1,1) = v[0];
    rData.A1(1,2) = 0;
    rData.A1(2,0) = H + eta;
    rData.A1(2,1) = 0;
    rData.A1(2,2) = v[0];

    /*
     * A_2 = {{ u_2    0      0  },
     *        {  0    u_2     g  },
     *        {  0  H + eta  u_2}}
     */
    rData.A2(0,0) = v[1];
    rData.A2(0,1) = 0;
    rData.A2(0,2) = 0;
    rData.A2(1,0) = 0;
    rData.A2(1,1) = v[1];
    rData.A2(1,2) = g;
    rData.A2(2,0) = 0;
    rData.A2(2,1) = H + eta;
    rData.A2(2,2) = v[1];

    /// b_1
    rData.b1[0] = 0;
    rData.b1[1] = 0;
    rData.b1[2] = -v[0];

    /// b_2
    rData.b2[0] = 0;
    rData.b2[1] = 0;
    rData.b2[2] = -v[1];
}


template<std::size_t TNumNodes>
double BoussinesqElement<TNumNodes>::StabilizationParameter(const ElementData& rData) const
{
    const double eigenvalue = norm_2(rData.velocity) + std::sqrt(rData.gravity * std::abs(rData.height));
    const double epsilon = 1e-6;
    const double w = WaveElementType::WetFraction(rData);
    return w * rData.length * rData.stab_factor / (eigenvalue + epsilon);
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::CalculateArtificialViscosity(
    BoundedMatrix<double,3,3>& rViscosity,
    BoundedMatrix<double,2,2>& rDiffusion,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX)
{
    // Computation of the residuals and gradients
    double residual;
    array_1d<double,2> gradient;
    AlgebraicResidual(residual, gradient, rData, rN, rDN_DX);

    // slope limits
    const double min_slope = 0.1;
    const double max_slope = 1.0;

    // Computation of the artificial viscosity/diffusion
    const double residual_norm = std::abs(residual);
    const double gradient_norm = std::min(std::max(norm_2(gradient), min_slope), max_slope);
    const double artificial_diffusion = 0.5 * rData.shock_stab_factor * rData.length * residual_norm / gradient_norm;

    // Assembly of the shock capturing tensors
    rDiffusion = artificial_diffusion * IdentityMatrix(2);
    rViscosity = artificial_diffusion * IdentityMatrix(3);
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::AlgebraicResidual(
    double& rMassResidual,
    array_1d<double,2>& rFreeSurfaceGradient,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX) const
{
    // constants
    const double beta = -0.531;
    const double C1 = 0.5 * std::pow(beta, 2) - 0.166666666666;
    const double C3 = beta + 0.5;
    const double H = rData.depth;
    const double H3 = std::pow(H, 3);
    const double H2 = std::pow(H, 2);

    // Spatial derivatives
    rFreeSurfaceGradient = prod(rData.nodal_f, rDN_DX);
    const double v_divergence = WaveElementType::VectorDivergence(rData.nodal_v, rDN_DX);

    // Mass conservation residual
    const double vertical_vel = inner_prod(rData.nodal_w, rN);
    const double wave_f = rData.height * v_divergence;
    const double convection_f = rData.velocity[0] * rFreeSurfaceGradient[0] + rData.velocity[1] * rFreeSurfaceGradient[1];
    double dispersion_f = C1 * H3 * WaveElementType::VectorDivergence(rData.nodal_v_lap, rDN_DX);
    dispersion_f       += C3 * H2 * WaveElementType::VectorDivergence(rData.nodal_q_lap, rDN_DX);
    rMassResidual = vertical_vel + wave_f + convection_f + dispersion_f;
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::AddDispersiveTerms(
    LocalVectorType& rVector,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    // Constants
    const double beta = -0.531;
    const double C1 = 0.5 * std::pow(beta, 2) - 0.166666666666;
    const double C3 = beta + 0.5;
    const double H = rData.depth;
    const double H3 = std::pow(H, 3);
    const double H2 = std::pow(H, 2);

    // Projected auxiliary field
    array_1d<array_1d<double,3>,TNumNodes> nodal_Jf;
    for (std::size_t i = 0; i < TNumNodes; ++i) {
        nodal_Jf[i] = C1 * H3 * rData.nodal_v_lap[i] + C3 * H2 * rData.nodal_q_lap[i];
    }

    // Stabilization constants
    const double l = StabilizationParameter(rData);
    const array_1d<double,3> A1i3 = column(rData.A1, 2);
    const array_1d<double,3> A2i3 = column(rData.A2, 2);

    // Adding the contribution of the dispersive fields
    for (std::size_t i = 0; i < TNumNodes; ++i)
    {
        for (std::size_t j = 0; j < TNumNodes; ++j)
        {
            double g1_ij;
            double g2_ij;
            double d_ij;
            if (rData.integrate_by_parts) {
                g1_ij = -rDN_DX(i,0) * rN[j];
                g2_ij = -rDN_DX(i,1) * rN[j];
            } else {
                g1_ij = rN[i] * rDN_DX(j,0);
                g2_ij = rN[i] * rDN_DX(j,1);
            }

            /// Gradient contribution
            rVector[3*i + 2] -= Weight*g1_ij*nodal_Jf[j][0];
            rVector[3*i + 2] -= Weight*g2_ij*nodal_Jf[j][1];

            /// Stabilization x-x
            d_ij = rDN_DX(i,0) * rDN_DX(j,0);
            MathUtils<double>::AddVector(rVector, -Weight*l*d_ij*A1i3*nodal_Jf[j][0], 3*i);

            /// Stabilization y-y
            d_ij = rDN_DX(i,1) * rDN_DX(j,1);
            MathUtils<double>::AddVector(rVector, -Weight*l*d_ij*A2i3*nodal_Jf[j][1], 3*i);

            /// Stabilization x-y
            d_ij = rDN_DX(i,0) * rDN_DX(j,1);
            MathUtils<double>::AddVector(rVector, -Weight*l*d_ij*A1i3*nodal_Jf[j][1], 3*i);

            /// Stabilization y-x
            d_ij = rDN_DX(i,1) * rDN_DX(j,0);
            MathUtils<double>::AddVector(rVector, -Weight*l*d_ij*A2i3*nodal_Jf[j][0], 3*i);
        }
    }
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::AddMassTerms(
    LocalMatrixType& rMatrix,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    // Constants
    const double beta = -0.531;
    const double C2 = 0.5 * std::pow(beta, 2);
    const double C4 = beta;
    const double H = rData.depth;
    const double H2 = std::pow(H, 2);
    const double l = StabilizationParameter(rData);
    BoundedMatrix<double,3,3> M = IdentityMatrix(3);
    BoundedMatrix<double,2,2> K;
    BoundedMatrix<double,2,2> Ju;
    array_1d<double,2> derivatives_i;
    array_1d<double,2> derivatives_j;

    for (std::size_t i = 0; i < TNumNodes; ++i)
    {
        derivatives_i[0] = rDN_DX(i,0);
        derivatives_i[1] = rDN_DX(i,1);
        for (std::size_t j = 0; j < TNumNodes; ++j)
        {
            /// Consistent mass matrix
            const double n_ij = WaveElementType::ShapeFunctionProduct(rN, i, j);
            MathUtils<double>::AddMatrix(rMatrix, Weight*n_ij*M, 3*i, 3*j);

            /// Stabilization x
            const double g1_ij = rDN_DX(i,0) * rN[j];
            MathUtils<double>::AddMatrix(rMatrix, Weight*l*g1_ij*rData.A1, 3*i, 3*j);

            /// Stabilization y
            const double g2_ij = rDN_DX(i,1) * rN[j];
            MathUtils<double>::AddMatrix(rMatrix, Weight*l*g2_ij*rData.A2, 3*i, 3*j);

            /// Dispersive term
            derivatives_j[0] = rDN_DX(j,0);
            derivatives_j[1] = rDN_DX(j,1);
            noalias(K) = -outer_prod(derivatives_i, derivatives_j);
            noalias(Ju) = C2 * H2 * K - C4 * H * rData.nodal_z[j] * K;
            MathUtils<double>::AddMatrix(rMatrix, Weight*Ju, 3*i, 3*j);
        }
    }
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::AddAuxiliaryLaplacian(
    LocalMatrixType& rVelocityLaplacian,
    LocalMatrixType& rMomentumLaplacian,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX,
    const double Weight)
{
    array_1d<double,3> gradients_vector_i;
    array_1d<double,3> gradients_vector_j;

    for (std::size_t i = 0; i < TNumNodes; ++i)
    {
        gradients_vector_i[0] = rDN_DX(i,0);
        gradients_vector_i[1] = rDN_DX(i,1);
        gradients_vector_i[2] = 0.0;

        for (std::size_t j = 0; j < TNumNodes; ++j)
        {
            gradients_vector_j[0] = rDN_DX(j,0);
            gradients_vector_j[1] = rDN_DX(j,1);
            gradients_vector_j[2] = 0.0;
            const double depth = std::max(0.0, -rData.nodal_z[j]);

            MathUtils<double>::AddMatrix(rVelocityLaplacian, -Weight*outer_prod(gradients_vector_i, gradients_vector_j), 3*i, 3*j);
            MathUtils<double>::AddMatrix(rMomentumLaplacian, -Weight*depth*outer_prod(gradients_vector_i, gradients_vector_j), 3*i, 3*j);
        }
    }
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    auto& r_geometry = this->GetGeometry();

    // Struct to pass around the data
    ElementData data;
    WaveElementType::InitializeData(data, rCurrentProcessInfo);
    GetNodalData(data, r_geometry);

    Vector weights;
    Matrix N_container;
    ShapeFunctionsGradientsType DN_DX_container;
    WaveElementType::CalculateGeometryData(r_geometry, weights, N_container, DN_DX_container);

    // Auxiliary field
    LocalMatrixType laplacian = ZeroMatrix(mLocalSize, mLocalSize);
    LocalMatrixType H_laplacian = ZeroMatrix(mLocalSize, mLocalSize);

    const std::size_t num_gauss_points = weights.size();
    for (std::size_t g = 0; g < num_gauss_points; ++g)
    {
        const double weight = weights[g];
        const array_1d<double,TNumNodes> N = row(N_container, g);
        const BoundedMatrix<double,TNumNodes,2> DN_DX = DN_DX_container[g];

        UpdateGaussPointData(data, N);

        AddAuxiliaryLaplacian(laplacian, H_laplacian, data, N, DN_DX, weight);
    }

    // Calculate the laplacian vector
    const LocalVectorType& vel_vector = this->GetUnknownVector(data);
    LocalVectorType vel_laplacian_vector = prod(laplacian, vel_vector);
    LocalVectorType vel_H_laplacian_vector = prod(H_laplacian, vel_vector);

    // Add the elemental contribution to the nodes
    array_1d<double,3> vel_laplacian = ZeroVector(3);
    array_1d<double,3> vel_H_laplacian = ZeroVector(3);
    for (std::size_t i = 0; i < TNumNodes; ++i)
    {
        std::size_t block = 3 * i;
        vel_laplacian[0] = vel_laplacian_vector[block];
        vel_laplacian[1] = vel_laplacian_vector[block + 1];
        vel_H_laplacian[0] = vel_H_laplacian_vector[block];
        vel_H_laplacian[1] = vel_H_laplacian_vector[block + 1];
        r_geometry[i].SetLock();
        r_geometry[i].FastGetSolutionStepValue(VELOCITY_LAPLACIAN) += vel_laplacian;
        r_geometry[i].FastGetSolutionStepValue(VELOCITY_H_LAPLACIAN) += vel_H_laplacian;
        r_geometry[i].UnSetLock();
    }
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::AddRightHandSide(
    LocalVectorType& rRHS,
    ElementData& rData,
    const Matrix& rNContainer,
    const ShapeFunctionsGradientsType& rDN_DXContainer,
    const Vector& rWeights)
{
    LocalMatrixType lhs = ZeroMatrix(mLocalSize, mLocalSize);
    LocalVectorType dummy; // since the free surface is a primary variable, the bottom artificial viscosity must not be substracted.

    const std::size_t num_gauss_points = rWeights.size();

    for (std::size_t g = 0; g < num_gauss_points; ++g)
    {
        const double weight = rWeights[g];
        const array_1d<double,TNumNodes> N = row(rNContainer, g);
        const BoundedMatrix<double,TNumNodes,2> DN_DX = rDN_DXContainer[g];

        UpdateGaussPointData(rData, N);

        this->AddWaveTerms(lhs, rRHS, rData, N, DN_DX, weight);
        this->AddFrictionTerms(lhs, rRHS, rData, N, DN_DX, weight);
        this->AddDispersiveTerms(rRHS, rData, N, DN_DX, weight);
        this->AddArtificialViscosityTerms(lhs, dummy, rData, N, DN_DX, weight);

        // Deactivating the dry domain
        const double w = WaveElementType::WetFraction(rData);
        lhs *= w;
        rRHS *= w;

        // Controlling the drainage mechanism
        const double factor = 2.0 * norm_2(rData.velocity);
        BoundedMatrix<double,3,3> S = ZeroMatrix(3,3);
        S(0,0) = factor;
        S(1,1) = factor;
        for (std::size_t i = 0; i < TNumNodes; ++i)
        {
            /// Lumped mass matrix
            const double n_ij = 1.0 / static_cast<double>(TNumNodes);
            MathUtils<double>::AddMatrix(lhs, weight*(1.0-w)*n_ij*S, 3*i, 3*i);
        }
    }
    noalias(rRHS) -= prod(lhs, this->GetUnknownVector(rData));
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if(rRightHandSideVector.size() != mLocalSize)
        rRightHandSideVector.resize(mLocalSize, false);

    const auto& r_geometry = this->GetGeometry();

    LocalVectorType f0 = ZeroVector(mLocalSize);
    LocalVectorType f1 = ZeroVector(mLocalSize);
    LocalVectorType f2 = ZeroVector(mLocalSize);
    LocalVectorType f3 = ZeroVector(mLocalSize);

    // Struct to pass around the data
    ElementData data;
    WaveElementType::InitializeData(data, rCurrentProcessInfo);

    Vector weights;
    Matrix N_container;
    ShapeFunctionsGradientsType DN_DX_container;
    WaveElementType::CalculateGeometryData(r_geometry, weights, N_container, DN_DX_container);

    GetNodalData(data, r_geometry, 0);
    AddRightHandSide(f0, data, N_container, DN_DX_container, weights);

    GetNodalData(data, r_geometry, 1);
    AddRightHandSide(f1, data, N_container, DN_DX_container, weights);

    GetNodalData(data, r_geometry, 2);
    AddRightHandSide(f2, data, N_container, DN_DX_container, weights);

    GetNodalData(data, r_geometry, 3);
    AddRightHandSide(f3, data, N_container, DN_DX_container, weights);

    noalias(rRightHandSideVector) = (9*f0 + 19*f1 - 5*f2 + f3) / 24.0;
}


template<std::size_t TNumNodes>
void BoussinesqElement<TNumNodes>::AddExplicitContribution(const ProcessInfo& rCurrentProcessInfo)
{
    auto& r_geometry = this->GetGeometry();

    LocalVectorType f1 = ZeroVector(mLocalSize);
    LocalVectorType f2 = ZeroVector(mLocalSize);
    LocalVectorType f3 = ZeroVector(mLocalSize);

    // Struct to pass around the data
    ElementData data;
    WaveElementType::InitializeData(data, rCurrentProcessInfo);

    Vector weights;
    Matrix N_container;
    ShapeFunctionsGradientsType DN_DX_container;
    WaveElementType::CalculateGeometryData(r_geometry, weights, N_container, DN_DX_container);

    GetNodalData(data, r_geometry, 1);
    AddRightHandSide(f1, data, N_container, DN_DX_container, weights);

    GetNodalData(data, r_geometry, 2);
    AddRightHandSide(f2, data, N_container, DN_DX_container, weights);

    GetNodalData(data, r_geometry, 3);
    AddRightHandSide(f3, data, N_container, DN_DX_container, weights);

    LocalVectorType increment = (23*f1 - 16*f2 + 5*f3) / 12.0;
    array_1d<double,3> nodal_increment;
    for (std::size_t i = 0; i < TNumNodes; ++i)
    {
        std::size_t block = 3*i;
        nodal_increment[0] = increment[block];
        nodal_increment[1] = increment[block + 1];
        nodal_increment[2] = increment[block + 2];

        r_geometry[i].SetLock();
        r_geometry[i].FastGetSolutionStepValue(RHS) += nodal_increment;
        r_geometry[i].UnSetLock();
    }
}


template class BoussinesqElement<3>;
template class BoussinesqElement<4>;

} // namespace Kratos

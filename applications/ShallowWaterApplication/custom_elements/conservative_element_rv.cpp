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
#include "conservative_element_rv.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/phase_function.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

template<std::size_t TNumNodes>
void ConservativeElementRV<TNumNodes>::GetNodalData(ElementData& rData, const GeometryType& rGeometry, int Step)
{
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rData.nodal_h[i] = rGeometry[i].FastGetSolutionStepValue(HEIGHT, Step);
        rData.nodal_z[i] = rGeometry[i].FastGetSolutionStepValue(TOPOGRAPHY, Step);
        rData.nodal_v[i] = rGeometry[i].FastGetSolutionStepValue(VELOCITY, Step);
        rData.nodal_q[i] = rGeometry[i].FastGetSolutionStepValue(MOMENTUM, Step);
        rData.nodal_a[i] = rGeometry[i].FastGetSolutionStepValue(ACCELERATION, Step);
        rData.nodal_w[i] = rGeometry[i].FastGetSolutionStepValue(VERTICAL_VELOCITY, Step);
    }
}

template<std::size_t TNumNodes>
void ConservativeElementRV<TNumNodes>::CalculateArtificialViscosity(
    BoundedMatrix<double,3,3>& rViscosity,
    BoundedMatrix<double,2,2>& rDiffusion,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX)
{
    double art_visc;
    double art_diff;
    ShockCapturingParameters(art_visc, art_diff, rData, rN, rDN_DX);

    // The viscosity previously added by the stabilization
    const double eigenvalue = norm_2(rData.velocity) + std::sqrt(rData.gravity * std::abs(rData.height));
    const double stab_viscosity = this->StabilizationParameter(rData) * std::pow(eigenvalue,2);

    // The orthogonal fourth order tensor
    BoundedMatrix<double,3,3> crosswind_tensor_4;
    CrossWindTensor(crosswind_tensor_4, rData.velocity);
    crosswind_tensor_4 *= art_visc;

    // The streamline fourth order tensor
    BoundedMatrix<double,3,3> streamline_tensor_4;
    StreamLineTensor(streamline_tensor_4, rData.velocity);
    streamline_tensor_4 *= std::max(0.0, art_visc - stab_viscosity);

    // The constitutive tensor
    BoundedMatrix<double,3,3> constitutive_tensor = IdentityMatrix(3,3);
    array_1d<double,3> m({1.0, 1.0, 0.0});
    constitutive_tensor -= outer_prod(m, m) / 3.0;
    rViscosity = prod(constitutive_tensor, crosswind_tensor_4 + streamline_tensor_4);

    // The diffusion previously added by the stabilization
    const double stab_diffusivity = this->StabilizationParameter(rData) * std::pow(eigenvalue,2);

    // The second order crosswind tensor
    BoundedMatrix<double,2,2> crosswind_tensor_2;
    CrossWindTensor(crosswind_tensor_2, rData.velocity);
    crosswind_tensor_2 *= art_diff;

    // The second order streamline tensor
    BoundedMatrix<double,2,2> streamline_tensor_2;
    StreamLineTensor(streamline_tensor_2, rData.velocity);
    streamline_tensor_2 *= std::max(0.0, art_diff - stab_diffusivity);

    // The constitutive matrix
    rDiffusion = crosswind_tensor_2 + streamline_tensor_2;
}

template<std::size_t TNumNodes>
void ConservativeElementRV<TNumNodes>::ShockCapturingParameters(
    double& rArtViscosity,
    double& rArtDiffusion,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    // Computation of the residuals and gradients
    array_1d<double,3> flow_residual;
    double height_residual;
    BoundedMatrix<double,3,3> flow_grad;
    array_1d<double,3> height_grad;
    AlgebraicResidual(flow_residual, height_residual, flow_grad, height_grad, rData, rN, rDN_DX);

    // slope limits
    const double min_slope = 0.1;
    const double max_slope = 1.0;
    const double eigenvalue = norm_2(rData.velocity) + std::sqrt(rData.gravity * std::abs(rData.height));
    const double min_q_slope = std::max(eigenvalue * 1.0, min_slope);
    const double max_q_slope = std::max(eigenvalue * 10.0, max_slope);

    // Final assembly of the parameters
    const double q_residual_norm = norm_2(flow_residual);
    const double q_grad_frobenius = norm_frobenius(flow_grad);
    const double q_gradient_norm = std::min(std::max(q_grad_frobenius, min_q_slope), max_q_slope);
    rArtViscosity = 0.5 * rData.shock_stab_factor * rData.length * q_residual_norm / q_gradient_norm;

    const double h_residual_norm = std::abs(height_residual);
    const double h_gradient_norm = std::min(std::max(norm_2(height_grad), min_slope), max_slope);
    rArtDiffusion = 0.5 * rData.shock_stab_factor * rData.length * h_residual_norm / h_gradient_norm;
}

template<std::size_t TNumNodes>
void ConservativeElementRV<TNumNodes>::AlgebraicResidual(
    array_1d<double,3>& rFlowResidual,
    double& rHeightResidual,
    BoundedMatrix<double,3,3>& rFlowGrad,
    array_1d<double,3>& rHeightGrad,
    const ElementData& rData,
    const array_1d<double,TNumNodes>& rN,
    const BoundedMatrix<double,3,2>& rDN_DX)
{
    const array_1d<double,3> flow_rate = WaveElementType::VectorProduct(rData.nodal_q, rN);
    const array_1d<double,3> flow_acc = WaveElementType::VectorProduct(rData.nodal_a, rN);
    const double height_acc = inner_prod(rData.nodal_w, rN);

    rHeightGrad = WaveElementType::ScalarGradient(rData.nodal_h, rDN_DX);
    rFlowGrad = WaveElementType::VectorGradient(rData.nodal_q, rDN_DX);
    const double flow_div = WaveElementType::VectorDivergence(rData.nodal_q, rDN_DX);
    const array_1d<double,3> topography_grad = WaveElementType::ScalarGradient(rData.nodal_z, rDN_DX);

    const double c2 = rData.gravity * rData.height;
    const array_1d<double,3> friction = rData.gravity * rData.height * rData.p_bottom_friction->CalculateRHS(rData.height, rData.velocity);
    array_1d<double,3> flux = ZeroVector(3);
    for (size_t i = 0; i < TNumNodes; ++i) {
        BoundedMatrix<double,3,3> aux = outer_prod(rData.nodal_q[i], rData.nodal_v[i]);
        flux[0] += aux(0,0) * rDN_DX(i,0);
        flux[0] += aux(0,1) * rDN_DX(i,1);
        flux[1] += aux(1,0) * rDN_DX(i,0);
        flux[1] += aux(1,1) * rDN_DX(i,1);
    }
    BoundedMatrix<double,3,3> art_s = ZeroMatrix(3,3);
    this->CalculateArtificialDamping(art_s, rData);
    array_1d<double,3> unknowns = flow_rate;
    unknowns[2] = rData.height;
    const array_1d<double,3> damping = prod(art_s, unknowns);

    rFlowResidual = flow_acc + flux + c2 * (rHeightGrad + topography_grad) + friction + damping; // + rData.height * mesh_acc;
    rHeightResidual = height_acc + flow_div;
}

template<std::size_t TNumNodes>
void ConservativeElementRV<TNumNodes>::StreamLineTensor(BoundedMatrix<double,2,2>& rTensor, const array_1d<double,3>& rVector)
{
    const double e = std::numeric_limits<double>::epsilon(); // small value to avoid division by zero
    array_1d<double,2> aux_vector;
    aux_vector[0] = rVector[0];
    aux_vector[1] = rVector[1];
    rTensor = outer_prod(aux_vector, aux_vector) / (inner_prod(rVector, rVector) + e);
}

template<std::size_t TNumNodes>
void ConservativeElementRV<TNumNodes>::CrossWindTensor(BoundedMatrix<double,2,2>& rTensor, const array_1d<double,3>& rVector)
{
    StreamLineTensor(rTensor, rVector);
    rTensor = IdentityMatrix(2) - rTensor;
}

template<std::size_t TNumNodes>
void ConservativeElementRV<TNumNodes>::StreamLineTensor(BoundedMatrix<double,3,3>& rTensor, const array_1d<double,3>& rVector)
{
    BoundedMatrix<double,2,2> stream_line_tensor;
    StreamLineTensor(stream_line_tensor, rVector);
    rTensor(0,0) = stream_line_tensor(0,0);
    rTensor(1,1) = stream_line_tensor(1,1);
    rTensor(0,1) = stream_line_tensor(0,1);
    rTensor(1,0) = stream_line_tensor(0,1);
    rTensor(2,2) = stream_line_tensor(0,1);
    rTensor(0,2) = 0.0;
    rTensor(1,2) = 0.0;
    rTensor(2,0) = 0.0;
    rTensor(2,1) = 0.0;
}

template<std::size_t TNumNodes>
void ConservativeElementRV<TNumNodes>::CrossWindTensor(BoundedMatrix<double,3,3>& rTensor, const array_1d<double,3>& rVector)
{
    StreamLineTensor(rTensor, rVector);
    rTensor = IdentityMatrix(3) - rTensor;
}

template class ConservativeElementRV<3>;

} // namespace Kratos

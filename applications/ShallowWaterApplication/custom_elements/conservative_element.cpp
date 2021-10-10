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
#include "conservative_element.h"
#include "utilities/geometry_utilities.h"
#include "shallow_water_application_variables.h"
#include "custom_utilities/shallow_water_utilities.h"

namespace Kratos
{

template<std::size_t TNumNodes>
const Variable<double>& ConservativeElement<TNumNodes>::GetUnknownComponent(int Index) const
{
    switch (Index) {
        case 0: return MOMENTUM_X;
        case 1: return MOMENTUM_Y;
        case 2: return HEIGHT;
        default: KRATOS_ERROR << "WaveElement::GetUnknownComponent index out of bounds." << std::endl;
    }
}

template<std::size_t TNumNodes>
typename ConservativeElement<TNumNodes>::LocalVectorType ConservativeElement<TNumNodes>::GetUnknownVector(ElementData& rData)
{
    std::size_t index = 0;
    array_1d<double,mLocalSize> unknown;
    for (std::size_t i = 0; i < TNumNodes; ++i) {
        unknown[index++] = rData.nodal_q[i][0];
        unknown[index++] = rData.nodal_q[i][1];
        unknown[index++] = rData.nodal_h[i];
    }
    return unknown;
}

template<std::size_t TNumNodes>
void ConservativeElement<TNumNodes>::CalculateGaussPointData(ElementData& rData, const array_1d<double,TNumNodes>& rN)
{
    const double h = inner_prod(rData.nodal_h, rN);
    const double c2 = rData.gravity * h;
    const array_1d<double,3> v = WaveElementType::VectorProduct(rData.nodal_v, rN);

    rData.height = h;
    rData.velocity = v;

    /**
     * A_1 = {{2 * u_1   0   -u_1^2 + gh},
     *        {  u_2    u_1   -u_1 * u_2},
     *        {   1      0        0     }}
     */
    rData.A1(0,0) = 2*v[0];
    rData.A1(0,1) = 0;
    rData.A1(0,2) = -v[0]*v[0] + c2;
    rData.A1(1,0) = v[1];
    rData.A1(1,1) = v[0];
    rData.A1(1,2) = -v[0]*v[1];
    rData.A1(2,0) = 1;
    rData.A1(2,1) = 0;
    rData.A1(2,2) = 0;

    /*
     * A_2 = {{u_2    u_1      -u_1 * u_2},
     *        { 0   2 * u_2   -u_2^2 + gh},
     *        { 0      1            0    }}
     */
    rData.A2(0,0) = v[1];
    rData.A2(0,1) = v[0];
    rData.A2(0,2) = -v[0]*v[1];
    rData.A2(1,0) = 0;
    rData.A2(1,1) = 2*v[1];
    rData.A2(1,2) = -v[1]*v[1] + c2;
    rData.A2(2,0) = 0;
    rData.A2(2,1) = 1;
    rData.A2(2,2) = 0;

    /// b_1
    rData.b1[0] = c2;
    rData.b1[1] = 0;
    rData.b1[2] = 0;

    /// b_2
    rData.b2[0] = 0;
    rData.b2[1] = c2;
    rData.b2[2] = 0;
}

template<std::size_t TNumNodes>
double ConservativeElement<TNumNodes>::StabilizationParameter(const ElementData& rData) const
{
    const double lambda = std::sqrt(rData.gravity * std::abs(rData.height)) + norm_2(rData.velocity);
    const double epsilon = 1e-6;
    const double threshold = rData.relative_dry_height * rData.length;
    const double w = ShallowWaterUtilities().WetFraction(rData.height, threshold);
    return w * rData.length * rData.stab_factor / (lambda + epsilon);
}

template<std::size_t TNumNodes>
void ConservativeElement<TNumNodes>::CalculateArtificialViscosity(
    BoundedMatrix<double,3,3>& rViscosity,
    BoundedMatrix<double,2,2>& rDiffusion,
    const ElementData& rData,
    const BoundedMatrix<double,TNumNodes,2>& rDN_DX)
{
    double jump = 0;
    array_1d<double,2> inner_grad_h = prod(rData.nodal_h, rDN_DX);
    for (auto& r_elem : this->GetValue(NEIGHBOUR_ELEMENTS)) {
        array_1d<double,2> outer_grad_h;
        array_1d<double,2> normal;
        CalculateGradient(outer_grad_h, r_elem.GetGeometry());
        CalculateEdgeUnitNormal(normal, r_elem.GetGeometry());

        const double gj_h = norm_2(inner_grad_h - outer_grad_h);
        const double gm_h = std::abs(inner_prod(inner_grad_h, normal)) + std::abs(inner_prod(outer_grad_h, normal)) + 1e-16;
        const double correction = std::abs(inner_prod(inner_grad_h, normal)) / (norm_2(inner_grad_h) + 1e-16);

        jump = std::max(jump, correction * gj_h / gm_h);
    }
    const double lambda = std::sqrt(rData.gravity * std::abs(rData.height)) + norm_2(rData.velocity);
    const double visc = rData.shock_stab_factor * rData.length * lambda * jump;
    rViscosity = visc * IdentityMatrix(3);
    rDiffusion = visc * IdentityMatrix(2);
}

template<std::size_t TNumNodes>
void ConservativeElement<TNumNodes>::CalculateArtificialDamping(
    BoundedMatrix<double,3,3>& rDamping,
    const ElementData& rData)
{
    double factor = 1e3 / rData.length;
    double threshold = rData.relative_dry_height * rData.length;
    factor *= 1.0 - ShallowWaterUtilities().WetFraction(rData.height, threshold);
    rDamping(0,0) = factor;
    rDamping(1,1) = factor;
}

template<std::size_t TNumNodes>
void ConservativeElement<TNumNodes>::CalculateGradient(
    array_1d<double,2>& rGradient,
    const GeometryType& rGeometry)
{
    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(rGeometry, DN_DX, N, area);
    array_1d<double,3> nodal_h;
    std::size_t i = 0;
    for (auto& r_node : rGeometry) {
        nodal_h[i++] = r_node.FastGetSolutionStepValue(HEIGHT);
    }
    rGradient = prod(nodal_h, DN_DX);
}

template<std::size_t TNumNodes>
void ConservativeElement<TNumNodes>::CalculateEdgeUnitNormal(
    array_1d<double,2>& rNormal,
    const GeometryType& rGeometry)
{
    array_1d<double,3> normal;
    normal = rGeometry.Center();
    normal -= this->GetGeometry().Center();
    normal /= norm_2(normal) + 1e-16;
    rNormal[0] = normal[0];
    rNormal[1] = normal[1];
}

template class ConservativeElement<3>;

} // namespace Kratos

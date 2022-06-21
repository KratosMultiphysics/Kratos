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
//                   Ruben Zorrilla
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "embedded_primitive_element.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/phase_function.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

template<std::size_t TNumNodes>
void EmbeddedPrimitiveElement<TNumNodes>::UpdateGaussPointData(ElementData& rData, const array_1d<double,TNumNodes>& rN)
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
double EmbeddedPrimitiveElement<TNumNodes>::StabilizationParameter(const ElementData& rData) const
{
    const double lambda = std::sqrt(rData.gravity * std::abs(rData.height)) + norm_2(rData.velocity);
    const double epsilon = 1e-6;
    const double threshold = rData.relative_dry_height * rData.length;
    const double w = PhaseFunction::WetFraction(rData.height, threshold);
    return w * rData.length * rData.stab_factor / (lambda + epsilon);
}

template<std::size_t TNumNodes>
void EmbeddedPrimitiveElement<TNumNodes>::CalculateGeometryData(
    const GeometryType& rGeometry,
    Vector &rGaussWeights,
    Matrix &rNContainer,
    ShapeFunctionsGradientsType &rDN_DX)
{

}

template class EmbeddedPrimitiveElement<3>;

} // namespace Kratos

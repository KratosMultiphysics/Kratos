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
#include "primitive_element.h"
#include "utilities/geometry_utilities.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

template<std::size_t TNumNodes>
void PrimitiveElement<TNumNodes>::UpdateGaussPointData(
    ElementData& rData,
    const array_1d<double,TNumNodes>& rN)
{
    const double h = inner_prod(rData.nodal_h, rN);
    const array_1d<double,3> v = BaseType::VectorProduct(rData.nodal_v, rN);

    rData.height = h;
    rData.velocity = v;

    /**
     * A_1 = {{ u_1   0   g },
     *        { 0   u_1   0 },
     *        { h   0   u_1 }}
     */
    rData.A1 = ZeroMatrix(3, 3);
    rData.A1(0,2) = rData.gravity;
    rData.A1(2,0) = h;
    rData.A1(0,0) = v[0];
    rData.A1(1,1) = v[0];
    rData.A1(2,2) = v[0];

    /*
     * A_2 = {{ u_2   0   0 },
     *        { 0   u_2   g },
     *        { 0   h   u_2 }}
     */
    rData.A2 = ZeroMatrix(3, 3);
    rData.A2(1,2) = rData.gravity;
    rData.A2(2,1) = h;
    rData.A2(0,0) = v[1];
    rData.A2(1,1) = v[1];
    rData.A2(2,2) = v[1];

    /// b_1
    rData.b1 = ZeroVector(3);
    rData.b1[0] = rData.gravity;

    /// b_2
    rData.b2 = ZeroVector(3);
    rData.b2[1] = rData.gravity;
}

template class PrimitiveElement<3>;
template class PrimitiveElement<4>;

} // namespace Kratos

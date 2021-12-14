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
#include "boussinesq_condition.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

template<std::size_t TNumNodes>
const Variable<double>& BoussinesqCondition<TNumNodes>::GetUnknownComponent(int Index) const
{
    switch (Index) {
        case 0: return VELOCITY_X;
        case 1: return VELOCITY_Y;
        case 2: return FREE_SURFACE_ELEVATION;
        default: KRATOS_ERROR << "BoussinesqCondition::GetUnknownComponent index out of bounds." << std::endl;
    }
}

template<std::size_t TNumNodes>
typename BoussinesqCondition<TNumNodes>::LocalVectorType BoussinesqCondition<TNumNodes>::GetUnknownVector(ConditionData& rData)
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
void BoussinesqCondition<TNumNodes>::CalculateGaussPointData(
    ConditionData& rData,
    const IndexType PointIndex,
    const array_1d<double,TNumNodes>& rN)
{
    const double eta = inner_prod(rData.nodal_f, rN);
    const double H = -inner_prod(rData.nodal_z, rN);
    const double g = rData.gravity;
    const double e = rData.amplitude / H;  // the non linearity ratio
    const array_1d<double,3> v = WaveConditionType::VectorProduct(rData.nodal_v, rN);

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

    // Calculate the normal vector
    auto integration_point = this->GetGeometry().IntegrationPoints()[PointIndex];
    rData.normal = this->GetGeometry().Normal(integration_point);
    rData.normal /= norm_2(rData.normal);
}

template class BoussinesqCondition<2>;

} // namespace Kratos

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
#include "conservative_condition.h"
#include "includes/checks.h"
#include "shallow_water_application_variables.h"
#include "custom_utilities/shallow_water_utilities.h"

namespace Kratos
{

template<std::size_t TNumNodes>
const Variable<double>& ConservativeCondition<TNumNodes>::GetUnknownComponent(int Index) const
{
    switch (Index) {
        case 0: return MOMENTUM_X;
        case 1: return MOMENTUM_Y;
        case 2: return HEIGHT;
        default: KRATOS_ERROR << "WaveCondition::GetUnknownComponent index out of bounds." << std::endl;
    }
}

template<std::size_t TNumNodes>
typename ConservativeCondition<TNumNodes>::LocalVectorType ConservativeCondition<TNumNodes>::GetUnknownVector(ConditionData& rData)
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
void ConservativeCondition<TNumNodes>::CalculateGaussPointData(
    ConditionData& rData,
    const IndexType PointIndex,
    const array_1d<double,TNumNodes>& rN)
{
    const double h = inner_prod(rData.nodal_h, rN);
    const double c2 = rData.gravity * h;
    const array_1d<double,3> v = WaveConditionType::VectorProduct(rData.nodal_v, rN);

    rData.height = h;
    rData.velocity = v;

    rData.A1(0,0) = 2*v[0];
    rData.A1(0,1) = 0;
    rData.A1(0,2) = -v[0]*v[0] + c2;
    rData.A1(1,0) = v[1];
    rData.A1(1,1) = v[0];
    rData.A1(1,2) = -v[0]*v[1];
    rData.A1(2,0) = 1;
    rData.A1(2,1) = 0;
    rData.A1(2,2) = 0;

    rData.A2(0,0) = v[1];
    rData.A2(0,1) = v[0];
    rData.A2(0,2) = -v[0]*v[1];
    rData.A2(1,0) = 0;
    rData.A2(1,1) = 2*v[1];
    rData.A2(1,2) = -v[1]*v[1] + c2;
    rData.A2(2,0) = 0;
    rData.A2(2,1) = 1;
    rData.A2(2,2) = 0;

    rData.b1[0] = c2;
    rData.b1[1] = 0;
    rData.b1[2] = 0;

    rData.b2[0] = 0;
    rData.b2[1] = c2;
    rData.b2[2] = 0;

    auto integration_point = this->GetGeometry().IntegrationPoints()[PointIndex];
    rData.normal = this->GetGeometry().Normal(integration_point);
    rData.normal /= norm_2(rData.normal);
}

template class ConservativeCondition<2>;

} // namespace Kratos

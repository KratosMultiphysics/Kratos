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
#include "includes/kratos_flags.h"
#include "includes/checks.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

template<std::size_t TNumNodes>
const Parameters ConservativeCondition<TNumNodes>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "required_variables"         : ["MOMENTUM","VELOCITY","HEIGHT","TOPOGRAPHY","ACCELERATION","VERTICAL_VELOCITY"],
        "required_dofs"              : ["MOMENTUM_X","MOMENTUM_Y","HEIGHT"],
        "compatible_geometries"      : ["Triangle2D3"],
        "element_integrates_in_time" : false
    })");
    return specifications;
}

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
    const double z = inner_prod(rData.nodal_z, rN);
    const array_1d<double,3> v = BaseType::VectorProduct(rData.nodal_v, rN);
    const bool is_supercritical = (norm_2(v) >= std::sqrt(rData.gravity*h));
    auto integration_point = this->GetGeometry().IntegrationPoints()[PointIndex];
    rData.normal = this->GetGeometry().UnitNormal(integration_point);

    rData.height = h;
    rData.velocity = v;

    if (this->Is(SLIP)) {
        rData.v_neumann = 0.0;
        rData.h_dirichlet = h;
    }
    else if (this->Is(INLET)) {
        rData.v_neumann = inner_prod(rData.normal, this->GetValue(VELOCITY));
        rData.h_dirichlet = (is_supercritical) ? this->GetValue(HEIGHT) : h;
    }
    else if (this->Is(OUTLET)) {
        rData.v_neumann = inner_prod(rData.normal, v);
        rData.h_dirichlet = (is_supercritical) ? h : this->GetValue(HEIGHT);
    }
    else {
        rData.v_neumann = inner_prod(rData.normal, v);
        rData.h_dirichlet = h;
    }

    /// Boundary traces
    const double vel_trace = rData.v_neumann * h;
    double press_trace = rData.gravity * std::pow(rData.h_dirichlet + z, 2);

    /// Convective flux
    array_1d<double,3> conv_flux = v;
    conv_flux[2] = 1;
    conv_flux *= vel_trace;

    /// Pressure flux
    array_1d<double,3> press_flux = rData.normal;
    press_flux[2] = 0.0;
    press_flux *= press_trace;

    /// Assembly of the flux
    rData.flux = conv_flux + press_flux;
}

template class ConservativeCondition<2>;

} // namespace Kratos

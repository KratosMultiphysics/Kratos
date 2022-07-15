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
#include "includes/kratos_flags.h"


// Application includes
#include "primitive_condition.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{

template<std::size_t TNumNodes>
void PrimitiveCondition<TNumNodes>::CalculateGaussPointData(
    ConditionData& rData,
    const IndexType PointIndex,
    const array_1d<double,TNumNodes>& rN)
{
    const double h = inner_prod(rData.nodal_h, rN);
    const double z = inner_prod(rData.nodal_z, rN);
    const array_1d<double,3> v = BaseType::VectorProduct(rData.nodal_v, rN);
    const bool is_supercritical = (norm_2(v) >= std::sqrt(rData.gravity*h));
    const auto integration_point = this->GetGeometry().IntegrationPoints()[PointIndex];
    rData.normal = this->GetGeometry().UnitNormal(integration_point);

    rData.height = h;
    rData.velocity = v;

    rData.A1 = ZeroMatrix(3, 3);
    rData.A2 = ZeroMatrix(3, 3);
    rData.b1 = ZeroVector(3);
    rData.b2 = ZeroVector(3);

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
    const double press_trace = rData.gravity * (rData.h_dirichlet + z);
    const double vel_trace = rData.v_neumann;

    /// Convective flux
    array_1d<double,3> conv_flux = v;
    conv_flux[2] = h;
    conv_flux *= vel_trace;

    /// Pressure flux
    array_1d<double,3> press_flux = rData.normal;
    press_flux[2] = 0.0;
    press_flux *= press_trace;

    /// Assembly of the flux
    rData.flux = conv_flux + press_flux;
}


template<std::size_t TNumNodes>
void PrimitiveCondition<TNumNodes>::AddFluxTerms(
    LocalVectorType& rVector,
    const ConditionData& rData,
    const array_1d<double,TNumNodes>& rN,
    const double Weight)
{
    const double penalty_factor = 1.0 * rData.length;
    const auto n = rData.normal;

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        double n_i;
        if (rData.integrate_by_parts) {
            n_i = rN[i];
        } else {
            n_i = 0.0;
        }

        /// Flux contribution
        MathUtils<double>::AddVector(rVector, -Weight * n_i * rData.flux, 3*i);

        /// Penalty stabilization
        double i_velocity_penalty = inner_prod(n, rData.nodal_v[i]) - rData.v_neumann;
        rVector[3*i    ] -= n[0] * Weight * penalty_factor * i_velocity_penalty;
        rVector[3*i + 1] -= n[1] * Weight * penalty_factor * i_velocity_penalty;
        rVector[3*i + 2] -=        Weight * penalty_factor * (rData.nodal_h[i] - rData.h_dirichlet);
    }
}


template class PrimitiveCondition<2>;

} // namespace Kratos

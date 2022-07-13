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
    const array_1d<double,3> v = BaseType::VectorProduct(rData.nodal_v, rN);

    rData.height = h;
    rData.velocity = v;

    rData.A1 = ZeroMatrix(3, 3);
    // rData.A1(0,2) = rData.gravity;
    // rData.A1(2,0) = h;
    // rData.A1(0,0) = v[0];
    // rData.A1(1,1) = v[0];
    // rData.A1(2,2) = v[0];

    rData.A2 = ZeroMatrix(3, 3);
    // rData.A2(1,2) = rData.gravity;
    // rData.A2(2,1) = h;
    // rData.A2(0,0) = v[1];
    // rData.A2(1,1) = v[1];
    // rData.A2(2,2) = v[1];

    rData.b1 = ZeroVector(3);
    // rData.b1[0] = rData.gravity;
    // rData.b1[1] = 0.0;
    // rData.b1[2] = 0.0;

    rData.b2 = ZeroVector(3);
    // rData.b2[0] = 0.0;
    // rData.b2[1] = rData.gravity;
    // rData.b2[2] = 0.0;

    auto integration_point = this->GetGeometry().IntegrationPoints()[PointIndex];
    rData.normal = this->GetGeometry().UnitNormal(integration_point);
}


template<std::size_t TNumNodes>
void PrimitiveCondition<TNumNodes>::AddFluxTerms(
    LocalVectorType& rVector,
    const ConditionData& rData,
    const array_1d<double,TNumNodes>& rN,
    const double Weight)
{
    const auto h = inner_prod(rData.nodal_h, rN);
    const auto z = inner_prod(rData.nodal_z, rN);
    const array_1d<double,3> v = BaseType::VectorProduct(rData.nodal_v, rN);
    const auto g = rData.gravity;
    const auto n = rData.normal;
    double vel_trace;
    double press_trace;
    double vel_penalty;
    double press_penalty;
    const double penalty_factor = 1e-2 * rData.length;
    const bool is_supercritical = (norm_2(v) >= std::sqrt(g*h));

    if (this->Is(SLIP)) {
        press_trace = g * (h + z);
        vel_trace = 0.0;
        press_penalty = 0.0;
        vel_penalty = inner_prod(n, v);
    }
    else if (this->Is(INLET)) {
        const double v_neumann = inner_prod(n, this->GetValue(VELOCITY));
        const double h_inlet = (is_supercritical) ? this->GetValue(HEIGHT) : h;
        press_trace = g * (h_inlet + z);
        vel_trace = v_neumann;
        press_penalty = h - h_inlet;
        vel_penalty = inner_prod(n, v) - v_neumann;
    }
    else if (this->Is(OUTLET)) {
        const double h_dirichlet = (is_supercritical) ? h : this->GetValue(HEIGHT);
        press_trace = g * (h_dirichlet + z);
        vel_trace = inner_prod(v, rN);
        press_penalty = h - h_dirichlet;
        vel_penalty = 0.0;
    }
    else {
        press_trace = g * (h + z);
        vel_trace = inner_prod(v, rN);
        press_penalty = 0.0;
        vel_penalty = 0.0;
    }

    /// Convective flux
    array_1d<double,3> conv_flux = v;
    conv_flux[2] = h;
    conv_flux *= vel_trace;

    /// Pressure flux
    array_1d<double,3> press_flux = rData.normal;
    press_flux[2] = 0.0;
    press_flux *= press_trace;

    /// Assembly of the flux
    const array_1d<double,3> flux = conv_flux + press_flux;

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        double n_i;
        if (rData.integrate_by_parts) {
            n_i = rN[i];
        } else {
            n_i = 0.0;
        }

        /// Flux contribution
        MathUtils<double>::AddVector(rVector, -Weight * n_i * flux, 3*i);

        /// Penalty stabilization
        rVector[3*i    ] -= n[0] * Weight * penalty_factor * vel_penalty;
        rVector[3*i + 1] -= n[1] * Weight * penalty_factor * vel_penalty;
        rVector[3*i + 2] -=        Weight * penalty_factor * press_penalty;
    }
}


template class PrimitiveCondition<2>;

} // namespace Kratos

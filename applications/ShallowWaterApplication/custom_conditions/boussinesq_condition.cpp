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
    rData.normal = this->GetGeometry().UnitNormal(integration_point);
}

template<std::size_t TNumNodes>
void BoussinesqCondition<TNumNodes>::AddAuxiliaryLaplacian(
    LocalVectorType& rLaplacian,
    const GeometryType& rParentGeometry,
    const ConditionData& rData,
    const array_1d<double,TNumNodes>& rN,
    const Matrix& rDN_DX,
    const double Weight)
{
    array_1d<double,TNumNodes> normal_i;
    double divergence_j;
    std::size_t elem_num_nodes = rParentGeometry.size();
    std::vector<array_1d<double,3>> nodal_v(elem_num_nodes);

    for (IndexType i = 0; i < elem_num_nodes; ++i) {
        nodal_v[i] = rParentGeometry[i].FastGetSolutionStepValue(VELOCITY);
    }

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        normal_i[0] = rData.normal[0] * rN[i];
        normal_i[1] = rData.normal[1] * rN[i];
        normal_i[2] = 0.0;

        for (IndexType j = 0; j < elem_num_nodes; ++j)
        {
            divergence_j  = rDN_DX(j,0) * nodal_v[j][0];
            divergence_j += rDN_DX(j,1) * nodal_v[j][1];

            MathUtils<double>::AddVector(rLaplacian, Weight*normal_i*divergence_j, 3*i);
        }
    }
}

template<std::size_t TNumNodes>
void BoussinesqCondition<TNumNodes>::CalculateShapeFunctionDerivaties(
    Matrix& rDN_DX,
    const GeometryType& rParentGeometry,
    const Point& rPoint)
{
    std::size_t num_nodes = rParentGeometry.size();
    Point local_coordinates;
    std::vector<array_1d<double,3>> global_derivatives;
    rParentGeometry.PointLocalCoordinates(local_coordinates, rPoint);
    rParentGeometry.GlobalSpaceDerivatives(global_derivatives, local_coordinates, 1);
    rDN_DX.resize(num_nodes, 2);
    for (std::size_t i = 0; i < num_nodes; ++i)
    {
        rDN_DX(i,0) = global_derivatives[1](i);
        rDN_DX(i,1) = global_derivatives[2](i);
    }
}

template<std::size_t TNumNodes>
void BoussinesqCondition<TNumNodes>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    auto& r_geom = this->GetGeometry();

    // Struct to pass around the data
    ConditionData data;
    WaveConditionType::InitializeData(data, rCurrentProcessInfo);

    // Geometrical data
    auto& parent_geom = this->GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry();
    Matrix DN_DX;           // Gradients of the parent element
    Vector weights;         // Line integration data
    Matrix N_container;     // Line integration data
    this->CalculateGeometryData(weights, N_container);
    const std::size_t num_gauss_points = weights.size();
    const auto& g_points = r_geom.IntegrationPoints();

    // Boundary term for the auxiliary fields
    LocalVectorType acc_laplacian_vector = ZeroVector(mLocalSize);
    LocalVectorType vel_laplacian_vector = ZeroVector(mLocalSize);

    // Gauss point contribution
    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const double weight = weights[g];
        const auto point = g_points[g];
        const array_1d<double,TNumNodes> N = row(N_container, g);

        CalculateGaussPointData(data, g, N);
        CalculateShapeFunctionDerivaties(DN_DX, parent_geom, point);
        AddAuxiliaryLaplacian(vel_laplacian_vector, parent_geom, data, N, DN_DX, weight);
    }

    array_1d<double,3> vel_laplacian = ZeroVector(3);
    array_1d<double,3> acc_laplacian = ZeroVector(3);
    for (std::size_t i = 0; i < TNumNodes; ++i)
    {
        std::size_t block = 3 * i;
        vel_laplacian[0] = vel_laplacian_vector[block];
        vel_laplacian[1] = vel_laplacian_vector[block + 1];
        r_geom[i].SetLock();
        r_geom[i].FastGetSolutionStepValue(VELOCITY_LAPLACIAN) += vel_laplacian;
        r_geom[i].UnSetLock();
    }
}

template class BoussinesqCondition<2>;

} // namespace Kratos

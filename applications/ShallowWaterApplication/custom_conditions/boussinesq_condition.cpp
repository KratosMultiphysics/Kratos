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
#include "utilities/geometry_utilities.h"

namespace Kratos
{

template<std::size_t TNumNodes>
const Parameters BoussinesqCondition<TNumNodes>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "required_variables"         : ["VELOCITY","FREE_SURFACE_ELEVATION","TOPOGRAPHY","ACCELERATION","VERTICAL_VELOCITY","VELOCITY_LAPLACIAN","VELOCITY_H_LAPLACIAN"],
        "required_dofs"              : ["VELOCITY_X","VELOCITY_Y","FREE_SURFACE_ELEVATION"],
        "compatible_geometries"      : ["Line2D4"],
        "element_integrates_in_time" : false
    })");
    return specifications;
}

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
    const array_1d<double,3> v = WaveConditionType::VectorProduct(rData.nodal_v, rN);

    rData.depth = H;
    rData.height = H + eta;
    rData.velocity = v;

    /**
     * A_1 = {{  u_1     0    g  },
     *        {   0     u_1   0  },
     *        {H + eta   0   u_1}}
     */
    rData.A1(0,0) = v[0];
    rData.A1(0,1) = 0;
    rData.A1(0,2) = g;
    rData.A1(1,0) = 0;
    rData.A1(1,1) = v[0];
    rData.A1(1,2) = 0;
    rData.A1(2,0) = H + eta;
    rData.A1(2,1) = 0;
    rData.A1(2,2) = v[0];

    /*
     * A_2 = {{ u_2     0      0  },
     *        {  0     u_2     g  },
     *        {  0   H + eta  u_2}}
     */
    rData.A2(0,0) = v[1];
    rData.A2(0,1) = 0;
    rData.A2(0,2) = 0;
    rData.A2(1,0) = 0;
    rData.A2(1,1) = v[1];
    rData.A2(1,2) = g;
    rData.A2(2,0) = 0;
    rData.A2(2,1) = H + eta;
    rData.A2(2,2) = v[1];

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
void BoussinesqCondition<TNumNodes>::AddLaplacianBoundary(
    LocalVectorType& rVelocityLaplacian,
    LocalVectorType& rMomentumLaplacian,
    const GeometryType& rParentGeometry,
    const ConditionData& rData,
    const array_1d<double,TNumNodes>& rN,
    const Matrix& rDN_DX,
    const double Weight)
{
    array_1d<double,3> normal_i;
    std::size_t elem_num_nodes = rParentGeometry.size();

    double divergence_u = 0.0;
    double divergence_Hu = 0.0;
    for (IndexType i = 0; i < elem_num_nodes; ++i) {
        const array_1d<double,3> velocity = rParentGeometry[i].FastGetSolutionStepValue(VELOCITY);
        const double depth = -rParentGeometry[i].FastGetSolutionStepValue(TOPOGRAPHY);
        divergence_u  +=  rDN_DX(i,0) * velocity[0] + rDN_DX(i,1) * velocity[1];
        divergence_Hu += (rDN_DX(i,0) * velocity[0] + rDN_DX(i,1) * velocity[1]) * depth;
    }

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        normal_i = rData.normal * rN[i];
        MathUtils<double>::AddVector(rVelocityLaplacian, Weight*normal_i*divergence_u, 3*i);
        MathUtils<double>::AddVector(rMomentumLaplacian, Weight*normal_i*divergence_Hu, 3*i);
    }
}


template<std::size_t TNumNodes>
void BoussinesqCondition<TNumNodes>::AddMomentumDispersionTerms(
    LocalVectorType& rLaplacianBoundary,
    const GeometryType& rParentGeometry,
    const ConditionData& rData,
    const array_1d<double,TNumNodes>& rN,
    const Matrix& rDN_DX,
    const double Weight)
{
    // Constants
    const double beta = -0.531;
    const double C2 = 0.5 * std::pow(beta, 2);
    const double C4 = beta;
    const double H = rData.depth;
    const double H2 = std::pow(H, 2);

    double Ju = 0.0;
    std::size_t elem_num_nodes = rParentGeometry.size();
    for (IndexType i = 0; i < elem_num_nodes; ++i) {
        const array_1d<double,3> acceleration = rParentGeometry[i].FastGetSolutionStepValue(ACCELERATION);
        const double depth  = -rParentGeometry[i].FastGetSolutionStepValue(TOPOGRAPHY);
        const double div_a  =  rDN_DX(i,0) * acceleration[0] + rDN_DX(i,1) * acceleration[1];
        const double div_Ha = (rDN_DX(i,0) * acceleration[0] + rDN_DX(i,1) * acceleration[1]) * depth;
        Ju += C2 * H2 * div_a + C4 * H * div_Ha;
    }

    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        const array_1d<double,3> normal_i = rData.normal * rN[i];
        MathUtils<double>::AddVector(rLaplacianBoundary, -Weight*normal_i*Ju, 3*i);
    }
}


template<std::size_t TNumNodes>
void BoussinesqCondition<TNumNodes>::CalculateShapeFunctionDerivaties(
    Matrix& rDN_DX,
    const GeometryType& rParentGeometry,
    const Point& rPoint)
{
    const std::size_t num_nodes = rParentGeometry.size();
    rDN_DX.resize(num_nodes, 2);
    if (num_nodes == 3) {
        BoundedMatrix<double, 3, 2> gradients;
        array_1d<double,3> shape_functions;
        double area;
        GeometryUtils::CalculateGeometryData(rParentGeometry, gradients, shape_functions, area);
        noalias(rDN_DX) = gradients;
    } else {
        const std::size_t working_space_dimension = 3;
        const std::size_t local_space_dimension = 2;
        Point local_coordinates;
        Matrix local_gradients(num_nodes, local_space_dimension);
        Matrix J(working_space_dimension, local_space_dimension);
        Matrix J_inv(local_space_dimension, working_space_dimension);
        double det_J;
        rParentGeometry.PointLocalCoordinates(local_coordinates, rPoint);
        rParentGeometry.ShapeFunctionsLocalGradients(local_gradients, local_coordinates);
        rParentGeometry.Jacobian(J, local_coordinates);
        MathUtils<double>::InvertMatrix(J, J_inv, det_J);
        noalias(rDN_DX) = prod(local_gradients, J_inv);
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
    this->CalculateGeometryData(r_geom, weights, N_container);
    const std::size_t num_gauss_points = weights.size();
    const auto& g_points = r_geom.IntegrationPoints();

    // Boundary term for the auxiliary fields
    LocalVectorType vel_laplacian_vector = ZeroVector(mLocalSize);
    LocalVectorType vel_H_laplacian_vector = ZeroVector(mLocalSize);

    // Gauss point contribution
    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const double weight = weights[g];
        const auto point = g_points[g];
        const array_1d<double,TNumNodes> N = row(N_container, g);

        CalculateGaussPointData(data, g, N);
        CalculateShapeFunctionDerivaties(DN_DX, parent_geom, point);
        AddLaplacianBoundary(vel_laplacian_vector, vel_H_laplacian_vector, parent_geom, data, N, DN_DX, weight);
    }

    array_1d<double,3> vel_laplacian = ZeroVector(3);
    array_1d<double,3> vel_H_laplacian = ZeroVector(3);
    for (std::size_t i = 0; i < TNumNodes; ++i)
    {
        std::size_t block = 3 * i;
        vel_laplacian[0] = vel_laplacian_vector[block];
        vel_laplacian[1] = vel_laplacian_vector[block + 1];
        vel_H_laplacian[0] = vel_H_laplacian_vector[block];
        vel_H_laplacian[1] = vel_H_laplacian_vector[block + 1];
        r_geom[i].SetLock();
        r_geom[i].FastGetSolutionStepValue(VELOCITY_LAPLACIAN) += vel_laplacian;
        r_geom[i].FastGetSolutionStepValue(VELOCITY_H_LAPLACIAN) += vel_H_laplacian;
        r_geom[i].UnSetLock();
    }
}


template<std::size_t TNumNodes>
void BoussinesqCondition<TNumNodes>::CalculateLocalSystem(Matrix& rLeftHandSideMatrix, Vector& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if(rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);

    if(rRightHandSideVector.size() != mLocalSize)
        rRightHandSideVector.resize(mLocalSize, false);

    // Geometries
    auto& r_geom = this->GetGeometry();
    auto& r_parent_geom = this->GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry();

    // Struct to pass around the data
    ConditionData data;
    WaveConditionType::InitializeData(data, rCurrentProcessInfo);

    // Geometrical data
    Matrix DN_DX;           // Gradients of the parent element
    Vector weights;         // Line integration data
    Matrix N_container;     // Line integration data
    this->CalculateGeometryData(r_geom, weights, N_container);
    const std::size_t num_gauss_points = weights.size();
    const auto& g_points = r_geom.IntegrationPoints();

    // Boundary term for the auxiliary fields
    LocalVectorType rhs = ZeroVector(mLocalSize);

    // Gauss point contribution
    for (IndexType g = 0; g < num_gauss_points; ++g)
    {
        const double weight = weights[g];
        const auto point = g_points[g];
        const array_1d<double,TNumNodes> N = row(N_container, g);

        CalculateGaussPointData(data, g, N);
        CalculateShapeFunctionDerivaties(DN_DX, r_parent_geom, point);
        AddMomentumDispersionTerms(rhs, r_parent_geom, data, N, DN_DX, weight);
        this->AddFluxTerms(rhs, data, N, weight);
    }

    noalias(rRightHandSideVector) = rhs;
}


template class BoussinesqCondition<2>;

} // namespace Kratos

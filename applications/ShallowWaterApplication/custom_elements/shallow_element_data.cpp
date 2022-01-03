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
#include "shallow_element_data.h"
#include "utilities/geometry_utilities.h"
#include "shallow_water_application_variables.h"
#include "custom_friction_laws/friction_laws_factory.h"

namespace Kratos
{

template<std::size_t TNumNodes>
void ShallowElementData<TNumNodes>::InitializeData(
    const GeometryType& rGeometry,
    const Properties& rProperties,
    const ProcessInfo& rCurrentProcessInfo)
{
    integrate_by_parts = rCurrentProcessInfo[INTEGRATE_BY_PARTS]; //since it is passed as const it will return false if it doesn't have INTEGRATE_BY_PARTS
    stab_factor = rCurrentProcessInfo[STABILIZATION_FACTOR];
    relative_dry_height = rCurrentProcessInfo[RELATIVE_DRY_HEIGHT];
    gravity = rCurrentProcessInfo[GRAVITY_Z];
    length = rGeometry.Length();
    p_bottom_friction = FrictionLawsFactory().CreateBottomFrictionLaw(
        rGeometry, rProperties, rCurrentProcessInfo);
}


template<std::size_t TNumNodes>
void ShallowElementData<TNumNodes>::SetNodalData(const GeometryType& rGeometry, int Step)
{
    KRATOS_ERROR << "Implement this method in the derived class." << std::endl;
}


template<std::size_t TNumNodes>
void ShallowElementData<TNumNodes>::UpdateGaussPointData(const array_1d<double,TNumNodes>& rN)
{
    KRATOS_ERROR << "Implement this method in the derived class." << std::endl;
}


template<std::size_t TNumNodes>
typename ShallowElementData<TNumNodes>::LocalVectorType& ShallowElementData<TNumNodes>::GetUnknownVector() const
{
    KRATOS_ERROR << "Implement this method in the derived class." << std::endl;
}


template<std::size_t TNumNodes>
const Variable<double>& ShallowElementData<TNumNodes>::UnknownComponent(int Index)
{
    KRATOS_ERROR << "Implement this method in the derived class." << std::endl;
}


template<>
void ShallowElementData<3>::CalculateGeometryData(const GeometryType& rGeometry)
{
    BoundedMatrix<double,3,2> DN_DX; // Gradients matrix
    array_1d<double,3> N;            // Position of the gauss point
    double area;
    GeometryUtils::CalculateGeometryData(rGeometry, DN_DX, N, area);

    if (mWeights.size() != 1) {
        mWeights.resize(1, false);
    }
    mWeights[0] = area;

    if (mN_Container.size1() != 1 && mN_Container.size2() != 3) {
        mN_Container.resize(1, 3, false);
    }
    for (std::size_t i = 0; i < 3; ++i) {
        mN_Container(0, i) = N[i];
    }

    if (mDN_DX_Container.size() != 1) {
        mDN_DX_Container.resize(1, false);
    }
    mDN_DX_Container[0] = DN_DX;
}


template<std::size_t TNumNodes>
void ShallowElementData<TNumNodes>::CalculateGeometryData(const GeometryType& rGeometry)
{
    Vector det_j_vector;
    const auto integration_method = rGeometry.GetDefaultIntegrationMethod();
    mN_Container = rGeometry.ShapeFunctionsValues(integration_method);
    rGeometry.ShapeFunctionsIntegrationPointsGradients(mDN_DX_Container, det_j_vector, integration_method);

    const unsigned int number_of_gauss_points = rGeometry.IntegrationPointsNumber(integration_method);
    const GeometryType::IntegrationPointsArrayType& integration_points = rGeometry.IntegrationPoints(integration_method);

    if (mWeights.size() != number_of_gauss_points)
        mWeights.resize(number_of_gauss_points, false);

    for (unsigned int g = 0; g < number_of_gauss_points; ++g)
        mWeights[g] = det_j_vector[g] * integration_points[g].Weight();
}


template<>
double ShallowElementData<3>::ShapeFunctionProduct(
    const array_1d<double,3>& rN,
    const std::size_t I,
    const std::size_t J)
{
    return (I == J) ? 1.0/6.0 : 1.0/12.0;
}


template<std::size_t TNumNodes>
double ShallowElementData<TNumNodes>::ShapeFunctionProduct(
    const array_1d<double,TNumNodes>& rN,
    const std::size_t I,
    const std::size_t J)
{
    return rN[I] * rN[J];
}


template<std::size_t TNumNodes>
array_1d<double,3> ShallowElementData<TNumNodes>::VectorInnerProduct(
    const array_1d<array_1d<double,3>,TNumNodes>& rV,
    const array_1d<double,TNumNodes>& rN)
{
    array_1d<double,3> result = ZeroVector(3);
    for (std::size_t i = 0; i < TNumNodes; ++i)
    {
        result += rV[i] * rN[i];
    }
    return result;
}


template<std::size_t TNumNodes>
std::size_t ShallowElementData<TNumNodes>::NumberOfGaussPoints()
{
    return mWeights.size();
}


template<std::size_t TNumNodes>
const double ShallowElementData<TNumNodes>::GetWeight(int GaussPointIndex)
{
    return mWeights[GaussPointIndex];
}


template<std::size_t TNumNodes>
const array_1d<double,TNumNodes> ShallowElementData<TNumNodes>::GetShapeFunctions(int GaussPointIndex)
{
    return row(mN_Container, GaussPointIndex);
}


template<std::size_t TNumNodes>
const BoundedMatrix<double,TNumNodes,2> ShallowElementData<TNumNodes>::GetShapeFunctionDerivatives(int GaussPointIndex)
{
    return mDN_DX_Container[GaussPointIndex];
}


template class ShallowElementData<3>;
template class ShallowElementData<4>;
template class ShallowElementData<6>;
template class ShallowElementData<8>;
template class ShallowElementData<9>;

} // namespace Kratos

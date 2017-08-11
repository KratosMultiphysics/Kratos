//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#include "integration_point_data.h"

namespace Kratos
{

template <class TElementData>
void IntegrationPointData<TElementData>::FillIntegrationPointData(
    IntegrationPointData& rIntegrationPointData,
    const TElementData& rElementData,
    int GaussPointIndex,
    const Vector& rGaussWeights,
    const Matrix& rNContainer,
    const Geometry<Node<3> >::ShapeFunctionsGradientsType& rDN_DX)
{
    rIntegrationPointData.IntegrationPointIndex = GaussPointIndex;
    noalias(rIntegrationPointData.N) = row(rNContainer,GaussPointIndex);
    noalias(rIntegrationPointData.DN_DX) = rDN_DX[GaussPointIndex];
    
    rIntegrationPointData.EvaluateInPoint(rIntegrationPointData.Pressure,rElementData.Pressure);
    rIntegrationPointData.EvaluateInPoint(rIntegrationPointData.Density, rElementData.Density);
    rIntegrationPointData.EvaluateInPoint(rIntegrationPointData.Viscosity, rElementData.Viscosity);
    rIntegrationPointData.EvaluateInPoint(rIntegrationPointData.MassProjection, rElementData.MassProjection);

    rIntegrationPointData.EvaluateInPoint(rIntegrationPointData.Velocity, rElementData.Velocity);
    rIntegrationPointData.EvaluateInPoint(rIntegrationPointData.MeshVelocity, rElementData.MeshVelocity);
    noalias(rIntegrationPointData.ConvectiveVelocity) = rIntegrationPointData.Velocity - rIntegrationPointData.MeshVelocity;

    rIntegrationPointData.EvaluateInPoint(rIntegrationPointData.BodyForce, rElementData.BodyForce);
    rIntegrationPointData.EvaluateInPoint(rIntegrationPointData.MomentumProjection, rElementData.MomentumProjection);
}

template <class TElementData>
void IntegrationPointData<TElementData>::EvaluateInPoint(
    double& rResult, const typename TElementData::ScalarDataType& rNodalValues)
{
    rResult = N[0] * rNodalValues[0];

    for (unsigned int i = 1; i < TElementData::NumNodes; i++)
    {
        rResult += N[i] * rNodalValues[i];
    }
}

template <class TElementData>
void IntegrationPointData<TElementData>::EvaluateInPoint(
    array_1d<double, TElementData::Dim>& rResult,
    const typename TElementData::VectorDataType& rNodalValues)
{
    for (unsigned int d = 0; d < TElementData::Dim; d++)
    {
        rResult[d] = N[0] * rNodalValues(0, d);
    }

    for (unsigned int i = 1; i < TElementData::NumNodes; i++)
    {
        for (unsigned int d = 0; d < TElementData::Dim; d++)
        {
            rResult[d] += N[i] * rNodalValues(i, d);
        }
    }
}

template class IntegrationPointData< FluidElementData<2,3> >;
template class IntegrationPointData< FluidElementData<3,4> >;

}
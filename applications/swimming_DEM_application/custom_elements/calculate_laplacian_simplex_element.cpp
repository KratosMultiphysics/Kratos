#include "swimming_DEM_application.h"
#include "calculate_laplacian_simplex_element.h"

namespace Kratos
{

template <>
void ComputeLaplacianSimplex<2,3>::CalculateWeights(ShapeFunctionDerivativesArrayType& rDN_DX,
        Matrix& rNContainer,
        Vector& rGaussWeights)
{

  const GeometryType& rGeom = this->GetGeometry();
  Vector DetJ;
  rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, GeometryData::GI_GAUSS_2);
  rNContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
  const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);

  rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2), false);

  for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2); g++)
      rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
}


template <>
void ComputeLaplacianSimplex<3,4>::CalculateWeights(ShapeFunctionDerivativesArrayType& rDN_DX,
        Matrix& rNContainer,
        Vector& rGaussWeights)
{
    const GeometryType& rGeom = this->GetGeometry();
    Vector DetJ;
    rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, GeometryData::GI_GAUSS_2);
    rNContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);

    rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2), false);

    for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2); g++)
        rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
}

template <>
void ComputeLaplacianSimplex<2,3>::AddRHSLaplacian(VectorType& F,
                             const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,
                             const double Weight)
{
    double Coef = Weight;
    array_1d<double, 3 > Velocity;

    int LocalIndex = 0;
    int LocalNodalIndex = 0;

    for (unsigned int iNode = 0; iNode < 3; ++iNode)
    {
        Velocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
        for (unsigned int d = 0; d < 2; ++d)
        {
            F[LocalIndex++] -= Coef * rShapeDeriv(LocalNodalIndex, d) * Velocity[d] * rShapeDeriv(iNode, d);
        }
        LocalNodalIndex++;
    }
}

template <>
void ComputeLaplacianSimplex<3,4>::AddRHSLaplacian(VectorType& F,
                             const boost::numeric::ublas::bounded_matrix<double, 4, 3>& rShapeDeriv,
                             const double Weight)
{
    double Coef = Weight;

    int LocalIndex = 0;
    for (unsigned int iNodeB = 0; iNodeB < 4; ++iNodeB){

        for (unsigned int dj = 0; dj < 3; ++dj){
            double value = 0.0;

            for (unsigned int iNodeA = 0; iNodeA < 4; ++iNodeA){
                const array_1d<double, 3 >& Velocity = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY);

                for (unsigned int di = 0; di < 3; ++di){
                    value -= rShapeDeriv(iNodeB, di) * Velocity[dj] * rShapeDeriv(iNodeA, di);
                }
            }

            F[LocalIndex++] += Coef * value;
        }
    }
}

} // namespace Kratos

#include "swimming_DEM_application.h"
#include "calculate_mat_deriv_simplex_element.h"

namespace Kratos
{

template <>
void ComputeMaterialDerivativeSimplex<2,3>::CalculateWeights(ShapeFunctionDerivativesArrayType& rDN_DX,
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
void ComputeMaterialDerivativeSimplex<3,4>::CalculateWeights(ShapeFunctionDerivativesArrayType& rDN_DX,
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
void ComputeMaterialDerivativeSimplex<2,3>::AddRHSMatAcceleration(VectorType& F,
                             const array_1d<double,3>& rShapeFunc,
                             const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,
                             const double Weight)
{
    double Coef = Weight;
    array_1d<double, 3 > Velocity;

    int LocalIndex = 0;
    for (unsigned int iNodeB = 0; iNodeB < 3; ++iNodeB){

        for (unsigned int dj = 0; dj < 2; ++dj){
            double value = 0.0;

            for (unsigned int iNodeA = 0; iNodeA < 3; ++iNodeA){
                Velocity = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY);

                for (unsigned int di = 0; di < 2; ++di){
                    value += rShapeFunc[iNodeB] * Velocity[di] * rShapeDeriv(iNodeA, di) * Velocity[dj];
                }
            }

            F[LocalIndex++] += Coef * value;
        }
    }
}

template <>
void ComputeMaterialDerivativeSimplex<3,4>::AddRHSMatAcceleration(VectorType& F,
                             const array_1d<double,4>& rShapeFunc,
                             const boost::numeric::ublas::bounded_matrix<double, 4, 3>& rShapeDeriv,
                             const double Weight)
{
    double Coef = Weight;
    array_1d<double, 3 > Velocity;

    int LocalIndex = 0;
    for (unsigned int iNodeB = 0; iNodeB < 4; ++iNodeB){

        for (unsigned int dj = 0; dj < 3; ++dj){
            double value = 0.0;

            for (unsigned int iNodeA = 0; iNodeA < 4; ++iNodeA){
                Velocity = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY);

                for (unsigned int di = 0; di < 3; ++di){
                    value += rShapeFunc[iNodeB] * Velocity[di] * rShapeDeriv(iNodeA, di) * Velocity[dj];
                }
            }

            F[LocalIndex++] += Coef * value;
        }
    }
}

} // namespace Kratos

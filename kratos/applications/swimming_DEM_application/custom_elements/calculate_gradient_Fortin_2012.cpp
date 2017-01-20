#include "swimming_DEM_application.h"
#include "calculate_gradient_Fortin_2012.h"

namespace Kratos
{

template <>
void ComputeGradientFortin2012<2,3>::CalculateWeights(ShapeFunctionDerivativesArrayType& rDN_DX,
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
void ComputeGradientFortin2012<3,4>::CalculateWeights(ShapeFunctionDerivativesArrayType& rDN_DX,
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
void ComputeGradientFortin2012<2,3>::EvaluateInPoint(array_1d< double, 3 > & rResult,
                             const Variable< array_1d< double, 3 > >& rVariable,
                             const array_1d< double, 3 >& rShapeFunc)
{
    // Compute the weighted value of the nodal variable in the (Gauss) Point
    noalias(rResult) = rShapeFunc[0] * this->GetGeometry()[0].FastGetSolutionStepValue(rVariable);
    for (unsigned int iNode = 1; iNode < 3; ++iNode)
        noalias(rResult) += rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable);
}

template <>
void ComputeGradientFortin2012<3,4>::EvaluateInPoint(array_1d< double, 3 > & rResult,
                             const Variable< array_1d< double, 3 > >& rVariable,
                             const array_1d< double, 4 >& rShapeFunc)
{
    // Compute the weighted value of the nodal variable in the (Gauss) Point
    rResult = rShapeFunc[0] * this->GetGeometry()[0].FastGetSolutionStepValue(rVariable);
    for (unsigned int iNode = 1; iNode < 4; ++iNode)
        rResult += rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable);
}

template <>
void ComputeGradientFortin2012<2,3>::AddRHSGradient(VectorType& F,
                             const array_1d<double,3>& rShapeFunc,
                             const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,
                             const double Weight)
{
    double Coef = Weight;
    int LocalIndex = 0;

    for (unsigned int iNodeB = 0; iNodeB < 3; ++iNodeB){

        for (unsigned int di = 0; di < 2; ++di){
            double value = 0.0;

            for (unsigned int iNodeA = 0; iNodeA < 3; ++iNodeA){
                double NodalComponent = 0.0;

                if (mCurrentComponent == 'X'){
                    NodalComponent = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY_X);
                }
                else if (mCurrentComponent == 'Y'){
                    NodalComponent = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY_Y);
                }

                value += rShapeFunc[iNodeB] * rShapeDeriv(iNodeA, di) * NodalComponent;
            }

            F[LocalIndex++] += Coef * value;
        }
    }
}

template <>
void ComputeGradientFortin2012<3,4>::AddRHSGradient(VectorType& F,
                             const array_1d<double,4>& rShapeFunc,
                             const boost::numeric::ublas::bounded_matrix<double, 4, 3>& rShapeDeriv,
                             const double Weight)
{
    double Coef = Weight;
    int LocalIndex = 0;

    for (unsigned int iNodeB = 0; iNodeB < 4; ++iNodeB){

        for (unsigned int di = 0; di < 3; ++di){
            double value = 0.0;

            for (unsigned int iNodeA = 0; iNodeA < 4; ++iNodeA){
                double NodalComponent = 0.0;

                if (mCurrentComponent == 'X'){
                    NodalComponent = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY_X);
                }
                else if (mCurrentComponent == 'Y'){
                    NodalComponent = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY_Y);
                }
                else if (mCurrentComponent == 'Z'){
                    NodalComponent = this->GetGeometry()[iNodeA].FastGetSolutionStepValue(VELOCITY_Z);
                }
                value += rShapeFunc[iNodeB] * rShapeDeriv(iNodeA, di) * NodalComponent;
            }

            F[LocalIndex++] += Coef * value;
        }
    }
}

} // namespace Kratos

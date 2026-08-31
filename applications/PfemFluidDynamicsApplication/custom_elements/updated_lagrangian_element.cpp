//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:              April 2018 $
//   Revision:            $Revision:                 0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_element.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

  /*
   * public UpdatedLagrangianElement<TDim> functions
   */

  template <unsigned int TDim>
  UpdatedLagrangianElement<TDim>::UpdatedLagrangianElement(UpdatedLagrangianElement const &rOther)
      : Element(rOther)
  {
    KRATOS_TRY;
    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  Element::Pointer UpdatedLagrangianElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
  {
    KRATOS_TRY;

    UpdatedLagrangianElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    return Element::Pointer(new UpdatedLagrangianElement(NewElement));

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  GeometryData::IntegrationMethod UpdatedLagrangianElement<TDim>::GetIntegrationMethod() const
  {
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
  }

  template <>
  void UpdatedLagrangianElement<2>::VelocityEquationIdVector(EquationIdVectorType &rResult,
                                                             const ProcessInfo &rCurrentProcessInfo) const
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = NumNodes * 2;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
      rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X, xpos).EquationId();
      rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y, xpos + 1).EquationId();
    }
  }

  template <>
  void UpdatedLagrangianElement<3>::VelocityEquationIdVector(EquationIdVectorType &rResult,
                                                             const ProcessInfo &rCurrentProcessInfo) const
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3 * NumNodes;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
      rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X, xpos).EquationId();
      rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y, xpos + 1).EquationId();
      rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z, xpos + 2).EquationId();
    }
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::PressureEquationIdVector(EquationIdVectorType &rResult,
                                                                const ProcessInfo &rCurrentProcessInfo) const
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rResult.size() != NumNodes)
      rResult.resize(NumNodes, false);

    const unsigned int pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (SizeType i = 0; i < NumNodes; ++i)
      rResult[i] = rGeom[i].GetDof(PRESSURE, pos).EquationId();
  }

  template <>
  void UpdatedLagrangianElement<2>::GetVelocityDofList(DofsVectorType &rElementalDofList,
                                                       const ProcessInfo &rCurrentProcessInfo) const
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2 * NumNodes;

    if (rElementalDofList.size() != LocalSize)
      rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
      rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
    }
  }

  template <>
  void UpdatedLagrangianElement<3>::GetVelocityDofList(DofsVectorType &rElementalDofList,
                                                       const ProcessInfo &rCurrentProcessInfo) const
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3 * NumNodes;

    if (rElementalDofList.size() != LocalSize)
      rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
      rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
      rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
    }
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::GetPressureDofList(DofsVectorType &rElementalDofList,
                                                          const ProcessInfo &rCurrentProcessInfo) const
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rElementalDofList.size() != NumNodes)
      rElementalDofList.resize(NumNodes);

    SizeType LocalIndex = 0;
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
    }
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::GetPressureValues(Vector &rValues,
                                                         const int Step)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rValues.size() != NumNodes)
      rValues.resize(NumNodes, false);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE, Step);
    }
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::GetFluidFractionRateValues(Vector &rValues)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rValues.size() != NumNodes)
      rValues.resize(NumNodes, false);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rValues[i] = rGeom[i].FastGetSolutionStepValue(FLUID_FRACTION_RATE);
    }
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::GetDensityValues(Vector &rValues,
                                                        const int Step)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rValues.size() != NumNodes)
      rValues.resize(NumNodes, false);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rValues[i] = rGeom[i].FastGetSolutionStepValue(DENSITY, Step);
    }
  }

  template <>
  void UpdatedLagrangianElement<2>::CalcMeanVelocityNorm(double &meanVelocity,
                                                         const int Step)
  {

    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    double velX = 0;
    double velY = 0;
    double coeff = 0.0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      velX += rGeom[i].FastGetSolutionStepValue(VELOCITY_X, Step);
      velY += rGeom[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
      coeff += 1.0;
    }

    meanVelocity = velX * velX + velY * velY;
    meanVelocity = sqrt(meanVelocity) / coeff;
  }

  template <>
  void UpdatedLagrangianElement<3>::CalcMeanVelocityNorm(double &meanVelocity,
                                                         const int Step)
  {

    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    double velX = 0;
    double velY = 0;
    double velZ = 0;
    double coeff = 0.0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      velX += rGeom[i].FastGetSolutionStepValue(VELOCITY_X, Step);
      velY += rGeom[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
      velZ += rGeom[i].FastGetSolutionStepValue(VELOCITY_Z, Step);
      coeff += 1.0;
    }

    meanVelocity = velX * velX + velY * velY + velZ * velZ;
    meanVelocity = sqrt(meanVelocity) / coeff;
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::CalculateDeltaPosition(Matrix &rDeltaPosition)
  {
    KRATOS_TRY

    const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    noalias(rDeltaPosition) = ZeroMatrix(number_of_nodes, dimension);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
      array_1d<double, 3> &CurrentDisplacement = this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
      array_1d<double, 3> &PreviousDisplacement = this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, 1);

      for (unsigned int j = 0; j < dimension; j++)
      {
        rDeltaPosition(i, j) = CurrentDisplacement[j] - PreviousDisplacement[j];
      }
    }

    KRATOS_CATCH("")
  }

  template <>
  void UpdatedLagrangianElement<2>::GetDisplacementValues(Vector &rValues,
                                                          const int Step)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2 * NumNodes;

    if (rValues.size() != LocalSize)
      rValues.resize(LocalSize, false);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_X, Step);
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
    }
  }

  template <>
  void UpdatedLagrangianElement<3>::GetDisplacementValues(Vector &rValues,
                                                          const int Step)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3 * NumNodes;

    if (rValues.size() != LocalSize)
      rValues.resize(LocalSize, false);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_X, Step);
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT_Z, Step);
    }
  }

  template <>
  void UpdatedLagrangianElement<2>::GetVelocityValues(Vector &rValues,
                                                      const int Step)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2 * NumNodes;

    if (rValues.size() != LocalSize)
      rValues.resize(LocalSize, false);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X, Step);
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
    }
  }

  template <>
  void UpdatedLagrangianElement<3>::GetVelocityValues(Vector &rValues,
                                                      const int Step)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3 * NumNodes;

    if (rValues.size() != LocalSize)
      rValues.resize(LocalSize, false);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X, Step);
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Z, Step);
    }
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::GetElementalAcceleration(Vector &meanAcceleration,
                                                                const int Step,
                                                                const double TimeStep)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    double count = 0;
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      meanAcceleration += 0.5 / TimeStep * (rGeom[i].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[i].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[i].FastGetSolutionStepValue(ACCELERATION, 1);
      count += 1.0;
    }
    meanAcceleration *= 1.0 / count;
  }

  template <>
  void UpdatedLagrangianElement<2>::GetAccelerationValues(Vector &rValues,
                                                          const int Step)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2 * NumNodes;

    if (rValues.size() != LocalSize)
      rValues.resize(LocalSize, false);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_X, Step);
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_Y, Step);
    }
  }

  template <>
  void UpdatedLagrangianElement<3>::GetAccelerationValues(Vector &rValues,
                                                          const int Step)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3 * NumNodes;

    if (rValues.size() != LocalSize)
      rValues.resize(LocalSize, false);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_X, Step);
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_Y, Step);
      rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_Z, Step);
    }
  }

  template <>
  void UpdatedLagrangianElement<2>::GetPositions(Vector &rValues, const ProcessInfo &rCurrentProcessInfo, const double theta)
  {

    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2 * NumNodes;

    if (rValues.size() != LocalSize)
      rValues.resize(LocalSize, false);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rValues[Index++] = rGeom[i].X();
      rValues[Index++] = rGeom[i].Y();
    }
  }

  template <>
  void UpdatedLagrangianElement<3>::GetPositions(Vector &rValues, const ProcessInfo &rCurrentProcessInfo, const double theta)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3 * NumNodes;

    if (rValues.size() != LocalSize)
      rValues.resize(LocalSize, false);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rValues[Index++] = rGeom[i].X();
      rValues[Index++] = rGeom[i].Y();
      rValues[Index++] = rGeom[i].Z();
    }
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::CalculateGeometryData(ShapeFunctionDerivativesArrayType &rDN_DX,
                                                             Matrix &NContainer,
                                                             Vector &rGaussWeights)
  {
    const GeometryType &rGeom = this->GetGeometry();
    Vector DetJ;
    rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, GeometryData::IntegrationMethod::GI_GAUSS_1);
    NContainer = rGeom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);
    const GeometryType::IntegrationPointsArrayType &IntegrationPoints = rGeom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_1);

    rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_1), false);

    for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_1); ++g)
    {
      rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
    }
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::CalculateGeometryData(Vector &rGaussWeights)
  {
    const GeometryType &rGeom = this->GetGeometry();
    Vector DetJ;
    rGeom.DeterminantOfJacobian(DetJ, GeometryData::IntegrationMethod::GI_GAUSS_1);
    const GeometryType::IntegrationPointsArrayType &IntegrationPoints = rGeom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_1);

    rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_1), false);

    for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::IntegrationMethod::GI_GAUSS_1); ++g)
    {
      rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
    }
  }

  template <unsigned int TDim>
  double UpdatedLagrangianElement<TDim>::ElementSize()
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    double ElemSize = 0.0;
    array_1d<double, 3> Edge = rGeom[1].Coordinates() - rGeom[0].Coordinates();
    ElemSize += norm_2(Edge);

    double count = 1.0;
    for (SizeType i = 2; i < NumNodes; i++)
    {
      for (SizeType j = 0; j < i; j++)
      {
        noalias(Edge) = rGeom[i].Coordinates() - rGeom[j].Coordinates();
        ElemSize += norm_2(Edge);
        count += 1.0;
      }
    }
    ElemSize *= 1.0 / count;
    return ElemSize;
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::CalcFGrad(const ShapeFunctionDerivativesType &rDN_DX,
                                                 MatrixType &Fgrad,
                                                 MatrixType &invFgrad,
                                                 double &FJacobian,
                                                 const ProcessInfo &rCurrentProcessInfo,
                                                 const double theta)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = TDim * NumNodes;
    VectorType NodePosition = ZeroVector(LocalSize);
    this->GetPositions(NodePosition, rCurrentProcessInfo, theta);

    Fgrad.resize(TDim, TDim, false);

    noalias(Fgrad) = ZeroMatrix(TDim, TDim);
    for (SizeType i = 0; i < TDim; i++)
    {
      for (SizeType j = 0; j < TDim; j++)
      {
        for (SizeType k = 0; k < NumNodes; k++)
        {
          Fgrad(i, j) += NodePosition[TDim * k + i] * rDN_DX(k, j);
        }
      }
    }

    // Inverse
    invFgrad.resize(TDim, TDim, false);
    noalias(invFgrad) = ZeroMatrix(TDim, TDim);
    FJacobian = 1;

    if constexpr (TDim == 2)
    {
      MathUtils<double>::InvertMatrix2(Fgrad, invFgrad, FJacobian);
    }
    else if constexpr (TDim == 3)
    {
      MathUtils<double>::InvertMatrix3(Fgrad, invFgrad, FJacobian);
    }
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::CalcVelDefGrad(const ShapeFunctionDerivativesType &rDN_DX,
                                                      MatrixType &FgradVel,
                                                      const double theta)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = TDim * NumNodes;
    VectorType VelocityValues = ZeroVector(LocalSize);
    VectorType RHSVelocities = ZeroVector(LocalSize);
    this->GetVelocityValues(RHSVelocities, 0);
    RHSVelocities *= theta;
    this->GetVelocityValues(VelocityValues, 1);
    RHSVelocities += VelocityValues * (1.0 - theta);

    FgradVel.resize(TDim, TDim, false);

    noalias(FgradVel) = ZeroMatrix(TDim, TDim);
    for (SizeType i = 0; i < TDim; i++)
    {
      for (SizeType j = 0; j < TDim; j++)
      {
        for (SizeType k = 0; k < NumNodes; k++)
        {
          FgradVel(i, j) += RHSVelocities[TDim * k + i] * rDN_DX(k, j);
        }
      }
    }
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::CalcVelDefGradAndInverse(const ShapeFunctionDerivativesType &rDN_DX,
                                                                MatrixType &FgradVel,
                                                                MatrixType &invFgradVel,
                                                                double &FVelJacobian,
                                                                const double theta)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = TDim * NumNodes;
    VectorType VelocityValues = ZeroVector(LocalSize);
    VectorType RHSVelocities = ZeroVector(LocalSize);
    this->GetVelocityValues(RHSVelocities, 0);
    RHSVelocities *= theta;
    this->GetVelocityValues(VelocityValues, 1);
    RHSVelocities += VelocityValues * (1.0 - theta);

    FgradVel.resize(TDim, TDim, false);

    noalias(FgradVel) = ZeroMatrix(TDim, TDim);
    for (SizeType i = 0; i < TDim; i++)
    {
      for (SizeType j = 0; j < TDim; j++)
      {
        for (SizeType k = 0; k < NumNodes; k++)
        {
          FgradVel(i, j) += RHSVelocities[TDim * k + i] * rDN_DX(k, j);
        }
      }
    }

    // Inverse
    invFgradVel.resize(TDim, TDim, false);
    noalias(invFgradVel) = ZeroMatrix(TDim, TDim);
    FVelJacobian = 1;

    if constexpr (TDim == 2)
    {
      MathUtils<double>::InvertMatrix2(FgradVel, invFgradVel, FVelJacobian);
    }
    else if constexpr (TDim == 3)
    {
      MathUtils<double>::InvertMatrix3(FgradVel, invFgradVel, FVelJacobian);
    }
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::CalcSpatialVelocityGrad(const MatrixType &invFgrad,
                                                               const MatrixType &VelDefgrad,
                                                               MatrixType &SpatialVelocityGrad)
  {
    SpatialVelocityGrad.resize(TDim, TDim, false);

    SpatialVelocityGrad = prod(VelDefgrad, invFgrad);
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::CalcVolDefRateFromSpatialVelGrad(double &volumetricDefRate,
                                                                        const MatrixType &SpatialVelocityGrad)
  {
    volumetricDefRate = 0;
    for (SizeType i = 0; i < TDim; i++)
    {
      volumetricDefRate += SpatialVelocityGrad(i, i);
    }
  }

  template <>
  void UpdatedLagrangianElement<2>::CalcMDGreenLagrangeMaterial(const MatrixType &Fgrad,
                                                                const MatrixType &VelDefgrad,
                                                                VectorType &MDGreenLagrangeMaterial)
  {

    // x-component
    MDGreenLagrangeMaterial[0] = VelDefgrad(0, 0) * Fgrad(0, 0) + VelDefgrad(1, 0) * Fgrad(1, 0);
    // y-component
    MDGreenLagrangeMaterial[1] = VelDefgrad(1, 1) * Fgrad(1, 1) + VelDefgrad(0, 1) * Fgrad(0, 1);
    // xy-component
    MDGreenLagrangeMaterial[2] = (VelDefgrad(0, 0) * Fgrad(0, 1) + VelDefgrad(1, 0) * Fgrad(1, 1) +
                                  VelDefgrad(0, 1) * Fgrad(0, 0) + VelDefgrad(1, 1) * Fgrad(1, 0)) *
                                 0.5;
  }

  template <>
  void UpdatedLagrangianElement<3>::CalcMDGreenLagrangeMaterial(const MatrixType &Fgrad,
                                                                const MatrixType &VelDefgrad,
                                                                VectorType &MDGreenLagrangeMaterial)
  {
    MatrixType FgradTransp = ZeroMatrix(3, 3);
    MatrixType VelDefgradTransp = ZeroMatrix(3, 3);
    MatrixType part1 = ZeroMatrix(3, 3);
    MatrixType part2 = ZeroMatrix(3, 3);

    FgradTransp = Fgrad;
    FgradTransp(0, 1) = Fgrad(1, 0);
    FgradTransp(0, 2) = Fgrad(2, 0);
    FgradTransp(1, 0) = Fgrad(0, 1);
    FgradTransp(1, 2) = Fgrad(2, 1);
    FgradTransp(2, 0) = Fgrad(0, 2);
    FgradTransp(2, 1) = Fgrad(1, 2);

    VelDefgradTransp = VelDefgrad;
    VelDefgradTransp(0, 1) = VelDefgrad(1, 0);
    VelDefgradTransp(0, 2) = VelDefgrad(2, 0);
    VelDefgradTransp(1, 0) = VelDefgrad(0, 1);
    VelDefgradTransp(1, 2) = VelDefgrad(2, 1);
    VelDefgradTransp(2, 0) = VelDefgrad(0, 2);
    VelDefgradTransp(2, 1) = VelDefgrad(1, 2);

    part1 = prod(VelDefgradTransp, Fgrad);
    part2 = prod(FgradTransp, VelDefgrad);

    MDGreenLagrangeMaterial[0] = (part1(0, 0) + part2(0, 0)) * 0.5; // xx-component
    MDGreenLagrangeMaterial[1] = (part1(1, 1) + part2(1, 1)) * 0.5; // yy-component
    MDGreenLagrangeMaterial[2] = (part1(2, 2) + part2(2, 2)) * 0.5; // zz-component
    MDGreenLagrangeMaterial[3] = (part1(0, 1) + part2(0, 1)) * 0.5; // xy-component
    MDGreenLagrangeMaterial[4] = (part1(0, 2) + part2(0, 2)) * 0.5; // xz-component
    MDGreenLagrangeMaterial[5] = (part1(1, 2) + part2(1, 2)) * 0.5; // yz-component
  }

  template <>
  void UpdatedLagrangianElement<2>::CalcSpatialDefRate(const VectorType &MDGreenLagrangeMaterial,
                                                       const MatrixType &invFgrad,
                                                       VectorType &SpatialDefRate)
  {
    // x-component
    SpatialDefRate[0] = invFgrad(0, 0) * MDGreenLagrangeMaterial[0] * invFgrad(0, 0) +
                        invFgrad(1, 0) * MDGreenLagrangeMaterial[2] * invFgrad(0, 0) * 2 +
                        invFgrad(1, 0) * MDGreenLagrangeMaterial[1] * invFgrad(1, 0);
    // y-component
    SpatialDefRate[1] = invFgrad(0, 1) * MDGreenLagrangeMaterial[0] * invFgrad(0, 1) +
                        invFgrad(0, 1) * MDGreenLagrangeMaterial[2] * invFgrad(1, 1) * 2 +
                        invFgrad(1, 1) * MDGreenLagrangeMaterial[1] * invFgrad(1, 1);
    // xy-component
    SpatialDefRate[2] = invFgrad(0, 0) * MDGreenLagrangeMaterial[0] * invFgrad(0, 1) +
                        invFgrad(0, 0) * MDGreenLagrangeMaterial[2] * invFgrad(1, 1) +
                        invFgrad(1, 0) * MDGreenLagrangeMaterial[2] * invFgrad(0, 1) +
                        invFgrad(1, 0) * MDGreenLagrangeMaterial[1] * invFgrad(1, 1);
  }

  template <>
  void UpdatedLagrangianElement<3>::CalcSpatialDefRate(const VectorType &MDGreenLagrangeMaterial,
                                                       const MatrixType &invFgrad,
                                                       VectorType &SpatialDefRate)
  {
    MatrixType MDGLM = ZeroMatrix(3, 3);
    MatrixType invFgradTransp = ZeroMatrix(3, 3);
    MatrixType part1 = ZeroMatrix(3, 3);
    MatrixType totalMatrix = ZeroMatrix(3, 3);

    invFgradTransp = invFgrad;
    invFgradTransp(0, 1) = invFgrad(1, 0);
    invFgradTransp(0, 2) = invFgrad(2, 0);
    invFgradTransp(1, 0) = invFgrad(0, 1);
    invFgradTransp(1, 2) = invFgrad(2, 1);
    invFgradTransp(2, 0) = invFgrad(0, 2);
    invFgradTransp(2, 1) = invFgrad(1, 2);

    MDGLM(0, 0) = MDGreenLagrangeMaterial[0]; // XX-component;
    MDGLM(1, 1) = MDGreenLagrangeMaterial[1]; // YY-component;
    MDGLM(2, 2) = MDGreenLagrangeMaterial[2]; // ZZ-component;
    MDGLM(0, 1) = MDGreenLagrangeMaterial[3]; // XY-component;
    MDGLM(1, 0) = MDGreenLagrangeMaterial[3]; // XY-component;
    MDGLM(0, 2) = MDGreenLagrangeMaterial[4]; // ZX-component;
    MDGLM(2, 0) = MDGreenLagrangeMaterial[4]; // ZX-component;
    MDGLM(1, 2) = MDGreenLagrangeMaterial[5]; // YZ-component;
    MDGLM(2, 1) = MDGreenLagrangeMaterial[5]; // YZ-component;

    part1 = prod(MDGLM, invFgrad);

    totalMatrix = prod(invFgradTransp, part1);

    SpatialDefRate[0] = totalMatrix(0, 0);
    SpatialDefRate[1] = totalMatrix(1, 1);
    SpatialDefRate[2] = totalMatrix(2, 2);
    SpatialDefRate[3] = totalMatrix(0, 1);
    SpatialDefRate[4] = totalMatrix(0, 2);
    SpatialDefRate[5] = totalMatrix(1, 2);
  }

  template <>
  void UpdatedLagrangianElement<2>::CalcDeviatoricInvariant(const VectorType &SpatialDefRate,
                                                            double &DeviatoricInvariant)
  {
    const double trace_d = SpatialDefRate[0] + SpatialDefRate[1];
    const double dev_X = SpatialDefRate[0] - trace_d / 3.0;
    const double dev_Y = SpatialDefRate[1] - trace_d / 3.0;
    DeviatoricInvariant = sqrt(2 * (dev_X * dev_X + SpatialDefRate[2] * SpatialDefRate[2] + dev_Y * dev_Y));
  }

  template <>
  void UpdatedLagrangianElement<2>::CalcEquivalentStrainRate(const VectorType &SpatialDefRate,
                                                             double &EquivalentStrainRate)
  {
    EquivalentStrainRate = sqrt(2.0 * (SpatialDefRate[0] * SpatialDefRate[0] +
                                       SpatialDefRate[1] * SpatialDefRate[1] +
                                       2.0 * SpatialDefRate[2] * SpatialDefRate[2]));
  }

  template <>
  void UpdatedLagrangianElement<3>::CalcDeviatoricInvariant(const VectorType &SpatialDefRate,
                                                            double &DeviatoricInvariant)
  {
    const double trace_d = SpatialDefRate[0] + SpatialDefRate[1] + SpatialDefRate[2];
    const double dev_X = SpatialDefRate[0] - trace_d / 3.0;
    const double dev_Y = SpatialDefRate[1] - trace_d / 3.0;
    const double dev_Z = SpatialDefRate[2] - trace_d / 3.0;
    DeviatoricInvariant = sqrt(2 * (dev_X * dev_X + dev_Y * dev_Y + dev_Z * dev_Z +
                                    SpatialDefRate[3] * SpatialDefRate[3] +
                                    SpatialDefRate[4] * SpatialDefRate[4] +
                                    SpatialDefRate[5] * SpatialDefRate[5]));
  }

  template <>
  void UpdatedLagrangianElement<3>::CalcEquivalentStrainRate(const VectorType &SpatialDefRate,
                                                             double &EquivalentStrainRate)
  {
    EquivalentStrainRate = sqrt(2.0 * (SpatialDefRate[0] * SpatialDefRate[0] +
                                       SpatialDefRate[1] * SpatialDefRate[1] +
                                       SpatialDefRate[2] * SpatialDefRate[2] +
                                       2.0 * SpatialDefRate[3] * SpatialDefRate[3] +
                                       2.0 * SpatialDefRate[4] * SpatialDefRate[4] +
                                       2.0 * SpatialDefRate[5] * SpatialDefRate[5]));
  }

  template <>
  double UpdatedLagrangianElement<2>::CalcNormalProjectionDefRate(const VectorType &SpatialDefRate,
                                                                  const array_1d<double, 3> NormalVector)
  {

    double NormalProjSpatialDefRate = NormalVector[0] * SpatialDefRate[0] * NormalVector[0] +
                                      NormalVector[1] * SpatialDefRate[1] * NormalVector[1] +
                                      2 * NormalVector[0] * SpatialDefRate[2] * NormalVector[1];

    return NormalProjSpatialDefRate;
  }

  template <>
  double UpdatedLagrangianElement<3>::CalcNormalProjectionDefRate(const VectorType &SpatialDefRate,
                                                                  const array_1d<double, 3> NormalVector)
  {

    double NormalProjSpatialDefRate = NormalVector[0] * SpatialDefRate[0] * NormalVector[0] +
                                      NormalVector[1] * SpatialDefRate[1] * NormalVector[1] +
                                      NormalVector[2] * SpatialDefRate[2] * NormalVector[2] +
                                      2 * NormalVector[0] * SpatialDefRate[3] * NormalVector[1] +
                                      2 * NormalVector[0] * SpatialDefRate[4] * NormalVector[2] +
                                      2 * NormalVector[1] * SpatialDefRate[5] * NormalVector[2];

    return NormalProjSpatialDefRate;
  }

  template <>
  double UpdatedLagrangianElement<2>::CalcNormalProjectionDefRate(const VectorType &SpatialDefRate)
  {

    double NormalProjSpatialDefRate = 0;
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    array_1d<double, 3> NormalMean(3, 0.0);

    for (SizeType i = 0; i < (NumNodes - 1); i++)
    {
      for (SizeType j = (i + 1); j < NumNodes; j++)
      {
        if (rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE))
        {
          noalias(NormalMean) += (rGeom[i].FastGetSolutionStepValue(NORMAL) +
                                  rGeom[j].FastGetSolutionStepValue(NORMAL)) *
                                 0.5;
        }
      }
    }

    NormalProjSpatialDefRate = NormalMean[0] * SpatialDefRate[0] * NormalMean[0] +
                               NormalMean[1] * SpatialDefRate[1] * NormalMean[1] +
                               2 * NormalMean[0] * SpatialDefRate[2] * NormalMean[1];

    return NormalProjSpatialDefRate;
  }

  template <>
  double UpdatedLagrangianElement<3>::CalcNormalProjectionDefRate(const VectorType &SpatialDefRate)
  {
    double NormalProjSpatialDefRate = 0;
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    array_1d<double, 3> NormalMean(3, 0.0);

    for (SizeType i = 0; i < (NumNodes - 2); i++)
    {
      for (SizeType j = (i + 1); j < (NumNodes - 1); j++)
      {
        for (SizeType k = (j + 1); k < NumNodes; k++)
        {
          if (rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE) && rGeom[k].Is(FREE_SURFACE))
          {
            noalias(NormalMean) += (rGeom[i].FastGetSolutionStepValue(NORMAL) +
                                    rGeom[j].FastGetSolutionStepValue(NORMAL) +
                                    rGeom[k].FastGetSolutionStepValue(NORMAL)) /
                                   3.0;
          }
        }
      }
    }

    NormalProjSpatialDefRate = NormalMean[0] * SpatialDefRate[0] * NormalMean[0] +
                               NormalMean[1] * SpatialDefRate[1] * NormalMean[1] +
                               NormalMean[2] * SpatialDefRate[2] * NormalMean[2] +
                               2 * NormalMean[0] * SpatialDefRate[3] * NormalMean[1] +
                               2 * NormalMean[0] * SpatialDefRate[4] * NormalMean[2] +
                               2 * NormalMean[1] * SpatialDefRate[5] * NormalMean[2];

    return NormalProjSpatialDefRate;
  }

  template <unsigned int TDim>
  void UpdatedLagrangianElement<TDim>::GetPressureVelocityValues(Vector &rValues,
                                                                 const int Step)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rValues.size() != NumNodes)
      rValues.resize(NumNodes, false);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE_VELOCITY, Step);
    }
  }

  template <unsigned int TDim>
  bool UpdatedLagrangianElement<TDim>::CalcMechanicsUpdated(ElementalVariables &rElementalVariables,
                                                            const ProcessInfo &rCurrentProcessInfo,
                                                            const ShapeFunctionDerivativesType &rDN_DX)
  {

    const double theta = this->GetThetaMomentum();
    bool computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta);
    return computeElement;
  }

  template <>
  bool UpdatedLagrangianElement<2>::CalcCompleteStrainRate(ElementalVariables &rElementalVariables,
                                                           const ProcessInfo &rCurrentProcessInfo,
                                                           const ShapeFunctionDerivativesType &rDN_DX,
                                                           const double theta)
  {
    bool computeElement = true;
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = dimension * NumNodes;
    VectorType NodePosition = ZeroVector(LocalSize);
    VectorType VelocityValues = ZeroVector(LocalSize);
    VectorType RHSVelocities = ZeroVector(LocalSize);
    this->GetPositions(NodePosition, rCurrentProcessInfo, theta);
    this->GetVelocityValues(RHSVelocities, 0);
    RHSVelocities *= theta;
    this->GetVelocityValues(VelocityValues, 1);
    RHSVelocities += VelocityValues * (1.0 - theta);

    rElementalVariables.Fgrad = ZeroMatrix(dimension, dimension);
    rElementalVariables.FgradVel = ZeroMatrix(dimension, dimension);
    for (SizeType i = 0; i < dimension; i++)
    {
      for (SizeType j = 0; j < dimension; j++)
      {
        for (SizeType k = 0; k < NumNodes; k++)
        {
          rElementalVariables.Fgrad(i, j) += NodePosition[dimension * k + i] * rDN_DX(k, j);
          rElementalVariables.FgradVel(i, j) += RHSVelocities[dimension * k + i] * rDN_DX(k, j);
        }
      }
    }

    // Inverse
    rElementalVariables.InvFgrad = ZeroMatrix(dimension, dimension);
    rElementalVariables.DetFgrad = 1;
    MathUtils<double>::InvertMatrix2(rElementalVariables.Fgrad,
                                     rElementalVariables.InvFgrad,
                                     rElementalVariables.DetFgrad);

    // it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
    rElementalVariables.SpatialVelocityGrad.resize(dimension, dimension, false);
    rElementalVariables.SpatialVelocityGrad = prod(rElementalVariables.FgradVel, rElementalVariables.InvFgrad);

    rElementalVariables.VolumetricDefRate = 0;
    for (SizeType i = 0; i < dimension; i++)
    {
      rElementalVariables.VolumetricDefRate += rElementalVariables.SpatialVelocityGrad(i, i);
    }

    rElementalVariables.SpatialDefRate[0] = rElementalVariables.SpatialVelocityGrad(0, 0);
    rElementalVariables.SpatialDefRate[1] = rElementalVariables.SpatialVelocityGrad(1, 1);
    rElementalVariables.SpatialDefRate[2] = 0.5 * (rElementalVariables.SpatialVelocityGrad(1, 0) + rElementalVariables.SpatialVelocityGrad(0, 1));

    constexpr double aThird = 1.0 / 3.0;
    const double dev_X = rElementalVariables.SpatialDefRate[0] -
                         (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1]) * aThird;
    const double dev_Y = rElementalVariables.SpatialDefRate[1] -
                         (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1]) * aThird;
    rElementalVariables.DeviatoricInvariant = sqrt(2 * (dev_X * dev_X + dev_Y * dev_Y +
                                                        rElementalVariables.SpatialDefRate[2] * rElementalVariables.SpatialDefRate[2]));

    rElementalVariables.EquivalentStrainRate = sqrt((2.0 * rElementalVariables.SpatialDefRate[0] * rElementalVariables.SpatialDefRate[0] +
                                                     2.0 * rElementalVariables.SpatialDefRate[1] * rElementalVariables.SpatialDefRate[1] +
                                                     4.0 * rElementalVariables.SpatialDefRate[2] * rElementalVariables.SpatialDefRate[2]));

    return computeElement;
  }

  template <>
  bool UpdatedLagrangianElement<3>::CalcCompleteStrainRate(ElementalVariables &rElementalVariables,
                                                           const ProcessInfo &rCurrentProcessInfo,
                                                           const ShapeFunctionDerivativesType &rDN_DX,
                                                           const double theta)
  {

    bool computeElement = true;
    const unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = dimension * NumNodes;
    VectorType NodePosition = ZeroVector(LocalSize);
    VectorType VelocityValues = ZeroVector(LocalSize);
    VectorType RHSVelocities = ZeroVector(LocalSize);
    this->GetPositions(NodePosition, rCurrentProcessInfo, theta);
    this->GetVelocityValues(RHSVelocities, 0);
    RHSVelocities *= theta;
    this->GetVelocityValues(VelocityValues, 1);
    RHSVelocities += VelocityValues * (1.0 - theta);

    rElementalVariables.Fgrad = ZeroMatrix(dimension, dimension);
    rElementalVariables.FgradVel = ZeroMatrix(dimension, dimension);
    for (SizeType i = 0; i < dimension; i++)
    {
      for (SizeType j = 0; j < dimension; j++)
      {
        for (SizeType k = 0; k < NumNodes; k++)
        {
          rElementalVariables.Fgrad(i, j) += NodePosition[dimension * k + i] * rDN_DX(k, j);
          rElementalVariables.FgradVel(i, j) += RHSVelocities[dimension * k + i] * rDN_DX(k, j);
        }
      }
    }

    // Inverse
    rElementalVariables.InvFgrad = ZeroMatrix(dimension, dimension);
    rElementalVariables.DetFgrad = 1.0;
    MathUtils<double>::InvertMatrix3(rElementalVariables.Fgrad,
                                     rElementalVariables.InvFgrad,
                                     rElementalVariables.DetFgrad);

    // it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
    rElementalVariables.SpatialVelocityGrad.resize(dimension, dimension, false);
    rElementalVariables.SpatialVelocityGrad = prod(rElementalVariables.FgradVel, rElementalVariables.InvFgrad);

    rElementalVariables.VolumetricDefRate = 0.0;
    for (SizeType i = 0; i < dimension; i++)
    {
      rElementalVariables.VolumetricDefRate += rElementalVariables.SpatialVelocityGrad(i, i);
    }

    rElementalVariables.SpatialDefRate[0] = rElementalVariables.SpatialVelocityGrad(0, 0);
    rElementalVariables.SpatialDefRate[1] = rElementalVariables.SpatialVelocityGrad(1, 1);
    rElementalVariables.SpatialDefRate[2] = rElementalVariables.SpatialVelocityGrad(2, 2);
    rElementalVariables.SpatialDefRate[3] = 0.5 * (rElementalVariables.SpatialVelocityGrad(1, 0) + rElementalVariables.SpatialVelocityGrad(0, 1));
    rElementalVariables.SpatialDefRate[4] = 0.5 * (rElementalVariables.SpatialVelocityGrad(2, 0) + rElementalVariables.SpatialVelocityGrad(0, 2));
    rElementalVariables.SpatialDefRate[5] = 0.5 * (rElementalVariables.SpatialVelocityGrad(2, 1) + rElementalVariables.SpatialVelocityGrad(1, 2));

    constexpr double aThird = 1.0 / 3.0;
    const double dev_X = rElementalVariables.SpatialDefRate[0] -
                         (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1] + rElementalVariables.SpatialDefRate[2]) * aThird;
    const double dev_Y = rElementalVariables.SpatialDefRate[1] -
                         (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1] + rElementalVariables.SpatialDefRate[2]) * aThird;
    const double dev_Z = rElementalVariables.SpatialDefRate[2] -
                         (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1] + rElementalVariables.SpatialDefRate[2]) * aThird;
    rElementalVariables.DeviatoricInvariant = sqrt(2 * (dev_X * dev_X + dev_Y * dev_Y + dev_Z * dev_Z +
                                                        rElementalVariables.SpatialDefRate[3] * rElementalVariables.SpatialDefRate[3] +
                                                        rElementalVariables.SpatialDefRate[4] * rElementalVariables.SpatialDefRate[4] +
                                                        rElementalVariables.SpatialDefRate[5] * rElementalVariables.SpatialDefRate[5]));

    rElementalVariables.EquivalentStrainRate = sqrt(2.0 * (rElementalVariables.SpatialDefRate[0] * rElementalVariables.SpatialDefRate[0] +
                                                           rElementalVariables.SpatialDefRate[1] * rElementalVariables.SpatialDefRate[1] +
                                                           rElementalVariables.SpatialDefRate[2] * rElementalVariables.SpatialDefRate[2] +
                                                           2.0 * rElementalVariables.SpatialDefRate[3] * rElementalVariables.SpatialDefRate[3] +
                                                           2.0 * rElementalVariables.SpatialDefRate[4] * rElementalVariables.SpatialDefRate[4] +
                                                           2.0 * rElementalVariables.SpatialDefRate[5] * rElementalVariables.SpatialDefRate[5]));

    return computeElement;
  }

  template <unsigned int TDim>
  bool UpdatedLagrangianElement<TDim>::CalcStrainRate(ElementalVariables &rElementalVariables,
                                                      const ProcessInfo &rCurrentProcessInfo,
                                                      const ShapeFunctionDerivativesType &rDN_DX,
                                                      const double theta)
  {

    bool computeElement = true;

    this->CalcFGrad(rDN_DX,
                    rElementalVariables.Fgrad,
                    rElementalVariables.InvFgrad,
                    rElementalVariables.DetFgrad,
                    rCurrentProcessInfo,
                    theta);

    // it computes the material time derivative of the deformation gradient and its jacobian and inverse
    this->CalcVelDefGrad(rDN_DX,
                         rElementalVariables.FgradVel,
                         theta);

    // it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
    this->CalcSpatialVelocityGrad(rElementalVariables.InvFgrad,
                                  rElementalVariables.FgradVel,
                                  rElementalVariables.SpatialVelocityGrad);

    this->CalcVolDefRateFromSpatialVelGrad(rElementalVariables.VolumetricDefRate,
                                           rElementalVariables.SpatialVelocityGrad);

    // it computes Material time Derivative of Green Lagrange strain tensor in MATERIAL configuration --> [D(E)/Dt]
    this->CalcMDGreenLagrangeMaterial(rElementalVariables.Fgrad,
                                      rElementalVariables.FgradVel,
                                      rElementalVariables.MDGreenLagrangeMaterial);

    // it computes Material time Derivative of Green Lagrange strain tensor in SPATIAL configuration  --> [d]
    this->CalcSpatialDefRate(rElementalVariables.MDGreenLagrangeMaterial,
                             rElementalVariables.InvFgrad,
                             rElementalVariables.SpatialDefRate);

    this->CalcDeviatoricInvariant(rElementalVariables.SpatialDefRate,
                                  rElementalVariables.DeviatoricInvariant);

    this->CalcEquivalentStrainRate(rElementalVariables.SpatialDefRate,
                                   rElementalVariables.EquivalentStrainRate);
    return computeElement;
  }

  template <unsigned int TDim>
  double UpdatedLagrangianElement<TDim>::EquivalentStrainRate(const ShapeFunctionDerivativesType &rDN_DX) const
  {
    const GeometryType &rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    // Calculate Symetric gradient
    MatrixType S = ZeroMatrix(TDim, TDim);
    for (unsigned int n = 0; n < NumNodes; ++n)
    {
      const array_1d<double, 3> &rVel = rGeom[n].FastGetSolutionStepValue(VELOCITY, 1); // OLD VELOCITY (which is incompressible, unlike the fractional step one)
      for (unsigned int i = 0; i < TDim; ++i)
        for (unsigned int j = 0; j < TDim; ++j)
          S(i, j) += 0.5 * (rDN_DX(n, j) * rVel[i] + rDN_DX(n, i) * rVel[j]);
    }

    // Norm of symetric gradient
    double NormS = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
      for (unsigned int j = 0; j < TDim; ++j)
        NormS += S(i, j) * S(i, j);

    return std::sqrt(2.0 * NormS);
  }

  template <>
  void UpdatedLagrangianElement<2>::ComputeMechanicalDissipation(ElementalVariables &rElementalVariables)
  {
    KRATOS_TRY;
    const double volumetric_strain = (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1]) * 0.5;

    const double mechanical_dissipation = rElementalVariables.UpdatedDeviatoricCauchyStress[0] * (rElementalVariables.SpatialDefRate[0] - volumetric_strain) +
                                          rElementalVariables.UpdatedDeviatoricCauchyStress[1] * (rElementalVariables.SpatialDefRate[1] - volumetric_strain) +
                                          2.0 * rElementalVariables.UpdatedDeviatoricCauchyStress[2] * rElementalVariables.SpatialDefRate[2];

    this->SetValue(MECHANICAL_DISSIPATION, mechanical_dissipation);
    KRATOS_CATCH("");
  }

  template <>
  void UpdatedLagrangianElement<3>::ComputeMechanicalDissipation(ElementalVariables &rElementalVariables)
  {
    KRATOS_TRY;
    const double volumetric_strain = (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1] + rElementalVariables.SpatialDefRate[2]) / 3.0;

    const double mechanical_dissipation = rElementalVariables.UpdatedDeviatoricCauchyStress[0] * (rElementalVariables.SpatialDefRate[0] - volumetric_strain) +
                                          rElementalVariables.UpdatedDeviatoricCauchyStress[1] * (rElementalVariables.SpatialDefRate[1] - volumetric_strain) +
                                          rElementalVariables.UpdatedDeviatoricCauchyStress[2] * (rElementalVariables.SpatialDefRate[2] - volumetric_strain) +
                                          2.0 * rElementalVariables.UpdatedDeviatoricCauchyStress[3] * rElementalVariables.SpatialDefRate[3] +
                                          2.0 * rElementalVariables.UpdatedDeviatoricCauchyStress[4] * rElementalVariables.SpatialDefRate[4] +
                                          2.0 * rElementalVariables.UpdatedDeviatoricCauchyStress[5] * rElementalVariables.SpatialDefRate[5];

    this->SetValue(MECHANICAL_DISSIPATION, mechanical_dissipation);
    KRATOS_CATCH("");
  }
  /*
   * Template class definition (this should allow us to compile the desired template instantiations)
   */

  template class UpdatedLagrangianElement<2>;
  template class UpdatedLagrangianElement<3>;

} // namespace Kratos

//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:               June 2021 $
//   Revision:            $Revision:                 0.0 $
//
//   Implementation of the Gauss-Seidel two step Updated Lagrangian Velocity-Pressure element
//     ( There is a ScalingConstant to multiply the mass balance equation for a number because i read it somewhere)
//

// System includes

// External includes

// Project includes
#include "custom_elements/three_step_updated_lagrangian_element.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

  template <unsigned int TDim>
  Element::Pointer ThreeStepUpdatedLagrangianElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
  {
    KRATOS_TRY;

    ThreeStepUpdatedLagrangianElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());
    return Element::Pointer(new ThreeStepUpdatedLagrangianElement(NewElement));

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                                                     VectorType &rRightHandSideVector,
                                                                     const ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;

    switch (rCurrentProcessInfo[FRACTIONAL_STEP])
    {
    case 1:
    {
      this->CalculateFirstVelocitySystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
      break;
    }
    case 5:
    {
      this->CalculateLocalPressureSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
      break;
    }
    case 6:
    {
      this->CalculateLastVelocitySystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
      break;
    }

    default:
    {
      KRATOS_THROW_ERROR(std::logic_error, "Unexpected value for THREE_STEP_UPDATED_LAGRANGIAN_V_P_ELEMENT index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
    }
    }

    KRATOS_CATCH("");
  }
  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::EquationIdVector(EquationIdVectorType &rResult,
                                                                 const ProcessInfo &rCurrentProcessInfo) const
  {
    KRATOS_TRY;
    switch (rCurrentProcessInfo[FRACTIONAL_STEP])
    {
    case 1:
    {
      this->VelocityEquationIdVector(rResult, rCurrentProcessInfo);
      break;
    }
    case 5:
    {
      this->PressureEquationIdVector(rResult, rCurrentProcessInfo);
      break;
    }
    case 6:
    {
      this->VelocityEquationIdVector(rResult, rCurrentProcessInfo);
      break;
    }
    default:
    {
      KRATOS_THROW_ERROR(std::logic_error, "Unexpected value for FRACTIONAL_STEP index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
    }
    }

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::GetDofList(DofsVectorType &rElementalDofList,
                                                           const ProcessInfo &rCurrentProcessInfo) const
  {
    KRATOS_TRY;
    switch (rCurrentProcessInfo[FRACTIONAL_STEP])
    {
    case 1:
    {
      this->GetVelocityDofList(rElementalDofList, rCurrentProcessInfo);
      break;
    }
    case 5:
    {
      this->GetPressureDofList(rElementalDofList, rCurrentProcessInfo);
      break;
    }
    case 6:
    {
      this->GetVelocityDofList(rElementalDofList, rCurrentProcessInfo);
      break;
    }
    default:
    {
      KRATOS_THROW_ERROR(std::logic_error, "Unexpected value for FRACTIONAL_STEP index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
    }
    }

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::CalculateFirstVelocitySystem(MatrixType &rLeftHandSideMatrix,
                                                                             VectorType &rRightHandSideVector,
                                                                             const ProcessInfo &rCurrentProcessInfo)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = TDim * NumNodes;

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != LocalSize)
      rLeftHandSideMatrix.resize(LocalSize, LocalSize);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);

    if (rRightHandSideVector.size() != LocalSize)
      rRightHandSideVector.resize(LocalSize);

    rRightHandSideVector = ZeroVector(LocalSize);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    MatrixType MassMatrix = ZeroMatrix(LocalSize, LocalSize);

    // const double eta = rCurrentProcessInfo[FS_PRESSURE_GRADIENT_RELAXATION_FACTOR];  //!!!!!!!!!!!!!!!!!!!
    const double eta = 0.0;

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
    {
      const double GaussWeight = GaussWeights[g];
      const ShapeFunctionsType &N = row(NContainer, g);
      const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

      // Evaluate required variables at the integration point
      double Density = 0.0;
      double Viscosity = 0.0;
      array_1d<double, 3> BodyForce = ZeroVector(3);

      this->EvaluateInPoint(Density, DENSITY, N);
      this->EvaluateInPoint(BodyForce, BODY_FORCE, N);

      // Evaluate the pressure and pressure gradient at this point (for the G * P_n term)
      double OldPressure = 0.0;
      this->EvaluateInPoint(OldPressure, PRESSURE, N, 0);

      this->EvaluateInPoint(Viscosity, DYNAMIC_VISCOSITY, N);

      // Add integration point contribution to the local mass matrix
      const double MassCoeff = Viscosity * GaussWeight;
      this->AddMomentumMassTerm(MassMatrix, N, MassCoeff);

      double etaOldPressure = OldPressure * eta;
      // Add RHS contributions to the local system equation
      this->AddMomentumSystemTerms(rLeftHandSideMatrix, rRightHandSideVector, Density, BodyForce, etaOldPressure, N, rDN_DX, GaussWeight);

      // Add viscous term
      const double ViscousCoeff = Viscosity * GaussWeight;
      this->AddViscousTerm(rLeftHandSideMatrix, rDN_DX, ViscousCoeff);
    }

    // Add residual of previous iteration to RHS
    VectorType LastValues = ZeroVector(LocalSize);
    this->GetVelocityValues(LastValues, 0);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, LastValues);

    // Add dynamic term
    const double TimeStep = rCurrentProcessInfo[DELTA_TIME];
    VectorType VelocityValues = ZeroVector(LocalSize);
    VectorType AccelerationValues = ZeroVector(LocalSize);
    this->GetAccelerationValues(AccelerationValues, 0);
    this->GetVelocityValues(VelocityValues, 0);
    noalias(AccelerationValues) += -2.0 * VelocityValues / TimeStep;
    this->GetVelocityValues(VelocityValues, 1);
    noalias(AccelerationValues) += 2.0 * VelocityValues / TimeStep; //these are negative accelerations
    noalias(rRightHandSideVector) += prod(MassMatrix, AccelerationValues);
    noalias(rLeftHandSideMatrix) += MassMatrix * 2 / TimeStep;
    // // using BDF coefficients
    // const Vector &rBDFCoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
    // noalias(rLeftHandSideMatrix) += rBDFCoeffs[0] * MassMatrix;
    // VectorType TimeTerm = rBDFCoeffs[0] * LastValues;
    // for (SizeType i = 1; i < rBDFCoeffs.size(); i++)
    // {
    //   this->GetVelocityValues(LastValues, i);
    //   noalias(TimeTerm) += rBDFCoeffs[i] * LastValues;
    // }
    // noalias(rRightHandSideVector) -= prod(MassMatrix, TimeTerm);
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::AddMomentumSystemTerms(Matrix &rLHSMatrix,
                                                                       Vector &rRHSVector,
                                                                       const double Density,
                                                                       const array_1d<double, 3> &rBodyForce,
                                                                       const double OldPressure,
                                                                       const ShapeFunctionsType &rN,
                                                                       const ShapeFunctionDerivativesType &rDN_DX,
                                                                       const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    SizeType FirstRow = 0;

 
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      // Build RHS

        for (SizeType d = 0; d < TDim; ++d)
        {
          // Body force
          double RHSi = Density * rN[i] * rBodyForce[d];
          // Pressure gradient (integrated by parts)
          RHSi += rDN_DX(i, d) * OldPressure;
          rRHSVector[FirstRow + d] += Weight * RHSi;
        }
      FirstRow += TDim;
    }
  }


  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::AddExternalForces(Vector &rRHSVector,
                                                                  const double Density,
                                                                  const ShapeFunctionsType &rN,
                                                                  const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    SizeType FirstRow = 0;

    array_1d<double, 3> VolumeAcceleration(3, 0.0);
    this->EvaluateInPoint(VolumeAcceleration, VOLUME_ACCELERATION, rN);
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      if (this->GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION))
      {
        for (SizeType d = 0; d < TDim; ++d)
        {
          // Volume Acceleration
          rRHSVector[FirstRow + d] += Weight * Density * rN[i] * VolumeAcceleration[d];
        }
      }
      FirstRow += TDim;
    }
  }

  template <>
  void ThreeStepUpdatedLagrangianElement<2>::AddViscousTerm(MatrixType &rDampingMatrix,
                                                            const ShapeFunctionDerivativesType &rShapeDeriv,
                                                            const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    SizeType FirstRow(0), FirstCol(0);

    for (SizeType j = 0; j < NumNodes; ++j)
    {
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        // First Row
        rDampingMatrix(FirstRow, FirstCol) += Weight * (FourThirds * rShapeDeriv(i, 0) * rShapeDeriv(j, 0) + rShapeDeriv(i, 1) * rShapeDeriv(j, 1));
        rDampingMatrix(FirstRow, FirstCol + 1) += Weight * (nTwoThirds * rShapeDeriv(i, 0) * rShapeDeriv(j, 1) + rShapeDeriv(i, 1) * rShapeDeriv(j, 0));

        // Second Row
        rDampingMatrix(FirstRow + 1, FirstCol) += Weight * (nTwoThirds * rShapeDeriv(i, 1) * rShapeDeriv(j, 0) + rShapeDeriv(i, 0) * rShapeDeriv(j, 1));
        rDampingMatrix(FirstRow + 1, FirstCol + 1) += Weight * (FourThirds * rShapeDeriv(i, 1) * rShapeDeriv(j, 1) + rShapeDeriv(i, 0) * rShapeDeriv(j, 0));

        // Update Counter
        FirstRow += 2;
      }
      FirstRow = 0;
      FirstCol += 2;
    }
  }

  template <>
  void ThreeStepUpdatedLagrangianElement<3>::AddViscousTerm(MatrixType &rDampingMatrix,
                                                            const ShapeFunctionDerivativesType &rShapeDeriv,
                                                            const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    const double OneThird = 1.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    unsigned int FirstRow(0), FirstCol(0);

    for (SizeType j = 0; j < NumNodes; ++j)
    {
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        // (dN_i/dx_k dN_j/dx_k)
        const double Diag = rShapeDeriv(i, 0) * rShapeDeriv(j, 0) + rShapeDeriv(i, 1) * rShapeDeriv(j, 1) + rShapeDeriv(i, 2) * rShapeDeriv(j, 2);

        // First Row
        rDampingMatrix(FirstRow, FirstCol) += Weight * (OneThird * rShapeDeriv(i, 0) * rShapeDeriv(j, 0) + Diag);
        rDampingMatrix(FirstRow, FirstCol + 1) += Weight * (nTwoThirds * rShapeDeriv(i, 0) * rShapeDeriv(j, 1) + rShapeDeriv(i, 1) * rShapeDeriv(j, 0));
        rDampingMatrix(FirstRow, FirstCol + 2) += Weight * (nTwoThirds * rShapeDeriv(i, 0) * rShapeDeriv(j, 2) + rShapeDeriv(i, 2) * rShapeDeriv(j, 0));

        // Second Row
        rDampingMatrix(FirstRow + 1, FirstCol) += Weight * (nTwoThirds * rShapeDeriv(i, 1) * rShapeDeriv(j, 0) + rShapeDeriv(i, 0) * rShapeDeriv(j, 1));
        rDampingMatrix(FirstRow + 1, FirstCol + 1) += Weight * (OneThird * rShapeDeriv(i, 1) * rShapeDeriv(j, 1) + Diag);
        rDampingMatrix(FirstRow + 1, FirstCol + 2) += Weight * (nTwoThirds * rShapeDeriv(i, 1) * rShapeDeriv(j, 2) + rShapeDeriv(i, 2) * rShapeDeriv(j, 1));

        // Third Row
        rDampingMatrix(FirstRow + 2, FirstCol) += Weight * (nTwoThirds * rShapeDeriv(i, 2) * rShapeDeriv(j, 0) + rShapeDeriv(i, 0) * rShapeDeriv(j, 2));
        rDampingMatrix(FirstRow + 2, FirstCol + 1) += Weight * (nTwoThirds * rShapeDeriv(i, 2) * rShapeDeriv(j, 1) + rShapeDeriv(i, 1) * rShapeDeriv(j, 2));
        rDampingMatrix(FirstRow + 2, FirstCol + 2) += Weight * (OneThird * rShapeDeriv(i, 2) * rShapeDeriv(j, 2) + Diag);

        // Update Counter
        FirstRow += 3;
      }
      FirstRow = 0;
      FirstCol += 3;
    }
  }

  template <unsigned int TDim>
  double ThreeStepUpdatedLagrangianElement<TDim>::ElementSize()
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    // calculate minimum element length (used in stabilization Tau)
    array_1d<double, 3> Edge = rGeom[1].Coordinates() - rGeom[0].Coordinates();
    double ElemSize = Edge[0] * Edge[0];
    for (SizeType d = 1; d < TDim; d++)
      ElemSize += Edge[d] * Edge[d];

    for (SizeType i = 2; i < NumNodes; i++)
      for (SizeType j = 0; j < i; j++)
      {
        Edge = rGeom[i].Coordinates() - rGeom[j].Coordinates();
        double Length = Edge[0] * Edge[0];
        for (SizeType d = 1; d < TDim; d++)
          Length += Edge[d] * Edge[d];
        if (Length < ElemSize)
          ElemSize = Length;
      }
    return sqrt(ElemSize);
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::ComputeLumpedMassMatrix(Matrix &rMassMatrix,
                                                                        const double Weight,
                                                                        double &MeanValue)
  {

    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    double Count = 0;

    if ((NumNodes == 3 && TDim == 2) || (NumNodes == 4 && TDim == 3))
    {
      double Coeff = 1.0 + TDim;
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        double Mij = Weight / Coeff;

        for (unsigned int j = 0; j < TDim; j++)
        {
          unsigned int index = i * TDim + j;
          rMassMatrix(index, index) += Mij;
          Count += 1.0;
          MeanValue += Mij;
        }
      }
    }
    else if (NumNodes == 6 && TDim == 2)
    {
      double Mij = Weight / 57.0;
      double consistent = 1.0;
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        if (i < 3)
        {
          consistent = 3.0;
        }
        else
        {
          consistent = 16.0;
        }
        for (unsigned int j = 0; j < TDim; j++)
        {
          unsigned int index = i * TDim + j;
          rMassMatrix(index, index) += Mij * consistent;
          Count += 1.0;
          MeanValue += Mij;
        }
      }
    }
    else
    {
      std::cout << "ComputeLumpedMassMatrix 3D quadratic not yet implemented!" << std::endl;
    }
    MeanValue *= 1.0 / Count;
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::AddMomentumMassTerm(Matrix &rMassMatrix,
                                                                    const ShapeFunctionsType &rN,
                                                                    const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    IndexType FirstRow = 0;
    IndexType FirstCol = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      for (SizeType j = 0; j < NumNodes; ++j)
      {
        const double Mij = Weight * rN[i] * rN[j];
        for (SizeType d = 0; d < TDim; ++d)
          rMassMatrix(FirstRow + d, FirstCol + d) += Mij;
        FirstCol += TDim;
      }
      FirstRow += TDim;
      FirstCol = 0;
    }
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::AddExternalForces(Vector &rRHSVector,
                                                                  const double Density,
                                                                  const ShapeFunctionsType &rN,
                                                                  const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    SizeType FirstRow = 0;

    array_1d<double, 3> VolumeAcceleration(3, 0.0);

    this->EvaluateInPoint(VolumeAcceleration, VOLUME_ACCELERATION, rN);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      if (this->GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION))
      {
        for (SizeType d = 0; d < TDim; ++d)
        {
          // Volume Acceleration
          rRHSVector[FirstRow + d] += Weight * Density * rN[i] * VolumeAcceleration[d];
        }
      }
      FirstRow += TDim;
    }
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::CalculateLocalPressureSystem(MatrixType &rLeftHandSideMatrix,
                                                                             VectorType &rRightHandSideVector,
                                                                             const ProcessInfo &rCurrentProcessInfo)
  {
    // GeometryType &rGeom = this->GetGeometry();
    // const SizeType NumNodes = rGeom.PointsNumber();

    // // Check sizes and initialize
    // if (rLeftHandSideMatrix.size1() != NumNodes)
    //   rLeftHandSideMatrix.resize(NumNodes, NumNodes);

    // rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);

    // if (rRightHandSideVector.size() != NumNodes)
    //   rRightHandSideVector.resize(NumNodes);

    // rRightHandSideVector = ZeroVector(NumNodes);

    // // Shape functions and integration points
    // ShapeFunctionDerivativesArrayType DN_DX;
    // Matrix NContainer;
    // VectorType GaussWeights;
    // this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    // const unsigned int NumGauss = GaussWeights.size();

    // // Stabilization parameters
    // double ElemSize = this->ElementSize();
    // double TauOne;
    // double TauTwo;

    // const double eta = rCurrentProcessInfo[FS_PRESSURE_GRADIENT_RELAXATION_FACTOR];

    // // Loop on integration points
    // for (unsigned int g = 0; g < NumGauss; g++)
    // {
    //   const double GaussWeight = GaussWeights[g];
    //   const ShapeFunctionsType &N = row(NContainer, g);
    //   const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

    //   // Evaluate required variables at the integration point
    //   double Density;
    //   array_1d<double, 3> BodyForce = ZeroVector(3);
    //   array_1d<double, 3> MomentumProjection = ZeroVector(3);

    //   this->EvaluateInPoint(Density, DENSITY, N);
    //   this->EvaluateInPoint(BodyForce, BODY_FORCE, N);
    //   this->EvaluateInPoint(MomentumProjection, PRESS_PROJ, N);

    //   //        // Evaluate the pressure and pressure gradient at this point (for the G * P_n term)
    //   array_1d<double, TDim> OldPressureGradient = ZeroVector(TDim);
    //   this->EvaluateGradientInPoint(OldPressureGradient, PRESSURE, rDN_DX);

    //   // Stabilization parameters
    //   array_1d<double, 3> ConvVel = ZeroVector(3);
    //   this->EvaluateConvVelocity(ConvVel, N);
    //   double Viscosity = this->EffectiveViscosity(Density, N, rDN_DX, ElemSize, rCurrentProcessInfo);
    //   this->CalculateTau(TauOne, TauTwo, ElemSize, ConvVel, Density, Viscosity, rCurrentProcessInfo);

    //   double DivU;
    //   this->EvaluateDivergenceInPoint(DivU, VELOCITY, rDN_DX);

    //   // constant coefficient multiplying the pressure Laplacian (See Codina, Badia 2006 paper for details in case of a BDF2 time scheme)
    //   const double LaplacianCoeff = 1.0 / (Density * rCurrentProcessInfo[BDF_COEFFICIENTS][0]);

    //   // Add convection, stabilization and RHS contributions to the local system equation
    //   for (SizeType i = 0; i < NumNodes; ++i)
    //   {
    //     // LHS contribution
    //     for (SizeType j = 0; j < NumNodes; ++j)
    //     {
    //       double Lij = 0.0;
    //       for (SizeType d = 0; d < TDim; ++d)
    //         Lij += rDN_DX(i, d) * rDN_DX(j, d);
    //       Lij *= (LaplacianCoeff + TauOne);

    //       rLeftHandSideMatrix(i, j) += GaussWeight * Lij;
    //     }

    //     // RHS contribution

    //     // Velocity divergence
    //     double RHSi = -N[i] * DivU;

    //     for (SizeType d = 0; d < TDim; ++d)
    //     {
    //       // Momentum stabilization
    //       RHSi += rDN_DX(i, d) * TauOne * (Density * (BodyForce[d] /* - Conv*/) - OldPressureGradient[d] - MomentumProjection[d]);
    //       RHSi += (eta - 1.0) * LaplacianCoeff * rDN_DX(i, d) * OldPressureGradient[d];
    //     }

    //     rRightHandSideVector[i] += GaussWeight * RHSi;
    //   }
    // }
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::CalculateLastVelocitySystem(MatrixType &rLeftHandSideMatrix,
                                                                            VectorType &rRightHandSideVector,
                                                                            const ProcessInfo &rCurrentProcessInfo) //useless, it is done in Calculate called from the strategy
  {
    // const GeometryType &rGeom = this->GetGeometry();
    // const SizeType NumNodes = rGeom.PointsNumber();
    // const SizeType LocalSize = TDim * NumNodes;

    // // Check sizes and initialize
    // if (rLeftHandSideMatrix.size1() != LocalSize)
    //   rLeftHandSideMatrix.resize(LocalSize, LocalSize);

    // rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);

    // if (rRightHandSideVector.size() != LocalSize)
    //   rRightHandSideVector.resize(LocalSize);

    // rRightHandSideVector = ZeroVector(LocalSize);

    // // Shape functions and integration points
    // ShapeFunctionDerivativesArrayType DN_DX;
    // Matrix NContainer;
    // VectorType GaussWeights;
    // this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    // const unsigned int NumGauss = GaussWeights.size();

    // // Stabilization parameters
    // double ElemSize = this->ElementSize();

    // double Bdf0 = rCurrentProcessInfo[BDF_COEFFICIENTS][0]; // for BDF2, 3/(2dt) /// !!!!!!!! play attention to this!

    // // Loop on integration points
    // for (unsigned int g = 0; g < NumGauss; g++)
    // {
    //   const double GaussWeight = GaussWeights[g];
    //   const ShapeFunctionsType &N = row(NContainer, g);
    //   const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

    //   // Evaluate required variables at the integration point
    //   double Density;
    //   this->EvaluateInPoint(Density, DENSITY, N);

    //   const double Coeff = GaussWeight * Density * Bdf0;

    //   // Calculate contribution to the gradient term (RHS)
    //   double DeltaPressure;
    //   this->EvaluateInPoint(DeltaPressure, PRESSURE_OLD_IT, N); /// !!!!!!!! play attention to this!

    //   SizeType RowIndex = 0;

    //   for (SizeType i = 0; i < NumNodes; ++i)
    //   {
    //     for (SizeType d = 0; d < TDim; ++d)
    //     {
    //       rRightHandSideVector[RowIndex++] += GaussWeight * rDN_DX(i, d) * DeltaPressure;
    //     }
    //   }

    //   // LHS matrix includes mass and viscous terms
    //   double Viscosity = 0;
    //   this->EvaluateInPoint(Viscosity, DYNAMIC_VISCOSITY, N);

    //   this->AddMomentumMassTerm(rLeftHandSideMatrix, N, Coeff);

    //   // Add viscous term
    //   const double ViscousCoeff = Viscosity * GaussWeight;
    //   this->AddViscousTerm(rLeftHandSideMatrix, rDN_DX, ViscousCoeff);
    // }
  }

  /*
   * Template class definition (this should allow us to compile the desired template instantiations)
   */

  template class ThreeStepUpdatedLagrangianElement<2>;
  template class ThreeStepUpdatedLagrangianElement<3>;

} // namespace Kratos

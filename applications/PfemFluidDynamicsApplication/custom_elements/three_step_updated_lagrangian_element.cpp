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
      //this->CalculateLocalPressureSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
      //CalculatePSPGLocalContinuityEqForPressure(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
      CalculateFICLocalContinuityEqForPressure(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
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
      this->EvaluateInPoint(BodyForce, VOLUME_ACCELERATION, N);

      // Evaluate the pressure and pressure gradient at this point (for the G * P_n term)
      double OldPressure = 0.0;
      this->EvaluateInPoint(OldPressure, PRESSURE, N, 0);

      this->EvaluateInPoint(Viscosity, DYNAMIC_VISCOSITY, N);

      // Add integration point contribution to the local mass matrix
      const double MassCoeff = Density * GaussWeight;
      AddMomentumMassTerm(MassMatrix, N, MassCoeff);
      //ComputeLumpedMassMatrix(MassMatrix, MassCoeff);

      double etaOldPressure = OldPressure * eta;
      // Add RHS contributions to the local system equation
      AddMomentumRHSTerms(rRightHandSideVector, Density, BodyForce, etaOldPressure, N, rDN_DX, GaussWeight);
      //AddExternalForces(rRightHandSideVector, Density, N, GaussWeight);

      // Add viscous term
      const double ViscousCoeff = Viscosity * GaussWeight;
      this->AddViscousTerm(rLeftHandSideMatrix, rDN_DX, ViscousCoeff);
    }

    // Add residual of previous iteration to RHS
    VectorType VelocityValues = ZeroVector(LocalSize);
    this->GetVelocityValues(VelocityValues, 0);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, VelocityValues);

    // Add dynamic term
    const double TimeStep = rCurrentProcessInfo[DELTA_TIME];
    VectorType AccelerationValues = ZeroVector(LocalSize);
    this->GetAccelerationValues(AccelerationValues, 0);
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
  void ThreeStepUpdatedLagrangianElement<TDim>::Calculate(const Variable<array_1d<double, 3>> &rVariable,
                                                          array_1d<double, 3> &rOutput,
                                                          const ProcessInfo &rCurrentProcessInfo)
  {
    if (rVariable == VELOCITY)
    {
      GeometryType &rGeom = this->GetGeometry();
      const SizeType NumNodes = rGeom.PointsNumber();
      const SizeType LocalSize = TDim * NumNodes;

      // Shape functions and integration points
      ShapeFunctionDerivativesArrayType DN_DX;
      Matrix NContainer;
      VectorType GaussWeights;
      this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
      const unsigned int NumGauss = GaussWeights.size();

      VectorType NodalVelCorrection = ZeroVector(LocalSize);

      // Loop on integration points
      for (unsigned int g = 0; g < NumGauss; ++g)
      {
        const ShapeFunctionsType &N = row(NContainer, g);
        const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

        double Density = 0;
        this->EvaluateInPoint(Density, DENSITY, N);

        const double TimeStep = rCurrentProcessInfo[DELTA_TIME];
        // double timeFactor = 1.0 / rCurrentProcessInfo[BDF_COEFFICIENTS][0];
        double timeFactor = 0.5 * TimeStep;
        const double Coeff = GaussWeights[g] * timeFactor / Density;

        // Calculate contribution to the gradient term (RHS)
        // double DeltaPressure;
        // this->EvaluateInPoint(DeltaPressure, PRESSURE_OLD_IT, N);
        double elementalPressure = 0;
        this->EvaluateInPoint(elementalPressure, PRESSURE, N);

        SizeType RowIndex = 0;

        for (SizeType i = 0; i < NumNodes; ++i)
        {
          for (SizeType d = 0; d < TDim; ++d)
          {
            // NodalVelCorrection[RowIndex++] += Coeff * rDN_DX(i, d) * DeltaPressure;
            NodalVelCorrection[RowIndex++] += Coeff * rDN_DX(i, d) * elementalPressure;
          }
        }
      }

      SizeType Index = 0;

      for (SizeType i = 0; i < NumNodes; ++i)
      {
        rGeom[i].SetLock(); // So it is safe to write in the node in OpenMP
        array_1d<double, 3> &rTemp = rGeom[i].FastGetSolutionStepValue(FRACT_VEL);
        for (SizeType d = 0; d < TDim; ++d)
        {
          rTemp[d] += NodalVelCorrection[Index++];
        }
        rGeom[i].UnSetLock(); // Free the node for other threads
      }
    }
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::AddMomentumRHSTerms(Vector &rRHSVector,
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
                                                                        const double Weight)
  {

    const SizeType NumNodes = this->GetGeometry().PointsNumber();

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
        }
      }
    }
    else
    {
      std::cout << "ComputeLumpedMassMatrix 3D quadratic not yet implemented!" << std::endl;
    }
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
  void ThreeStepUpdatedLagrangianElement<TDim>::CalculateLocalPressureSystem(MatrixType &rLeftHandSideMatrix,
                                                                             VectorType &rRightHandSideVector,
                                                                             const ProcessInfo &rCurrentProcessInfo)
  {
    GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const double TimeStep = rCurrentProcessInfo[DELTA_TIME];

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != NumNodes)
      rLeftHandSideMatrix.resize(NumNodes, NumNodes);

    rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);

    if (rRightHandSideVector.size() != NumNodes)
      rRightHandSideVector.resize(NumNodes);

    rRightHandSideVector = ZeroVector(NumNodes);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);
    // Stabilization parameters
    // double ElemSize = this->ElementSize();
    // double TauOne = 0;
    // double TauTwo = 0;

    //const double eta = rCurrentProcessInfo[FS_PRESSURE_GRADIENT_RELAXATION_FACTOR]; // !!!!!!!!!!!!!!!!!!!!!!!!!!!
    const double eta = 0;
    double VolumetricCoeff = 0;
    double lumpedBulkCoeff = 0;
    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
    {
      const double GaussWeight = GaussWeights[g];
      const ShapeFunctionsType &N = row(NContainer, g);
      const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

      // double theta = 1.0;
      // bool computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta);

      // Evaluate required variables at the integration point
      double Density;
      array_1d<double, 3> BodyForce = ZeroVector(3);

      this->EvaluateInPoint(Density, DENSITY, N);
      this->EvaluateInPoint(BodyForce, VOLUME_ACCELERATION, N);

      // array_1d<double, 3> MomentumProjection = ZeroVector(3);
      // this->EvaluateInPoint(MomentumProjection, PRESS_PROJ, N);

      // // Evaluate the pressure and pressure gradient at this point (for the G * P_n term)
      array_1d<double, TDim> OldPressureGradient = ZeroVector(TDim);
      this->EvaluateGradientInPoint(OldPressureGradient, PRESSURE, rDN_DX);

      // Stabilization parameters
      // array_1d<double, 3> ConvVel = ZeroVector(3);

      double Viscosity = 0;
      this->EvaluateInPoint(Viscosity, DYNAMIC_VISCOSITY, N);

      this->EvaluateInPoint(VolumetricCoeff, BULK_MODULUS, N);
      lumpedBulkCoeff = GaussWeight / (VolumetricCoeff * TimeStep);
      // this->CalculateTau(TauOne, TauTwo, ElemSize, ConvVel, Density, Viscosity, rCurrentProcessInfo);

      double DivU = 0;
      this->EvaluateDivergenceInPoint(DivU, VELOCITY, rDN_DX);

      // // constant coefficient multiplying the pressure Laplacian (See Codina, Badia 2006 paper for details in case of a BDF2 time scheme)
      // const double LaplacianCoeff = 1.0 / (Density * rCurrentProcessInfo[BDF_COEFFICIENTS][0]);
      const double LaplacianCoeff = 2.0 * TimeStep / Density;

      // Add convection, stabilization and RHS contributions to the local system equation
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        // LHS contribution
        for (SizeType j = 0; j < NumNodes; ++j)
        {
          double Lij = 0.0;
          for (SizeType d = 0; d < TDim; ++d)
            Lij += rDN_DX(i, d) * rDN_DX(j, d);
          // Lij *= (LaplacianCoeff + TauOne);
          Lij *= (LaplacianCoeff);
          rLeftHandSideMatrix(i, j) += GaussWeight * Lij;
        }

        // RHS contribution

        // Velocity divergence
        double RHSi = -N[i] * DivU;
        // RHSi = -N[i] * rElementalVariables.VolumetricDefRate;
        // double difference = DivU - rElementalVariables.VolumetricDefRate;
        // std::cout << "difference" << difference<< DivU<<std::endl;

        for (SizeType d = 0; d < TDim; ++d)
        {
          // Momentum stabilization
          // RHSi += rDN_DX(i, d) * TauOne * (Density * (BodyForce[d] /* - Conv*/) - OldPressureGradient[d] - MomentumProjection[d]);
          RHSi += (eta - 1.0) * LaplacianCoeff * rDN_DX(i, d) * OldPressureGradient[d];
        }

        rRightHandSideVector[i] += GaussWeight * RHSi;
      }
    }

    ///////////////////////////////// bulk matrix ////////////////////////////////
    VectorType PressureValues = ZeroVector(NumNodes);
    VectorType PressureValuesForRHS = ZeroVector(NumNodes);
    this->GetPressureValues(PressureValuesForRHS, 0);

    this->GetPressureValues(PressureValues, 1);
    noalias(PressureValuesForRHS) += -PressureValues;
    MatrixType BulkMatrix = ZeroMatrix(NumNodes, NumNodes);

    this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkCoeff);
    noalias(rLeftHandSideMatrix) += BulkMatrix;
    noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);
    ///////////////////////////////// bulk matrix ////////////////////////////////
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::CalculatePSPGLocalContinuityEqForPressure(MatrixType &rLeftHandSideMatrix,
                                                                                          VectorType &rRightHandSideVector,
                                                                                          const ProcessInfo &rCurrentProcessInfo)
  {

    GeometryType &rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != NumNodes)
      rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

    rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);

    if (rRightHandSideVector.size() != NumNodes)
      rRightHandSideVector.resize(NumNodes);

    rRightHandSideVector = ZeroVector(NumNodes);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    double TimeStep = rCurrentProcessInfo[DELTA_TIME];
    double theta = 1.0;
    double ElemSize = this->ElementSize();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    double DeviatoricCoeff = 0;
    double VolumetricCoeff = 0;
    double Density = 0;
    double totalVolume = 0;
    bool computeElement = false;

    MatrixType DynamicStabilizationMatrix = ZeroMatrix(NumNodes, NumNodes);

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; ++g)
    {
      const double GaussWeight = GaussWeights[g];
      totalVolume += GaussWeight;
      const ShapeFunctionsType &N = row(NContainer, g);
      const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];
      computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta);
      computeElement = true;

      this->EvaluateInPoint(DeviatoricCoeff, DYNAMIC_VISCOSITY, N);
      this->EvaluateInPoint(Density, DENSITY, N);
      this->EvaluateInPoint(VolumetricCoeff, BULK_MODULUS, N);

      VolumetricCoeff *= TimeStep;

      double Tau = 0;
      this->CalculateTauPSPG(Tau, ElemSize, Density, DeviatoricCoeff, rCurrentProcessInfo);

      if (computeElement == true && this->IsNot(BLOCKED) && this->IsNot(ISOLATED))
      {

        double StabilizedWeight = Tau * GaussWeight;
        this->ComputeStabLaplacianMatrix(rLeftHandSideMatrix, rDN_DX, StabilizedWeight);

        array_1d<double, TDim> OldPressureGradient = ZeroVector(TDim);
        this->EvaluateGradientInPoint(OldPressureGradient, PRESSURE, rDN_DX);

        for (SizeType i = 0; i < NumNodes; ++i)
        {
          // RHS contribution
          // Velocity divergence
          rRightHandSideVector[i] += GaussWeight * N[i] * rElementalVariables.VolumetricDefRate;

          this->AddPspgDynamicPartStabilization(rRightHandSideVector, Tau, Density, GaussWeight, TimeStep, rDN_DX, N, i);

          double laplacianRHSi = 0;
          double bodyForceStabilizedRHSi = 0;
          array_1d<double, 3> &VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
          for (SizeType d = 0; d < TDim; ++d)
          {
            laplacianRHSi += StabilizedWeight * rDN_DX(i, d) * OldPressureGradient[d];

            bodyForceStabilizedRHSi += StabilizedWeight * rDN_DX(i, d) * (Density * VolumeAcceleration[d]);
          }
          rRightHandSideVector[i] += -laplacianRHSi - bodyForceStabilizedRHSi;
        }
      }
    }
    VectorType PressureValues = ZeroVector(NumNodes);
    VectorType PressureValuesForRHS = ZeroVector(NumNodes);
    this->GetPressureValues(PressureValuesForRHS, 0);

    // VectorType AccelerationValues = ZeroVector(NumNodes);
    // this->GetAccelerationValues(AccelerationValues, 0);

    // noalias(rRightHandSideVector) += prod(DynamicStabilizationMatrix, AccelerationValues);

    //the LHS matrix up to now just contains the laplacian term and the bound term
    // noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, PressureValuesForRHS);

    this->GetPressureValues(PressureValues, 1);
    noalias(PressureValuesForRHS) += -PressureValues;
    MatrixType BulkMatrix = ZeroMatrix(NumNodes, NumNodes);
    double lumpedBulkCoeff = totalVolume / VolumetricCoeff;

    this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkCoeff);
    noalias(rLeftHandSideMatrix) += BulkMatrix;
    noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::CalculateFICLocalContinuityEqForPressure(MatrixType &rLeftHandSideMatrix,
                                                                                         VectorType &rRightHandSideVector,
                                                                                         const ProcessInfo &rCurrentProcessInfo)
  {

    GeometryType &rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != NumNodes)
      rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

    rLeftHandSideMatrix = ZeroMatrix(NumNodes, NumNodes);
    MatrixType LaplacianMatrix = ZeroMatrix(NumNodes, NumNodes);

    if (rRightHandSideVector.size() != NumNodes)
      rRightHandSideVector.resize(NumNodes);

    rRightHandSideVector = ZeroVector(NumNodes);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    double TimeStep = rCurrentProcessInfo[DELTA_TIME];
    double theta = 1.0;
    double ElemSize = this->ElementSize();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    double maxViscousValueForStabilization = 0.1;

    double Density = 0;
    double VolumetricCoeff = 0;
    double DeviatoricCoeff = 0;

    VectorType NewRhsLaplacian = ZeroVector(NumNodes);

    double totalVolume = 0;
    bool computeElement = false;
    double Tau = 0;
    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; ++g)
    {
      const double GaussWeight = GaussWeights[g];
      totalVolume += GaussWeight;
      const ShapeFunctionsType &N = row(NContainer, g);
      const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

      this->EvaluateInPoint(DeviatoricCoeff, DYNAMIC_VISCOSITY, N);
      this->EvaluateInPoint(Density, DENSITY, N);
      this->EvaluateInPoint(VolumetricCoeff, BULK_MODULUS, N);

      VolumetricCoeff *= TimeStep;

      if (DeviatoricCoeff > maxViscousValueForStabilization)
      {
        DeviatoricCoeff = maxViscousValueForStabilization;
      }

            this->CalculateTauFIC(Tau, ElemSize, Density, DeviatoricCoeff, rCurrentProcessInfo);
      computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta);

      if (computeElement == true && this->IsNot(BLOCKED) && this->IsNot(ISOLATED))
      {
        double BoundLHSCoeff = Tau * 4.0 * GaussWeight / (ElemSize * ElemSize);

        this->ComputeBoundLHSMatrix(rLeftHandSideMatrix, N, BoundLHSCoeff);

        double BoundRHSCoeffAcc = Tau * Density * 2 * GaussWeight / ElemSize;
        double BoundRHSCoeffDev = Tau * 8.0 * DeviatoricCoeff * GaussWeight / (ElemSize * ElemSize);
        this->ComputeBoundRHSVectorComplete(rRightHandSideVector, TimeStep, BoundRHSCoeffAcc, BoundRHSCoeffDev, rElementalVariables.SpatialDefRate);

        double StabLaplacianWeight = Tau * GaussWeight;
        this->ComputeStabLaplacianMatrix(LaplacianMatrix, rDN_DX, StabLaplacianWeight);

        array_1d<double, TDim> OldPressureGradient = ZeroVector(TDim);
        this->EvaluateGradientInPoint(OldPressureGradient, PRESSURE, rDN_DX);

        for (SizeType i = 0; i < NumNodes; ++i)
        {
          // RHS contribution
          // Velocity divergence
          rRightHandSideVector[i] += GaussWeight * N[i] * rElementalVariables.VolumetricDefRate;
          this->AddStabilizationNodalTermsRHS(rRightHandSideVector, Tau, Density, GaussWeight, rDN_DX, i);
          double laplacianRHSi = 0;
          for (SizeType d = 0; d < TDim; ++d)
          {
            laplacianRHSi += StabLaplacianWeight * rDN_DX(i, d) * OldPressureGradient[d];
          }
          rRightHandSideVector[i] += -laplacianRHSi;
        }
      }
    }

    if (computeElement == true && this->IsNot(BLOCKED) && this->IsNot(ISOLATED))
    {

      VectorType PressureValues = ZeroVector(NumNodes);
      VectorType PressureValuesForRHS = ZeroVector(NumNodes);
      this->GetPressureValues(PressureValuesForRHS, 0);
      //the LHS matrix up to now just contains the laplacian term and the bound term
      noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, PressureValuesForRHS);
      rLeftHandSideMatrix += LaplacianMatrix;

      this->GetPressureValues(PressureValues, 1);
      noalias(PressureValuesForRHS) += -PressureValues;
      MatrixType BulkMatrix = ZeroMatrix(NumNodes, NumNodes);
      MatrixType BulkMatrixConsistent = ZeroMatrix(NumNodes, NumNodes);
      double lumpedBulkCoeff = totalVolume / (VolumetricCoeff);
      double lumpedBulkStabCoeff = lumpedBulkCoeff * Tau * Density / TimeStep;

      this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkCoeff);
      noalias(rLeftHandSideMatrix) += BulkMatrix;
      noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);

      this->GetPressureVelocityValues(PressureValues, 0);
      noalias(PressureValuesForRHS) += -PressureValues * TimeStep;
      noalias(BulkMatrix) = ZeroMatrix(NumNodes, NumNodes);
      this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkStabCoeff);
      noalias(rLeftHandSideMatrix) += BulkMatrix;
      noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);
    }
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::AddStabilizationNodalTermsRHS(VectorType &rRightHandSideVector,
                                                                              const double Tau,
                                                                              const double Density,
                                                                              const double Weight,
                                                                              const ShapeFunctionDerivativesType &rDN_DX,
                                                                              const SizeType i)
  {

    double RHSi = 0;
    if (this->GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION))
    { // it must be checked once at the begining only
      array_1d<double, 3> &VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

      for (SizeType d = 0; d < TDim; ++d)
      {
        RHSi += -rDN_DX(i, d) * Tau * (Density * VolumeAcceleration[d]);
      }
    }
    rRightHandSideVector[i] += Weight * RHSi;
  }

  template <>
  void ThreeStepUpdatedLagrangianElement<2>::ComputeBoundRHSVectorComplete(VectorType &BoundRHSVector,
                                                                           const double TimeStep,
                                                                           const double BoundRHSCoeffAcc,
                                                                           const double BoundRHSCoeffDev,
                                                                           const VectorType SpatialDefRate)
  {
    GeometryType &rGeom = this->GetGeometry();
    const double coeff = 1.0 / 3.0;
    const double timeFactor = 0.5 / TimeStep;

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE))
    {
      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);

      this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 1, 2);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(MeanAcc) = 0.5 * (AccA + AccB);

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

      if (rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
        BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[1].IsNot(INLET))
        BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    {

      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);

      this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 2, 1);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(MeanAcc) = 0.5 * (AccA + AccB);

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

      if (rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
        BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[2].IsNot(INLET))
        BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }

    if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    {

      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);

      this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 1, 2, 0);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(MeanAcc) = 0.5 * (AccA + AccB);

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

      if (rGeom[1].IsNot(INLET))
        BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[2].IsNot(INLET))
        BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }
  }

  template <>
  void ThreeStepUpdatedLagrangianElement<3>::ComputeBoundRHSVectorComplete(VectorType &BoundRHSVector,
                                                                           const double TimeStep,
                                                                           const double BoundRHSCoeffAcc,
                                                                           const double BoundRHSCoeffDev,
                                                                           const VectorType SpatialDefRate)
  {
    GeometryType &rGeom = this->GetGeometry();
    const double coeff = 0.25;
    const double timeFactor = 0.5 / TimeStep;
    const double one_third = 1.0 / 3.0;

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    {

      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> AccC(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);

      this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 1, 2, 3);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccC) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(MeanAcc) = (AccA + AccB + AccC) * one_third;

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

      if (rGeom[0].IsNot(INLET))
        BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[1].IsNot(INLET))
        BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[2].IsNot(INLET))
        BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    {

      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> AccC(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);
      this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 1, 3, 2);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccC) = timeFactor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(MeanAcc) = (AccA + AccB + AccC) * one_third;

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

      if (rGeom[0].IsNot(INLET))
        BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[1].IsNot(INLET))
        BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[3].IsNot(INLET))
        BoundRHSVector[3] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    {

      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> AccC(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);

      this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 2, 3, 1);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccC) = timeFactor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(MeanAcc) = (AccA + AccB + AccC) * one_third;

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

      if (rGeom[0].IsNot(INLET))
        BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[2].IsNot(INLET))
        BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[3].IsNot(INLET))
        BoundRHSVector[3] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }

    if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    {

      array_1d<double, 3> AccA(3, 0.0);
      array_1d<double, 3> AccB(3, 0.0);
      array_1d<double, 3> AccC(3, 0.0);
      array_1d<double, 3> MeanAcc(3, 0.0);
      array_1d<double, 3> NormalVector(3, 0.0);

      this->GetOutwardsUnitNormalForThreePoints(NormalVector, 1, 2, 3, 0);

      double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

      noalias(AccA) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
      noalias(AccC) = timeFactor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(MeanAcc) = (AccA + AccB + AccC) * one_third;

      const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

      if (rGeom[1].IsNot(INLET))
        BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[2].IsNot(INLET))
        BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

      if (rGeom[3].IsNot(INLET))
        BoundRHSVector[3] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    }
  }

  template <>
  void ThreeStepUpdatedLagrangianElement<2>::ComputeBoundLHSMatrix(Matrix &BoundLHSMatrix,
                                                                   const ShapeFunctionsType &rN,
                                                                   const double Weight)
  {
    GeometryType &rGeom = this->GetGeometry();
    double coeff = 1.0 / 3.0;

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE))
    {
      if (rGeom[0].IsNot(INLET))
        BoundLHSMatrix(0, 0) += Weight * coeff;
      if (rGeom[1].IsNot(INLET))
        BoundLHSMatrix(1, 1) += Weight * coeff;
    }
    if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    {
      if (rGeom[0].IsNot(INLET))
        BoundLHSMatrix(0, 0) += Weight * coeff;
      if (rGeom[2].IsNot(INLET))
        BoundLHSMatrix(2, 2) += Weight * coeff;
    }
    if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    {
      if (rGeom[1].IsNot(INLET))
        BoundLHSMatrix(1, 1) += Weight * coeff;
      if (rGeom[2].IsNot(INLET))
        BoundLHSMatrix(2, 2) += Weight * coeff;
    }
    // }
  }

  template <>
  void ThreeStepUpdatedLagrangianElement<3>::ComputeBoundLHSMatrix(Matrix &BoundLHSMatrix,
                                                                   const ShapeFunctionsType &rN,
                                                                   const double Weight)
  {
    GeometryType &rGeom = this->GetGeometry();
    double coeff = 0.25;

    if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    {
      if (rGeom[0].IsNot(INLET))
        BoundLHSMatrix(0, 0) += Weight * coeff;
      if (rGeom[1].IsNot(INLET))
        BoundLHSMatrix(1, 1) += Weight * coeff;
      if (rGeom[2].IsNot(INLET))
        BoundLHSMatrix(2, 2) += Weight * coeff;
    }
    if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    {
      if (rGeom[0].IsNot(INLET))
        BoundLHSMatrix(0, 0) += Weight * coeff;
      if (rGeom[1].IsNot(INLET))
        BoundLHSMatrix(1, 1) += Weight * coeff;
      if (rGeom[3].IsNot(INLET))
        BoundLHSMatrix(3, 3) += Weight * coeff;
    }
    if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    {
      if (rGeom[0].IsNot(INLET))
        BoundLHSMatrix(0, 0) += Weight * coeff;
      if (rGeom[2].IsNot(INLET))
        BoundLHSMatrix(2, 2) += Weight * coeff;
      if (rGeom[3].IsNot(INLET))
        BoundLHSMatrix(3, 3) += Weight * coeff;
    }
    if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    {
      if (rGeom[1].IsNot(INLET))
        BoundLHSMatrix(1, 1) += Weight * coeff;
      if (rGeom[2].IsNot(INLET))
        BoundLHSMatrix(2, 2) += Weight * coeff;
      if (rGeom[3].IsNot(INLET))
        BoundLHSMatrix(3, 3) += Weight * coeff;
    }
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::CalculateTauFIC(double &Tau,
                                                                double ElemSize,
                                                                const double Density,
                                                                const double Viscosity,
                                                                const ProcessInfo &rCurrentProcessInfo)
  {
    double DeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
    if (rCurrentProcessInfo.GetValue(DELTA_TIME) < rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME))
    {
      DeltaTime = 0.5 * rCurrentProcessInfo.GetValue(DELTA_TIME) + 0.5 * rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME);
    }

    double MeanVelocity = 0;
    this->CalcMeanVelocityNorm(MeanVelocity, 0);

    // Tau = 1.0 / (2.0 * Density *(0.5 * MeanVelocity / ElemSize + 0.5/DeltaTime) +  8.0 * Viscosity / (ElemSize * ElemSize) );
    Tau = (ElemSize * ElemSize * DeltaTime) / (Density * MeanVelocity * DeltaTime * ElemSize + Density * ElemSize * ElemSize + 8.0 * Viscosity * DeltaTime);

    const double tolerance = 1.0e-13;
    if (MeanVelocity < tolerance)
    {
      Tau = 0;
    }
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::ComputeBulkMatrixLump(Matrix &BulkMatrix,
                                                                      const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    if ((NumNodes == 3 && TDim == 2) || (NumNodes == 4 && TDim == 3))
    {
      double coeff = 1.0 + TDim;
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        // LHS contribution
        double Mij = Weight / coeff;
        BulkMatrix(i, i) += Mij;
      }
    }
    else
    {
      std::cout << "... ComputeBulkMatrixLump TO IMPLEMENT" << std::endl;
    }
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::CalculateTauPSPG(double &Tau,
                                                                 double ElemSize,
                                                                 const double Density,
                                                                 const double Viscosity,
                                                                 const ProcessInfo &rCurrentProcessInfo)
  {
    double DeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
    if (rCurrentProcessInfo.GetValue(DELTA_TIME) < rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME))
    {
      DeltaTime = 0.5 * rCurrentProcessInfo.GetValue(DELTA_TIME) + 0.5 * rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME);
    }

    double MeanVelocity = 0;
    this->CalcMeanVelocityNorm(MeanVelocity, 0);

    // // Tau Fic
    //Tau = (ElemSize * ElemSize * DeltaTime) / (Density * MeanVelocity * DeltaTime * ElemSize + Density * ElemSize * ElemSize + 8.0 * Viscosity * DeltaTime);

    // // Tau Tezduyar first proposal
    // double timeScale = 0.5 * (ElemSize / (2 * MeanVelocity) + DeltaTime);
    // double Reynolds = MeanVelocity * ElemSize * Density / (2.0 * Viscosity);
    // double ReynoldsFactor = 1.0;
    // if (Reynolds < 3)
    // {
    //   ReynoldsFactor = (ElemSize * ReynoldsFactor / (2 * MeanVelocity)) / 3.0;
    // }
    // Tau = timeScale * ReynoldsFactor / Density;

    // // Tau Tezduyar second proposal
    // double timeStepTerm = 2.0 / DeltaTime;
    // double timeScaleTerm = 2.0 * MeanVelocity / ElemSize;
    // double viscousTerm = 3.0 * (4 * Viscosity / (Density * pow(ElemSize, 2)));
    // // std::cout << "timeStepTerm= "<<timeStepTerm<< " timeScaleTerm= "<<timeScaleTerm<< " viscousTerm= "<<viscousTerm << std::endl;
    // double alternativeTau = pow(timeStepTerm, 2) + pow(timeScaleTerm, 2) + pow(viscousTerm, 2);
    // Tau = 1.0 / (std::sqrt(alternativeTau) * Density);

    // // Simplest Tau
    Tau = DeltaTime / Density;

    // std::cout << "Tau= " << Tau << std::endl;

    const double tolerance = 1.0e-13;
    if (MeanVelocity < tolerance)
    {
      Tau = 0;
    }
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::ComputeStabLaplacianMatrix(MatrixType &StabLaplacianMatrix,
                                                                           const ShapeFunctionDerivativesType &rDN_DX,
                                                                           const double Weight)

  {
    // LHS contribution
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      for (SizeType j = 0; j < NumNodes; ++j)
      {
        double Lij = 0.0;
        for (SizeType d = 0; d < TDim; ++d)
        {
          Lij += rDN_DX(i, d) * rDN_DX(j, d);
        }
        StabLaplacianMatrix(i, j) += Weight * Lij;
      }
    }
  }

  template <>
  void ThreeStepUpdatedLagrangianElement<2>::AddPspgDynamicPartStabilization(VectorType &rRightHandSideVector,
                                                                             const double Tau,
                                                                             const double Density,
                                                                             const double Weight,
                                                                             const double TimeStep,
                                                                             const ShapeFunctionDerivativesType &rDN_DX,
                                                                             const ShapeFunctionsType &rN,
                                                                             const SizeType i)
  {

    double RHSi = 0;
    // LHS contribution
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType j = 0; j < NumNodes; ++j)
    {
      RHSi += rDN_DX(i, 0) * rN[j] * this->GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION_X, 0) + rDN_DX(i, 1) * rN[j] * this->GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION_Y, 0);
    }
    rRightHandSideVector[i] += Weight * Tau * Density * RHSi;
  }

  template <>
  void ThreeStepUpdatedLagrangianElement<3>::AddPspgDynamicPartStabilization(VectorType &rRightHandSideVector,
                                                                             const double Tau,
                                                                             const double Density,
                                                                             const double Weight,
                                                                             const double TimeStep,
                                                                             const ShapeFunctionDerivativesType &rDN_DX,
                                                                             const ShapeFunctionsType &rN,
                                                                             const SizeType i)
  {

    double RHSi = 0;
    //////////////////////////// TO CHEEEEEEEEEEEEEEEEEEEEEEECCCCCCCCCCCCCCCCCCCCKKKKKKKKKKKKKKK///////////////////////

    // LHS contribution
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType j = 0; j < NumNodes; ++j)
    {
      double termX = rDN_DX(i, 0) * rN[j] * this->GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION_X, 0);
      double termY = rDN_DX(i, 1) * rN[j] * this->GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION_Y, 0);
      double termZ = rDN_DX(i, 2) * rN[j] * this->GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION_Z, 0);
      RHSi += Tau * Density * (termX + termY + termZ);
    }
    rRightHandSideVector[i] += Weight * RHSi;
  }

  template <unsigned int TDim>
  void ThreeStepUpdatedLagrangianElement<TDim>::InitializeElementalVariables(ElementalVariables &rElementalVariables)
  {
    KRATOS_TRY;

    unsigned int voigtsize = 3;
    if (TDim == 3)
    {
      voigtsize = 6;
    }
    rElementalVariables.voigtsize = voigtsize;

    rElementalVariables.ConstitutiveMatrix = ZeroMatrix(voigtsize, voigtsize);

    rElementalVariables.DetFgrad = 1.0;

    rElementalVariables.DetFgradVel = 1.0;

    rElementalVariables.DeviatoricInvariant = 1.0;

    rElementalVariables.EquivalentStrainRate = 1.0;

    rElementalVariables.VolumetricDefRate = 1.0;

    rElementalVariables.SpatialDefRate = ZeroVector(voigtsize);

    rElementalVariables.MDGreenLagrangeMaterial.resize(voigtsize, false);

    noalias(rElementalVariables.MDGreenLagrangeMaterial) = ZeroVector(voigtsize);

    rElementalVariables.Fgrad = ZeroMatrix(TDim, TDim);

    rElementalVariables.InvFgrad = ZeroMatrix(TDim, TDim);

    rElementalVariables.FgradVel = ZeroMatrix(TDim, TDim);

    rElementalVariables.InvFgradVel = ZeroMatrix(TDim, TDim);

    rElementalVariables.SpatialVelocityGrad = ZeroMatrix(TDim, TDim);

    rElementalVariables.MeanPressure = 0;

    rElementalVariables.CurrentTotalCauchyStress = ZeroVector(voigtsize);

    rElementalVariables.UpdatedTotalCauchyStress = ZeroVector(voigtsize);

    rElementalVariables.CurrentDeviatoricCauchyStress = ZeroVector(voigtsize);

    rElementalVariables.UpdatedDeviatoricCauchyStress = ZeroVector(voigtsize);

    KRATOS_CATCH("");
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

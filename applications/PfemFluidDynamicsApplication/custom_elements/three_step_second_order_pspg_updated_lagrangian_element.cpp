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
#include "custom_elements/three_step_second_order_pspg_updated_lagrangian_element.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

  template <unsigned int TDim>
  Element::Pointer ThreeStepSecondOrderPspgUpdatedLagrangianElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
  {
    KRATOS_TRY;

    ThreeStepSecondOrderPspgUpdatedLagrangianElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());
    
    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());
    
    return Element::Pointer(new ThreeStepSecondOrderPspgUpdatedLagrangianElement(NewElement));

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void ThreeStepSecondOrderPspgUpdatedLagrangianElement<TDim>::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                                                                    VectorType &rRightHandSideVector,
                                                                                    const ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;

    switch (rCurrentProcessInfo[FRACTIONAL_STEP])
    {
    case 1:
    {
      this->CalculateSecondVelocitySystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
      break;
    }
    case 5:
    {
      this->CalculateFSplusPSPGPressureSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
      break;
    }

    default:
    {
      KRATOS_THROW_ERROR(std::logic_error, "Unexpected value for THREE_STEP_SECOND_ORDER_PSPG_UPDATED_LAGRANGIAN_ELEMENT index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
    }
    }

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void ThreeStepSecondOrderPspgUpdatedLagrangianElement<TDim>::CalculateFSplusPSPGPressureSystem(MatrixType &rLeftHandSideMatrix,
                                                                                                 VectorType &rRightHandSideVector,
                                                                                                 const ProcessInfo &rCurrentProcessInfo)
  {

    GeometryType &rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != NumNodes)
      rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(NumNodes, NumNodes);

    if (rRightHandSideVector.size() != NumNodes)
      rRightHandSideVector.resize(NumNodes, false);

    noalias(rRightHandSideVector) = ZeroVector(NumNodes);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    double TimeStep = rCurrentProcessInfo[DELTA_TIME];
    double ElemSize = this->ElementSize();

    double Viscosity = 0;
    double Density = 0;
    double totalVolume = 0;

    MatrixType DynamicStabilizationMatrix = ZeroMatrix(NumNodes, NumNodes);

    // ElementalVariables rElementalVariables;
    // this->InitializeElementalVariables(rElementalVariables);

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; ++g)
    {
      const double GaussWeight = GaussWeights[g];
      totalVolume += GaussWeight;
      const ShapeFunctionsType &rN = row(NContainer, g);
      const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

      this->EvaluateInPoint(Viscosity, DYNAMIC_VISCOSITY, rN);
      this->EvaluateInPoint(Density, DENSITY, rN);
      double PSPGweight = 0.1;
      double Tau = 0;
      this->CalculateTauPSPG(Tau, ElemSize, Density, Viscosity, rCurrentProcessInfo);

      double StabilizedWeight = Tau * GaussWeight;

      array_1d<double, TDim> OldPressureGradient = ZeroVector(TDim);
      this->EvaluateGradientDifferenceInPoint(OldPressureGradient, PRESSURE, rDN_DX, PSPGweight);

      double DivU = 0;
      this->EvaluateDivergenceInPoint(DivU, VELOCITY, rDN_DX);

      for (SizeType i = 0; i < NumNodes; ++i)
      {
        // LHS contribution
        double dynamicRHSi = 0;
        unsigned int freeSurfaceNodes = 0;
        for (SizeType j = 0; j < NumNodes; ++j)
        {
          double Lij = 0.0;
          for (SizeType d = 0; d < TDim; ++d)
            Lij += (1.0 + PSPGweight) * rDN_DX(i, d) * rDN_DX(j, d);

          Lij *= StabilizedWeight;
          rLeftHandSideMatrix(i, j) += Lij;

          dynamicRHSi += rDN_DX(i, 0) * rN[j] * (this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_X, 0) - this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_X, 1)) / TimeStep;
          dynamicRHSi += rDN_DX(i, 1) * rN[j] * (this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Y, 0) - this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Y, 1)) / TimeStep;
          if (TDim == 3)
          {
            dynamicRHSi += rDN_DX(i, 2) * rN[j] * (this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Z, 0) - this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY_Z, 1)) / TimeStep;
          }
          if (rGeom[j].Is(FREE_SURFACE) && rGeom[j].IsNot(INLET))
          {
            freeSurfaceNodes++;
          }
        }

        if (freeSurfaceNodes >= TDim)
        {
          // bool computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, 1.0);
          // VectorType deviatoricSpatialDefRate = rElementalVariables.SpatialDefRate;
          // deviatoricSpatialDefRate[0] -= DivU / 3.0;
          // deviatoricSpatialDefRate[1] -= DivU / 3.0;
          const double lagMultiplier = TimeStep / (ElemSize * ElemSize * Density);
          this->ComputeBoundaryTermsForPressureSystem(rLeftHandSideMatrix, rRightHandSideVector, rN, lagMultiplier);
        }

        // RHS contribution
        rRightHandSideVector[i] += -PSPGweight * GaussWeight * Tau * Density * dynamicRHSi;

        rRightHandSideVector[i] += -(1.0 + PSPGweight) * GaussWeight * rN[i] * DivU;

        double laplacianRHSi = 0;
        double bodyForceStabilizedRHSi = 0;
        array_1d<double, 3> VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

        // // stored body force information for the closed domain test with analytical solution
        // array_1d<double, 3> VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
        // double posX = (this->GetGeometry()[0].X() + this->GetGeometry()[1].X() + this->GetGeometry()[2].X()) / 3.0;
        // double posY = (this->GetGeometry()[0].Y() + this->GetGeometry()[1].Y() + this->GetGeometry()[2].Y()) / 3.0;
        // double coeffX = (12.0 - 24.0 * posY) * pow(posX, 4);
        // coeffX += (-24.0 + 48.0 * posY) * pow(posX, 3);
        // coeffX += (-48.0 * posY + 72.0 * pow(posY, 2) - 48.0 * pow(posY, 3) + 12.0) * pow(posX, 2);
        // coeffX += (-2.0 + 24.0 * posY - 72.0 * pow(posY, 2) + 48.0 * pow(posY, 3)) * posX;
        // coeffX += 1.0 - 4.0 * posY + 12.0 * pow(posY, 2) - 8.0 * pow(posY, 3);
        // double coeffY = (8.0 - 48.0 * posY + 48.0 * pow(posY, 2)) * pow(posX, 3);
        // coeffY += (-12.0 + 72.0 * posY - 72.0 * pow(posY, 2)) * pow(posX, 2);
        // coeffY += (4.0 - 24.0 * posY + 48.0 * pow(posY, 2) - 48.0 * pow(posY, 3) + 24.0 * pow(posY, 4)) * posX;
        // coeffY += -12.0 * pow(posY, 2) + 24.0 * pow(posY, 3) - 12.0 * pow(posY, 4);
        // VolumeAcceleration[0] *= coeffX;
        // VolumeAcceleration[1] *= coeffY;
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        for (SizeType d = 0; d < TDim; ++d)
        {
          laplacianRHSi += -StabilizedWeight * rDN_DX(i, d) * OldPressureGradient[d];

          bodyForceStabilizedRHSi += PSPGweight * StabilizedWeight * rDN_DX(i, d) * (Density * VolumeAcceleration[d]);
        }
        rRightHandSideVector[i] += laplacianRHSi + bodyForceStabilizedRHSi;
      }
    }
  }

  /*
   * Template class definition (this should allow us to compile the desired template instantiations)
   */

  template class ThreeStepSecondOrderPspgUpdatedLagrangianElement<2>;
  template class ThreeStepSecondOrderPspgUpdatedLagrangianElement<3>;

} // namespace Kratos

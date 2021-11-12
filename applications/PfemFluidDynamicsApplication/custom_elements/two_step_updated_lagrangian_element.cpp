//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:              April 2018 $
//   Revision:            $Revision:                 0.0 $
//
//   Implementation of the Gauss-Seidel two step Updated Lagrangian Velocity-Pressure element
//     ( There is a ScalingConstant to multiply the mass balance equation for a number because i read it somewhere)
//

// System includes

// External includes

// Project includes
#include "custom_elements/two_step_updated_lagrangian_element.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

  template <unsigned int TDim>
  Element::Pointer TwoStepUpdatedLagrangianElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
  {
    KRATOS_TRY;

    TwoStepUpdatedLagrangianElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());
    return Element::Pointer(new TwoStepUpdatedLagrangianElement(NewElement));

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianElement<TDim>::EquationIdVector(EquationIdVectorType &rResult,
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
    default:
    {
      KRATOS_THROW_ERROR(std::logic_error, "Unexpected value for FRACTIONAL_STEP index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
    }
    }

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianElement<TDim>::GetDofList(DofsVectorType &rElementalDofList,
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
    default:
    {
      KRATOS_THROW_ERROR(std::logic_error, "Unexpected value for FRACTIONAL_STEP index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
    }
    }

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianElement<TDim>::ComputeLumpedMassMatrix(Matrix &rMassMatrix,
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
  void TwoStepUpdatedLagrangianElement<TDim>::ComputeMassMatrix(Matrix &rMassMatrix,
                                                                const ShapeFunctionsType &rN,
                                                                const double Weight,
                                                                double &MeanValue)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    IndexType FirstRow = 0;
    IndexType FirstCol = 0;
    double Count = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      for (SizeType j = 0; j < NumNodes; ++j)
      {
        const double Mij = Weight * rN[i] * rN[j];

        for (SizeType d = 0; d < TDim; ++d)
        {
          rMassMatrix(FirstRow + d, FirstCol + d) += Mij;
          MeanValue += fabs(Mij);
          Count += 1.0;
        }
        FirstCol += TDim;
      }
      FirstRow += TDim;
      FirstCol = 0;
    }
    MeanValue *= 1.0 / Count;
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianElement<TDim>::AddExternalForces(Vector &rRHSVector,
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
      { // it must be checked once at the begining only
        // Build RHS
        // double posX=this->GetGeometry()[i].X();
        // double posY=this->GetGeometry()[i].Y();

        // double posX=(this->GetGeometry()[0].X() + this->GetGeometry()[1].X() + this->GetGeometry()[2].X())/3.0;

        // double posY=(this->GetGeometry()[0].Y() + this->GetGeometry()[1].Y() + this->GetGeometry()[2].Y())/3.0;

        // double coeffX =(12.0-24.0*posY)*pow(posX,4);

        // coeffX += (-24.0+48.0*posY)*pow(posX,3);

        // coeffX += (-48.0*posY+72.0*pow(posY,2)-48.0*pow(posY,3)+12.0)*pow(posX,2);

        // coeffX += (-2.0+24.0*posY-72.0*pow(posY,2)+48.0*pow(posY,3))*posX;

        // coeffX += 1.0-4.0*posY+12.0*pow(posY,2)-8.0*pow(posY,3);

        // double coeffY =(8.0-48.0*posY+48.0*pow(posY,2))*pow(posX,3);

        // coeffY += (-12.0+72.0*posY-72.0*pow(posY,2))*pow(posX,2);

        // coeffY += (4.0-24.0*posY+48.0*pow(posY,2)-48.0*pow(posY,3)+24.0*pow(posY,4))*posX;

        // coeffY += -12.0*pow(posY,2)+24.0*pow(posY,3)-12.0*pow(posY,4);

        // rRHSVector[FirstRow] += Weight * Density * rN[i] * VolumeAcceleration[0]*coeffX;

        // rRHSVector[FirstRow+1] += Weight * Density * rN[i] * VolumeAcceleration[1]*coeffY;

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
  void TwoStepUpdatedLagrangianElement<2>::AddInternalForces(Vector &rRHSVector,
                                                             const ShapeFunctionDerivativesType &rDN_DX,
                                                             ElementalVariables &rElementalVariables,
                                                             const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    SizeType FirstRow = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      const double lagDNXi = rDN_DX(i, 0) * rElementalVariables.InvFgrad(0, 0) + rDN_DX(i, 1) * rElementalVariables.InvFgrad(1, 0);
      const double lagDNYi = rDN_DX(i, 0) * rElementalVariables.InvFgrad(0, 1) + rDN_DX(i, 1) * rElementalVariables.InvFgrad(1, 1);
      // lagDNXi=rDN_DX(i,0);
      // lagDNYi=rDN_DX(i,1);

      rRHSVector[FirstRow] += -Weight * (lagDNXi * rElementalVariables.UpdatedTotalCauchyStress[0] +
                                         lagDNYi * rElementalVariables.UpdatedTotalCauchyStress[2]);

      rRHSVector[FirstRow + 1] += -Weight * (lagDNYi * rElementalVariables.UpdatedTotalCauchyStress[1] +
                                             lagDNXi * rElementalVariables.UpdatedTotalCauchyStress[2]);

      FirstRow += 2;
    }
  }

  template <>
  void TwoStepUpdatedLagrangianElement<3>::AddInternalForces(Vector &rRHSVector,
                                                             const ShapeFunctionDerivativesType &rDN_DX,
                                                             ElementalVariables &rElementalVariables,
                                                             const double Weight)
  {

    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    SizeType FirstRow = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      double lagDNXi = rDN_DX(i, 0) * rElementalVariables.InvFgrad(0, 0) + rDN_DX(i, 1) * rElementalVariables.InvFgrad(1, 0) + rDN_DX(i, 2) * rElementalVariables.InvFgrad(2, 0);
      double lagDNYi = rDN_DX(i, 0) * rElementalVariables.InvFgrad(0, 1) + rDN_DX(i, 1) * rElementalVariables.InvFgrad(1, 1) + rDN_DX(i, 2) * rElementalVariables.InvFgrad(2, 1);
      double lagDNZi = rDN_DX(i, 0) * rElementalVariables.InvFgrad(0, 2) + rDN_DX(i, 1) * rElementalVariables.InvFgrad(1, 2) + rDN_DX(i, 2) * rElementalVariables.InvFgrad(2, 2);
      // lagDNXi=rDN_DX(i,0);
      // lagDNYi=rDN_DX(i,1);
      // lagDNZi=rDN_DX(i,2);

      rRHSVector[FirstRow] += -Weight * (lagDNXi * rElementalVariables.UpdatedTotalCauchyStress[0] +
                                         lagDNYi * rElementalVariables.UpdatedTotalCauchyStress[3] +
                                         lagDNZi * rElementalVariables.UpdatedTotalCauchyStress[4]);

      rRHSVector[FirstRow + 1] += -Weight * (lagDNYi * rElementalVariables.UpdatedTotalCauchyStress[1] +
                                             lagDNXi * rElementalVariables.UpdatedTotalCauchyStress[3] +
                                             lagDNZi * rElementalVariables.UpdatedTotalCauchyStress[5]);

      rRHSVector[FirstRow + 2] += -Weight * (lagDNZi * rElementalVariables.UpdatedTotalCauchyStress[2] +
                                             lagDNXi * rElementalVariables.UpdatedTotalCauchyStress[4] +
                                             lagDNYi * rElementalVariables.UpdatedTotalCauchyStress[5]);

      FirstRow += 3;
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianElement<TDim>::ComputeBulkMatrix(Matrix &BulkMatrix,
                                                                const ShapeFunctionsType &rN,
                                                                const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    for (SizeType i = 0; i < NumNodes; ++i)
    {
      for (SizeType j = 0; j < NumNodes; ++j)
      {
        // LHS contribution
        double Mij = Weight * rN[i] * rN[j];
        BulkMatrix(i, j) += Mij;
      }
    }
  }

  template <>
  void TwoStepUpdatedLagrangianElement<2>::ComputeBulkMatrixConsistent(Matrix &BulkMatrix,
                                                                       const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      for (SizeType j = 0; j < NumNodes; ++j)
      {
        // LHS contribution
        double Mij = Weight / 12.0;
        if (i == j)
          Mij *= 2.0;
        BulkMatrix(i, j) += Mij;
      }
    }
  }

  template <>
  void TwoStepUpdatedLagrangianElement<3>::ComputeBulkMatrixConsistent(Matrix &BulkMatrix,
                                                                       const double Weight)
  {
    KRATOS_ERROR << "TO IMPLEMENT AND CHECK " << std::endl;
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < NumNodes; ++i)
    {
      for (SizeType j = 0; j < NumNodes; ++j)
      {
        // LHS contribution
        double Mij = Weight / 12;
        if (i == j)
          Mij *= 2.0;
        BulkMatrix(i, j) += Mij;
      }
    }
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianElement<TDim>::ComputeBulkMatrixLump(Matrix &BulkMatrix,
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
      KRATOS_ERROR << "... ComputeBulkMatrixLump TO IMPLEMENT" << std::endl;
    }
  }

  /*
   * Template class definition (this should allow us to compile the desired template instantiations)
   */

  template class TwoStepUpdatedLagrangianElement<2>;
  template class TwoStepUpdatedLagrangianElement<3>;

} // namespace Kratos

#include "vms_adjoint_element.h"

namespace Kratos {

///@name Specialized implementation for functions that depend on TDim
///@{

template<>
void VMSAdjointElement<2>::GetDofList(DofsVectorType& rElementalDofList,
                                      ProcessInfo& /*rCurrentProcessInfo*/)
{
  const unsigned int NumNodes(3), LocalSize(9);

  if (rElementalDofList.size() != LocalSize)
    rElementalDofList.resize(LocalSize);

  unsigned int LocalIndex = 0;

  for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
  {
    rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(
        ADJOINT_FLUID_VECTOR_1_X);
    rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(
        ADJOINT_FLUID_VECTOR_1_Y);
    rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(
        ADJOINT_FLUID_SCALAR_1);
  }
}

template<>
void VMSAdjointElement<3>::GetDofList(DofsVectorType& rElementalDofList,
                                      ProcessInfo& /*rCurrentProcessInfo*/)
{
  const SizeType NumNodes(4), LocalSize(16);

  if (rElementalDofList.size() != LocalSize)
    rElementalDofList.resize(LocalSize);

  IndexType LocalIndex = 0;

  for (IndexType iNode = 0; iNode < NumNodes; ++iNode)
  {
    rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(
        ADJOINT_FLUID_VECTOR_1_X);
    rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(
        ADJOINT_FLUID_VECTOR_1_Y);
    rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(
        ADJOINT_FLUID_VECTOR_1_Z);
    rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(
        ADJOINT_FLUID_SCALAR_1);
  }
}

template<>
void VMSAdjointElement<2>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& /*rCurrentProcessInfo*/)
{
  const SizeType NumNodes(3), LocalSize(9);

  if (rResult.size() != LocalSize)
    rResult.resize(LocalSize, false);

  IndexType LocalIndex = 0;

  for (IndexType iNode = 0; iNode < NumNodes; ++iNode)
  {
    rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_VECTOR_1_X)
        .EquationId();
    rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_VECTOR_1_Y)
        .EquationId();
    rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_SCALAR_1)
        .EquationId();
  }
}

template<>
void VMSAdjointElement<3>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& /*rCurrentProcessInfo*/)
{
  const SizeType NumNodes(4), LocalSize(16);

  if (rResult.size() != LocalSize)
    rResult.resize(LocalSize, false);

  IndexType LocalIndex = 0;

  for (IndexType iNode = 0; iNode < NumNodes; ++iNode)
  {
    rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_VECTOR_1_X)
        .EquationId();
    rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_VECTOR_1_Y)
        .EquationId();
    rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_VECTOR_1_Z)
        .EquationId();
    rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(ADJOINT_FLUID_SCALAR_1)
        .EquationId();
  }
}

template<>
double VMSAdjointElement<2>::CalculateElementSize(const double Area)
{
  // Diameter of a circle of given Area.
  return 1.128379167 * std::sqrt(Area);
}

template<>
double VMSAdjointElement<3>::CalculateElementSize(const double Volume)
{
  // Diameter of the sphere circumscribed to a regular tetrahedron
  // with the same volume.
  return 0.60046878 * std::pow(Volume,0.333333333333333333333);
}

template<>
void VMSAdjointElement<2>::CalculateDeterminantOfJacobianDerivatives(
    array_1d< double, 6 >& rDetJDerivatives)
{
  const SizeType LocalSize(6);

  if (rDetJDerivatives.size() != LocalSize)
    rDetJDerivatives.resize(LocalSize,false);

  rDetJDerivatives[0] = this->GetGeometry()[1].Y() - this->GetGeometry()[2].Y();
  rDetJDerivatives[1] = this->GetGeometry()[2].X() - this->GetGeometry()[1].X();
  rDetJDerivatives[2] = this->GetGeometry()[2].Y() - this->GetGeometry()[0].Y();
  rDetJDerivatives[3] = this->GetGeometry()[0].X() - this->GetGeometry()[2].X();
  rDetJDerivatives[4] = this->GetGeometry()[0].Y() - this->GetGeometry()[1].Y();
  rDetJDerivatives[5] = this->GetGeometry()[1].X() - this->GetGeometry()[0].X();
}

template<>
void VMSAdjointElement<3>::CalculateDeterminantOfJacobianDerivatives(
    array_1d< double, 12 >& rDetJDerivatives)
{
  const SizeType LocalSize(12);

  if (rDetJDerivatives.size() != LocalSize)
    rDetJDerivatives.resize(LocalSize,false);

  const double x0 = this->GetGeometry()[0].X();
  const double x1 = this->GetGeometry()[1].X();
  const double x2 = this->GetGeometry()[2].X();
  const double x3 = this->GetGeometry()[3].X();
  
  const double y0 = this->GetGeometry()[0].Y();
  const double y1 = this->GetGeometry()[1].Y();
  const double y2 = this->GetGeometry()[2].Y();
  const double y3 = this->GetGeometry()[3].Y();
  
  const double z0 = this->GetGeometry()[0].Z();
  const double z1 = this->GetGeometry()[1].Z();
  const double z2 = this->GetGeometry()[2].Z();
  const double z3 = this->GetGeometry()[3].Z();
        
  rDetJDerivatives[0] =-z1*y3 + z1*y2 + y3*z2 - y2*z3 + y1*z3 - y1*z2;
  rDetJDerivatives[1] =-z1*x2 + z1*x3 + x1*z2 - z2*x3 + x2*z3 - x1*z3;
  rDetJDerivatives[2] =-x2*y3 + y2*x3 - y1*x3 + x1*y3 + y1*x2 - x1*y2;
  rDetJDerivatives[3] = y2*z3 - y2*z0 - y0*z3 - y3*z2 + y3*z0 + y0*z2;
  rDetJDerivatives[4] = z2*x3 - z2*x0 - z0*x3 - x2*z3 + x2*z0 + x0*z3;
  rDetJDerivatives[5] = x2*y3 - x2*y0 - x0*y3 - y2*x3 + y2*x0 + y0*x3;
  rDetJDerivatives[6] = z1*y3 - z1*y0 - y3*z0 - y1*z3 + y0*z3 + y1*z0;
  rDetJDerivatives[7] =-z1*x3 + z1*x0 + z0*x3 - x0*z3 + x1*z3 - x1*z0;
  rDetJDerivatives[8] = x1*y0 + x0*y3 - y0*x3 - x1*y3 + y1*x3 - y1*x0;
  rDetJDerivatives[9] =-z1*y2 + z1*y0 + y2*z0 - y0*z2 + y1*z2 - y1*z0;
  rDetJDerivatives[10] = z1*x2 - z1*x0 - x2*z0 + z2*x0 - x1*z2 + x1*z0;
  rDetJDerivatives[11] =-y2*x0 + x1*y2 + y1*x0 - y1*x2 + x2*y0 - x1*y0;
}

template<>
void VMSAdjointElement<2>::CalculateStabilizationParametersDerivative(
    double& rTauOneDeriv,
    double& rTauTwoDeriv,
    const double TauOne,
    const double TauTwo,
    const double VelNorm,
    const double ElemSize,
    const double Density,
    const double Viscosity,
    const double DetJDeriv)
{
  const double TwoPi = 6.283185307179586;
  const double TwoOverPi = 0.636619772367581;
  rTauOneDeriv = TwoOverPi * TauOne * TauOne / std::pow(ElemSize,3)
      * (Density * VelNorm + 4.0 * Viscosity / ElemSize) * DetJDeriv;

  rTauTwoDeriv = Density * VelNorm / (TwoPi * ElemSize) * DetJDeriv;
}

template<>
void VMSAdjointElement<3>::CalculateStabilizationParametersDerivative(
    double& rTauOneDeriv,
    double& rTauTwoDeriv,
    const double TauOne,
    const double TauTwo,
    const double VelNorm,
    const double ElemSize,
    const double Density,
    const double Viscosity,
    const double DetJDeriv)
{
  const double CoefOne = 0.0240562975623840;
  const double CoefTwo = 0.00601407439059599;

  rTauOneDeriv = CoefOne * TauOne * TauOne / std::pow(ElemSize,4)
      * (Density * VelNorm + 4.0 * Viscosity / ElemSize) * DetJDeriv;

  rTauTwoDeriv = CoefTwo * Density * VelNorm / (ElemSize * ElemSize) * DetJDeriv;
}

template<>
void VMSAdjointElement<2>::AddViscousTerm(
    MatrixType& rResult,
    const VMSAdjointElement<2>::ShapeFunctionDerivativesType& rDN_DX,
    const double Weight)
{
  const SizeType NumNodes = 3;

  const double FourThirds = 4.0 / 3.0;
  const double nTwoThirds = -2.0 / 3.0;

  IndexType FirstRow(0), FirstCol(0);

  for (IndexType j = 0; j < NumNodes; ++j)
  {
    for (IndexType i = 0; i < NumNodes; ++i)
    {
      // First Row
      rResult(FirstRow,FirstCol) += Weight
          * (FourThirds * rDN_DX(i,0) * rDN_DX(j,0) + rDN_DX(i,1) * rDN_DX(j,1));
      rResult(FirstRow,FirstCol+1) += Weight
          * (nTwoThirds * rDN_DX(i,0) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,0));

      // Second Row
      rResult(FirstRow+1,FirstCol) += Weight
          * (nTwoThirds * rDN_DX(i,1) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,1));
      rResult(FirstRow+1,FirstCol+1) += Weight
          * (FourThirds * rDN_DX(i,1) * rDN_DX(j,1) + rDN_DX(i,0) * rDN_DX(j,0));

      // Update Counter
      FirstRow += 3;
    }
    FirstRow = 0;
    FirstCol += 3;
  }
}

template<>
void VMSAdjointElement<3>::AddViscousTerm(
    MatrixType& rResult,
    const VMSAdjointElement<3>::ShapeFunctionDerivativesType& rDN_DX,
    const double Weight)
{
  const unsigned int NumNodes = 4;

  const double OneThird = 1.0 / 3.0;
  const double nTwoThirds = -2.0 / 3.0;

  unsigned int FirstRow(0), FirstCol(0);

  for (unsigned int j = 0; j < NumNodes; ++j)
  {
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
      // (dN_i/dx_k dN_j/dx_k)
      const double Diag = rDN_DX(i,0) * rDN_DX(j,0)
          + rDN_DX(i,1) * rDN_DX(j,1) + rDN_DX(i,2) * rDN_DX(j,2);

      // First Row
      rResult(FirstRow,FirstCol) += Weight
          * (OneThird * rDN_DX(i,0) * rDN_DX(j,0) + Diag);
      rResult(FirstRow,FirstCol+1) += Weight
          * (nTwoThirds * rDN_DX(i,0) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,0));
      rResult(FirstRow,FirstCol+2) += Weight
          * (nTwoThirds * rDN_DX(i,0) * rDN_DX(j,2) + rDN_DX(i,2) * rDN_DX(j,0));

      // Second Row
      rResult(FirstRow+1,FirstCol) += Weight
          * (nTwoThirds * rDN_DX(i,1) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,1));
      rResult(FirstRow + 1, FirstCol + 1) += Weight
          * (OneThird * rDN_DX(i,1) * rDN_DX(j,1) + Diag);
      rResult(FirstRow + 1, FirstCol + 2) += Weight
          * (nTwoThirds * rDN_DX(i,1) * rDN_DX(j,2) + rDN_DX(i,2) * rDN_DX(j,1));

      // Third Row
      rResult(FirstRow+2,FirstCol) += Weight
          * (nTwoThirds * rDN_DX(i,2) * rDN_DX(j,0) + rDN_DX(i,0) * rDN_DX(j,2));
      rResult(FirstRow+2,FirstCol+1) += Weight
          * (nTwoThirds * rDN_DX(i,2) * rDN_DX(j,1) + rDN_DX(i,1) * rDN_DX(j,2));
      rResult(FirstRow+2,FirstCol+2) += Weight
          * (OneThird * rDN_DX(i,2) * rDN_DX(j,2) + Diag);

      // Update Counter
      FirstRow += 4;
    }
    FirstRow = 0;
    FirstCol += 4;
  }
}

template<>
void VMSAdjointElement<2>::AddViscousTermDerivative(
    BoundedMatrix< double, 9, 9 >& rResult,
    const VMSAdjointElement<2>::ShapeFunctionDerivativesType& rDN_DX,
    const VMSAdjointElement<2>::ShapeFunctionDerivativesType& rDN_DX_Deriv,
    const double Weight,
    const double WeightDeriv)
{
  const SizeType NumNodes = 3;

  const double FourThirds = 4.0 / 3.0;
  const double nTwoThirds = -2.0 / 3.0;

  IndexType FirstRow(0), FirstCol(0);

  for (IndexType j = 0; j < NumNodes; ++j)
  {
    for (IndexType i = 0; i < NumNodes; ++i)
    {
      // First Row
      rResult(FirstRow,FirstCol) +=
          WeightDeriv * (FourThirds * rDN_DX(i,0) * rDN_DX(j,0)
                         + rDN_DX(i,1) * rDN_DX(j,1))
          + Weight * (FourThirds * rDN_DX_Deriv(i,0) * rDN_DX(j,0)
                      + rDN_DX_Deriv(i,1) * rDN_DX(j,1))
          + Weight * (FourThirds * rDN_DX(i,0) * rDN_DX_Deriv(j,0)
                      + rDN_DX(i,1) * rDN_DX_Deriv(j,1));
      rResult(FirstRow,FirstCol+1) +=
          WeightDeriv * (nTwoThirds * rDN_DX(i,0) * rDN_DX(j,1)
                         + rDN_DX(i,1) * rDN_DX(j,0))
          + Weight * (nTwoThirds * rDN_DX_Deriv(i,0) * rDN_DX(j,1)
                      + rDN_DX_Deriv(i,1) * rDN_DX(j,0))
          + Weight * (nTwoThirds * rDN_DX(i,0) * rDN_DX_Deriv(j,1)
                      + rDN_DX(i,1) * rDN_DX_Deriv(j,0));

      // Second Row
      rResult(FirstRow+1,FirstCol) +=
          WeightDeriv * (nTwoThirds * rDN_DX(i,1) * rDN_DX(j,0)
                         + rDN_DX(i,0) * rDN_DX(j,1))
          + Weight * (nTwoThirds * rDN_DX_Deriv(i,1) * rDN_DX(j,0)
                      + rDN_DX_Deriv(i,0) * rDN_DX(j,1))
          + Weight * (nTwoThirds * rDN_DX(i,1) * rDN_DX_Deriv(j,0)
                      + rDN_DX(i,0) * rDN_DX_Deriv(j,1));
      rResult(FirstRow+1,FirstCol+1) +=
          WeightDeriv * (FourThirds * rDN_DX(i,1) * rDN_DX(j,1)
                         + rDN_DX(i,0) * rDN_DX(j,0))
          + Weight * (FourThirds * rDN_DX_Deriv(i,1) * rDN_DX(j,1)
                      + rDN_DX_Deriv(i,0) * rDN_DX(j,0))
          + Weight * (FourThirds * rDN_DX(i,1) * rDN_DX_Deriv(j,1)
                      + rDN_DX(i,0) * rDN_DX_Deriv(j,0));

      // Update Counter
      FirstRow += 3;
    }
    FirstRow = 0;
    FirstCol += 3;
  }
}

template<>
void VMSAdjointElement<3>::AddViscousTermDerivative(
    BoundedMatrix< double, 16, 16 >& rResult,
    const VMSAdjointElement<3>::ShapeFunctionDerivativesType& rDN_DX,
    const VMSAdjointElement<3>::ShapeFunctionDerivativesType& rDN_DX_Deriv,
    const double Weight,
    const double WeightDeriv)
{
  const unsigned int NumNodes = 4;

  const double OneThird = 1.0 / 3.0;
  const double nTwoThirds = -2.0 / 3.0;

  unsigned int FirstRow(0), FirstCol(0);

  for (unsigned int j = 0; j < NumNodes; ++j)
  {
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
      // (dN_i/dx_k dN_j/dx_k)
      const double Diag = rDN_DX(i,0) * rDN_DX(j,0) + rDN_DX(i,1) * rDN_DX(j,1)
          + rDN_DX(i,2) * rDN_DX(j,2);
      const double DiagDerivNi = rDN_DX_Deriv(i,0) * rDN_DX(j,0)
          + rDN_DX_Deriv(i,1) * rDN_DX(j,1) + rDN_DX_Deriv(i,2) * rDN_DX(j,2);
      const double DiagDerivNj = rDN_DX(i,0) * rDN_DX_Deriv(j,0)
          + rDN_DX(i,1) * rDN_DX_Deriv(j,1) + rDN_DX(i,2) * rDN_DX_Deriv(j,2);

      // First Row
      rResult(FirstRow,FirstCol) +=
          WeightDeriv * (OneThird * rDN_DX(i,0) * rDN_DX(j,0) + Diag)
          + Weight * (OneThird * rDN_DX_Deriv(i,0) * rDN_DX(j,0) + DiagDerivNi)
          + Weight * (OneThird * rDN_DX(i,0) * rDN_DX_Deriv(j,0) + DiagDerivNj);
      rResult(FirstRow,FirstCol+1) +=
          WeightDeriv * (nTwoThirds * rDN_DX(i,0) * rDN_DX(j,1)
                         + rDN_DX(i,1) * rDN_DX(j,0))
          + Weight * (nTwoThirds * rDN_DX_Deriv(i,0) * rDN_DX(j,1)
                         + rDN_DX_Deriv(i,1) * rDN_DX(j,0))
          + Weight * (nTwoThirds * rDN_DX(i,0) * rDN_DX_Deriv(j,1)
                         + rDN_DX(i,1) * rDN_DX_Deriv(j,0));
      rResult(FirstRow,FirstCol+2) +=
          WeightDeriv * (nTwoThirds * rDN_DX(i,0) * rDN_DX(j,2)
                         + rDN_DX(i,2) * rDN_DX(j,0))
          + Weight * (nTwoThirds * rDN_DX_Deriv(i,0) * rDN_DX(j,2)
                         + rDN_DX_Deriv(i,2) * rDN_DX(j,0))
          + Weight * (nTwoThirds * rDN_DX(i,0) * rDN_DX_Deriv(j,2)
                         + rDN_DX(i,2) * rDN_DX_Deriv(j,0));

      // Second Row
      rResult(FirstRow+1,FirstCol) +=
          WeightDeriv * (nTwoThirds * rDN_DX(i,1) * rDN_DX(j,0)
                         + rDN_DX(i,0) * rDN_DX(j,1))
          + Weight * (nTwoThirds * rDN_DX_Deriv(i,1) * rDN_DX(j,0)
                         + rDN_DX_Deriv(i,0) * rDN_DX(j,1))
          + Weight * (nTwoThirds * rDN_DX(i,1) * rDN_DX_Deriv(j,0)
                         + rDN_DX(i,0) * rDN_DX_Deriv(j,1));
      rResult(FirstRow+1,FirstCol+1) +=
          WeightDeriv * (OneThird * rDN_DX(i,1) * rDN_DX(j,1) + Diag)
          + Weight * (OneThird * rDN_DX_Deriv(i,1) * rDN_DX(j,1) + DiagDerivNi)
          + Weight * (OneThird * rDN_DX(i,1) * rDN_DX_Deriv(j,1) + DiagDerivNj);
      rResult(FirstRow+1,FirstCol+2) +=
          WeightDeriv * (nTwoThirds * rDN_DX(i,1) * rDN_DX(j,2)
                         + rDN_DX(i,2) * rDN_DX(j,1))
          + Weight * (nTwoThirds * rDN_DX_Deriv(i,1) * rDN_DX(j,2)
                         + rDN_DX_Deriv(i,2) * rDN_DX(j,1))
          + Weight * (nTwoThirds * rDN_DX(i,1) * rDN_DX_Deriv(j,2)
                         + rDN_DX(i,2) * rDN_DX_Deriv(j,1));

      // Third Row
      rResult(FirstRow+2,FirstCol) +=
          WeightDeriv * (nTwoThirds * rDN_DX(i,2) * rDN_DX(j,0)
                         + rDN_DX(i,0) * rDN_DX(j,2))
          + Weight * (nTwoThirds * rDN_DX_Deriv(i,2) * rDN_DX(j,0)
                         + rDN_DX_Deriv(i,0) * rDN_DX(j,2))
          + Weight * (nTwoThirds * rDN_DX(i,2) * rDN_DX_Deriv(j,0)
                         + rDN_DX(i,0) * rDN_DX_Deriv(j,2));
      rResult(FirstRow+2,FirstCol+1) +=
          WeightDeriv * (nTwoThirds * rDN_DX(i,2) * rDN_DX(j,1)
                         + rDN_DX(i,1) * rDN_DX(j,2))
          + Weight * (nTwoThirds * rDN_DX_Deriv(i,2) * rDN_DX(j,1)
                         + rDN_DX_Deriv(i,1) * rDN_DX(j,2))
          + Weight * (nTwoThirds * rDN_DX(i,2) * rDN_DX_Deriv(j,1)
                         + rDN_DX(i,1) * rDN_DX_Deriv(j,2));
      rResult(FirstRow+2,FirstCol+2) +=
          WeightDeriv * (OneThird * rDN_DX(i,2) * rDN_DX(j,2) + Diag)
          + Weight * (OneThird * rDN_DX_Deriv(i,2) * rDN_DX(j,2) + DiagDerivNi)
          + Weight * (OneThird * rDN_DX(i,2) * rDN_DX_Deriv(j,2) + DiagDerivNj);

      // Update Counter
      FirstRow += 4;
    }
    FirstRow = 0;
    FirstCol += 4;
  }
}

///@} // Specialized implementations

template class VMSAdjointElement<2> ;
template class VMSAdjointElement<3> ;

}  // namespace Kratos

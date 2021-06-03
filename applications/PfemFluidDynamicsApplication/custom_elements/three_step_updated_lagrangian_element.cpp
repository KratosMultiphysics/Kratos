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

  /*
   * public ThreeStepUpdatedLagrangianElement<TDim> functions
   */
  // template <unsigned int TDim>
  // ThreeStepUpdatedLagrangianElement<TDim>::ThreeStepUpdatedLagrangianElement(ThreeStepUpdatedLagrangianElement const &rOther)
  // {
  //   KRATOS_TRY;
  //   KRATOS_CATCH("");
  // }

  template <unsigned int TDim>
  Element::Pointer ThreeStepUpdatedLagrangianElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
  {
    KRATOS_TRY;

    ThreeStepUpdatedLagrangianElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());
    return Element::Pointer(new ThreeStepUpdatedLagrangianElement(NewElement));

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
  bool ThreeStepUpdatedLagrangianElement<TDim>::CalcMechanicsUpdated(ElementalVariables &rElementalVariables,
                                                                   const ProcessInfo &rCurrentProcessInfo,
                                                                   const ShapeFunctionDerivativesType &rDN_DX,
                                                                   unsigned int g)
  {

    double theta = this->GetThetaMomentum();
    // bool computeElement=this->CalcStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
    bool computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta);
    return computeElement;
  }

  template <>
  bool ThreeStepUpdatedLagrangianElement<2>::CalcCompleteStrainRate(ElementalVariables &rElementalVariables,
                                                                  const ProcessInfo &rCurrentProcessInfo,
                                                                  const ShapeFunctionDerivativesType &rDN_DX,
                                                                  const double theta)
  {
    bool computeElement = true;
    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
    GeometryType &rGeom = this->GetGeometry();
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

    //Inverse
    rElementalVariables.InvFgrad = ZeroMatrix(dimension, dimension);
    rElementalVariables.DetFgrad = 1;
    MathUtils<double>::InvertMatrix2(rElementalVariables.Fgrad,
                                     rElementalVariables.InvFgrad,
                                     rElementalVariables.DetFgrad);

    // rElementalVariables.InvFgradVel=ZeroMatrix(dimension,dimension);
    // rElementalVariables.DetFgradVel=1;
    // MathUtils<double>::InvertMatrix2(rElementalVariables.FgradVel,
    // 		    rElementalVariables.InvFgradVel,
    // 		    rElementalVariables.DetFgradVel);

    //it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
    rElementalVariables.SpatialVelocityGrad.resize(dimension, dimension, false);
    rElementalVariables.SpatialVelocityGrad = prod(rElementalVariables.FgradVel, rElementalVariables.InvFgrad);

    rElementalVariables.VolumetricDefRate = 0;
    for (SizeType i = 0; i < dimension; i++)
    {
      rElementalVariables.VolumetricDefRate += rElementalVariables.SpatialVelocityGrad(i, i);
    }

    // //it checks whether tr(l) == div(v)
    // CheckStrain1(rElementalVariables.VolumetricDefRate,rElementalVariables.SpatialVelocityGrad);
    // CheckStrain2(rElementalVariables.SpatialVelocityGrad,rElementalVariables.Fgrad,ElementalVariables.FgradVel);
    // //it computes Material time Derivative of Green Lagrange strain tensor in MATERIAL configuration --> [D(E)/Dt]
    // // x-component
    // rElementalVariables.MDGreenLagrangeMaterial[0]=rElementalVariables.FgradVel(0,0)*rElementalVariables.Fgrad(0,0) +
    //   rElementalVariables.FgradVel(1,0)*rElementalVariables.Fgrad(1,0);
    // // y-component
    // rElementalVariables.MDGreenLagrangeMaterial[1]=rElementalVariables.FgradVel(1,1)*rElementalVariables.Fgrad(1,1) +
    //   rElementalVariables.FgradVel(0,1)*rElementalVariables.Fgrad(0,1);
    // // xy-component
    // rElementalVariables.MDGreenLagrangeMaterial[2]=(rElementalVariables.FgradVel(0,0)*rElementalVariables.Fgrad(0,1) +
    // 						  rElementalVariables.FgradVel(1,0)*rElementalVariables.Fgrad(1,1) +
    // 						  rElementalVariables.FgradVel(0,1)*rElementalVariables.Fgrad(0,0) +
    // 						  rElementalVariables.FgradVel(1,1)*rElementalVariables.Fgrad(1,0))*0.5;
    // //it computes Material time Derivative of Green Lagrange strain tensor in SPATIAL configuration  --> [d]
    // // x-component
    // rElementalVariables.SpatialDefRate[0]= rElementalVariables.InvFgrad(0,0)*rElementalVariables.MDGreenLagrangeMaterial[0]*rElementalVariables.InvFgrad(0,0) +
    //   rElementalVariables.InvFgrad(1,0)*rElementalVariables.MDGreenLagrangeMaterial[2]*rElementalVariables.InvFgrad(0,0)*2 +
    //   rElementalVariables.InvFgrad(1,0)*rElementalVariables.MDGreenLagrangeMaterial[1]*rElementalVariables.InvFgrad(1,0);
    // // y-component
    // rElementalVariables.SpatialDefRate[1]= rElementalVariables.InvFgrad(0,1)*rElementalVariables.MDGreenLagrangeMaterial[0]*rElementalVariables.InvFgrad(0,1) +
    //   rElementalVariables.InvFgrad(0,1)*rElementalVariables.MDGreenLagrangeMaterial[2]*rElementalVariables.InvFgrad(1,1)*2 +
    //   rElementalVariables.InvFgrad(1,1)*rElementalVariables.MDGreenLagrangeMaterial[1]*rElementalVariables.InvFgrad(1,1);
    // // xy-component
    // rElementalVariables.SpatialDefRate[2]=rElementalVariables.InvFgrad(0,0)*rElementalVariables.MDGreenLagrangeMaterial[0]*rElementalVariables.InvFgrad(0,1) +
    //   rElementalVariables.InvFgrad(0,0)*rElementalVariables.MDGreenLagrangeMaterial[2]*rElementalVariables.InvFgrad(1,1) +
    //   rElementalVariables.InvFgrad(1,0)*rElementalVariables.MDGreenLagrangeMaterial[2]*rElementalVariables.InvFgrad(0,1) +
    //   rElementalVariables.InvFgrad(1,0)*rElementalVariables.MDGreenLagrangeMaterial[1]*rElementalVariables.InvFgrad(1,1);

    // computeElement=CheckStrain3(rElementalVariables.SpatialDefRate,rElementalVariables.SpatialVelocityGrad);

    rElementalVariables.SpatialDefRate[0] = rElementalVariables.SpatialVelocityGrad(0, 0);
    rElementalVariables.SpatialDefRate[1] = rElementalVariables.SpatialVelocityGrad(1, 1);
    rElementalVariables.SpatialDefRate[2] = 0.5 * (rElementalVariables.SpatialVelocityGrad(1, 0) + rElementalVariables.SpatialVelocityGrad(0, 1));

    double aThird = 1.0 / 3.0;
    double dev_X = rElementalVariables.SpatialDefRate[0] -
                   (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1]) * aThird;
    double dev_Y = rElementalVariables.SpatialDefRate[1] -
                   (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1]) * aThird;
    rElementalVariables.DeviatoricInvariant = sqrt(2 * (dev_X * dev_X + dev_Y * dev_Y +
                                                        rElementalVariables.SpatialDefRate[2] * rElementalVariables.SpatialDefRate[2]));

    rElementalVariables.EquivalentStrainRate = sqrt((2.0 * rElementalVariables.SpatialDefRate[0] * rElementalVariables.SpatialDefRate[0] +
                                                     2.0 * rElementalVariables.SpatialDefRate[1] * rElementalVariables.SpatialDefRate[1] +
                                                     4.0 * rElementalVariables.SpatialDefRate[2] * rElementalVariables.SpatialDefRate[2]));

    return computeElement;
  }

  template <>
  bool ThreeStepUpdatedLagrangianElement<3>::CalcCompleteStrainRate(ElementalVariables &rElementalVariables,
                                                                  const ProcessInfo &rCurrentProcessInfo,
                                                                  const ShapeFunctionDerivativesType &rDN_DX,
                                                                  const double theta)
  {

    bool computeElement = true;
    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
    GeometryType &rGeom = this->GetGeometry();
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

    //Inverse
    rElementalVariables.InvFgrad = ZeroMatrix(dimension, dimension);
    rElementalVariables.DetFgrad = 1;
    MathUtils<double>::InvertMatrix3(rElementalVariables.Fgrad,
                                     rElementalVariables.InvFgrad,
                                     rElementalVariables.DetFgrad);

    // rElementalVariables.InvFgradVel=ZeroMatrix(dimension,dimension);
    // rElementalVariables.DetFgradVel=1;
    // MathUtils<double>::InvertMatrix3(rElementalVariables.FgradVel,
    // 		    rElementalVariables.InvFgradVel,
    // 		    rElementalVariables.DetFgradVel);

    //it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
    rElementalVariables.SpatialVelocityGrad.resize(dimension, dimension, false);
    rElementalVariables.SpatialVelocityGrad = prod(rElementalVariables.FgradVel, rElementalVariables.InvFgrad);

    rElementalVariables.VolumetricDefRate = 0;
    for (SizeType i = 0; i < dimension; i++)
    {
      rElementalVariables.VolumetricDefRate += rElementalVariables.SpatialVelocityGrad(i, i);
    }

    // //it checks whether tr(l) == div(v)
    // CheckStrain1(rElementalVariables.VolumetricDefRate,rElementalVariables.SpatialVelocityGrad);
    // CheckStrain2(rElementalVariables.SpatialVelocityGrad,rElementalVariables.Fgrad,ElementalVariables.FgradVel);

    // //it computes Material time Derivative of Green Lagrange strain tensor in MATERIAL configuration --> [D(E)/Dt]
    // MatrixType MatrixA= ZeroMatrix(3,3);
    // MatrixType MatrixB= ZeroMatrix(3,3);
    // MatrixType Matrix1= ZeroMatrix(3,3);
    // MatrixType Matrix2= ZeroMatrix(3,3);

    // MatrixA=rElementalVariables.Fgrad;
    // MatrixA(0,1)=rElementalVariables.Fgrad(1,0);
    // MatrixA(0,2)=rElementalVariables.Fgrad(2,0);
    // MatrixA(1,0)=rElementalVariables.Fgrad(0,1);
    // MatrixA(1,2)=rElementalVariables.Fgrad(2,1);
    // MatrixA(2,0)=rElementalVariables.Fgrad(0,2);
    // MatrixA(2,1)=rElementalVariables.Fgrad(1,2);

    // MatrixB=rElementalVariables.FgradVel;
    // MatrixB(0,1)=rElementalVariables.FgradVel(1,0);
    // MatrixB(0,2)=rElementalVariables.FgradVel(2,0);
    // MatrixB(1,0)=rElementalVariables.FgradVel(0,1);
    // MatrixB(1,2)=rElementalVariables.FgradVel(2,1);
    // MatrixB(2,0)=rElementalVariables.FgradVel(0,2);
    // MatrixB(2,1)=rElementalVariables.FgradVel(1,2);

    // noalias(Matrix1)=prod(MatrixB,rElementalVariables.Fgrad);
    // noalias(Matrix2)=prod(MatrixA,rElementalVariables.FgradVel);

    // rElementalVariables.MDGreenLagrangeMaterial[0]= ( Matrix1(0,0) + Matrix2(0,0) ) * 0.5;  //xx-component
    // rElementalVariables.MDGreenLagrangeMaterial[1]= ( Matrix1(1,1) + Matrix2(1,1) ) * 0.5;  //yy-component
    // rElementalVariables.MDGreenLagrangeMaterial[2]= ( Matrix1(2,2) + Matrix2(2,2) ) * 0.5;  //zz-component
    // rElementalVariables.MDGreenLagrangeMaterial[3]= ( Matrix1(0,1) + Matrix2(0,1) ) * 0.5;  //xy-component
    // rElementalVariables.MDGreenLagrangeMaterial[4]= ( Matrix1(0,2) + Matrix2(0,2) ) * 0.5;  //xz-component
    // rElementalVariables.MDGreenLagrangeMaterial[5]= ( Matrix1(1,2) + Matrix2(1,2) ) * 0.5;  //yz-component

    // //it computes Material time Derivative of Green Lagrange strain tensor in SPATIAL configuration  --> [d]
    // MatrixA=rElementalVariables.InvFgrad;
    // MatrixA(0,1)=rElementalVariables.InvFgrad(1,0);
    // MatrixA(0,2)=rElementalVariables.InvFgrad(2,0);
    // MatrixA(1,0)=rElementalVariables.InvFgrad(0,1);
    // MatrixA(1,2)=rElementalVariables.InvFgrad(2,1);
    // MatrixA(2,0)=rElementalVariables.InvFgrad(0,2);
    // MatrixA(2,1)=rElementalVariables.InvFgrad(1,2);

    // MatrixB(0,0)=rElementalVariables.MDGreenLagrangeMaterial[0];  //XX-component;
    // MatrixB(1,1)=rElementalVariables.MDGreenLagrangeMaterial[1];  //YY-component;
    // MatrixB(2,2)=rElementalVariables.MDGreenLagrangeMaterial[2];  //ZZ-component;
    // MatrixB(0,1)=rElementalVariables.MDGreenLagrangeMaterial[3];  //XY-component;
    // MatrixB(1,0)=rElementalVariables.MDGreenLagrangeMaterial[3];  //XY-component;
    // MatrixB(0,2)=rElementalVariables.MDGreenLagrangeMaterial[4];  //ZX-component;
    // MatrixB(2,0)=rElementalVariables.MDGreenLagrangeMaterial[4];  //ZX-component;
    // MatrixB(1,2)=rElementalVariables.MDGreenLagrangeMaterial[5];  //YZ-component;
    // MatrixB(2,1)=rElementalVariables.MDGreenLagrangeMaterial[5];  //YZ-component;

    // noalias(Matrix1)=prod(MatrixB,rElementalVariables.InvFgrad);
    // noalias(Matrix2)=prod(MatrixA,Matrix1);

    // rElementalVariables.SpatialDefRate[0]=Matrix2(0,0);
    // rElementalVariables.SpatialDefRate[1]=Matrix2(1,1);
    // rElementalVariables.SpatialDefRate[2]=Matrix2(2,2);
    // rElementalVariables.SpatialDefRate[3]=Matrix2(0,1);
    // rElementalVariables.SpatialDefRate[4]=Matrix2(0,2);
    // rElementalVariables.SpatialDefRate[5]=Matrix2(1,2);

    rElementalVariables.SpatialDefRate[0] = rElementalVariables.SpatialVelocityGrad(0, 0);
    rElementalVariables.SpatialDefRate[1] = rElementalVariables.SpatialVelocityGrad(1, 1);
    rElementalVariables.SpatialDefRate[2] = rElementalVariables.SpatialVelocityGrad(2, 2);
    rElementalVariables.SpatialDefRate[3] = 0.5 * (rElementalVariables.SpatialVelocityGrad(1, 0) + rElementalVariables.SpatialVelocityGrad(0, 1));
    rElementalVariables.SpatialDefRate[4] = 0.5 * (rElementalVariables.SpatialVelocityGrad(2, 0) + rElementalVariables.SpatialVelocityGrad(0, 2));
    rElementalVariables.SpatialDefRate[5] = 0.5 * (rElementalVariables.SpatialVelocityGrad(2, 1) + rElementalVariables.SpatialVelocityGrad(1, 2));
    // computeElement=CheckStrain3(rElementalVariables.SpatialDefRate,rElementalVariables.SpatialVelocityGrad);

    double aThird = 1.0 / 3.0;
    double dev_X = rElementalVariables.SpatialDefRate[0] -
                   (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1] + rElementalVariables.SpatialDefRate[2]) * aThird;
    double dev_Y = rElementalVariables.SpatialDefRate[1] -
                   (rElementalVariables.SpatialDefRate[0] + rElementalVariables.SpatialDefRate[1] + rElementalVariables.SpatialDefRate[2]) * aThird;
    double dev_Z = rElementalVariables.SpatialDefRate[2] -
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
  bool ThreeStepUpdatedLagrangianElement<TDim>::CalcStrainRate(ElementalVariables &rElementalVariables,
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

    //it computes the material time derivative of the deformation gradient and its jacobian and inverse
    this->CalcVelDefGrad(rDN_DX,
                         rElementalVariables.FgradVel,
                         theta);

    //it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
    this->CalcSpatialVelocityGrad(rElementalVariables.InvFgrad,
                                  rElementalVariables.FgradVel,
                                  rElementalVariables.SpatialVelocityGrad);

    this->CalcVolDefRateFromSpatialVelGrad(rElementalVariables.VolumetricDefRate,
                                           rElementalVariables.SpatialVelocityGrad);

    // this->CalcVolumetricDefRate(rDN_DX,
    // 			      rElementalVariables.VolumetricDefRate,
    // 			      rElementalVariables.InvFgrad,
    //                          theta );

    // //it checks whether tr(l) == div(v)
    // CheckStrain1(rElementalVariables.VolumetricDefRate,
    // 	       rElementalVariables.SpatialVelocityGrad);

    // CheckStrain2(rElementalVariables.SpatialVelocityGrad,
    // 	       rElementalVariables.Fgrad,
    // 	       rElementalVariables.FgradVel);

    //it computes Material time Derivative of Green Lagrange strain tensor in MATERIAL configuration --> [D(E)/Dt]
    this->CalcMDGreenLagrangeMaterial(rElementalVariables.Fgrad,
                                      rElementalVariables.FgradVel,
                                      rElementalVariables.MDGreenLagrangeMaterial);

    //it computes Material time Derivative of Green Lagrange strain tensor in SPATIAL configuration  --> [d]
    this->CalcSpatialDefRate(rElementalVariables.MDGreenLagrangeMaterial,
                             rElementalVariables.InvFgrad,
                             rElementalVariables.SpatialDefRate);

    // computeElement=CheckStrain3(rElementalVariables.SpatialDefRate,
    // 			      rElementalVariables.SpatialVelocityGrad);

    this->CalcDeviatoricInvariant(rElementalVariables.SpatialDefRate,
                                  rElementalVariables.DeviatoricInvariant);

    this->CalcEquivalentStrainRate(rElementalVariables.SpatialDefRate,
                                   rElementalVariables.EquivalentStrainRate);
    return computeElement;
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
  void ThreeStepUpdatedLagrangianElement<TDim>::ComputeMassMatrix(Matrix &rMassMatrix,
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
  void ThreeStepUpdatedLagrangianElement<2>::AddInternalForces(Vector &rRHSVector,
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
  void ThreeStepUpdatedLagrangianElement<3>::AddInternalForces(Vector &rRHSVector,
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
  void ThreeStepUpdatedLagrangianElement<TDim>::ComputeBulkMatrix(Matrix &BulkMatrix,
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
  void ThreeStepUpdatedLagrangianElement<2>::ComputeBulkMatrixConsistent(Matrix &BulkMatrix,
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
  void ThreeStepUpdatedLagrangianElement<3>::ComputeBulkMatrixConsistent(Matrix &BulkMatrix,
                                                                       const double Weight)
  {
    std::cout << "TO IMPLEMENT AND CHECK " << std::endl;
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

  /*
   * Template class definition (this should allow us to compile the desired template instantiations)
   */

  template class ThreeStepUpdatedLagrangianElement<2>;
  template class ThreeStepUpdatedLagrangianElement<3>;

} // namespace Kratos

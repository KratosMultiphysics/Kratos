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
#include "custom_elements/three_step_first_order_updated_lagrangian_element.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

  template <unsigned int TDim>
  Element::Pointer ThreeStepFirstOrderUpdatedLagrangianElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
  {
    KRATOS_TRY;
    
    ThreeStepFirstOrderUpdatedLagrangianElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer(new ThreeStepFirstOrderUpdatedLagrangianElement(NewElement));

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void ThreeStepFirstOrderUpdatedLagrangianElement<TDim>::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
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
      CalculateStandardFSPressureSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
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
  void ThreeStepFirstOrderUpdatedLagrangianElement<TDim>::CalculateFirstVelocitySystem(MatrixType &rLeftHandSideMatrix,
                                                                                       VectorType &rRightHandSideVector,
                                                                                       const ProcessInfo &rCurrentProcessInfo)
  {
    const GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = TDim * NumNodes;

    // Check sizes and initialize
    if (rLeftHandSideMatrix.size1() != LocalSize)
      rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

    if (rRightHandSideVector.size() != LocalSize)
      rRightHandSideVector.resize(LocalSize, false);

    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    MatrixType MassMatrix = ZeroMatrix(LocalSize, LocalSize);
    MatrixType ViscousMatrix = ZeroMatrix(LocalSize, LocalSize);

    double theta_velocity = 1.0;
    double gamma_pressure = 0.0;

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
    {
      const double GaussWeight = GaussWeights[g];
      const ShapeFunctionsType &rN = row(NContainer, g);
      const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

      // Evaluate required variables at the integration point
      double Density = 0.0;
      double Viscosity = 0.0;
      array_1d<double, 3> BodyForce = ZeroVector(3);

      this->EvaluateInPoint(Density, DENSITY, rN);
      this->EvaluateInPoint(BodyForce, VOLUME_ACCELERATION, rN);

      this->EvaluateInPoint(Viscosity, DYNAMIC_VISCOSITY, rN);

      // Add integration point contribution to the local mass matrix
      const double MassCoeff = Density * GaussWeight;
      //AddMomentumMassTerm(MassMatrix, rN, MassCoeff);
      this->ComputeLumpedMassMatrix(MassMatrix, MassCoeff);

      // // Evaluate the pressure and pressure gradient at this point (for the G * P_n term)
      double OldPressure = 0.0;
      this->EvaluateInPoint(OldPressure, PRESSURE, rN, 1);
      double oldPressureContribution = OldPressure * gamma_pressure;

      // Add RHS contributions to the local system equation
      this->AddMomentumRHSTerms(rRightHandSideVector, Density, BodyForce, oldPressureContribution, rN, rDN_DX, GaussWeight);
      // AddExternalForces(rRightHandSideVector, Density, rN, GaussWeight);

      // Add viscous term
      const double ViscousCoeff = Viscosity * GaussWeight;
      this->AddViscousTerm(rLeftHandSideMatrix, rDN_DX, ViscousCoeff, theta_velocity);

      this->AddViscousTerm(ViscousMatrix, rDN_DX, ViscousCoeff, 1.0);
    }

    // Add residual of previous iteration to RHS
    VectorType VelocityValues = ZeroVector(LocalSize);
    this->GetVelocityValues(VelocityValues, 0);
    // noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, VelocityValues);

    VectorType VelocityValuesTheta = ZeroVector(LocalSize);
    this->GetVelocityValues(VelocityValuesTheta, 1);
    VelocityValuesTheta *= (1.0 - theta_velocity);
    noalias(VelocityValuesTheta) += VelocityValues * theta_velocity;
    noalias(rRightHandSideVector) -= prod(ViscousMatrix, VelocityValuesTheta);

    // Add dynamic term
    const double TimeStep = rCurrentProcessInfo[DELTA_TIME];

    // // first order
    VectorType AccelerationValues = ZeroVector(LocalSize);
    noalias(AccelerationValues) += -VelocityValues / TimeStep;
    this->GetVelocityValues(VelocityValues, 1);
    noalias(AccelerationValues) += VelocityValues / TimeStep; //these are negative accelerations
    noalias(rRightHandSideVector) += prod(MassMatrix, AccelerationValues);
    noalias(rLeftHandSideMatrix) += MassMatrix / TimeStep;
  }

  // this is a part of the last iteration of the fractional step
  template <unsigned int TDim>
  void ThreeStepFirstOrderUpdatedLagrangianElement<TDim>::Calculate(const Variable<array_1d<double, 3>> &rVariable,
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
        const ShapeFunctionsType &rN = row(NContainer, g);
        const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

        double Density = 0;
        this->EvaluateInPoint(Density, DENSITY, rN);

        const double TimeStep = rCurrentProcessInfo[DELTA_TIME];
        double timeFactor = TimeStep; //first order accelerations
        const double Coeff = GaussWeights[g] * timeFactor / Density;

        double elementalPressure = 0;
        this->EvaluateInPoint(elementalPressure, PRESSURE, rN);

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
  void ThreeStepFirstOrderUpdatedLagrangianElement<TDim>::CalculateStandardFSPressureSystem(MatrixType &rLeftHandSideMatrix,
                                                                                            VectorType &rRightHandSideVector,
                                                                                            const ProcessInfo &rCurrentProcessInfo)
  {
    GeometryType &rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const double TimeStep = rCurrentProcessInfo[DELTA_TIME];
    double ElemSize = this->ElementSize();

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

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
    {
      const double GaussWeight = GaussWeights[g];
      const ShapeFunctionsType &rN = row(NContainer, g);
      const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];

      // Evaluate required variables at the integration point
      double Density;
      array_1d<double, 3> BodyForce = ZeroVector(3);

      this->EvaluateInPoint(Density, DENSITY, rN);
      this->EvaluateInPoint(BodyForce, VOLUME_ACCELERATION, rN);

      array_1d<double, TDim> OldPressureGradient = ZeroVector(TDim);
      this->EvaluateGradientInPoint(OldPressureGradient, PRESSURE, rDN_DX);

      double Viscosity = 0;
      this->EvaluateInPoint(Viscosity, DYNAMIC_VISCOSITY, rN);

      double DivU = 0;
      this->EvaluateDivergenceInPoint(DivU, VELOCITY, rDN_DX);

      const double LaplacianCoeff = TimeStep / Density;

      // Add convection, stabilization and RHS contributions to the local system equation
      for (SizeType i = 0; i < NumNodes; ++i)
      {
        // LHS contribution
        unsigned int freeSurfaceNodes = 0;
        for (SizeType j = 0; j < NumNodes; ++j)
        {
          double Lij = 0.0;
          for (SizeType d = 0; d < TDim; ++d)
            Lij += rDN_DX(i, d) * rDN_DX(j, d);

          Lij *= LaplacianCoeff;
          rLeftHandSideMatrix(i, j) += GaussWeight * Lij;

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
        // Velocity divergence
        double RHSi = -rN[i] * DivU;

        for (SizeType d = 0; d < TDim; ++d)
        {
          RHSi += -LaplacianCoeff * rDN_DX(i, d) * OldPressureGradient[d];
        }
        rRightHandSideVector[i] += GaussWeight * RHSi;
      }
    }
  }

  /*
   * Template class definition (this should allow us to compile the desired template instantiations)
   */

  template class ThreeStepFirstOrderUpdatedLagrangianElement<2>;
  template class ThreeStepFirstOrderUpdatedLagrangianElement<3>;

} // namespace Kratos

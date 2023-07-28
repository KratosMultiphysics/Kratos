//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_FIC_cut_fem_element.h"
#include "includes/cfd_variables.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include <cmath>

namespace Kratos
{

    /*
     * public TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim> functions
     */

    template <unsigned int TDim>
    Element::Pointer TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
    {

        TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

        NewElement.SetData(this->GetData());
        NewElement.SetFlags(this->GetFlags());

        return Element::Pointer(new TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement(NewElement));
    }

    // template <unsigned int TDim>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::Initialize(const ProcessInfo &rCurrentProcessInfo)
    // {
    //   KRATOS_TRY;

    //   // If we are restarting, the constitutive law will be already defined
    //   if (mpConstitutiveLaw == nullptr)
    //   {
    //     const Properties &r_properties = this->GetProperties();
    //     KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
    //         << "In initialization of Element " << this->Info() << ": No CONSTITUTIVE_LAW defined for property "
    //         << r_properties.Id() << "." << std::endl;
    //     mpConstitutiveLaw = r_properties[CONSTITUTIVE_LAW]->Clone();
    //   }

    //   KRATOS_CATCH("");
    // }

    template <unsigned int TDim>
    int TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        unsigned int ierr = BaseType::Check(rCurrentProcessInfo);
        // TODO: Check distance

        return ierr;

        KRATOS_CATCH("");
    }

    // template <>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<2>::ComputeBoundLHSMatrix(Matrix &BoundLHSMatrix,
    //                                                                                  const ShapeFunctionsType &rN,
    //                                                                                  const double Weight)
    // {
    //   GeometryType &rGeom = this->GetGeometry();
    //   double coeff = 1.0 / 3.0;

    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE))
    //   {
    //     if (rGeom[0].IsNot(INLET))
    //       BoundLHSMatrix(0, 0) += Weight * coeff;
    //     if (rGeom[1].IsNot(INLET))
    //       BoundLHSMatrix(1, 1) += Weight * coeff;
    //   }
    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    //   {
    //     if (rGeom[0].IsNot(INLET))
    //       BoundLHSMatrix(0, 0) += Weight * coeff;
    //     if (rGeom[2].IsNot(INLET))
    //       BoundLHSMatrix(2, 2) += Weight * coeff;
    //   }
    //   if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    //   {
    //     if (rGeom[1].IsNot(INLET))
    //       BoundLHSMatrix(1, 1) += Weight * coeff;
    //     if (rGeom[2].IsNot(INLET))
    //       BoundLHSMatrix(2, 2) += Weight * coeff;
    //   }
    //   // }
    // }

    // template <>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<3>::ComputeBoundLHSMatrix(Matrix &BoundLHSMatrix,
    //                                                                                  const ShapeFunctionsType &rN,
    //                                                                                  const double Weight)
    // {
    //   GeometryType &rGeom = this->GetGeometry();
    //   double coeff = 0.25;

    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    //   {
    //     if (rGeom[0].IsNot(INLET))
    //       BoundLHSMatrix(0, 0) += Weight * coeff;
    //     if (rGeom[1].IsNot(INLET))
    //       BoundLHSMatrix(1, 1) += Weight * coeff;
    //     if (rGeom[2].IsNot(INLET))
    //       BoundLHSMatrix(2, 2) += Weight * coeff;
    //   }
    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    //   {
    //     if (rGeom[0].IsNot(INLET))
    //       BoundLHSMatrix(0, 0) += Weight * coeff;
    //     if (rGeom[1].IsNot(INLET))
    //       BoundLHSMatrix(1, 1) += Weight * coeff;
    //     if (rGeom[3].IsNot(INLET))
    //       BoundLHSMatrix(3, 3) += Weight * coeff;
    //   }
    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    //   {
    //     if (rGeom[0].IsNot(INLET))
    //       BoundLHSMatrix(0, 0) += Weight * coeff;
    //     if (rGeom[2].IsNot(INLET))
    //       BoundLHSMatrix(2, 2) += Weight * coeff;
    //     if (rGeom[3].IsNot(INLET))
    //       BoundLHSMatrix(3, 3) += Weight * coeff;
    //   }
    //   if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    //   {
    //     if (rGeom[1].IsNot(INLET))
    //       BoundLHSMatrix(1, 1) += Weight * coeff;
    //     if (rGeom[2].IsNot(INLET))
    //       BoundLHSMatrix(2, 2) += Weight * coeff;
    //     if (rGeom[3].IsNot(INLET))
    //       BoundLHSMatrix(3, 3) += Weight * coeff;
    //   }
    // }

    // template <>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<2>::ComputeBoundRHSVectorComplete(VectorType &BoundRHSVector,
    //                                                                                          const double TimeStep,
    //                                                                                          const double BoundRHSCoeffAcc,
    //                                                                                          const double BoundRHSCoeffDev,
    //                                                                                          const VectorType SpatialDefRate)
    // {
    //   GeometryType &rGeom = this->GetGeometry();
    //   const double coeff = 1.0 / 3.0;
    //   const double timeFactor = 0.5 / TimeStep;

    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE))
    //   {
    //     array_1d<double, 3> AccA(3, 0.0);
    //     array_1d<double, 3> AccB(3, 0.0);
    //     array_1d<double, 3> MeanAcc(3, 0.0);
    //     array_1d<double, 3> NormalVector(3, 0.0);

    //     this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 1, 2);

    //     double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    //     noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(MeanAcc) = 0.5 * (AccA + AccB);

    //     const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

    //     if (rGeom[0].IsNot(INLET)) // to change into moving wall!!!!!
    //       BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    //     if (rGeom[1].IsNot(INLET))
    //       BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    //   }

    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    //   {

    //     array_1d<double, 3> AccA(3, 0.0);
    //     array_1d<double, 3> AccB(3, 0.0);
    //     array_1d<double, 3> MeanAcc(3, 0.0);
    //     array_1d<double, 3> NormalVector(3, 0.0);

    //     this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 0, 2, 1);

    //     double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    //     noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(MeanAcc) = 0.5 * (AccA + AccB);

    //     const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

    //     if (rGeom[0].IsNot(INLET)) // to change into moving wall!!!!!
    //       BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    //     if (rGeom[2].IsNot(INLET))
    //       BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    //   }

    //   if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    //   {

    //     array_1d<double, 3> AccA(3, 0.0);
    //     array_1d<double, 3> AccB(3, 0.0);
    //     array_1d<double, 3> MeanAcc(3, 0.0);
    //     array_1d<double, 3> NormalVector(3, 0.0);

    //     this->GetOutwardsUnitNormalForTwoPoints(NormalVector, 1, 2, 0);

    //     double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    //     noalias(AccA) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(MeanAcc) = 0.5 * (AccA + AccB);

    //     const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1];

    //     if (rGeom[1].IsNot(INLET))
    //       BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    //     if (rGeom[2].IsNot(INLET))
    //       BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    //   }
    // }

    // template <>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<3>::ComputeBoundRHSVectorComplete(VectorType &BoundRHSVector,
    //                                                                                          const double TimeStep,
    //                                                                                          const double BoundRHSCoeffAcc,
    //                                                                                          const double BoundRHSCoeffDev,
    //                                                                                          const VectorType SpatialDefRate)
    // {
    //   GeometryType &rGeom = this->GetGeometry();
    //   const double coeff = 0.25;
    //   const double timeFactor = 0.5 / TimeStep;
    //   const double one_third = 1.0 / 3.0;

    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    //   {

    //     array_1d<double, 3> AccA(3, 0.0);
    //     array_1d<double, 3> AccB(3, 0.0);
    //     array_1d<double, 3> AccC(3, 0.0);
    //     array_1d<double, 3> MeanAcc(3, 0.0);
    //     array_1d<double, 3> NormalVector(3, 0.0);

    //     this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 1, 2, 3);

    //     double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    //     noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccC) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);

    //     noalias(MeanAcc) = (AccA + AccB + AccC) * one_third;

    //     const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    //     if (rGeom[0].IsNot(INLET))
    //       BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    //     if (rGeom[1].IsNot(INLET))
    //       BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    //     if (rGeom[2].IsNot(INLET))
    //       BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    //   }

    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    //   {

    //     array_1d<double, 3> AccA(3, 0.0);
    //     array_1d<double, 3> AccB(3, 0.0);
    //     array_1d<double, 3> AccC(3, 0.0);
    //     array_1d<double, 3> MeanAcc(3, 0.0);
    //     array_1d<double, 3> NormalVector(3, 0.0);
    //     this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 1, 3, 2);

    //     double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    //     noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccC) = timeFactor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

    //     noalias(MeanAcc) = (AccA + AccB + AccC) * one_third;

    //     const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    //     if (rGeom[0].IsNot(INLET))
    //       BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    //     if (rGeom[1].IsNot(INLET))
    //       BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    //     if (rGeom[3].IsNot(INLET))
    //       BoundRHSVector[3] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    //   }

    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    //   {

    //     array_1d<double, 3> AccA(3, 0.0);
    //     array_1d<double, 3> AccB(3, 0.0);
    //     array_1d<double, 3> AccC(3, 0.0);
    //     array_1d<double, 3> MeanAcc(3, 0.0);
    //     array_1d<double, 3> NormalVector(3, 0.0);

    //     this->GetOutwardsUnitNormalForThreePoints(NormalVector, 0, 2, 3, 1);

    //     double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    //     noalias(AccA) = timeFactor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccC) = timeFactor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

    //     noalias(MeanAcc) = (AccA + AccB + AccC) * one_third;

    //     const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    //     if (rGeom[0].IsNot(INLET))
    //       BoundRHSVector[0] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    //     if (rGeom[2].IsNot(INLET))
    //       BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    //     if (rGeom[3].IsNot(INLET))
    //       BoundRHSVector[3] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    //   }

    //   if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    //   {

    //     array_1d<double, 3> AccA(3, 0.0);
    //     array_1d<double, 3> AccB(3, 0.0);
    //     array_1d<double, 3> AccC(3, 0.0);
    //     array_1d<double, 3> MeanAcc(3, 0.0);
    //     array_1d<double, 3> NormalVector(3, 0.0);

    //     this->GetOutwardsUnitNormalForThreePoints(NormalVector, 1, 2, 3, 0);

    //     double SpatialDefRateNormalProjection = this->CalcNormalProjectionDefRate(SpatialDefRate, NormalVector);

    //     noalias(AccA) = timeFactor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = timeFactor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccC) = timeFactor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);

    //     noalias(MeanAcc) = (AccA + AccB + AccC) * one_third;

    //     const double accelerationsNormalProjection = MeanAcc[0] * NormalVector[0] + MeanAcc[1] * NormalVector[1] + MeanAcc[2] * NormalVector[2];

    //     if (rGeom[1].IsNot(INLET))
    //       BoundRHSVector[1] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    //     if (rGeom[2].IsNot(INLET))
    //       BoundRHSVector[2] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);

    //     if (rGeom[3].IsNot(INLET))
    //       BoundRHSVector[3] += coeff * (BoundRHSCoeffAcc * accelerationsNormalProjection + BoundRHSCoeffDev * SpatialDefRateNormalProjection);
    //   }
    // }

    // template <>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<2>::ComputeBoundRHSVector(VectorType &BoundRHSVector,
    //                                                                                  const ShapeFunctionsType &rN,
    //                                                                                  const double TimeStep,
    //                                                                                  const double BoundRHSCoeffAcc,
    //                                                                                  const double BoundRHSCoeffDev)
    // {
    //   GeometryType &rGeom = this->GetGeometry();
    //   // const SizeType NumNodes = rGeom.PointsNumber();
    //   array_1d<double, 3> AccA(3, 0.0);
    //   array_1d<double, 3> AccB(3, 0.0);

    //   // for (SizeType i = 0; i < (NumNodes-1); i++)
    //   //   {
    //   // 	for (SizeType j = (i+1); j < NumNodes; j++)
    //   // 	  {
    //   // 	    if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE)){
    //   // 	      AccA= 0.5/TimeStep*(rGeom[i].FastGetSolutionStepValue(VELOCITY,0)-rGeom[i].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[i].FastGetSolutionStepValue(ACCELERATION,1);
    //   // 	      AccB= 0.5/TimeStep*(rGeom[j].FastGetSolutionStepValue(VELOCITY,0)-rGeom[j].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[j].FastGetSolutionStepValue(ACCELERATION,1);
    //   // 	      const array_1d<double, 3> &NormalA    = rGeom[i].FastGetSolutionStepValue(NORMAL);
    //   // 	      const array_1d<double, 3> &NormalB    = rGeom[j].FastGetSolutionStepValue(NORMAL);
    //   // 	      double coeff=3.0;
    //   // 	      if(rGeom[i].IsNot(INLET)) //to change into moving wall!!!!!
    //   // 		BoundRHSVector[i] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) +
    //   // 				      BoundRHSCoeffDev)/coeff ;
    //   // 	      if(rGeom[j].IsNot(INLET))
    //   // 		BoundRHSVector[j] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) +
    //   // 				      BoundRHSCoeffDev)/coeff ;
    //   // 	    }
    //   // 	  }

    //   //   }

    //   // for (SizeType i = 0; i < (NumNodes-1); i++)
    //   //   {
    //   // 	for (SizeType j = (i+1); j < NumNodes; j++)
    //   // 	  {
    //   // 	    if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE)){
    //   // 	      AccA= 0.5/TimeStep*(rGeom[i].FastGetSolutionStepValue(VELOCITY,0)-rGeom[i].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[i].FastGetSolutionStepValue(ACCELERATION,1);
    //   // 	      AccB= 0.5/TimeStep*(rGeom[j].FastGetSolutionStepValue(VELOCITY,0)-rGeom[j].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[j].FastGetSolutionStepValue(ACCELERATION,1);
    //   // 	      const array_1d<double, 3> &NormalA    = rGeom[i].FastGetSolutionStepValue(NORMAL);
    //   // 	      const array_1d<double, 3> &NormalB    = rGeom[j].FastGetSolutionStepValue(NORMAL);
    //   // 	      if(rGeom[i].IsNot(INLET))
    //   // 		BoundRHSVector[i] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) +
    //   // 				      BoundRHSCoeffDev) * rN[i];
    //   // 	      if(rGeom[j].IsNot(INLET))
    //   // 		BoundRHSVector[j] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) +
    //   // 				      BoundRHSCoeffDev) * rN[j] ;
    //   // 	    }
    //   // 	  }

    //   //   }

    //   const double factor = 0.5 / TimeStep;
    //   const double one_third = 1.0 / 3.0;

    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE))
    //   {
    //     noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    //     // noalias(AccA)=rGeom[0].FastGetSolutionStepValue(ACCELERATION,0);
    //     // noalias(AccB)=rGeom[1].FastGetSolutionStepValue(ACCELERATION,0);
    //     const array_1d<double, 3> &NormalA = rGeom[0].FastGetSolutionStepValue(NORMAL);
    //     const array_1d<double, 3> &NormalB = rGeom[1].FastGetSolutionStepValue(NORMAL);
    //     if (rGeom[0].IsNot(INLET)) // to change into moving wall!!!!!
    //       BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1]) + BoundRHSCoeffDev);
    //     if (rGeom[1].IsNot(INLET))
    //       BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1]) + BoundRHSCoeffDev);
    //   }
    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    //   {
    //     noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    //     const array_1d<double, 3> &NormalA = rGeom[0].FastGetSolutionStepValue(NORMAL);
    //     const array_1d<double, 3> &NormalB = rGeom[2].FastGetSolutionStepValue(NORMAL);
    //     if (rGeom[0].IsNot(INLET))
    //       BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1]) + BoundRHSCoeffDev);
    //     if (rGeom[2].IsNot(INLET))
    //       BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1]) + BoundRHSCoeffDev);
    //   }
    //   if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    //   {
    //     noalias(AccA) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    //     const array_1d<double, 3> &NormalA = rGeom[1].FastGetSolutionStepValue(NORMAL);
    //     const array_1d<double, 3> &NormalB = rGeom[2].FastGetSolutionStepValue(NORMAL);
    //     if (rGeom[1].IsNot(INLET))
    //       BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1]) + BoundRHSCoeffDev);
    //     if (rGeom[2].IsNot(INLET))
    //       BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1]) + BoundRHSCoeffDev);
    //   }
    // }

    // template <>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<3>::ComputeBoundRHSVector(VectorType &BoundRHSVector,
    //                                                                                  const ShapeFunctionsType &rN,
    //                                                                                  const double TimeStep,
    //                                                                                  const double BoundRHSCoeffAcc,
    //                                                                                  const double BoundRHSCoeffDev)
    // {
    //   GeometryType &rGeom = this->GetGeometry();
    //   // const SizeType NumNodes = rGeom.PointsNumber();
    //   array_1d<double, 3> AccA(3, 0.0);
    //   array_1d<double, 3> AccB(3, 0.0);
    //   array_1d<double, 3> AccC(3, 0.0);

    //   // for (SizeType i = 0; i < (NumNodes-2); i++)
    //   //   {
    //   // 	for (SizeType j = (i+1); j < (NumNodes-1); j++)
    //   // 	  {
    //   // 	    for (SizeType k = (j+1); k < NumNodes; k++)
    //   // 	      {
    //   // 		if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE) && rGeom[k].Is(FREE_SURFACE)){
    //   // 		  AccA= 0.5/TimeStep*(rGeom[i].FastGetSolutionStepValue(VELOCITY,0)-rGeom[i].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[i].FastGetSolutionStepValue(ACCELERATION,1);
    //   // 		  AccB= 0.5/TimeStep*(rGeom[j].FastGetSolutionStepValue(VELOCITY,0)-rGeom[j].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[j].FastGetSolutionStepValue(ACCELERATION,1);
    //   // 		  AccC= 0.5/TimeStep*(rGeom[k].FastGetSolutionStepValue(VELOCITY,0)-rGeom[k].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[k].FastGetSolutionStepValue(ACCELERATION,1);

    //   // 		  const array_1d<double, 3> &NormalA    = rGeom[i].FastGetSolutionStepValue(NORMAL);
    //   // 		  const array_1d<double, 3> &NormalB    = rGeom[j].FastGetSolutionStepValue(NORMAL);
    //   // 		  const array_1d<double, 3> &NormalC    = rGeom[k].FastGetSolutionStepValue(NORMAL);
    //   // 		  if(rGeom[i].IsNot(INLET))
    //   // 		    BoundRHSVector[i] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
    //   // 					  BoundRHSCoeffDev) * rN[i];
    //   // 		  if(rGeom[j].IsNot(INLET))
    //   // 		    BoundRHSVector[j] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
    //   // 					  BoundRHSCoeffDev) * rN[j] ;
    //   // 		  if(rGeom[k].IsNot(INLET))
    //   // 		    BoundRHSVector[k] += (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
    //   // 					  BoundRHSCoeffDev) * rN[k] ;
    //   // 		}
    //   // 	      }
    //   // 	  }

    //   //   }
    //   const double factor = 0.5 / TimeStep;

    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE))
    //   {
    //     noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccC) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    //     const array_1d<double, 3> &NormalA = rGeom[0].FastGetSolutionStepValue(NORMAL);
    //     const array_1d<double, 3> &NormalB = rGeom[1].FastGetSolutionStepValue(NORMAL);
    //     const array_1d<double, 3> &NormalC = rGeom[2].FastGetSolutionStepValue(NORMAL);
    //     if (rGeom[0].IsNot(INLET))
    //       BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1] + AccA[2] * NormalA[2]) +
    //                                    BoundRHSCoeffDev);
    //     if (rGeom[1].IsNot(INLET))
    //       BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1] + AccB[2] * NormalB[2]) +
    //                                    BoundRHSCoeffDev);
    //     if (rGeom[2].IsNot(INLET))
    //       BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc * (AccC[0] * NormalC[0] + AccC[1] * NormalC[1] + AccC[2] * NormalC[2]) +
    //                                    BoundRHSCoeffDev);
    //   }
    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    //   {
    //     noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccC) = factor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);
    //     const array_1d<double, 3> &NormalA = rGeom[0].FastGetSolutionStepValue(NORMAL);
    //     const array_1d<double, 3> &NormalB = rGeom[1].FastGetSolutionStepValue(NORMAL);
    //     const array_1d<double, 3> &NormalC = rGeom[3].FastGetSolutionStepValue(NORMAL);
    //     if (rGeom[0].IsNot(INLET))
    //       BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1] + AccA[2] * NormalA[2]) +
    //                                    BoundRHSCoeffDev);
    //     if (rGeom[1].IsNot(INLET))
    //       BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1] + AccB[2] * NormalB[2]) +
    //                                    BoundRHSCoeffDev);
    //     if (rGeom[3].IsNot(INLET))
    //       BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc * (AccC[0] * NormalC[0] + AccC[1] * NormalC[1] + AccC[2] * NormalC[2]) +
    //                                    BoundRHSCoeffDev);
    //   }
    //   if (rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    //   {
    //     noalias(AccA) = factor * (rGeom[0].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[0].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccC) = factor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);
    //     const array_1d<double, 3> &NormalA = rGeom[0].FastGetSolutionStepValue(NORMAL);
    //     const array_1d<double, 3> &NormalB = rGeom[2].FastGetSolutionStepValue(NORMAL);
    //     const array_1d<double, 3> &NormalC = rGeom[3].FastGetSolutionStepValue(NORMAL);
    //     if (rGeom[0].IsNot(INLET))
    //       BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1] + AccA[2] * NormalA[2]) +
    //                                    BoundRHSCoeffDev);
    //     if (rGeom[2].IsNot(INLET))
    //       BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1] + AccB[2] * NormalB[2]) +
    //                                    BoundRHSCoeffDev);
    //     if (rGeom[3].IsNot(INLET))
    //       BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc * (AccC[0] * NormalC[0] + AccC[1] * NormalC[1] + AccC[2] * NormalC[2]) +
    //                                    BoundRHSCoeffDev);
    //   }
    //   if (rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE) && rGeom[3].Is(FREE_SURFACE))
    //   {
    //     noalias(AccA) = factor * (rGeom[1].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[1].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccB) = factor * (rGeom[2].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[2].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION, 1);
    //     noalias(AccC) = factor * (rGeom[3].FastGetSolutionStepValue(VELOCITY, 0) - rGeom[3].FastGetSolutionStepValue(VELOCITY, 1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION, 1);
    //     const array_1d<double, 3> &NormalA = rGeom[1].FastGetSolutionStepValue(NORMAL);
    //     const array_1d<double, 3> &NormalB = rGeom[2].FastGetSolutionStepValue(NORMAL);
    //     const array_1d<double, 3> &NormalC = rGeom[3].FastGetSolutionStepValue(NORMAL);
    //     if (rGeom[1].IsNot(INLET))
    //       BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc * (AccA[0] * NormalA[0] + AccA[1] * NormalA[1] + AccA[2] * NormalA[2]) +
    //                                    BoundRHSCoeffDev);
    //     if (rGeom[2].IsNot(INLET))
    //       BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc * (AccB[0] * NormalB[0] + AccB[1] * NormalB[1] + AccB[2] * NormalB[2]) +
    //                                    BoundRHSCoeffDev);
    //     if (rGeom[3].IsNot(INLET))
    //       BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc * (AccC[0] * NormalC[0] + AccC[1] * NormalC[1] + AccC[2] * NormalC[2]) +
    //                                    BoundRHSCoeffDev);
    //   }
    // }

    // template <unsigned int TDim>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateTauFIC(double &Tau,
    //                                                                               double ElemSize,
    //                                                                               const double Density,
    //                                                                               const double Viscosity,
    //                                                                               const ProcessInfo &rCurrentProcessInfo)
    // {
    //   double DeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
    //   if (rCurrentProcessInfo.GetValue(DELTA_TIME) < rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME))
    //   {
    //     DeltaTime = 0.5 * rCurrentProcessInfo.GetValue(DELTA_TIME) + 0.5 * rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME);
    //   }

    //   double MeanVelocity = 0;
    //   this->CalcMeanVelocityNorm(MeanVelocity, 0);

    //   // Tau = 1.0 / (2.0 * Density *(0.5 * MeanVelocity / ElemSize + 0.5/DeltaTime) +  8.0 * Viscosity / (ElemSize * ElemSize) );
    //   Tau = (ElemSize * ElemSize * DeltaTime) / (Density * MeanVelocity * DeltaTime * ElemSize + Density * ElemSize * ElemSize + 8.0 * Viscosity * DeltaTime);

    //   const double tolerance = 1.0e-13;
    //   if (MeanVelocity < tolerance)
    //   {
    //     Tau = 0;
    //   }
    // }

    // template <unsigned int TDim>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::AddStabilizationMatrixLHS(MatrixType &rLeftHandSideMatrix,
    //                                                                                         Matrix &BulkAccMatrix,
    //                                                                                         const ShapeFunctionsType &rN,
    //                                                                                         const double Weight)
    // {
    //   const SizeType NumNodes = this->GetGeometry().PointsNumber();

    //   if (BulkAccMatrix.size1() != NumNodes)
    //     BulkAccMatrix.resize(NumNodes, NumNodes, false);

    //   noalias(BulkAccMatrix) = ZeroMatrix(NumNodes, NumNodes);
    //   for (SizeType i = 0; i < NumNodes; ++i)
    //   {
    //     // LHS contribution
    //     for (SizeType j = 0; j < NumNodes; ++j)
    //     {
    //       double Mij = 0.0;
    //       Mij = Weight * rN[i] * rN[j];
    //       BulkAccMatrix(i, j) += Mij;
    //     }
    //   }
    //   rLeftHandSideMatrix += BulkAccMatrix;
    // }

    // template <unsigned int TDim>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::ComputeStabLaplacianMatrix(MatrixType &StabLaplacianMatrix,
    //                                                                                          const ShapeFunctionDerivativesType &rDN_DX,
    //                                                                                          const double Weight)

    // {
    //   // LHS contribution
    //   const SizeType NumNodes = this->GetGeometry().PointsNumber();
    //   for (SizeType i = 0; i < NumNodes; ++i)
    //   {
    //     for (SizeType j = 0; j < NumNodes; ++j)
    //     {
    //       double Lij = 0.0;
    //       for (SizeType d = 0; d < TDim; ++d)
    //       {
    //         Lij += rDN_DX(i, d) * rDN_DX(j, d);
    //       }
    //       StabLaplacianMatrix(i, j) += Weight * Lij;
    //     }
    //   }
    // }

    // template <unsigned int TDim>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::AddStabilizationNodalTermsLHS(MatrixType &rLeftHandSideMatrix,
    //                                                                                             const double Tau,
    //                                                                                             const double Weight,
    //                                                                                             const ShapeFunctionDerivativesType &rDN_DX,
    //                                                                                             const SizeType i)
    // {
    //   // LHS contribution
    //   const SizeType NumNodes = this->GetGeometry().PointsNumber();
    //   for (SizeType j = 0; j < NumNodes; ++j)
    //   {
    //     double Lij = 0.0;
    //     for (SizeType d = 0; d < TDim; ++d)
    //     {
    //       Lij += rDN_DX(i, d) * rDN_DX(j, d);
    //     }
    //     Lij *= Tau;

    //     rLeftHandSideMatrix(i, j) += Weight * Lij;
    //   }
    // }

    // template <unsigned int TDim>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::AddStabilizationNodalTermsRHS(VectorType &rRightHandSideVector,
    //                                                                                             const double Tau,
    //                                                                                             const double Density,
    //                                                                                             const double Weight,
    //                                                                                             const ShapeFunctionDerivativesType &rDN_DX,
    //                                                                                             const SizeType i)
    // {

    //   double RHSi = 0;
    //   if (this->GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION))
    //   { // it must be checked once at the begining only
    //     array_1d<double, 3> &VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

    //     // double posX=(this->GetGeometry()[0].X() + this->GetGeometry()[1].X() + this->GetGeometry()[2].X())/3.0;

    //     // double posY=(this->GetGeometry()[0].Y() + this->GetGeometry()[1].Y() + this->GetGeometry()[2].Y())/3.0;

    //     // double coeffX =(12.0-24.0*posY)*pow(posX,4);

    //     // coeffX += (-24.0+48.0*posY)*pow(posX,3);

    //     // coeffX += (-48.0*posY+72.0*pow(posY,2)-48.0*pow(posY,3)+12.0)*pow(posX,2);

    //     // coeffX += (-2.0+24.0*posY-72.0*pow(posY,2)+48.0*pow(posY,3))*posX;

    //     // coeffX += 1.0-4.0*posY+12.0*pow(posY,2)-8.0*pow(posY,3);

    //     // double coeffY =(8.0-48.0*posY+48.0*pow(posY,2))*pow(posX,3);

    //     // coeffY += (-12.0+72.0*posY-72.0*pow(posY,2))*pow(posX,2);

    //     // coeffY += (4.0-24.0*posY+48.0*pow(posY,2)-48.0*pow(posY,3)+24.0*pow(posY,4))*posX;

    //     // coeffY += -12.0*pow(posY,2)+24.0*pow(posY,3)-12.0*pow(posY,4);

    //     // RHSi += - rDN_DX(i,0) * Tau * ( Density * VolumeAcceleration[0]*coeffX );

    //     // RHSi += - rDN_DX(i,1) * Tau * ( Density * VolumeAcceleration[1]*coeffY );

    //     for (SizeType d = 0; d < TDim; ++d)
    //     {
    //       RHSi += -rDN_DX(i, d) * Tau * (Density * VolumeAcceleration[d]);
    //     }
    //   }
    //   rRightHandSideVector[i] += Weight * RHSi;
    // }

    // template <unsigned int TDim>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType &rLeftHandSideMatrix,
    //                                                                                                     VectorType &rRightHandSideVector,
    //                                                                                                     const ProcessInfo &rCurrentProcessInfo)
    // {

    //   GeometryType &rGeom = this->GetGeometry();
    //   const unsigned int NumNodes = rGeom.PointsNumber();

    //   // Check sizes and initialize
    //   if (rLeftHandSideMatrix.size1() != NumNodes)
    //     rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

    //   noalias(rLeftHandSideMatrix) = ZeroMatrix(NumNodes, NumNodes);
    //   MatrixType LaplacianMatrix = ZeroMatrix(NumNodes, NumNodes);

    //   if (rRightHandSideVector.size() != NumNodes)
    //     rRightHandSideVector.resize(NumNodes, false);

    //   noalias(rRightHandSideVector) = ZeroVector(NumNodes);

    //   // Shape functions and integration points
    //   ShapeFunctionDerivativesArrayType DN_DX;
    //   Matrix NContainer;
    //   VectorType GaussWeights;
    //   this->CalculateGeometryData(DN_DX, NContainer, GaussWeights);
    //   const unsigned int NumGauss = GaussWeights.size();

    //   double TimeStep = rCurrentProcessInfo[DELTA_TIME];
    //   double theta = this->GetThetaContinuity();
    //   double ElemSize = this->ElementSize();

    //   ElementalVariables rElementalVariables;
    //   this->InitializeElementalVariables(rElementalVariables);

    //   double maxViscousValueForStabilization = 0.1;
    //   double Density = this->mMaterialDensity;
    //   double VolumetricCoeff = this->mMaterialVolumetricCoefficient;
    //   double DeviatoricCoeff = this->mMaterialDeviatoricCoefficient;

    //   if (DeviatoricCoeff > maxViscousValueForStabilization)
    //   {
    //     DeviatoricCoeff = maxViscousValueForStabilization;
    //   }

    //   VectorType NewRhsLaplacian = ZeroVector(NumNodes);

    //   double Tau = 0;
    //   this->CalculateTauFIC(Tau, ElemSize, Density, DeviatoricCoeff, rCurrentProcessInfo);

    //   double totalVolume = 0;
    //   bool computeElement = false;
    //   // Loop on integration points
    //   for (unsigned int g = 0; g < NumGauss; ++g)
    //   {
    //     const double GaussWeight = GaussWeights[g];
    //     totalVolume += GaussWeight;
    //     const ShapeFunctionsType &N = row(NContainer, g);
    //     const ShapeFunctionDerivativesType &rDN_DX = DN_DX[g];
    //     computeElement = this->CalcCompleteStrainRate(rElementalVariables, rCurrentProcessInfo, rDN_DX, theta);

    //     if (computeElement == true && this->IsNot(BLOCKED) && this->IsNot(ISOLATED))
    //     {

    //       double BoundLHSCoeff = Tau * 4.0 * GaussWeight / (ElemSize * ElemSize);
    //       // if constexpr (TDim==3){
    //       //   BoundLHSCoeff=Tau*2*GaussWeight/(0.81649658*ElemSize*ElemSize);
    //       // }

    //       this->ComputeBoundLHSMatrix(rLeftHandSideMatrix, N, BoundLHSCoeff);

    //       double BoundRHSCoeffAcc = Tau * Density * 2 * GaussWeight / ElemSize;
    //       double BoundRHSCoeffDev = Tau * 8.0 * DeviatoricCoeff * GaussWeight / (ElemSize * ElemSize);
    //       // double NProjSpatialDefRate=this->CalcNormalProjectionDefRate(rElementalVariables.SpatialDefRate);
    //       // double BoundRHSCoeffDev=Tau*8.0*NProjSpatialDefRate*DeviatoricCoeff*GaussWeight/(ElemSize*ElemSize);
    //       // this->ComputeBoundRHSVector(rRightHandSideVector,N,TimeStep,BoundRHSCoeffAcc,BoundRHSCoeffDev);
    //       this->ComputeBoundRHSVectorComplete(rRightHandSideVector, TimeStep, BoundRHSCoeffAcc, BoundRHSCoeffDev, rElementalVariables.SpatialDefRate);

    //       double StabLaplacianWeight = Tau * GaussWeight;
    //       this->ComputeStabLaplacianMatrix(LaplacianMatrix, rDN_DX, StabLaplacianWeight);

    //       array_1d<double, TDim> OldPressureGradient = ZeroVector(TDim);
    //       this->EvaluateGradientInPoint(OldPressureGradient, PRESSURE, rDN_DX);
    //       // KRATOS_WATCH(OldPressureGradient);

    //       for (SizeType i = 0; i < NumNodes; ++i)
    //       {
    //         // RHS contribution
    //         // Velocity divergence
    //         rRightHandSideVector[i] += GaussWeight * N[i] * rElementalVariables.VolumetricDefRate;
    //         this->AddStabilizationNodalTermsRHS(rRightHandSideVector, Tau, Density, GaussWeight, rDN_DX, i);
    //         double laplacianRHSi = 0;
    //         for (SizeType d = 0; d < TDim; ++d)
    //         {
    //           laplacianRHSi += StabLaplacianWeight * rDN_DX(i, d) * OldPressureGradient[d];
    //         }
    //         rRightHandSideVector[i] += -laplacianRHSi;
    //         // NewRhsLaplacian[i] += -laplacianRHSi;
    //       }
    //     }
    //   }

    //   if (computeElement == true && this->IsNot(BLOCKED) && this->IsNot(ISOLATED))
    //   {

    //     VectorType PressureValues = ZeroVector(NumNodes);
    //     VectorType PressureValuesForRHS = ZeroVector(NumNodes);
    //     this->GetPressureValues(PressureValuesForRHS, 0);
    //     // the LHS matrix up to now just contains the laplacian term and the bound term
    //     noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, PressureValuesForRHS);
    //     rLeftHandSideMatrix += LaplacianMatrix;
    //     // noalias(rRightHandSideVector) -= prod(LaplacianMatrix, PressureValuesForRHS);

    //     // VectorType RhsLaplacian = ZeroVector(NumNodes);

    //     // RhsLaplacian = -prod(LaplacianMatrix, PressureValuesForRHS);

    //     // VectorType differenceRhsLaplacian = ZeroVector(NumNodes);
    //     // differenceRhsLaplacian = RhsLaplacian - NewRhsLaplacian;
    //     // // KRATOS_WATCH(RhsLaplacian);
    //     // // KRATOS_WATCH(NewRhsLaplacian);
    //     // KRATOS_WATCH(differenceRhsLaplacian);

    //     this->GetPressureValues(PressureValues, 1);
    //     noalias(PressureValuesForRHS) += -PressureValues;
    //     MatrixType BulkMatrix = ZeroMatrix(NumNodes, NumNodes);
    //     MatrixType BulkMatrixConsistent = ZeroMatrix(NumNodes, NumNodes);
    //     double lumpedBulkCoeff = totalVolume / (VolumetricCoeff);
    //     double lumpedBulkStabCoeff = lumpedBulkCoeff * Tau * Density / TimeStep;

    //     this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkCoeff);
    //     noalias(rLeftHandSideMatrix) += BulkMatrix;
    //     // noalias(rLeftHandSideMatrix)+=BulkMatrixConsistent;
    //     noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);
    //     // noalias(rRightHandSideVector) -=prod(BulkMatrixConsistent,PressureValuesForRHS);

    //     this->GetPressureVelocityValues(PressureValues, 0);
    //     noalias(PressureValuesForRHS) += -PressureValues * TimeStep;
    //     noalias(BulkMatrix) = ZeroMatrix(NumNodes, NumNodes);
    //     this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkStabCoeff);
    //     noalias(rLeftHandSideMatrix) += BulkMatrix;
    //     // noalias(rLeftHandSideMatrix)+=BulkMatrixConsistent;
    //     noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);
    //     // noalias(rRightHandSideVector) -=prod(BulkMatrixConsistent,PressureValuesForRHS);
    //   }
    //   else if (this->IsNot(BLOCKED) && this->IsNot(ISOLATED))
    //   {
    //     double lumpedBulkCoeff = totalVolume * Tau * Density / (TimeStep * VolumetricCoeff);
    //     MatrixType BulkVelMatrixLump = ZeroMatrix(NumNodes, NumNodes);
    //     this->ComputeBulkMatrixLump(BulkVelMatrixLump, lumpedBulkCoeff);
    //     noalias(rLeftHandSideMatrix) += BulkVelMatrixLump;
    //     VectorType PressureValues = ZeroVector(NumNodes);
    //     VectorType PressureValuesForRHS = ZeroVector(NumNodes);
    //     this->GetPressureValues(PressureValuesForRHS, 0);
    //     this->GetPressureValues(PressureValues, 1);
    //     noalias(PressureValuesForRHS) += -PressureValues;
    //     noalias(rRightHandSideVector) -= prod(BulkVelMatrixLump, PressureValuesForRHS);
    //   }
    //   else if (this->Is(BLOCKED) && this->IsNot(ISOLATED))
    //   {
    //     VectorType PressureValues = ZeroVector(NumNodes);
    //     VectorType PressureValuesForRHS = ZeroVector(NumNodes);
    //     this->GetPressureValues(PressureValuesForRHS, 0);
    //     // the LHS matrix up to now is void

    //     this->GetPressureValues(PressureValues, 1);
    //     noalias(PressureValuesForRHS) += -PressureValues;
    //     MatrixType BulkMatrix = ZeroMatrix(NumNodes, NumNodes);
    //     double lumpedBulkCoeff = totalVolume / (VolumetricCoeff);

    //     this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkCoeff);
    //     noalias(rLeftHandSideMatrix) += BulkMatrix;
    //     noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);
    //   }
    //   else if (this->Is(ISOLATED))
    //   {
    //     // VectorType PressureValuesForRHS = ZeroVector(NumNodes);
    //     MatrixType BulkMatrix = ZeroMatrix(NumNodes, NumNodes);
    //     double lumpedBulkCoeff = totalVolume / (VolumetricCoeff);

    //     this->ComputeBulkMatrixLump(BulkMatrix, lumpedBulkCoeff);
    //     noalias(rLeftHandSideMatrix) += BulkMatrix;
    //     // noalias(rRightHandSideVector) -= prod(BulkMatrix, PressureValuesForRHS);
    //   }
    // }

    // template <unsigned int TDim>
    // void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::GetPressureAccelerationValues(Vector &rValues,
    //                                                                                             const int Step)
    // {
    //   GeometryType &rGeom = this->GetGeometry();
    //   const SizeType NumNodes = rGeom.PointsNumber();

    //   if (rValues.size() != NumNodes)
    //     rValues.resize(NumNodes, false);

    //   for (SizeType i = 0; i < NumNodes; ++i)
    //   {
    //     rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE_ACCELERATION, Step);
    //   }
    // }

    template <unsigned int TDim>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateLocalMomentumEquations(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo)
    {
        // Volume Navier-Stokes contribution
        // Note that this uses the CalculateGeometryData below, meaning that if it is cut, it already does the subintegration
        BaseType::CalculateLocalMomentumEquations(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

        // If intersected, add the boundary contribution
        if (IsCut())
        {
            // Calculate intersection Gauss points geometry data
            Matrix interface_N;
            ShapeFunctionDerivativesArrayType interface_DN_DX;
            Vector interface_gauss_weights;
            ModifiedShapeFunctions::AreaNormalsContainerType interface_unit_normals;
            CalculateIntersectionGeometryData(interface_DN_DX, interface_N, interface_gauss_weights, interface_unit_normals);

            // Create an auxiliary elemental data container for the boundary terms integration
            ElementalVariables elemental_variables;
            this->InitializeElementalVariables(elemental_variables);

            // Get other data
            const double kappa = rCurrentProcessInfo[PENALTY_COEFFICIENT];
            KRATOS_ERROR_IF(kappa < 1.0e-12) << "'PENALTY_COEFFICIENT' is zero." << std::endl;
            const double h = this->ElementSize();

            // Check if the element is "slip"
            // Note that if one node is flagged as SLIP the entire element is considered slip
            bool is_slip = false;
            const auto &r_geom = this->GetGeometry();
            for (const auto& r_node : r_geom) {
                if (r_node.Is(SLIP)) {
                    is_slip = true;
                    break;
                }
            }

            // Interface Gauss points loop
            array_1d<double,3> vel_gauss;
            const double rho = this->mMaterialDensity;
            const double dt = rCurrentProcessInfo[DELTA_TIME];
            const double theta = this->GetThetaMomentum();
            array_1d<double, TDim> proj_dev_stress;
            const std::size_t n_nodes = r_geom.PointsNumber();
            const unsigned int n_int_gauss_pts = interface_gauss_weights.size();
            for (unsigned int g = 0; g < n_int_gauss_pts; g++)
            {
                // Get interface Gauss point data
                const double g_weight = interface_gauss_weights[g];
                const auto g_DN_DX = interface_DN_DX[g];
                const auto g_shape_functions = row(interface_N, g);
                const auto &r_g_unit_normal = interface_unit_normals[g];

                // Calculate the mechanical response at the interface Gauss point to get the viscous stress
                this->CalcMechanicsUpdated(elemental_variables, rCurrentProcessInfo, g_DN_DX);

                auto &r_strain_vector = elemental_variables.SpatialDefRate;
                auto &r_dev_stress_vector = elemental_variables.UpdatedDeviatoricCauchyStress;
                auto &r_constitutive_matrix = elemental_variables.ConstitutiveMatrix;

                auto p_cons_law = this->GetProperties().GetValue(CONSTITUTIVE_LAW);
                auto constitutive_law_values = ConstitutiveLaw::Parameters(
                    r_geom,
                    this->GetProperties(),
                    rCurrentProcessInfo);

                auto &r_constitutive_law_options = constitutive_law_values.GetOptions();
                r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
                r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

                constitutive_law_values.SetShapeFunctionsValues(g_shape_functions);
                constitutive_law_values.SetStrainVector(r_strain_vector);
                constitutive_law_values.SetStressVector(r_dev_stress_vector);
                constitutive_law_values.SetConstitutiveMatrix(r_constitutive_matrix);

                p_cons_law->CalculateMaterialResponseCauchy(constitutive_law_values);
                this->UpdateStressTensor(elemental_variables);

                // Take dynamic viscosity from the bottom right corner of the constitutive matrix
                double mu = r_constitutive_matrix(StrainSize-1, StrainSize-1);
                const double max_mu_value = 10000; // TODO: check it 
                if (mu > max_mu_value)
                {
                    mu = max_mu_value;
                }
                // Interpolate the pressure and velocity at the interface Gauss point to calculate the isochoric stress
                double pres_gauss = 0.0;
                vel_gauss = ZeroVector(3);
                for (std::size_t j = 0; j < n_nodes; ++j)
                {
                    const double p_0 = r_geom[j].FastGetSolutionStepValue(PRESSURE);
                    const double p_1 = r_geom[j].FastGetSolutionStepValue(PRESSURE,1);
                    pres_gauss += g_shape_functions[j] * (theta * p_0 + (1 - theta) * p_1);
                    const auto &r_vel_j_0 = r_geom[j].FastGetSolutionStepValue(VELOCITY);
                    const auto &r_vel_j_1 = r_geom[j].FastGetSolutionStepValue(VELOCITY,1);
                    noalias(vel_gauss) += g_shape_functions[j] * (theta * r_vel_j_0 + (1.0 - theta) * r_vel_j_1) ;
                }

                // Navier-Stokes traction boundary term
                VoigtStressNormalProjection(r_dev_stress_vector, r_g_unit_normal, proj_dev_stress);
                for (std::size_t i = 0; i < n_nodes; ++i)
                {
                    const double aux = g_weight * g_shape_functions[i];
                    for (std::size_t d = 0; d < TDim; ++d)
                    {
                        // Add Right Hand Side contribution
                        // TODO: Add the LHS terms
                        rRightHandSideVector(i * TDim + d) += aux * (pres_gauss * r_g_unit_normal[d] - proj_dev_stress[d]);
                    }
                }

                // Allocate and calculate auxiliary arrays
                array_1d<double,3> vel_j;
                array_1d<double,3> bc_vel;
                array_1d<double,3> wall_vel;
                BoundedMatrix<double, 3, 3> norm_proj_mat;
                BoundedMatrix<double, StrainSize, TDim * NumNodes> B;
                BoundedMatrix<double, TDim, StrainSize> voigt_normal;
                CalculateBMatrix(g_DN_DX, B);
                VoigtTransformForProduct(r_g_unit_normal, voigt_normal);
                BoundedMatrix<double, StrainSize, TDim *NumNodes> aux_BC = prod(trans(B), trans(r_constitutive_matrix));
                BoundedMatrix<double, TDim * NumNodes, TDim> aux_BC_proj = prod(aux_BC, trans(voigt_normal));

                // Cut-FEM boundary condition Nitsche imposition
                const double penalty_parameter = kappa * (mu / h + rho * (norm_2(vel_gauss) + h / dt));
                for (IndexType i = 0; i < n_nodes; ++i)
                {
                    for (IndexType j = 0; j < n_nodes; ++j)
                    {
                        // j-node data
                        const auto &r_vel_j_0 = r_geom[j].FastGetSolutionStepValue(VELOCITY);
                        const auto &r_vel_j_1 = r_geom[j].FastGetSolutionStepValue(VELOCITY,1);
                        vel_j = theta * r_vel_j_0 + (1.0 - theta) * r_vel_j_1;
                        noalias(wall_vel) = ZeroVector(3); // TODO: This should be the interpolation of the "structure" velocity in the future

                        // Check the boundary condition to be imposed (no-slip or pure slip)
                        if (is_slip) {
                            noalias(norm_proj_mat) = outer_prod(r_g_unit_normal, r_g_unit_normal);
                            noalias(bc_vel) = prod(norm_proj_mat, vel_j - wall_vel);
                        } else {
                            noalias(bc_vel) = vel_j - wall_vel;
                        }

                        // Assemble boundary condition RHS and LHS contributions
                        const double aux_1 = g_weight * penalty_parameter * g_shape_functions[i] * g_shape_functions[j];
                        const double aux_2 = g_weight * g_shape_functions[j];
                        for (IndexType d1 = 0; d1 < TDim; ++d1)
                        {
                            // Penalty term
                            rLeftHandSideMatrix(i * TDim + d1, j * TDim + d1) += aux_1;
                            rRightHandSideVector(i * TDim + d1) -= aux_1 * bc_vel[d1];
                            // Nitsche term (only viscous component)
                            for (IndexType d2 = 0; d2 < TDim; ++d2)
                            {
                                rLeftHandSideMatrix(i * TDim + d1, j * TDim + d2) -= aux_BC_proj(i * TDim + d1,  d2) * aux_2;
                                rRightHandSideVector(i * TDim + d1) += aux_BC_proj(i * TDim + d1, d2) * aux_2 * bc_vel[d1];
                            }
                        }
                    }
                }
            }
        }
    }

    template <unsigned int TDim>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateGeometryData(
        ShapeFunctionDerivativesArrayType &rDN_DX,
        Matrix &NContainer,
        Vector &rGaussWeights)
    {
        if (IsCut())
        {
            // Calculate cut element Gauss point values
            CalculateCutGeometryData(rDN_DX, NContainer, rGaussWeights);
        }
        else
        {
            // If not cut, we use the standard shape functions data calculator from the parent
            BaseType::CalculateGeometryData(rDN_DX, NContainer, rGaussWeights);
        }
    }

    template <unsigned int TDim>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateGeometryData(Vector &rGaussWeights)
    {
        if (IsCut())
        {
            CalculateCutGeometryData(rGaussWeights);
        }
        else
        {
            BaseType::CalculateGeometryData(rGaussWeights);
        }
    }

    template <unsigned int TDim>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateIntersectionGeometryData(
        ShapeFunctionDerivativesArrayType &rInterfaceDNDX,
        Matrix &rInterfaceN,
        Vector &rInterfaceGaussWeights,
        ModifiedShapeFunctions::AreaNormalsContainerType &rInterfaceUnitNormals)
    {
        const auto &r_geom = this->GetGeometry();

        // Auxiliary distance vector for the element subdivision utility
        Vector distances_vector(NumNodes);
        for (std::size_t i = 0; i < NumNodes; ++i)
        {
            distances_vector[i] = r_geom[i].FastGetSolutionStepValue(DISTANCE);
        }

        // Get the subintegration utility
        ModifiedShapeFunctions::Pointer p_mod_sh_func = nullptr;
        if constexpr (TDim == 2)
        {
            p_mod_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), distances_vector);
        }
        else
        {
            p_mod_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(this->pGetGeometry(), distances_vector);
        }

        // Fluid side interface
        p_mod_sh_func->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
            rInterfaceN,
            rInterfaceDNDX,
            rInterfaceGaussWeights,
            GeometryData::IntegrationMethod::GI_GAUSS_1);

        // Fluid side interface normals
        p_mod_sh_func->ComputePositiveSideInterfaceAreaNormals(
            rInterfaceUnitNormals,
            GeometryData::IntegrationMethod::GI_GAUSS_1);

        for (unsigned int i = 0; i < rInterfaceUnitNormals.size(); ++i)
        {
            const double norm = norm_2(rInterfaceUnitNormals[i]);
            KRATOS_WARNING_IF("CalculateIntersectionGeometryData", norm < 1.0e-12) << "Normal is close to zero in element " << this->Id() << " cut interface." << std::endl;
            rInterfaceUnitNormals[i] /= norm;
        }
    }

    template <unsigned int TDim>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateCutGeometryData(
        ShapeFunctionDerivativesArrayType &rDNDX,
        Matrix &rN,
        Vector &rGaussWeights)
    {
        const auto &r_geom = this->GetGeometry();

        // Auxiliary distance vector for the element subdivision utility
        Vector distances_vector(NumNodes);
        for (std::size_t i = 0; i < NumNodes; ++i)
        {
            distances_vector[i] = r_geom[i].FastGetSolutionStepValue(DISTANCE);
        }

        // Get the subintegration utility
        ModifiedShapeFunctions::Pointer p_mod_sh_func = nullptr;
        if constexpr (TDim == 2)
        {
            p_mod_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), distances_vector);
        }
        else
        {
            p_mod_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(this->pGetGeometry(), distances_vector);
        }

        // Fluid side
        p_mod_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
            rN,
            rDNDX,
            rGaussWeights,
            GeometryData::IntegrationMethod::GI_GAUSS_1);

        // TODO: Remove this after developing
        //  Matrix pos_rN;
        //  Vector pos_rGaussWeights;
        //  ShapeFunctionDerivativesArrayType pos_rDNDX;
        //  p_mod_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
        //      pos_rN,
        //      pos_rDNDX,
        //      pos_rGaussWeights,
        //      GeometryData::IntegrationMethod::GI_GAUSS_1);

        // Matrix neg_rN;
        // Vector neg_rGaussWeights;
        // ShapeFunctionDerivativesArrayType neg_rDNDX;
        // p_mod_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        //     neg_rN,
        //     neg_rDNDX,
        //     neg_rGaussWeights,
        //     GeometryData::IntegrationMethod::GI_GAUSS_1);

        // std::size_t n_pos_gauss = pos_rGaussWeights.size();
        // std::size_t n_neg_gauss = neg_rGaussWeights.size();
        // std::size_t n_tot_gauss = n_pos_gauss + n_neg_gauss;
        // KRATOS_WATCH(this->Id())
        // KRATOS_WATCH(n_tot_gauss)
        // rDNDX.resize(n_tot_gauss);
        // rN = ZeroMatrix(n_tot_gauss, NumNodes);
        // rGaussWeights = ZeroVector(n_tot_gauss);
        // std::size_t i_gauss = 0;
        // for (std::size_t gpos = 0; gpos < n_pos_gauss; ++gpos) {
        //     rDNDX[i_gauss] = pos_rDNDX[gpos];
        //     for (std::size_t i = 0; i < NumNodes; ++i) {
        //         rN(i_gauss, i) = pos_rN(gpos, i);
        //     }
        //     rGaussWeights(i_gauss) = pos_rGaussWeights(gpos);
        //     i_gauss++;
        // }
        // for (std::size_t gneg = 0; gneg < n_neg_gauss; ++gneg) {
        //     rDNDX[i_gauss] = neg_rDNDX[gneg];
        //     for (std::size_t i = 0; i < NumNodes; ++i) {
        //         rN(i_gauss, i) = neg_rN(gneg, i);
        //     }
        //     rGaussWeights(i_gauss) = neg_rGaussWeights(gneg);
        //     i_gauss++;
        // }
    }

    template <unsigned int TDim>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::CalculateCutGeometryData(Vector &rGaussWeights)
    {
        const auto &r_geom = this->GetGeometry();

        // Auxiliary distance vector for the element subdivision utility
        Vector distances_vector(NumNodes);
        for (std::size_t i = 0; i < NumNodes; ++i)
        {
            distances_vector[i] = r_geom[i].FastGetSolutionStepValue(DISTANCE);
        }

        // Get the subintegration utility
        ModifiedShapeFunctions::Pointer p_mod_sh_func = nullptr;
        if constexpr (TDim == 2)
        {
            p_mod_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), distances_vector);
        }
        else
        {
            p_mod_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(this->pGetGeometry(), distances_vector);
        }

        // Fluid side
        Matrix aux_N_container;
        p_mod_sh_func->ComputePositiveSideShapeFunctionsAndWeights(
            aux_N_container,
            rGaussWeights,
            GeometryData::IntegrationMethod::GI_GAUSS_1);
    }

    template <unsigned int TDim>
    bool TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::IsCut() const
    {
        const auto &r_geom = this->GetGeometry();
        SizeType n_pos = 0;
        SizeType n_neg = 0;
        for (const auto &r_node : r_geom)
        {
            if (r_node.FastGetSolutionStepValue(DISTANCE) > 0.0)
            {
                n_pos++;
            }
            else
            {
                n_neg++;
            }
        }

        return n_pos != 0 && n_neg != 0 ? true : false;
    }

    template <unsigned int TDim>
    bool TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<TDim>::IsPositive() const
    {
        const auto &r_geom = this->GetGeometry();
        SizeType n_pos = 0;
        for (const auto &r_node : r_geom)
        {
            if (r_node.FastGetSolutionStepValue(DISTANCE) > 0.0)
            {
                n_pos++;
            }
        }

        return n_pos == NumNodes ? true : false;
    }

    template <>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<2>::VoigtStressNormalProjection(
        const Vector &rVoigtStress,
        const array_1d<double, 3> &rUnitNormal,
        array_1d<double, 2> &rProjectedStress)
    {
        rProjectedStress[0] = rVoigtStress[0] * rUnitNormal[0] + rVoigtStress[2] * rUnitNormal[1];
        rProjectedStress[1] = rVoigtStress[2] * rUnitNormal[0] + rVoigtStress[1] * rUnitNormal[1];
    }

    template <>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<3>::VoigtStressNormalProjection(
        const Vector &rVoigtStress,
        const array_1d<double, 3> &rUnitNormal,
        array_1d<double, 3> &rProjectedStress)
    {
        rProjectedStress[0] = rVoigtStress[0] * rUnitNormal[0] + rVoigtStress[3] * rUnitNormal[1] + rVoigtStress[5] * rUnitNormal[2];
        rProjectedStress[1] = rVoigtStress[3] * rUnitNormal[0] + rVoigtStress[1] * rUnitNormal[1] + rVoigtStress[4] * rUnitNormal[2];
        rProjectedStress[2] = rVoigtStress[5] * rUnitNormal[0] + rVoigtStress[4] * rUnitNormal[1] + rVoigtStress[2] * rUnitNormal[2];
    }

    template <>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<2>::CalculateBMatrix(
        const Matrix &rDNDX,
        BoundedMatrix<double, StrainSize, 2 * NumNodes> &rB)
    {
        IndexType index;
        rB = ZeroMatrix(StrainSize, 2 * NumNodes);
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            index = 2 * i;
            rB(0, index + 0) = rDNDX(i, 0);
            rB(1, index + 1) = rDNDX(i, 1);
            rB(2, index + 0) = rDNDX(i, 1);
            rB(2, index + 1) = rDNDX(i, 0);
        }
    }

    template <>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<3>::CalculateBMatrix(
        const Matrix &rDNDX,
        BoundedMatrix<double, StrainSize, 3 * NumNodes> &rB)
    {
        IndexType index;
        rB = ZeroMatrix(StrainSize, 3 * NumNodes);
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            index = 3 * i;
            rB(0, index + 0) = rDNDX(i, 0);
            rB(1, index + 1) = rDNDX(i, 1);
            rB(2, index + 2) = rDNDX(i, 2);
            rB(3, index + 0) = rDNDX(i, 1);
            rB(3, index + 1) = rDNDX(i, 0);
            rB(4, index + 1) = rDNDX(i, 2);
            rB(4, index + 2) = rDNDX(i, 1);
            rB(5, index + 0) = rDNDX(i, 2);
            rB(5, index + 2) = rDNDX(i, 0);
        }
    }

    template <>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<2>::VoigtTransformForProduct(
        const array_1d<double, 3> &rVector,
        BoundedMatrix<double, 2, StrainSize> &rVoigtMatrix)
    {

        rVoigtMatrix.clear();

        rVoigtMatrix(0, 0) = rVector(0);
        rVoigtMatrix(0, 2) = rVector(1);
        rVoigtMatrix(1, 1) = rVector(1);
        rVoigtMatrix(1, 2) = rVector(0);
    }

    template <>
    void TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<3>::VoigtTransformForProduct(
        const array_1d<double, 3> &rVector,
        BoundedMatrix<double, 3, StrainSize> &rVoigtMatrix)
    {

        rVoigtMatrix.clear();

        rVoigtMatrix(0, 0) = rVector(0);
        rVoigtMatrix(0, 3) = rVector(1);
        rVoigtMatrix(0, 5) = rVector(2);
        rVoigtMatrix(1, 1) = rVector(1);
        rVoigtMatrix(1, 3) = rVector(0);
        rVoigtMatrix(1, 4) = rVector(2);
        rVoigtMatrix(2, 2) = rVector(2);
        rVoigtMatrix(2, 4) = rVector(1);
        rVoigtMatrix(2, 5) = rVector(0);
    }

    template class TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<2>;
    template class TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<3>;

} // namespace Kratos

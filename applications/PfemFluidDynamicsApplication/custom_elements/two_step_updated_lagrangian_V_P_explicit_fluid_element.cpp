//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 0.0 $
//
//   Implementation of the Gauss-Seidel two step Updated Lagrangian Velocity-Pressure element
//     ( There is a ScalingConstant to multiply the mass balance equation for a number because i read it somewhere)
//

// System includes

// External includes

// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_explicit_fluid_element.h"
#include "includes/cfd_variables.h"
#include <math.h>

namespace Kratos {

  template< unsigned int TDim >
  Element::Pointer TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {

    TwoStepUpdatedLagrangianVPExplicitFluidElement NewElement(NewId, this->GetGeometry().Create( rThisNodes ), this->pGetProperties() );

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer( new TwoStepUpdatedLagrangianVPExplicitFluidElement(NewElement) );

  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::Initialize()
  {
    KRATOS_TRY;
    KRATOS_CATCH( "" );
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
  {

  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;
    KRATOS_CATCH( "" );
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::ComputeMaterialParameters(double& Density,
										       double& DeviatoricCoeff,
										       double& VolumetricCoeff,
										       ProcessInfo &currentProcessInfo,
										       ElementalVariables &rElementalVariables)
  {
    double FluidBulkModulus=0;
    double FluidYieldShear=0;
    double staticFrictionCoefficient=0;
    double regularizationCoefficient=0;
    // double inertialNumberThreshold=0;
    double timeStep=currentProcessInfo[DELTA_TIME];

    this->EvaluatePropertyFromANotRigidNode(Density,DENSITY);
    this->EvaluatePropertyFromANotRigidNode(FluidBulkModulus,BULK_MODULUS);
    this->EvaluatePropertyFromANotRigidNode(FluidYieldShear,YIELD_SHEAR);
    this->EvaluatePropertyFromANotRigidNode(staticFrictionCoefficient,STATIC_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(regularizationCoefficient,REGULARIZATION_COEFFICIENT);
    // this->EvaluatePropertyFromANotRigidNode(inertialNumberThreshold,INERTIAL_NUMBER_ONE);

    if(FluidBulkModulus==0){
      FluidBulkModulus = 1000000000.0;
    }
    VolumetricCoeff = FluidBulkModulus*timeStep;


    if(FluidYieldShear!=0){
      // std::cout<<"For a Newtonian fluid I should not enter here"<<std::endl;
      DeviatoricCoeff=this->ComputeNonLinearViscosity(rElementalVariables.EquivalentStrainRate);
    }else if(staticFrictionCoefficient!=0){
      DeviatoricCoeff=this->ComputePapanastasiouMuIrheologyViscosity(rElementalVariables);
      // if(regularizationCoefficient!=0 && inertialNumberThreshold==0){
      // 	// DeviatoricCoeff=this->ComputeBercovierMuIrheologyViscosity(rElementalVariables);
      // 	DeviatoricCoeff=this->ComputePapanastasiouMuIrheologyViscosity(rElementalVariables);
      // }
      // else if(regularizationCoefficient==0 && inertialNumberThreshold!=0){
      // 	DeviatoricCoeff=this->ComputeBarkerMuIrheologyViscosity(rElementalVariables);
      // }else if(regularizationCoefficient!=0 && inertialNumberThreshold!=0){
      // 	DeviatoricCoeff=this->ComputeBarkerBercovierMuIrheologyViscosity(rElementalVariables);
      // }else{
      // 	DeviatoricCoeff=this->ComputeJopMuIrheologyViscosity(rElementalVariables);
      // }
    }else{
      // std::cout<<"For a Newtonian fluid I should  enter here"<<std::endl;
      this->EvaluatePropertyFromANotRigidNode(DeviatoricCoeff,DYNAMIC_VISCOSITY);
    }

    // this->ComputeMaterialParametersGranularGas(rElementalVariables,VolumetricCoeff,DeviatoricCoeff);
    // std::cout<<"Density "<<Density<<std::endl;
    // std::cout<<"FluidBulkModulus "<<FluidBulkModulus<<std::endl;
    // std::cout<<"staticFrictionCoefficient "<<staticFrictionCoefficient<<std::endl;
    // std::cout<<"DeviatoricCoeff "<<DeviatoricCoeff<<std::endl;


    this->mMaterialDeviatoricCoefficient=DeviatoricCoeff;
    this->mMaterialVolumetricCoefficient=VolumetricCoeff;
    this->mMaterialDensity=Density;

    // const SizeType NumNodes = this->GetGeometry().PointsNumber();
    // for (SizeType i = 0; i < NumNodes; ++i)
    //   {
    // 	this->GetGeometry()[i].FastGetSolutionStepValue(ADAPTIVE_EXPONENT)=VolumetricCoeff;
    // 	this->GetGeometry()[i].FastGetSolutionStepValue(ALPHA_PARAMETER)=DeviatoricCoeff;
    // 	this->GetGeometry()[i].FastGetSolutionStepValue(FLOW_INDEX)=rElementalVariables.EquivalentStrainRate;
    //   }

  }



  template< unsigned int TDim>
  double TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::ComputeNonLinearViscosity(double & equivalentStrainRate)
  {
    double FluidViscosity=0;

    double FluidFlowIndex=0;
    double FluidYieldShear=0;
    double FluidAdaptiveExponent=0;
    this->EvaluatePropertyFromANotRigidNode(FluidViscosity,DYNAMIC_VISCOSITY);
    this->EvaluatePropertyFromANotRigidNode(FluidFlowIndex,FLOW_INDEX);
    this->EvaluatePropertyFromANotRigidNode(FluidYieldShear,YIELD_SHEAR);
    this->EvaluatePropertyFromANotRigidNode(FluidAdaptiveExponent,ADAPTIVE_EXPONENT);
    double exponent=-FluidAdaptiveExponent*equivalentStrainRate;
    if(equivalentStrainRate!=0){
      FluidViscosity+=(FluidYieldShear/equivalentStrainRate)*(1-exp(exponent));
    }
    if(equivalentStrainRate<0.00001 && FluidYieldShear!=0 && FluidAdaptiveExponent!=0){
      // for gamma_dot very small the limit of the Papanastasiou viscosity is mu=m*tau_yield
      FluidViscosity=FluidAdaptiveExponent*FluidYieldShear;
    }
    return FluidViscosity;
  }


  template< unsigned int TDim>
  void TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::ComputeMaterialParametersGranularGas(double& Density,
												  double& DeviatoricCoeff,
												  double& VolumetricCoeff,
												  ProcessInfo &currentProcessInfo,
												  ElementalVariables &rElementalVariables)
  {

    this->EvaluatePropertyFromANotRigidNode(Density,DENSITY);

    double staticFrictionCoefficient=0;
    double dynamicFrictionCoefficient=0;
    double inertialNumberZero=0;
    double grainDiameter=0;
    double grainDensity=0;

    this->EvaluatePropertyFromANotRigidNode(staticFrictionCoefficient,STATIC_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(dynamicFrictionCoefficient,DYNAMIC_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(inertialNumberZero,INERTIAL_NUMBER_ZERO);
    this->EvaluatePropertyFromANotRigidNode(grainDiameter,GRAIN_DIAMETER);
    this->EvaluatePropertyFromANotRigidNode(grainDensity,GRAIN_DENSITY);
    double temperature=5.0;
    double gZero=0;
    double voidRatioS=0.35;
    double epsilonR=0.6;
    double nuDot=rElementalVariables.VolumetricDefRate;//????
    double voidRatio=1.45;
    double concentration=1.0/(1.0+voidRatio);
    double concentrationS=1.0/(1.0+voidRatioS);
    if(concentration>0.49){
      gZero=5.69*(concentrationS-0.49)/(concentrationS-concentration);
    }else{
      gZero=(2-concentration)/(2*(1-concentration)*(1-concentration)*(1-concentration));
    }
    // if(voidRatio<1.04){
    //   gZero=(2.9019-2.7881*voidRatioS)*(1.0+voidRatio)/((voidRatio-voidRatioS));
    // }else{
    //   double lC=1.0;// ???
    //   gZero=lC*(2.0*voidRatio+1.0)*(1.0+voidRatio)*(1.0+voidRatio)/(2*pow(voidRatio,3));
    // }
    const double PI  =3.141592653589793238463;
    double G=gZero*concentration;
    double F=0.5*(1.0+epsilonR)+0.25*G;
    double J=0.5*(1.0+epsilonR) + 0.09817477*(5.0+2.0*(1.0+epsilonR)*(-1.0+3.0*epsilonR)*G)*(5.0+4.0*(1.0+epsilonR)*G)/
      ((24.0-6.0*(1.0+epsilonR)*(1.0+epsilonR)-5.0*(1.0-pow(epsilonR,2)))*pow(G,2));
    double cStar=32*(1.0-epsilonR)*(1.0-pow(epsilonR,2))/(81.0-17.0*epsilonR+30.0*pow(epsilonR,2)*(1.0-epsilonR));
    double Gamma=8.148733086*G*(1.0+epsilonR)*(1-0.03125*cStar);
    double nuK=(1.0-0.4*(1.0+epsilonR)*(1.0-3.0*epsilonR)*gZero/(1.0+voidRatio))/
      (gZero*(1-0.25*(1.0-epsilonR)*(1.0-epsilonR))*(1-0.015625*cStar)-0.20833333333*gZero*(1-pow(epsilonR,2))*(1+0.09375*cStar));
    double f1Coeff=4*G*F/(1.0+voidRatio);
    double f2Coeff=0.90270333333*G*J*concentration;
    double psiStar0=0.4166666667*gZero*(1+pow(epsilonR,2))*(1+0.09375*cStar);
    double pStar=1+2*(1+epsilonR)*G;
    double nuGammaStar=(1+epsilonR)*0.020833333333*gZero*(128-96*epsilonR+15*pow(epsilonR,2)-15*pow(epsilonR,3)+cStar*0.015625*(15*pow(epsilonR,3)-15*pow(epsilonR,2)+498*epsilonR-434));
    double lambda=0.375*((1.0-epsilonR)*(5*pow(epsilonR,2)+4*epsilonR-1.0)+0.083333333333*cStar*(-15*pow(epsilonR,3)+3*pow(epsilonR,2)-19*epsilonR-159));
    double cD=0.01188634*(1.0+voidRatio)*(0.26666666667*lambda*gZero*concentration+(pStar-1.0)*(0.6666666667-epsilonR)*cStar)/
      (0.5*psiStar0+nuGammaStar+0.078125*cStar*(1+0.046875*cStar)*gZero*(1-pow(epsilonR,2)));
    double psiStar1=(0.0520833333*sqrt(PI/temperature)*grainDiameter*(1.0+voidRatio)*(pStar-1.0)-0.15625*(1-pow(epsilonR,2))*(1+0.046875*cStar)*gZero*cD)*nuDot;
    double psiStar=psiStar0+psiStar1;
    double f3Coeff=20.2646753*psiStar*pow(concentration,2);
    double f4Coeff=0.092315304*(0.6*Gamma-0.6666666667*nuK*(1+0.8*gZero*(1.0+epsilonR)*concentration));
    double correlationLength=grainDiameter;
    double cMaterialParameter=0.5;
    double value=pow(cMaterialParameter,2)*pow(G,0.666666666667)*f3Coeff/(4.0*f2Coeff);
    if(value>1.0){
      correlationLength*=value;
    }
    double f5Coeff=correlationLength*f2Coeff/(grainDiameter*f3Coeff);
    temperature=pow(grainDiameter,2)*f5Coeff*rElementalVariables.EquivalentStrainRate;
    // temperature=pow(grainDiameter,2)*f5Coeff*ElementalVariables.SpatialDefRate[2];


    VolumetricCoeff=grainDensity*grainDiameter*f4Coeff*sqrt(fabs(temperature));
    DeviatoricCoeff=grainDensity*grainDiameter*f2Coeff*sqrt(fabs(temperature));

    // std::cout<<"  temperature "<<temperature<<"  grainDensity "<<grainDensity<<"  f1Coeff "<<f1Coeff<<"  f2Coeff "<<f2Coeff<<"  f3Coeff "<<f3Coeff<<"  f4Coeff "<<f4Coeff<<"  f5Coeff "<<f5Coeff<<std::endl;
    // std::cout<<"Density "<<Density<<"   VolumetricCoeff "<<VolumetricCoeff<<"   DeviatoricCoeff "<<DeviatoricCoeff<<std::endl;
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < NumNodes; ++i)
      {
    	this->GetGeometry()[i].FastGetSolutionStepValue(ADAPTIVE_EXPONENT)=VolumetricCoeff;
    	this->GetGeometry()[i].FastGetSolutionStepValue(ALPHA_PARAMETER)=DeviatoricCoeff;
	this->GetGeometry()[i].FastGetSolutionStepValue(FLOW_INDEX)=rElementalVariables.EquivalentStrainRate;
      }
    rElementalVariables.UpdatedTotalCauchyStress[0]=grainDensity*grainDiameter*f1Coeff*f5Coeff*rElementalVariables.SpatialDefRate[2]*rElementalVariables.SpatialDefRate[2];
    rElementalVariables.UpdatedTotalCauchyStress[1]=grainDensity*grainDiameter*f1Coeff*f5Coeff*rElementalVariables.SpatialDefRate[2]*rElementalVariables.SpatialDefRate[2];
    rElementalVariables.UpdatedTotalCauchyStress[2]=DeviatoricCoeff*rElementalVariables.SpatialDefRate[2];
  }


  template< unsigned int TDim>
  double TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::ComputeJopMuIrheologyViscosity(ElementalVariables & rElementalVariables)
  {
    double FluidViscosity=0;

    double staticFrictionCoefficient=0;
    double dynamicFrictionCoefficient=0;
    double inertialNumberZero=0;
    double grainDiameter=0;
    double grainDensity=0;

    this->EvaluatePropertyFromANotRigidNode(staticFrictionCoefficient,STATIC_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(dynamicFrictionCoefficient,DYNAMIC_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(inertialNumberZero,INERTIAL_NUMBER_ZERO);
    this->EvaluatePropertyFromANotRigidNode(grainDiameter,GRAIN_DIAMETER);
    this->EvaluatePropertyFromANotRigidNode(grainDensity,GRAIN_DENSITY);

    double meanPressure=rElementalVariables.MeanPressure;
    if(meanPressure>0){
      meanPressure=0.0000001;
    }

    double deltaFrictionCoefficient=dynamicFrictionCoefficient-staticFrictionCoefficient;
    double inertialNumber=0;
    if(meanPressure!=0){
      inertialNumber=rElementalVariables.EquivalentStrainRate*grainDiameter/sqrt(fabs(meanPressure)/grainDensity);
    }
    // double inertialNumber=rElementalVariables.EquivalentStrainRate*grainDiameter/sqrt(fabs(rElementalVariables.MeanPressure)/grainDensity);

    if(rElementalVariables.EquivalentStrainRate!=0 && fabs(meanPressure)!=0){
      double firstViscousTerm=staticFrictionCoefficient/rElementalVariables.EquivalentStrainRate;
      double secondViscousTerm=deltaFrictionCoefficient*inertialNumber/ ((inertialNumberZero+inertialNumber)*rElementalVariables.EquivalentStrainRate);
      FluidViscosity=(firstViscousTerm+secondViscousTerm)*fabs(meanPressure);
    }else{
      FluidViscosity=1.0;
    }
    return FluidViscosity;
  }


  template< unsigned int TDim>
  double TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::ComputeBercovierMuIrheologyViscosity(ElementalVariables & rElementalVariables)
  {
    double FluidViscosity=0;
    double staticFrictionCoefficient=0;
    double dynamicFrictionCoefficient=0;
    double inertialNumberZero=0;
    double grainDiameter=0;
    double grainDensity=0;
    double regularizationCoefficient=0;

    this->EvaluatePropertyFromANotRigidNode(staticFrictionCoefficient,STATIC_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(dynamicFrictionCoefficient,DYNAMIC_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(inertialNumberZero,INERTIAL_NUMBER_ZERO);
    this->EvaluatePropertyFromANotRigidNode(grainDiameter,GRAIN_DIAMETER);
    this->EvaluatePropertyFromANotRigidNode(grainDensity,GRAIN_DENSITY);
    this->EvaluatePropertyFromANotRigidNode(regularizationCoefficient,REGULARIZATION_COEFFICIENT);

    double meanPressure=rElementalVariables.MeanPressure;
    if(meanPressure>0){
      meanPressure=0.0000001;
    }

    double deltaFrictionCoefficient=dynamicFrictionCoefficient-staticFrictionCoefficient;
    double inertialNumber=0;
    if(meanPressure!=0){
      inertialNumber=rElementalVariables.EquivalentStrainRate*grainDiameter/sqrt(fabs(meanPressure)/grainDensity);
    }

    if(rElementalVariables.EquivalentStrainRate!=0 && fabs(meanPressure)!=0){
      double firstViscousTerm=staticFrictionCoefficient / sqrt(pow(rElementalVariables.EquivalentStrainRate,2)+pow(regularizationCoefficient,2));
      double secondViscousTerm=deltaFrictionCoefficient*inertialNumber/ ((inertialNumberZero+inertialNumber)*rElementalVariables.EquivalentStrainRate);
      FluidViscosity=(firstViscousTerm+secondViscousTerm)*fabs(meanPressure);
    }else{
      FluidViscosity=1.0;
    }

    return FluidViscosity;
  }


  template< unsigned int TDim>
  double TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::ComputePapanastasiouMuIrheologyViscosity(ElementalVariables & rElementalVariables)
  {
    double FluidViscosity=0;
    double staticFrictionCoefficient=0;
    double dynamicFrictionCoefficient=0;
    double inertialNumberZero=0;
    double grainDiameter=0;
    double grainDensity=0;
    double regularizationCoefficient=0;

    this->EvaluatePropertyFromANotRigidNode(staticFrictionCoefficient,STATIC_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(dynamicFrictionCoefficient,DYNAMIC_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(inertialNumberZero,INERTIAL_NUMBER_ZERO);
    this->EvaluatePropertyFromANotRigidNode(grainDiameter,GRAIN_DIAMETER);
    this->EvaluatePropertyFromANotRigidNode(grainDensity,GRAIN_DENSITY);
    this->EvaluatePropertyFromANotRigidNode(regularizationCoefficient,REGULARIZATION_COEFFICIENT);

    double pressure=rElementalVariables.MeanPressure;
    if(pressure>0){
      pressure=0.0000001;
    }

    double deltaFrictionCoefficient=dynamicFrictionCoefficient-staticFrictionCoefficient;
    double inertialNumber=0;
    if(rElementalVariables.MeanPressure!=0){
      inertialNumber=rElementalVariables.EquivalentStrainRate*grainDiameter/sqrt(fabs(pressure)/grainDensity);
    }

    double exponent=-rElementalVariables.EquivalentStrainRate/regularizationCoefficient;

    if(rElementalVariables.EquivalentStrainRate!=0 && fabs(pressure)!=0){
      double firstViscousTerm=staticFrictionCoefficient*(1-exp(exponent))/rElementalVariables.EquivalentStrainRate;
      double secondViscousTerm=deltaFrictionCoefficient*inertialNumber/ ((inertialNumberZero+inertialNumber)*rElementalVariables.EquivalentStrainRate);
      FluidViscosity=(firstViscousTerm+secondViscousTerm)*fabs(pressure);
    }else{
      FluidViscosity=1.0;
    }

    // const SizeType NumNodes = this->GetGeometry().PointsNumber();
    // for (SizeType i = 0; i < NumNodes; ++i)
    //   {
    // 	this->GetGeometry()[i].FastGetSolutionStepValue(FLOW_INDEX)=FluidViscosity;
    // 	this->GetGeometry()[i].FastGetSolutionStepValue(ADAPTIVE_EXPONENT)=rElementalVariables.EquivalentStrainRate;
    // 	this->GetGeometry()[i].FastGetSolutionStepValue(ALPHA_PARAMETER)=inertialNumber;
    // 	// std::cout<<"FluidViscosity "<<FluidViscosity<<"  StrainRate "<<rElementalVariables.EquivalentStrainRate<<"  inertialNumber "<<inertialNumber<<"  pressure "<<rElementalVariables.MeanPressure<<std::endl;
    //   }

    return FluidViscosity;
  }

  template< unsigned int TDim>
  double TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::ComputeBarkerMuIrheologyViscosity(ElementalVariables & rElementalVariables)
  {
    double FluidViscosity=0;
    double staticFrictionCoefficient=0;
    double dynamicFrictionCoefficient=0;
    double inertialNumberZero=0;
    double grainDiameter=0;
    double grainDensity=0;
    double inertialNumberThreshold=0;
    double infiniteFrictionCoefficient=0;
    double alphaParameter=0;

    this->EvaluatePropertyFromANotRigidNode(staticFrictionCoefficient,STATIC_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(dynamicFrictionCoefficient,DYNAMIC_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(inertialNumberZero,INERTIAL_NUMBER_ZERO);
    this->EvaluatePropertyFromANotRigidNode(grainDiameter,GRAIN_DIAMETER);
    this->EvaluatePropertyFromANotRigidNode(grainDensity,GRAIN_DENSITY);
    this->EvaluatePropertyFromANotRigidNode(inertialNumberThreshold,INERTIAL_NUMBER_ONE);
    this->EvaluatePropertyFromANotRigidNode(infiniteFrictionCoefficient,INFINITE_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(alphaParameter,ALPHA_PARAMETER);

    double meanPressure=rElementalVariables.MeanPressure;
    if(meanPressure>0){
      meanPressure=0.0000001;
    }

    double inertialNumber=0;
    if(meanPressure!=0){
      inertialNumber=rElementalVariables.EquivalentStrainRate*grainDiameter/sqrt(fabs(meanPressure)/grainDensity);
    }

    if(inertialNumber>inertialNumberThreshold){
      FluidViscosity=(staticFrictionCoefficient*inertialNumberZero+dynamicFrictionCoefficient*inertialNumber+infiniteFrictionCoefficient*pow(inertialNumber,2))/(inertialNumberZero+inertialNumber);
    }else{
      double denominator=staticFrictionCoefficient*inertialNumberZero+dynamicFrictionCoefficient*inertialNumberThreshold+infiniteFrictionCoefficient*pow(inertialNumberThreshold,2);
      double exponent=alphaParameter*(inertialNumberZero+inertialNumberThreshold)*(inertialNumberZero+inertialNumberThreshold)/pow(denominator,2);
      double firstAconstant=inertialNumberThreshold*exp(exponent);
      FluidViscosity=sqrt(alphaParameter/log(firstAconstant/inertialNumber));
    }

    if(rElementalVariables.EquivalentStrainRate!=0 && fabs(meanPressure)!=0){
      FluidViscosity*=fabs(meanPressure)/rElementalVariables.EquivalentStrainRate;
    }else{
      FluidViscosity=1.0;
    }
    return FluidViscosity;
  }


  template< unsigned int TDim>
  double TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::ComputeBarkerBercovierMuIrheologyViscosity(ElementalVariables & rElementalVariables)
  {

    double FluidViscosity=0;
    double staticFrictionCoefficient=0;
    double dynamicFrictionCoefficient=0;
    double inertialNumberZero=0;
    double grainDiameter=0;
    double grainDensity=0;
    double inertialNumberThreshold=0;
    double infiniteFrictionCoefficient=0;
    double alphaParameter=0;
    double regularizationCoefficient=0;

    this->EvaluatePropertyFromANotRigidNode(staticFrictionCoefficient,STATIC_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(dynamicFrictionCoefficient,DYNAMIC_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(inertialNumberZero,INERTIAL_NUMBER_ZERO);
    this->EvaluatePropertyFromANotRigidNode(grainDiameter,GRAIN_DIAMETER);
    this->EvaluatePropertyFromANotRigidNode(grainDensity,GRAIN_DENSITY);
    this->EvaluatePropertyFromANotRigidNode(inertialNumberThreshold,INERTIAL_NUMBER_ONE);
    this->EvaluatePropertyFromANotRigidNode(infiniteFrictionCoefficient,INFINITE_FRICTION);
    this->EvaluatePropertyFromANotRigidNode(alphaParameter,ALPHA_PARAMETER);
    this->EvaluatePropertyFromANotRigidNode(regularizationCoefficient,REGULARIZATION_COEFFICIENT);

    double meanPressure=rElementalVariables.MeanPressure;
    if(meanPressure>0){
      meanPressure=0.0000001;
    }

    double inertialNumber=0;
    if(meanPressure!=0){
      inertialNumber=rElementalVariables.EquivalentStrainRate*grainDiameter/sqrt(fabs(meanPressure)/grainDensity);
    }

    if(inertialNumber>inertialNumberThreshold){
      double deltaFrictionCoefficient=dynamicFrictionCoefficient-staticFrictionCoefficient;
      double firstViscousTerm=staticFrictionCoefficient;
      double secondViscousTerm=deltaFrictionCoefficient*inertialNumber / (inertialNumberZero+inertialNumber);
      FluidViscosity=(firstViscousTerm+secondViscousTerm);
    }else{
      double denominator=staticFrictionCoefficient*inertialNumberZero+dynamicFrictionCoefficient*inertialNumberThreshold+infiniteFrictionCoefficient*pow(inertialNumberThreshold,2);
      double exponent=alphaParameter*(inertialNumberZero+inertialNumberThreshold)*(inertialNumberZero+inertialNumberThreshold)/pow(denominator,2);
      double firstAconstant=inertialNumberThreshold*exp(exponent);
      FluidViscosity=sqrt(alphaParameter/log(firstAconstant/inertialNumber));
    }

    if(rElementalVariables.EquivalentStrainRate!=0 && fabs(meanPressure)!=0){
      double exponent=-rElementalVariables.EquivalentStrainRate/regularizationCoefficient;
      FluidViscosity*=fabs(meanPressure)*(1-exp(exponent))/rElementalVariables.EquivalentStrainRate;
    }else{
      if(meanPressure==0 && rElementalVariables.EquivalentStrainRate!=0){
    	FluidViscosity*=1.0 / sqrt(pow(rElementalVariables.EquivalentStrainRate,2)+pow(regularizationCoefficient,2));
      }else if(meanPressure!=0 && rElementalVariables.EquivalentStrainRate==0){
    	FluidViscosity*=fabs(meanPressure) / sqrt(0.001+pow(regularizationCoefficient,2));
      }else{
    	FluidViscosity=1.0;
      }
    }

    return FluidViscosity;
  }



  template< unsigned int TDim >
  int TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    if(VELOCITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check that the application was correctly registered.","");
    if(ACCELERATION.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION Key is 0. Check that the application was correctly registered.","");
    if(PRESSURE.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check that the application was correctly registered.","");
    if(BODY_FORCE.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"BODY_FORCE Key is 0. Check that the application was correctly registered.","");
    if(DENSITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check that the application was correctly registered.","");
    if(DYNAMIC_VISCOSITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"DYNAMIC_VISCOSITY Key is 0. Check that the application was correctly registered.","");
    if(DELTA_TIME.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check that the application was correctly registered.","");

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
      {
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing BODY_FORCE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing DENSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(DYNAMIC_VISCOSITY) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing DYNAMIC_VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
           this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
           this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
      }

    // If this is a 2D problem, check that nodes are in XY plane
    if (this->GetGeometry().WorkingSpaceDimension() == 2)
      {
        for (unsigned int i=0; i<this->GetGeometry().size(); ++i)
	  {
            if (this->GetGeometry()[i].Z() != 0.0)
	      KRATOS_THROW_ERROR(std::invalid_argument,"Node with non-zero Z coordinate found. Id: ",this->GetGeometry()[i].Id());
	  }
      }

    return ierr;

    KRATOS_CATCH("");
  }





  template<>
  void TwoStepUpdatedLagrangianVPExplicitFluidElement<2>::GetPositions(Vector& rValues,const ProcessInfo& rCurrentProcessInfo,const double theta)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
	rValues[Index++] = rGeom[i].X();
	rValues[Index++] = rGeom[i].Y();
      }
  }



  template<>
  void TwoStepUpdatedLagrangianVPExplicitFluidElement<3>::GetPositions(Vector& rValues,const ProcessInfo& rCurrentProcessInfo,const double theta)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
 	rValues[Index++] = rGeom[i].X();
        rValues[Index++] = rGeom[i].Y();
        rValues[Index++] = rGeom[i].Z();
      }
  }



  template <  unsigned int TDim>
  void TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>:: InitializeElementalVariables(ElementalVariables & rElementalVariables)
  {
    KRATOS_TRY;

    unsigned int voigtsize  = 3;
    if( TDim == 3 )
      {
	voigtsize  = 6;
      }
    rElementalVariables.voigtsize=voigtsize;


    rElementalVariables.DetFgrad=1.0;

    rElementalVariables.DetFgradVel=1.0;

    rElementalVariables.DeviatoricInvariant=1.0;

    rElementalVariables.EquivalentStrainRate=1.0;

    rElementalVariables.VolumetricDefRate=1.0;

    rElementalVariables.SpatialDefRate= ZeroVector(voigtsize);

    rElementalVariables.MDGreenLagrangeMaterial.resize(voigtsize,false);

    noalias(rElementalVariables.MDGreenLagrangeMaterial) = ZeroVector(voigtsize);

    rElementalVariables.Fgrad = ZeroMatrix(TDim,TDim);

    rElementalVariables.InvFgrad= ZeroMatrix(TDim,TDim);

    rElementalVariables.FgradVel= ZeroMatrix(TDim,TDim);

    rElementalVariables.InvFgradVel= ZeroMatrix(TDim,TDim);

    rElementalVariables.SpatialVelocityGrad= ZeroMatrix(TDim,TDim);

    rElementalVariables.MeanPressure=0;

    rElementalVariables.CurrentTotalCauchyStress= ZeroVector(voigtsize);

    rElementalVariables.UpdatedTotalCauchyStress=  ZeroVector(voigtsize);

    rElementalVariables.CurrentDeviatoricCauchyStress=  ZeroVector(voigtsize);

    rElementalVariables.UpdatedDeviatoricCauchyStress= ZeroVector(voigtsize);

    KRATOS_CATCH("");

  }


  template < >
  void TwoStepUpdatedLagrangianVPExplicitFluidElement<2>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,double TimeStep, unsigned int g)
  {

    double CurrSecondLame  = this->mMaterialDeviatoricCoefficient;
    // double CurrBulkModulus = this->mMaterialVolumetricCoefficient;
    // double CurrFirstLame  = CurrBulkModulus - 2.0*CurrSecondLame/3.0;

    double DefX=rElementalVariables.SpatialDefRate[0];
    double DefY=rElementalVariables.SpatialDefRate[1];
    double DefXY=rElementalVariables.SpatialDefRate[2];

    double DefVol=rElementalVariables.VolumetricDefRate;

    double sigmaDev_xx= 2*CurrSecondLame*(DefX - DefVol/3.0);
    double sigmaDev_yy= 2*CurrSecondLame*(DefY - DefVol/3.0);
    double sigmaDev_xy= 2*CurrSecondLame*DefXY;

    // double sigmaTot_xx= CurrFirstLame*DefVol + 2.0*CurrSecondLame*DefX;
    // double sigmaTot_yy= CurrFirstLame*DefVol + 2.0*CurrSecondLame*DefY;
    // double sigmaTot_xy= 2.0*CurrSecondLame*DefXY;

    // sigmaDev_xx=rElementalVariables.CurrentDeviatoricCauchyStress[0];
    // sigmaDev_yy=rElementalVariables.CurrentDeviatoricCauchyStress[1];
    // sigmaDev_xy=rElementalVariables.CurrentDeviatoricCauchyStress[2];

    // sigmaTot_xx+=rElementalVariables.CurrentTotalCauchyStress[0];
    // sigmaTot_yy+=rElementalVariables.CurrentTotalCauchyStress[1];
    // sigmaTot_xy+=rElementalVariables.CurrentTotalCauchyStress[2];

    double sigmaTot_xx= sigmaDev_xx + rElementalVariables.MeanPressure;
    double sigmaTot_yy= sigmaDev_yy + rElementalVariables.MeanPressure;
    double sigmaTot_xy= sigmaDev_xy;

    // sigmaDev_xx= sigmaTot_xx - rElementalVariables.MeanPressure;
    // sigmaDev_yy= sigmaTot_yy - rElementalVariables.MeanPressure;
    // sigmaDev_xy= sigmaTot_xy;

    rElementalVariables.UpdatedDeviatoricCauchyStress[0]=sigmaDev_xx;
    rElementalVariables.UpdatedDeviatoricCauchyStress[1]=sigmaDev_yy;
    rElementalVariables.UpdatedDeviatoricCauchyStress[2]=sigmaDev_xy;

    rElementalVariables.UpdatedTotalCauchyStress[0]=sigmaTot_xx;
    rElementalVariables.UpdatedTotalCauchyStress[1]=sigmaTot_yy;
    rElementalVariables.UpdatedTotalCauchyStress[2]=sigmaTot_xy;

    // double TauNorm=sqrt((0.5*sigmaDev_xx*sigmaDev_xx + 0.5*sigmaDev_yy*sigmaDev_yy + sigmaDev_xy*sigmaDev_xy));
    // double FluidYieldShear=0;
    // this->EvaluatePropertyFromANotRigidNode(FluidYieldShear,YIELD_SHEAR);

    // const SizeType NumNodes = this->GetGeometry().PointsNumber();
    // for (SizeType i = 0; i < NumNodes; ++i)
    //   {
    // 	this->GetGeometry()[i].FastGetSolutionStepValue(ADAPTIVE_EXPONENT)=VolumetricCoeff;
    // 	this->GetGeometry()[i].FastGetSolutionStepValue(ALPHA_PARAMETER)=DeviatoricCoeff;
    // 	this->GetGeometry()[i].FastGetSolutionStepValue(FLOW_INDEX)=rElementalVariables.EquivalentStrainRate;
    //   }
    // GeometryType& rGeom = this->GetGeometry();
    // const SizeType NumNodes = rGeom.PointsNumber();
    // if(TauNorm>FluidYieldShear){
    //   this->SetValue(YIELDED,1.0);

    //   for (SizeType i = 0; i < NumNodes; ++i){

    // 	// rGeom[i].FastGetSolutionStepValue(FLOW_INDEX) = 1;
    // 	rGeom[i].FastGetSolutionStepValue(FREESURFACE) = 1;
    //   }
    // }else{

    //   this->SetValue(YIELDED,0.0);
    //   // std::vector<double> rOutput;
    //   // this->GetElementalValueForOutput(YIELDED,rOutput);

    //   for (SizeType i = 0; i < NumNodes; ++i){
    // 	// rGeom[i].FastGetSolutionStepValue(FLOW_INDEX) = 0;
    // 	rGeom[i].FastGetSolutionStepValue(FREESURFACE) = 0;
    //   }
    // }

  }

  //    template< unsigned int TDim >
  //   void TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
  // 									      std::vector<ConstitutiveLaw::Pointer>& rValues,
  // 									      const ProcessInfo& rCurrentProcessInfo )
  //   {

  //     if(rVariable == YIELDED)
  //     {

  //       rValues[0] = 1.5
  //         // if ( rValues.size() != mConstitutiveLawVector.size() )
  //         // {
  //         //     rValues.resize(mConstitutiveLawVector.size());
  //         // }

  //         // for(unsigned int i=0; i<rValues.size(); i++)
  //         // {
  //         //     rValues[i] = 1.5;
  //         // }
  //     }

  // }






  template < >
  void TwoStepUpdatedLagrangianVPExplicitFluidElement<3>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables, double TimeStep, unsigned int g)
  {

    double CurrSecondLame  = this->mMaterialDeviatoricCoefficient;
    // double CurrBulkModulus = this->mMaterialVolumetricCoefficient;
    // double CurrFirstLame  = CurrBulkModulus - 2.0*CurrSecondLame/3.0;

    double DefX=rElementalVariables.SpatialDefRate[0];
    double DefY=rElementalVariables.SpatialDefRate[1];
    double DefZ=rElementalVariables.SpatialDefRate[2];
    double DefXY=rElementalVariables.SpatialDefRate[3];
    double DefXZ=rElementalVariables.SpatialDefRate[4];
    double DefYZ=rElementalVariables.SpatialDefRate[5];

    double DefVol=rElementalVariables.VolumetricDefRate;

    double sigmaDev_xx= 2*CurrSecondLame*(DefX - DefVol/3.0);
    double sigmaDev_yy= 2*CurrSecondLame*(DefY - DefVol/3.0);
    double sigmaDev_zz= 2*CurrSecondLame*(DefZ - DefVol/3.0);
    double sigmaDev_xy= 2*CurrSecondLame*DefXY;
    double sigmaDev_xz= 2*CurrSecondLame*DefXZ;
    double sigmaDev_yz= 2*CurrSecondLame*DefYZ;

    // double sigmaTot_xx= CurrFirstLame*DefVol + 2*CurrSecondLame*DefX;
    // double sigmaTot_yy= CurrFirstLame*DefVol + 2*CurrSecondLame*DefY;
    // double sigmaTot_zz= CurrFirstLame*DefVol + 2*CurrSecondLame*DefZ;
    // double sigmaTot_xy= 2*CurrSecondLame*DefXY;
    // double sigmaTot_xz= 2*CurrSecondLame*DefXZ;
    // double sigmaTot_yz= 2*CurrSecondLame*DefYZ;

    // sigmaDev_xx+=rElementalVariables.CurrentDeviatoricCauchyStress[0];
    // sigmaDev_yy+=rElementalVariables.CurrentDeviatoricCauchyStress[1];
    // sigmaDev_zz+=rElementalVariables.CurrentDeviatoricCauchyStress[2];
    // sigmaDev_xy+=rElementalVariables.CurrentDeviatoricCauchyStress[3];
    // sigmaDev_xz+=rElementalVariables.CurrentDeviatoricCauchyStress[4];
    // sigmaDev_yz+=rElementalVariables.CurrentDeviatoricCauchyStress[5];

    double sigmaTot_xx= sigmaDev_xx + rElementalVariables.MeanPressure;
    double sigmaTot_yy= sigmaDev_yy + rElementalVariables.MeanPressure;
    double sigmaTot_zz= sigmaDev_zz + rElementalVariables.MeanPressure;
    double sigmaTot_xy= sigmaDev_xy;
    double sigmaTot_xz= sigmaDev_xz;
    double sigmaTot_yz= sigmaDev_yz;

    // sigmaTot_xx+=rElementalVariables.CurrentTotalCauchyStress[0];
    // sigmaTot_yy+=rElementalVariables.CurrentTotalCauchyStress[1];
    // sigmaTot_zz+=rElementalVariables.CurrentTotalCauchyStress[2];
    // sigmaTot_xy+=rElementalVariables.CurrentTotalCauchyStress[3];
    // sigmaTot_xz+=rElementalVariables.CurrentTotalCauchyStress[4];
    // sigmaTot_yz+=rElementalVariables.CurrentTotalCauchyStress[5];

    // sigmaDev_xx= sigmaTot_xx - rElementalVariables.MeanPressure;
    // sigmaDev_yy= sigmaTot_yy - rElementalVariables.MeanPressure;
    // sigmaDev_zz= sigmaTot_zz - rElementalVariables.MeanPressure;
    // sigmaDev_xy= sigmaTot_xy;
    // sigmaDev_xz= sigmaTot_xz;
    // sigmaDev_yz= sigmaTot_yz;


    rElementalVariables.UpdatedDeviatoricCauchyStress[0]=sigmaDev_xx;
    rElementalVariables.UpdatedDeviatoricCauchyStress[1]=sigmaDev_yy;
    rElementalVariables.UpdatedDeviatoricCauchyStress[2]=sigmaDev_zz;
    rElementalVariables.UpdatedDeviatoricCauchyStress[3]=sigmaDev_xy;
    rElementalVariables.UpdatedDeviatoricCauchyStress[4]=sigmaDev_xz;
    rElementalVariables.UpdatedDeviatoricCauchyStress[5]=sigmaDev_yz;

    rElementalVariables.UpdatedTotalCauchyStress[0]=sigmaTot_xx;
    rElementalVariables.UpdatedTotalCauchyStress[1]=sigmaTot_yy;
    rElementalVariables.UpdatedTotalCauchyStress[2]=sigmaTot_zz;
    rElementalVariables.UpdatedTotalCauchyStress[3]=sigmaTot_xy;
    rElementalVariables.UpdatedTotalCauchyStress[4]=sigmaTot_xz;
    rElementalVariables.UpdatedTotalCauchyStress[5]=sigmaTot_yz;

  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitFluidElement<TDim>::ComputeBulkMatrixRHS(Matrix& BulkMatrix,
										  const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < NumNodes; ++i)
      {
	for (SizeType j = 0; j < NumNodes; ++j)
	  {
	    // LHS contribution
	    double Mij  = Weight/12.0;
	    if(i==j)
	      Mij  *= 2.0;
	    BulkMatrix(i,j) +=  Mij;
	  }
      }
  }



  template class TwoStepUpdatedLagrangianVPExplicitFluidElement<2>;
  template class TwoStepUpdatedLagrangianVPExplicitFluidElement<3>;

}

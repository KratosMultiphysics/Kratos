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
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_element.h"
#include "includes/cfd_variables.h"
#include <math.h>

namespace Kratos {

  /*
   * public TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim> functions
   */

  template< unsigned int TDim >
  Element::Pointer TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {

    TwoStepUpdatedLagrangianVPImplicitFluidElement NewElement(NewId, this->GetGeometry().Create( rThisNodes ), this->pGetProperties() );

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer( new TwoStepUpdatedLagrangianVPImplicitFluidElement(NewElement) );

  }




  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::Initialize()
  {
    KRATOS_TRY;
    KRATOS_CATCH( "" );
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
  {

  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;
    KRATOS_CATCH( "" );
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::ComputeMaterialParameters(double& Density,
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
  double TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::ComputeNonLinearViscosity(double & equivalentStrainRate)
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
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::ComputeMaterialParametersGranularGas(double& Density,
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
  double TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::ComputeJopMuIrheologyViscosity(ElementalVariables & rElementalVariables)
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
  double TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::ComputeBercovierMuIrheologyViscosity(ElementalVariables & rElementalVariables)
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
  double TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::ComputePapanastasiouMuIrheologyViscosity(ElementalVariables & rElementalVariables)
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
  double TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::ComputeBarkerMuIrheologyViscosity(ElementalVariables & rElementalVariables)
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
  double TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::ComputeBarkerBercovierMuIrheologyViscosity(ElementalVariables & rElementalVariables)
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
  int TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
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
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<2>::ComputeMeanValueMaterialTangentMatrix(ElementalVariables & rElementalVariables,double& MeanValue,const ShapeFunctionDerivativesType& rDN_DX,const double secondLame,double & bulkModulus,const double Weight,double& MeanValueMass,const double TimeStep)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    double theta=0.5;
    double Count=0;
    for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    double lagDNXi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,0);
	    double lagDNYi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,1);
	    double lagDNXj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,0);
	    double lagDNYj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,1);
	    // lagDNXi=rDN_DX(i,0);
	    // lagDNYi=rDN_DX(i,1);
	    // lagDNXj=rDN_DX(j,0);
	    // lagDNYj=rDN_DX(j,1);


            // First Row
            MeanValue += fabs(Weight * ( (FourThirds * secondLame + bulkModulus)*  lagDNXi * lagDNXj + lagDNYi * lagDNYj * secondLame )*theta);
            MeanValue += fabs(Weight * ( (nTwoThirds* secondLame + bulkModulus) *  lagDNXi * lagDNYj + lagDNYi * lagDNXj * secondLame )*theta);

            // Second Row
            MeanValue += fabs(Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNXj + lagDNXi * lagDNYj *  secondLame )*theta);
            MeanValue += fabs(Weight * ( (FourThirds * secondLame + bulkModulus) * lagDNYi * lagDNYj + lagDNXi * lagDNXj * secondLame )*theta);

	    Count+=4.0;

	  }

      }

    MeanValue*=1.0/Count;

    if(MeanValueMass!=0 && MeanValue!=0){
      bulkModulus*=MeanValueMass*2.0/TimeStep/MeanValue;
    }else{
      std::cout<<" DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      std::cout<<" MeanValueMass="<<MeanValueMass;
      std::cout<<"\t MeanValueMaterial= "<<MeanValue;
      std::cout<<"\t VolumetricCoeff= "<<bulkModulus<<std::endl;
      bulkModulus*=TimeStep;
    }
  }

  template<>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<3>::ComputeMeanValueMaterialTangentMatrix(ElementalVariables & rElementalVariables,double& MeanValue,const ShapeFunctionDerivativesType& rDN_DX,const double secondLame,double & bulkModulus,const double Weight,double& MeanValueMass,const double TimeStep)
  {

    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    double theta=0.5;
    double Count=0;
    for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    double lagDNXi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,0)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,0);
	    double lagDNYi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,1)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,1);
	    double lagDNZi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,2)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,2)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,2);
	    double lagDNXj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,0)+rDN_DX(j,2)*rElementalVariables.InvFgrad(2,0);
	    double lagDNYj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,1)+rDN_DX(j,2)*rElementalVariables.InvFgrad(2,1);
	    double lagDNZj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,2)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,2)+rDN_DX(j,2)*rElementalVariables.InvFgrad(2,2);
	    // lagDNXi=rDN_DX(i,0);
	    // lagDNYi=rDN_DX(i,1);
	    // lagDNZi=rDN_DX(i,2);
	    // lagDNXj=rDN_DX(j,0);
	    // lagDNYj=rDN_DX(j,1);
	    // lagDNZj=rDN_DX(j,2);

            // First Row
            MeanValue += fabs(Weight * ( (FourThirds * secondLame + bulkModulus)*  lagDNXi * lagDNXj + (lagDNYi * lagDNYj +lagDNZi * lagDNZj) * secondLame ) *theta);
            MeanValue += fabs(Weight * ( (nTwoThirds* secondLame + bulkModulus) *  lagDNXi * lagDNYj + lagDNYi * lagDNXj * secondLame )*theta);
            MeanValue += fabs(Weight * ( (nTwoThirds* secondLame + bulkModulus) *  lagDNXi * lagDNZj + lagDNZi * lagDNXj * secondLame )*theta);

            // Second Row
            MeanValue += fabs(Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNXj + lagDNXi * lagDNYj *  secondLame )*theta);
            MeanValue += fabs(Weight * ( (FourThirds * secondLame + bulkModulus) * lagDNYi * lagDNYj + (lagDNXi * lagDNXj + lagDNZi * lagDNZj) * secondLame )*theta);
            MeanValue += fabs(Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNZj + lagDNZi * lagDNYj *  secondLame )*theta);

            // Third Row
            MeanValue += fabs(Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNZi * lagDNXj + lagDNXi * lagDNZj *  secondLame )*theta);
            MeanValue += fabs(Weight * ( (nTwoThirds* secondLame + bulkModulus)  * lagDNZi * lagDNYj + lagDNYi * lagDNZj *  secondLame )*theta);
	    MeanValue += fabs(Weight * ( (FourThirds * secondLame + bulkModulus) * lagDNZi * lagDNZj + (lagDNXi * lagDNXj + lagDNYi * lagDNYj) * secondLame )*theta);
	    Count+=9.0;
	  }

      }
    MeanValue*=1.0/Count;

    if(MeanValueMass!=0 && MeanValue!=0){
      bulkModulus*=MeanValueMass*2.0/TimeStep/MeanValue;
    }else{
      std::cout<<" DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      std::cout<<" MeanValueMass="<<MeanValueMass;
      std::cout<<"\t MeanValueMaterial= "<<MeanValue;
      std::cout<<"\t VolumetricCoeff= "<<bulkModulus<<std::endl;
      bulkModulus*=TimeStep;
    }
  }


  template<>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<2>::ComputeBulkReductionCoefficient(MatrixType MassMatrix,
											  MatrixType StiffnessMatrix,
											  double& meanValueStiff,
											  double& bulkCoefficient,
											  double timeStep)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    IndexType FirstRow = 0;
    IndexType FirstCol = 0;
    double meanValueMass=0;
    double countStiff=0;
    double countMass=0;
    for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow,FirstCol));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow,FirstCol+1));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+1,FirstCol));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+1,FirstCol+1));
            // Update Counter
	    countStiff+=4.0;

	    meanValueMass+= fabs(MassMatrix(FirstRow,FirstCol));
	    meanValueMass+= fabs(MassMatrix(FirstRow+1,FirstCol+1));
            // Update Counter
	    countMass+=2.0;

            FirstRow += 2;

	  }
	FirstRow = 0;
        FirstCol += 2;
      }

    meanValueStiff*=1.0/countStiff;
    meanValueMass*=1.0/countMass;

    if(meanValueMass!=0 && meanValueStiff!=0){
      bulkCoefficient=meanValueMass*4/(timeStep*meanValueStiff);
    }else{
      std::cout<<" DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      std::cout<<" coordinates "<<this->GetGeometry()[0].X()<<" "<<this->GetGeometry()[0].Y()<<std::endl;
      std::cout<<" MeanValueMass="<<meanValueMass;
      std::cout<<"\t MeanValueMaterial= "<<meanValueStiff;
      bulkCoefficient=timeStep;
    }
  }


  template<>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<3>::ComputeBulkReductionCoefficient(MatrixType MassMatrix,
											  MatrixType StiffnessMatrix,
											  double& meanValueStiff,
											  double& bulkCoefficient,
											  double timeStep)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    IndexType FirstRow = 0;
    IndexType FirstCol = 0;
    double meanValueMass=0;
    double countStiff=0;
    double countMass=0;
    for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow,FirstCol));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow,FirstCol+1));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow,FirstCol+2));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+1,FirstCol));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+1,FirstCol+1));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+1,FirstCol+2));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+2,FirstCol));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+2,FirstCol+1));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+2,FirstCol+2));
	    countStiff+=9.0;

	    meanValueMass+= fabs(MassMatrix(FirstRow,FirstCol));
	    meanValueMass+= fabs(MassMatrix(FirstRow+1,FirstCol+1));
	    meanValueMass+= fabs(MassMatrix(FirstRow+2,FirstCol+2));
	    countMass+=3.0;

	    // Update Counter
            FirstRow += 3;
	  }
	FirstRow = 0;
        FirstCol += 3;
      }

    meanValueStiff*=1.0/countStiff;
    meanValueMass*=1.0/countMass;

    if(meanValueMass!=0 && meanValueStiff!=0){
      bulkCoefficient=meanValueMass*2.0/timeStep/meanValueStiff;
    }else{
      std::cout<<" DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      std::cout<<" MeanValueMass="<<meanValueMass;
      std::cout<<"\t MeanValueMaterial= "<<meanValueStiff;
      bulkCoefficient=timeStep;
    }
  }



  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::ComputeBulkMatrix(Matrix& BulkMatrix,
									       const ShapeFunctionsType& rN,
									       const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    for (SizeType i = 0; i < NumNodes; ++i)
      {
	for (SizeType j = 0; j < NumNodes; ++j)
	  {
	    // LHS contribution
	    double Mij  = Weight*rN[i]*rN[j];
	    BulkMatrix(i,j) +=  Mij;
	  }

      }
  }

  template< >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<2>::ComputeBulkMatrixConsistent(Matrix& BulkMatrix,
										      const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < NumNodes; ++i)
      {
	for (SizeType j = 0; j < NumNodes; ++j)
	  {
	    // LHS contribution
	    double Mij  = Weight/12;
	    if(i==j)
	      Mij  *= 2.0;
	    BulkMatrix(i,j) +=  Mij;
	  }

      }

  }

  template< >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<3>::ComputeBulkMatrixConsistent(Matrix& BulkMatrix,
										      const double Weight)
  {
    std::cout<<"TO IMPLEMENT AND CHECK "<<std::endl;
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < NumNodes; ++i)
      {
	for (SizeType j = 0; j < NumNodes; ++j)
	  {
	    // LHS contribution
	    double Mij  = Weight/12;
	    if(i==j)
	      Mij  *= 2.0;
	    BulkMatrix(i,j) +=  Mij;
	  }

      }

  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::ComputeBulkMatrixLump(Matrix& BulkMatrix,
										   const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    double coeff=1.0+TDim;
    if((NumNodes==3 && TDim==2) || (NumNodes==4 && TDim==3)){
      for (SizeType i = 0; i < NumNodes; ++i)
	{
	  // LHS contribution
	  double Mij  = Weight /coeff;
	  BulkMatrix(i,i) +=  Mij;
	}
    }else{
      std::cout<<"... ComputeBulkMatrixLump TO IMPLEMENT"<<std::endl;
    }
  }



  template< >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<2>::ComputeBoundLHSMatrix(Matrix& BoundLHSMatrix,
										const ShapeFunctionsType& rN,
										const double Weight)
  {
    GeometryType& rGeom = this->GetGeometry();
    //const SizeType NumNodes = rGeom.PointsNumber();

    // for (SizeType i = 0; i < (NumNodes-1); i++)
    //   {
    //     for (SizeType j = (i+1); j < NumNodes; j++)
    // 	  {
    // 	    if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE)){
    // 	      if(rGeom[i].IsNot(INLET))
    // 		BoundLHSMatrix(i,i) +=  Weight / 3.0;
    // 	      if(rGeom[j].IsNot(INLET))
    // 		BoundLHSMatrix(j,j) +=  Weight / 3.0;
    // 	    }
    // 	  }

    //   }

    // for (SizeType i = 0; i < (NumNodes-1); i++)
    //   {
    //     for (SizeType j = (i+1); j < NumNodes; j++)
    // 	  {
    // 	    if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE)){
    // 	      if(rGeom[i].IsNot(INLET)){
    // 		BoundLHSMatrix(i,i) +=  Weight*rN[i]*rN[i];
    // 		BoundLHSMatrix(i,j) +=  Weight*rN[i]*rN[j];
    // 	      }
    // 	      if(rGeom[j].IsNot(INLET)){
    // 		BoundLHSMatrix(j,j) +=  Weight*rN[j]*rN[j];
    // 		BoundLHSMatrix(j,i) +=  Weight*rN[j]*rN[i];
    // 	      }
    // 	    }
    // 	  }

    //   }

    // if(this->Is(TO_ERASE)){
    //   WeakPointerVector<Node<3> >& rN0 = rGeom[0].GetValue(NEIGHBOUR_NODES);
    //   WeakPointerVector<Node<3> >& rN1 = rGeom[1].GetValue(NEIGHBOUR_NODES);
    //   WeakPointerVector<Node<3> >& rN2 = rGeom[2].GetValue(NEIGHBOUR_NODES);
    //   if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE) &&
    // 	 (rN0.size()<NumNodes || rN1.size()<NumNodes)){
    // 	if(rGeom[0].IsNot(INLET))
    // 	  BoundLHSMatrix(0,0) +=  Weight / 3.0;
    // 	if(rGeom[1].IsNot(INLET))
    // 	  BoundLHSMatrix(1,1) +=  Weight / 3.0;
    //   }else if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE) &&
    // 	 (rN0.size()<NumNodes || rN2.size()<NumNodes)){
    // 	if(rGeom[0].IsNot(INLET))
    // 	  BoundLHSMatrix(0,0) +=  Weight / 3.0;
    // 	if(rGeom[2].IsNot(INLET))
    // 	  BoundLHSMatrix(2,2) +=  Weight / 3.0;
    //   }else if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE) &&
    // 	 (rN2.size()<NumNodes || rN1.size()<NumNodes)){
    // 	if(rGeom[1].IsNot(INLET))
    // 	  BoundLHSMatrix(1,1) +=  Weight / 3.0;
    // 	if(rGeom[2].IsNot(INLET))
    // 	  BoundLHSMatrix(2,2) +=  Weight / 3.0;
    //   }
    // }else{

    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)){
      if(rGeom[0].IsNot(INLET))
	BoundLHSMatrix(0,0) +=  Weight / 3.0;
      if(rGeom[1].IsNot(INLET))
	BoundLHSMatrix(1,1) +=  Weight / 3.0;
    }
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
      if(rGeom[0].IsNot(INLET))
	BoundLHSMatrix(0,0) +=  Weight / 3.0;
      if(rGeom[2].IsNot(INLET))
	BoundLHSMatrix(2,2) +=  Weight / 3.0;
    }
    if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
      if(rGeom[1].IsNot(INLET))
	BoundLHSMatrix(1,1) +=  Weight / 3.0;
      if(rGeom[2].IsNot(INLET))
	BoundLHSMatrix(2,2) +=  Weight / 3.0;
    }
    // }

  }

  template<  >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<3>::ComputeBoundLHSMatrix(Matrix& BoundLHSMatrix,
										const ShapeFunctionsType& rN,
										const double Weight)
  {
    GeometryType& rGeom = this->GetGeometry();
    //const SizeType NumNodes = rGeom.PointsNumber();

    // for (SizeType i = 0; i < (NumNodes-2); i++)
    //   {
    //     for (SizeType j = (i+1); j < (NumNodes-1); j++)
    // 	  {
    // 	    for (SizeType k = (j+1); k < NumNodes; k++)
    // 	      {
    // 		if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE) && rGeom[k].Is(FREE_SURFACE)){
    // 		  if(rGeom[i].IsNot(INLET)){
    // 		    BoundLHSMatrix(i,i) +=  Weight*rN[i]*rN[i];
    // 		    BoundLHSMatrix(i,j) +=  Weight*rN[i]*rN[j];
    // 		    BoundLHSMatrix(i,k) +=  Weight*rN[i]*rN[k];
    // 		  }
    // 		  if(rGeom[j].IsNot(INLET)){
    // 		    BoundLHSMatrix(j,i) +=  Weight*rN[j]*rN[i];
    // 		    BoundLHSMatrix(j,j) +=  Weight*rN[j]*rN[j];
    // 		    BoundLHSMatrix(j,k) +=  Weight*rN[j]*rN[k];
    // 		  }
    // 		  if(rGeom[k].IsNot(INLET)){
    // 		    BoundLHSMatrix(k,i) +=  Weight*rN[k]*rN[i];
    // 		    BoundLHSMatrix(k,j) +=  Weight*rN[k]*rN[j];
    // 		    BoundLHSMatrix(k,k) +=  Weight*rN[k]*rN[k];
    // 		  }
    // 		}
    // 	      }
    // 	  }

    //   }


    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
      if(rGeom[0].IsNot(INLET))
    	BoundLHSMatrix(0,0) +=  Weight / 4.0;
      if(rGeom[1].IsNot(INLET))
    	BoundLHSMatrix(1,1) +=  Weight / 4.0;
      if(rGeom[2].IsNot(INLET))
    	BoundLHSMatrix(2,2) +=  Weight / 4.0;
    }
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
      if(rGeom[0].IsNot(INLET))
    	BoundLHSMatrix(0,0) +=  Weight / 4.0;
      if(rGeom[1].IsNot(INLET))
    	BoundLHSMatrix(1,1) +=  Weight / 4.0;
      if(rGeom[3].IsNot(INLET))
    	BoundLHSMatrix(3,3) +=  Weight / 4.0;
    }
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
      if(rGeom[0].IsNot(INLET))
    	BoundLHSMatrix(0,0) +=  Weight / 4.0;
      if(rGeom[2].IsNot(INLET))
    	BoundLHSMatrix(2,2) +=  Weight / 4.0;
      if(rGeom[3].IsNot(INLET))
    	BoundLHSMatrix(3,3) +=  Weight / 4.0;
    }
    if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
      if(rGeom[1].IsNot(INLET))
    	BoundLHSMatrix(1,1) +=  Weight / 4.0;
      if(rGeom[2].IsNot(INLET))
    	BoundLHSMatrix(2,2) +=  Weight / 4.0;
      if(rGeom[3].IsNot(INLET))
    	BoundLHSMatrix(3,3) +=  Weight / 4.0;
    }

  }


  template< >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<2>::ComputeBoundRHSVectorComplete(VectorType& BoundRHSVector, 
											const double TimeStep, 
											const double BoundRHSCoeffAcc, 
											const double BoundRHSCoeffDev, 
											const VectorType SpatialDefRate) 
  { 
    GeometryType& rGeom = this->GetGeometry(); 
 
 
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE) ){ 
      array_1d<double, 3>  AccA(3,0.0); 
      array_1d<double, 3>  AccB(3,0.0); 
      array_1d<double, 3>  MeanAcc(3,0.0); 
      array_1d<double, 3> NormalVector(3,0.0);       
      const double factor = 0.5/TimeStep; 
      const double one_third = 1.0/3.0; 
 
      this->GetOutwardsUnitNormalForTwoPoints(NormalVector,0,1,2); 
 
      double SpatialDefRateNormalProjection=this->CalcNormalProjectionDefRate(SpatialDefRate,NormalVector); 
       
      noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);  
      noalias(AccB)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);  
      noalias(MeanAcc)= 0.5*AccA+0.5*AccB; 
 
      const double accelerationsNormalProjection=MeanAcc[0]*NormalVector[0]+MeanAcc[1]*NormalVector[1]; 
 
      if(rGeom[0].IsNot(INLET)) //to change into moving wall!!!!! 
	BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
 
      if(rGeom[1].IsNot(INLET)) 
	BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
    } 
     
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE) ){ 
       
      array_1d<double, 3>  AccA(3,0.0); 
      array_1d<double, 3>  AccB(3,0.0); 
      array_1d<double, 3>  MeanAcc(3,0.0); 
      array_1d<double, 3> NormalVector(3,0.0);       
      const double factor = 0.5/TimeStep; 
      const double one_third = 1.0/3.0; 

      this->GetOutwardsUnitNormalForTwoPoints(NormalVector,0,2,1); 
 
      double SpatialDefRateNormalProjection=this->CalcNormalProjectionDefRate(SpatialDefRate,NormalVector); 
       
      noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);  
      noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
      noalias(MeanAcc)= 0.5*AccA+0.5*AccB; 
 
      const double accelerationsNormalProjection=MeanAcc[0]*NormalVector[0]+MeanAcc[1]*NormalVector[1]; 
 
      if(rGeom[0].IsNot(INLET)) //to change into moving wall!!!!! 
	BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
 
      if(rGeom[2].IsNot(INLET)) 
	BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
    } 
     
    if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE) ){ 
       
      array_1d<double, 3>  AccA(3,0.0); 
      array_1d<double, 3>  AccB(3,0.0); 
      array_1d<double, 3>  MeanAcc(3,0.0); 
      array_1d<double, 3> NormalVector(3,0.0);       
      const double factor = 0.5/TimeStep; 
      const double one_third = 1.0/3.0; 
       
      this->GetOutwardsUnitNormalForTwoPoints(NormalVector,1,2,0); 
 
      double SpatialDefRateNormalProjection=this->CalcNormalProjectionDefRate(SpatialDefRate,NormalVector); 
 
      noalias(AccA)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);  
      noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1);       
      noalias(MeanAcc)= 0.5*AccA+0.5*AccB; 
       
      const double accelerationsNormalProjection=MeanAcc[0]*NormalVector[0]+MeanAcc[1]*NormalVector[1]; 
       
      if(rGeom[1].IsNot(INLET)) 
	BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
       
      if(rGeom[2].IsNot(INLET)) 
	BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
    } 
   
  } 





  template< > 
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<3>::ComputeBoundRHSVectorComplete(VectorType& BoundRHSVector, 
											const double TimeStep, 
											const double BoundRHSCoeffAcc, 
											const double BoundRHSCoeffDev, 
											const VectorType SpatialDefRate) 
  { 
    GeometryType& rGeom = this->GetGeometry(); 
  
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){ 
 
      array_1d<double, 3> AccA(3,0.0); 
      array_1d<double, 3> AccB(3,0.0); 
      array_1d<double, 3> AccC(3,0.0); 
      array_1d<double, 3> MeanAcc(3,0.0); 
      array_1d<double, 3> NormalVector(3,0.0); 
      const double factor = 0.5/TimeStep; 
      const double one_third = 1.0/3.0; 
 
      this->GetOutwardsUnitNormalForThreePoints(NormalVector,0,1,2,3); 
 
      double SpatialDefRateNormalProjection=this->CalcNormalProjectionDefRate(SpatialDefRate,NormalVector); 
 
      noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);  
      noalias(AccB)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);  
      noalias(AccC)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
       
      noalias(MeanAcc)= one_third*AccA + one_third*AccB + one_third*AccC; 
 
      const double accelerationsNormalProjection=MeanAcc[0]*NormalVector[0] + MeanAcc[1]*NormalVector[1] + MeanAcc[2]*NormalVector[2]; 
       
      if(rGeom[0].IsNot(INLET)) 
	BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
 
      if(rGeom[1].IsNot(INLET)) 
	BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
 
      if(rGeom[2].IsNot(INLET)) 
	BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
 
    }

    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){ 
 
      array_1d<double, 3> AccA(3,0.0); 
      array_1d<double, 3> AccB(3,0.0); 
      array_1d<double, 3> AccC(3,0.0); 
      array_1d<double, 3> MeanAcc(3,0.0); 
      array_1d<double, 3> NormalVector(3,0.0); 
      const double factor = 0.5/TimeStep; 
      const double one_third = 1.0/3.0; 
 
      this->GetOutwardsUnitNormalForThreePoints(NormalVector,0,1,3,2); 
 
      double SpatialDefRateNormalProjection=this->CalcNormalProjectionDefRate(SpatialDefRate,NormalVector); 
 
      noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);  
      noalias(AccB)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);  
      noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
       
      noalias(MeanAcc)= one_third*AccA + one_third*AccB + one_third*AccC; 
 
      const double accelerationsNormalProjection=MeanAcc[0]*NormalVector[0] + MeanAcc[1]*NormalVector[1] + MeanAcc[2]*NormalVector[2]; 
       
      if(rGeom[0].IsNot(INLET)) 
	BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
 
      if(rGeom[1].IsNot(INLET)) 
	BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
 
      if(rGeom[3].IsNot(INLET)) 
	BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
       
    }

    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){ 
 
      array_1d<double, 3> AccA(3,0.0); 
      array_1d<double, 3> AccB(3,0.0); 
      array_1d<double, 3> AccC(3,0.0); 
      array_1d<double, 3> MeanAcc(3,0.0); 
      array_1d<double, 3> NormalVector(3,0.0); 
      const double factor = 0.5/TimeStep; 
      const double one_third = 1.0/3.0; 
 
      this->GetOutwardsUnitNormalForThreePoints(NormalVector,0,2,3,1); 
 
      double SpatialDefRateNormalProjection=this->CalcNormalProjectionDefRate(SpatialDefRate,NormalVector); 
 
      noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);  
      noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1);  
      noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
       
      noalias(MeanAcc)= one_third*AccA + one_third*AccB + one_third*AccC; 
 
      const double accelerationsNormalProjection=MeanAcc[0]*NormalVector[0] + MeanAcc[1]*NormalVector[1] + MeanAcc[2]*NormalVector[2]; 
       
      if(rGeom[0].IsNot(INLET)) 
	BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
 
      if(rGeom[2].IsNot(INLET)) 
	BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
 
      if(rGeom[3].IsNot(INLET)) 
	BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
       
    } 

    if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){ 
 
      array_1d<double, 3> AccA(3,0.0); 
      array_1d<double, 3> AccB(3,0.0); 
      array_1d<double, 3> AccC(3,0.0); 
      array_1d<double, 3> MeanAcc(3,0.0); 
      array_1d<double, 3> NormalVector(3,0.0); 
      const double factor = 0.5/TimeStep; 
      const double one_third = 1.0/3.0; 
 
      this->GetOutwardsUnitNormalForThreePoints(NormalVector,1,2,3,0); 
 
      double SpatialDefRateNormalProjection=this->CalcNormalProjectionDefRate(SpatialDefRate,NormalVector); 
       
      noalias(AccA)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);  
      noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
      noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1);  
       
      noalias(MeanAcc)= one_third*AccA + one_third*AccB + one_third*AccC; 
 
      const double accelerationsNormalProjection=MeanAcc[0]*NormalVector[0] + MeanAcc[1]*NormalVector[1] + MeanAcc[2]*NormalVector[2]; 
   
      if(rGeom[1].IsNot(INLET)) 
	BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
 
      if(rGeom[2].IsNot(INLET)) 
	BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
       
      if(rGeom[3].IsNot(INLET)) 
	BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc*accelerationsNormalProjection + BoundRHSCoeffDev*SpatialDefRateNormalProjection); 
 
    } 
 
 
  } 


  

  template< >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<2>::ComputeBoundRHSVector(VectorType& BoundRHSVector,
										const ShapeFunctionsType& rN,
										const double TimeStep,
										const double BoundRHSCoeffAcc,
										const double BoundRHSCoeffDev)
  {
    GeometryType& rGeom = this->GetGeometry();
    //const SizeType NumNodes = rGeom.PointsNumber();
    array_1d<double, 3>  AccA(3,0.0);
    array_1d<double, 3>  AccB(3,0.0);

    // for (SizeType i = 0; i < (NumNodes-1); i++)
    //   {
    // 	for (SizeType j = (i+1); j < NumNodes; j++)
    // 	  {
    // 	    if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE)){
    // 	      AccA= 0.5/TimeStep*(rGeom[i].FastGetSolutionStepValue(VELOCITY,0)-rGeom[i].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[i].FastGetSolutionStepValue(ACCELERATION,1);
    // 	      AccB= 0.5/TimeStep*(rGeom[j].FastGetSolutionStepValue(VELOCITY,0)-rGeom[j].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[j].FastGetSolutionStepValue(ACCELERATION,1);
    // 	      const array_1d<double, 3> &NormalA    = rGeom[i].FastGetSolutionStepValue(NORMAL);
    // 	      const array_1d<double, 3> &NormalB    = rGeom[j].FastGetSolutionStepValue(NORMAL);
    // 	      double coeff=3.0;
    // 	      if(rGeom[i].IsNot(INLET)) //to change into moving wall!!!!!
    // 		BoundRHSVector[i] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) +
    // 				      BoundRHSCoeffDev)/coeff ;
    // 	      if(rGeom[j].IsNot(INLET))
    // 		BoundRHSVector[j] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) +
    // 				      BoundRHSCoeffDev)/coeff ;
    // 	    }
    // 	  }

    //   }

    // for (SizeType i = 0; i < (NumNodes-1); i++)
    //   {
    // 	for (SizeType j = (i+1); j < NumNodes; j++)
    // 	  {
    // 	    if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE)){
    // 	      AccA= 0.5/TimeStep*(rGeom[i].FastGetSolutionStepValue(VELOCITY,0)-rGeom[i].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[i].FastGetSolutionStepValue(ACCELERATION,1);
    // 	      AccB= 0.5/TimeStep*(rGeom[j].FastGetSolutionStepValue(VELOCITY,0)-rGeom[j].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[j].FastGetSolutionStepValue(ACCELERATION,1);
    // 	      const array_1d<double, 3> &NormalA    = rGeom[i].FastGetSolutionStepValue(NORMAL);
    // 	      const array_1d<double, 3> &NormalB    = rGeom[j].FastGetSolutionStepValue(NORMAL);
    // 	      if(rGeom[i].IsNot(INLET))
    // 		BoundRHSVector[i] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) +
    // 				      BoundRHSCoeffDev) * rN[i];
    // 	      if(rGeom[j].IsNot(INLET))
    // 		BoundRHSVector[j] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) +
    // 				      BoundRHSCoeffDev) * rN[j] ;
    // 	    }
    // 	  }

    //   }

    const double factor = 0.5/TimeStep;
    const double one_third = 1.0/3.0;


    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE) ){
      noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);
      noalias(AccB)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);
      // noalias(AccA)=rGeom[0].FastGetSolutionStepValue(ACCELERATION,0);
      // noalias(AccB)=rGeom[1].FastGetSolutionStepValue(ACCELERATION,0);
      const array_1d<double, 3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
      const array_1d<double, 3> &NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
      if(rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
	BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + BoundRHSCoeffDev);
      if(rGeom[1].IsNot(INLET))
	BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + BoundRHSCoeffDev);
    }
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE) ){
      noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);
      noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1);
      const array_1d<double, 3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
      const array_1d<double, 3> &NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
      if(rGeom[0].IsNot(INLET))
	BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + BoundRHSCoeffDev);
      if(rGeom[2].IsNot(INLET))
	BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + BoundRHSCoeffDev);
    }
    if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE) ){
      noalias(AccA)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);
      noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1);
      const array_1d<double, 3> &NormalA    = rGeom[1].FastGetSolutionStepValue(NORMAL);
      const array_1d<double, 3> &NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
      if(rGeom[1].IsNot(INLET))
	BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + BoundRHSCoeffDev);
      if(rGeom[2].IsNot(INLET))
	BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + BoundRHSCoeffDev);
    }

  }



  template< >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<3>::ComputeBoundRHSVector(VectorType& BoundRHSVector,
										const ShapeFunctionsType& rN,
										const double TimeStep,
										const double BoundRHSCoeffAcc,
										const double BoundRHSCoeffDev)
  {
    GeometryType& rGeom = this->GetGeometry();
    //const SizeType NumNodes = rGeom.PointsNumber();
    array_1d<double, 3>  AccA(3,0.0);
    array_1d<double, 3>  AccB(3,0.0);
    array_1d<double, 3>  AccC(3,0.0);


    // for (SizeType i = 0; i < (NumNodes-2); i++)
    //   {
    // 	for (SizeType j = (i+1); j < (NumNodes-1); j++)
    // 	  {
    // 	    for (SizeType k = (j+1); k < NumNodes; k++)
    // 	      {
    // 		if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE) && rGeom[k].Is(FREE_SURFACE)){
    // 		  AccA= 0.5/TimeStep*(rGeom[i].FastGetSolutionStepValue(VELOCITY,0)-rGeom[i].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[i].FastGetSolutionStepValue(ACCELERATION,1);
    // 		  AccB= 0.5/TimeStep*(rGeom[j].FastGetSolutionStepValue(VELOCITY,0)-rGeom[j].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[j].FastGetSolutionStepValue(ACCELERATION,1);
    // 		  AccC= 0.5/TimeStep*(rGeom[k].FastGetSolutionStepValue(VELOCITY,0)-rGeom[k].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[k].FastGetSolutionStepValue(ACCELERATION,1);

    // 		  const array_1d<double, 3> &NormalA    = rGeom[i].FastGetSolutionStepValue(NORMAL);
    // 		  const array_1d<double, 3> &NormalB    = rGeom[j].FastGetSolutionStepValue(NORMAL);
    // 		  const array_1d<double, 3> &NormalC    = rGeom[k].FastGetSolutionStepValue(NORMAL);
    // 		  if(rGeom[i].IsNot(INLET))
    // 		    BoundRHSVector[i] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
    // 					  BoundRHSCoeffDev) * rN[i];
    // 		  if(rGeom[j].IsNot(INLET))
    // 		    BoundRHSVector[j] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
    // 					  BoundRHSCoeffDev) * rN[j] ;
    // 		  if(rGeom[k].IsNot(INLET))
    // 		    BoundRHSVector[k] += (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
    // 					  BoundRHSCoeffDev) * rN[k] ;
    // 		}
    // 	      }
    // 	  }

    //   }
    const double factor = 0.5/TimeStep;

    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
      noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);
      noalias(AccB)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);
      noalias(AccC)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1);
      const array_1d<double,3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalC    = rGeom[2].FastGetSolutionStepValue(NORMAL);
      if(rGeom[0].IsNot(INLET))
	BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
				     BoundRHSCoeffDev);
      if(rGeom[1].IsNot(INLET))
	BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
				     BoundRHSCoeffDev);
      if(rGeom[2].IsNot(INLET))
	BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
				     BoundRHSCoeffDev);
    }
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
      noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);
      noalias(AccB)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);
      noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1);
      const array_1d<double,3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
      if(rGeom[0].IsNot(INLET))
	BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
				     BoundRHSCoeffDev);
      if(rGeom[1].IsNot(INLET))
	BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
				     BoundRHSCoeffDev);
      if(rGeom[3].IsNot(INLET))
	BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
				     BoundRHSCoeffDev);
    }
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
      noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1);
      noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1);
      noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1);
      const array_1d<double,3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
      if(rGeom[0].IsNot(INLET))
	BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
				     BoundRHSCoeffDev);
      if(rGeom[2].IsNot(INLET))
	BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
				     BoundRHSCoeffDev);
      if(rGeom[3].IsNot(INLET))
	BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
				     BoundRHSCoeffDev);
    }
    if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
      noalias(AccA)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1);
      noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1);
      noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1);
      const array_1d<double,3> &NormalA    = rGeom[1].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
      if(rGeom[1].IsNot(INLET))
	BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
				     BoundRHSCoeffDev);
      if(rGeom[2].IsNot(INLET))
	BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
				     BoundRHSCoeffDev);
      if(rGeom[3].IsNot(INLET))
	BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
				     BoundRHSCoeffDev);
    }


  }



  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::CalculateTauFIC(double& Tau,
									     double ElemSize,
									     const double Density,
									     const double Viscosity,
									     const ProcessInfo& rCurrentProcessInfo)
  {
    double DeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
    if(rCurrentProcessInfo.GetValue(DELTA_TIME)<rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME)){
      DeltaTime = 0.5*rCurrentProcessInfo.GetValue(DELTA_TIME)+0.5*rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME);
    }

    double MeanVelocity=0;
    this->CalcMeanVelocity(MeanVelocity,0);

    // Tau = 1.0 / (2.0 * Density *(0.5 * MeanVelocity / ElemSize + 0.5/DeltaTime) +  8.0 * Viscosity / (ElemSize * ElemSize) );
    Tau = (ElemSize * ElemSize * DeltaTime) / ( Density * MeanVelocity * DeltaTime * ElemSize + Density * ElemSize * ElemSize +  8.0 * Viscosity * DeltaTime  );
    // if(Tau<0.0000001){
    //   Tau=0.0000001;
    // }
    // if(Tau>0.0001){
    //   Tau=0.0001;
    // }

    if(MeanVelocity==0){
      Tau=0;
    }
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::AddStabilizationMatrixLHS(MatrixType& rLeftHandSideMatrix,
										       Matrix& BulkAccMatrix,
										       const ShapeFunctionsType& rN,
										       const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    if( BulkAccMatrix.size1() != NumNodes )
      BulkAccMatrix.resize(NumNodes,NumNodes);

    BulkAccMatrix = ZeroMatrix(NumNodes,NumNodes);
    for (SizeType i = 0; i < NumNodes; ++i)
      {
	// LHS contribution
	for (SizeType j = 0; j < NumNodes; ++j)
	  {
	    double Mij = 0.0;
	    Mij = Weight * rN[i] * rN[j];
	    BulkAccMatrix(i,j) +=  Mij;
	  }
      }
    rLeftHandSideMatrix+=BulkAccMatrix;

  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::ComputeStabLaplacianMatrix(MatrixType& StabLaplacianMatrix,
											const ShapeFunctionDerivativesType& rDN_DX,
											const double Weight)

  {
    // LHS contribution
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < NumNodes; ++i)
      {
	for (SizeType j = 0; j < NumNodes; ++j)
	  {
	    double Lij = 0.0;
	    for (SizeType d = 0; d < TDim; ++d){
	      Lij += rDN_DX(i,d) * rDN_DX(j,d);
	    }
	    StabLaplacianMatrix(i,j) += Weight * Lij ;
	  }
      }
  }



  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::AddStabilizationNodalTermsLHS(MatrixType& rLeftHandSideMatrix,
											   const double Tau,
											   const double Weight,
											   const ShapeFunctionDerivativesType& rDN_DX,
											   const SizeType i)
  {
    // LHS contribution
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType j = 0; j < NumNodes; ++j)
      {
        double Lij = 0.0;
        for (SizeType d = 0; d < TDim; ++d){
    	  Lij += rDN_DX(i,d) * rDN_DX(j,d);
    	}
        Lij *= Tau;

        rLeftHandSideMatrix(i,j) += Weight * Lij;
      }
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::AddStabilizationNodalTermsRHS(VectorType& rRightHandSideVector,
											   const double Tau,
											   const double Density,
											   const double Weight,
											   const ShapeFunctionDerivativesType& rDN_DX,
											   const SizeType i)
  {

    double RHSi = 0;
    if( this->GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION) ){ // it must be checked once at the begining only
      array_1d<double, 3 >& VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
      for (SizeType d = 0; d < TDim; ++d)
	{
	  RHSi += - rDN_DX(i,d) * Tau * ( Density * VolumeAcceleration[d] );
	}
    }
    rRightHandSideVector[i] += Weight * RHSi;

  }


  template<>
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<2>::GetPositions(Vector& rValues,const ProcessInfo& rCurrentProcessInfo,const double theta)
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
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<3>::GetPositions(Vector& rValues,const ProcessInfo& rCurrentProcessInfo,const double theta)
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
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>:: InitializeElementalVariables(ElementalVariables & rElementalVariables)
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
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<2>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,double TimeStep, unsigned int g)
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


  template < >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<3>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables, double TimeStep, unsigned int g)
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
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {


    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != NumNodes )
      rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);

    rLeftHandSideMatrix = ZeroMatrix(NumNodes,NumNodes);

    if( rRightHandSideVector.size() != NumNodes )
      rRightHandSideVector.resize(NumNodes);

    rRightHandSideVector = ZeroVector(NumNodes);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    double TimeStep=rCurrentProcessInfo[DELTA_TIME];
    double theta=this->GetThetaContinuity();
    double ElemSize = this->ElementSize();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    double maxViscousValueForStabilization=0.1;
    double Density = this->mMaterialDensity;
    double VolumetricCoeff = this->mMaterialVolumetricCoefficient;
    double DeviatoricCoeff = this->mMaterialDeviatoricCoefficient;

    if(DeviatoricCoeff>maxViscousValueForStabilization){ 
      DeviatoricCoeff=maxViscousValueForStabilization;
    }

    double Tau=0;
    this->CalculateTauFIC(Tau,ElemSize,Density,DeviatoricCoeff,rCurrentProcessInfo);

    double totalVolume=0;
    bool computeElement=false;
    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
      {
	const double GaussWeight = GaussWeights[g];
	totalVolume+=GaussWeight;
	const ShapeFunctionsType& N = row(NContainer,g);
	const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];
	computeElement=this->CalcCompleteStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);

	if(computeElement==true){
	  // double BulkCoeff =GaussWeight/(VolumetricCoeff);
	  // this->ComputeBulkMatrix(BulkVelMatrix,N,BulkCoeff);
	  // double BulkStabCoeff=BulkCoeff*Tau*Density/TimeStep;
	  // this->ComputeBulkMatrix(BulkAccMatrix,N,BulkStabCoeff);

	  double BoundLHSCoeff=Tau*4.0*GaussWeight/(ElemSize*ElemSize);
 	  // if(TDim==3){
	  //   BoundLHSCoeff=Tau*2*GaussWeight/(0.81649658*ElemSize*ElemSize);
	  // }

	  this->ComputeBoundLHSMatrix(rLeftHandSideMatrix,N,BoundLHSCoeff);

	  double BoundRHSCoeffAcc=Tau*Density*2*GaussWeight/ElemSize;
	  double BoundRHSCoeffDev=Tau*8.0*DeviatoricCoeff*GaussWeight/(ElemSize*ElemSize);
	  // double NProjSpatialDefRate=this->CalcNormalProjectionDefRate(rElementalVariables.SpatialDefRate);
	  // double BoundRHSCoeffDev=Tau*8.0*NProjSpatialDefRate*DeviatoricCoeff*GaussWeight/(ElemSize*ElemSize);

	  // this->ComputeBoundRHSVector(rRightHandSideVector,N,TimeStep,BoundRHSCoeffAcc,BoundRHSCoeffDev);
	  this->ComputeBoundRHSVectorComplete(rRightHandSideVector,TimeStep,BoundRHSCoeffAcc,BoundRHSCoeffDev,rElementalVariables.SpatialDefRate); 

	  double StabLaplacianWeight=Tau*GaussWeight;
	  this->ComputeStabLaplacianMatrix(rLeftHandSideMatrix,rDN_DX,StabLaplacianWeight);

	  for (SizeType i = 0; i < NumNodes; ++i)
	    {
	      // RHS contribution
	      // Velocity divergence
	      rRightHandSideVector[i] += GaussWeight * N[i] * rElementalVariables.VolumetricDefRate;
	      this->AddStabilizationNodalTermsRHS(rRightHandSideVector,Tau,Density,GaussWeight,rDN_DX,i);
	    }
	}

      }


    if(computeElement==true){

      VectorType PressureValues= ZeroVector(NumNodes);
      VectorType PressureValuesForRHS= ZeroVector(NumNodes);
      this->GetPressureValues(PressureValuesForRHS,0);
      //the LHS matrix up to now just contains the laplacian term and the bound term
      noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,PressureValuesForRHS);

      this->GetPressureValues(PressureValues,1);
      noalias(PressureValuesForRHS)+=-PressureValues;
      MatrixType BulkMatrix = ZeroMatrix(NumNodes,NumNodes);
      MatrixType BulkMatrixConsistent = ZeroMatrix(NumNodes,NumNodes);
      double lumpedBulkCoeff =totalVolume/(VolumetricCoeff);
      double lumpedBulkStabCoeff=lumpedBulkCoeff*Tau*Density/TimeStep;

      this->ComputeBulkMatrixLump(BulkMatrix,lumpedBulkCoeff);
      // this->ComputeBulkMatrixConsistent(BulkMatrixConsistent,lumpedBulkCoeff);
      noalias(rLeftHandSideMatrix)+=BulkMatrix;
      // noalias(rLeftHandSideMatrix)+=BulkMatrixConsistent;
      noalias(rRightHandSideVector) -= prod(BulkMatrix,PressureValuesForRHS);
      // noalias(rRightHandSideVector) -=prod(BulkMatrixConsistent,PressureValuesForRHS);

      this->GetPressureVelocityValues(PressureValues,0);
      noalias(PressureValuesForRHS)+=-PressureValues*TimeStep;
      noalias(BulkMatrix) = ZeroMatrix(NumNodes,NumNodes);
      this->ComputeBulkMatrixLump(BulkMatrix,lumpedBulkStabCoeff);
      // this->ComputeBulkMatrixConsistent(BulkMatrixConsistent,lumpedBulkStabCoeff);
      noalias(rLeftHandSideMatrix)+=BulkMatrix;
      // noalias(rLeftHandSideMatrix)+=BulkMatrixConsistent;
      noalias(rRightHandSideVector) -= prod(BulkMatrix,PressureValuesForRHS);
      // noalias(rRightHandSideVector) -=prod(BulkMatrixConsistent,PressureValuesForRHS);

    }else{
      double lumpedBulkCoeff =totalVolume*Tau*Density/(TimeStep*VolumetricCoeff);
      MatrixType BulkVelMatrixLump = ZeroMatrix(NumNodes,NumNodes);
      this->ComputeBulkMatrixLump(BulkVelMatrixLump,lumpedBulkCoeff);
      noalias(rLeftHandSideMatrix)+=BulkVelMatrixLump;
      VectorType PressureValues= ZeroVector(NumNodes);
      VectorType PressureValuesForRHS= ZeroVector(NumNodes);
      this->GetPressureValues(PressureValuesForRHS,0);
      this->GetPressureValues(PressureValues,1);
      noalias(PressureValuesForRHS)+=-PressureValues;
      noalias(rRightHandSideVector) -= prod(BulkVelMatrixLump,PressureValuesForRHS);
    }



  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::GetPressureVelocityValues(Vector& rValues,
										       const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rValues.size() != NumNodes) rValues.resize(NumNodes);

    for (SizeType i = 0; i < NumNodes; ++i){
      rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE_VELOCITY,Step);

    }
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitFluidElement<TDim>::GetPressureAccelerationValues(Vector& rValues,
											   const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rValues.size() != NumNodes) rValues.resize(NumNodes);

    for (SizeType i = 0; i < NumNodes; ++i){
      rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE_ACCELERATION,Step);

    }
  }




  template class TwoStepUpdatedLagrangianVPImplicitFluidElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitFluidElement<3>;

}

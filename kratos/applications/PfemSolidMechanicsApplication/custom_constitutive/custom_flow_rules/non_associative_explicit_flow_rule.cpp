//

#include <iostream>
#include<cmath>

#include "../PfemSolidMechanicsApplication/custom_constitutive/custom_flow_rules/non_associative_explicit_flow_rule.hpp"
#include "utilities/math_utils.h"

#include "includes/ublas_interface.h"
#include "solid_mechanics_application.h"
namespace Kratos
{


NonAssociativeExplicitPlasticFlowRule::NonAssociativeExplicitPlasticFlowRule()
  :FlowRule()
{

}

NonAssociativeExplicitPlasticFlowRule::NonAssociativeExplicitPlasticFlowRule(const NonAssociativeExplicitPlasticFlowRule& rOther)
  :FlowRule(rOther)
{

}

FlowRule::Pointer NonAssociativeExplicitPlasticFlowRule::Clone() const
{
  FlowRule::Pointer p_clone(new NonAssociativeExplicitPlasticFlowRule(*this));
  return p_clone;
}


NonAssociativeExplicitPlasticFlowRule::~NonAssociativeExplicitPlasticFlowRule()
{
}




//NonAssociativeExplicitPlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables,  Matrix & rStressMatrix)
//NonAssociatievExplicitPlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables, const Matrix& DeltaHenckyStrain, ) 
bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables, const Matrix& rDeltaDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{


   bool PlasticityActive = false;
   rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);

   //1- Initialize some variables
   Vector StressVector = ZeroVector(6);
   InternalVariables PlasticVariables = mInternalVariables;
   PlasticVariables.DeltaPlasticStrain = 0.0;
   
   //mElasticLeftCauchyGreen
   Vector PreviousElasticHenckyStrain = ZeroVector(6);
   Vector DeltaHenckyStrain = ZeroVector(6);

   Matrix DeltaHenckyStrainMatrix;
   DeltaHenckyStrainMatrix = prod(rDeltaDeformationGradient, trans(rDeltaDeformationGradient));

   Matrix InitialLeftCauchyGreen = rNewElasticLeftCauchyGreen;
   PreviousElasticHenckyStrain = this->ConvertCauchyGreenTensorToHenckyStrain(InitialLeftCauchyGreen);
   DeltaHenckyStrain = this->ConvertCauchyGreenTensorToHenckyStrain(DeltaHenckyStrainMatrix);

   //computed with the actual 
   Vector NewElasticHenckyStrain = Vector(6);
   double Tolerance = 1E-5;


   Matrix AuxPrevious;
   AuxPrevious = rNewElasticLeftCauchyGreen;
   
   //2- Check for yield Condition (at the beggining)
   Vector NewStressVector;
   Vector PreviousStressVector;
   //TO DO SOMETHING WITH THE STR VEC

   this->CalculateKirchhoffStressVector(PreviousElasticHenckyStrain, PreviousStressVector);

   rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, PreviousStressVector, PlasticVariables.EquivalentPlasticStrain);

    //1. Compute an Elastic Trial State
    double StressErrorMeasure;
    double NewEquivalentPlasticStrain;

bool DebuggingMode = false;
    if (DebuggingMode ) {
    
		  PlasticityActive = true;

       if (rReturnMappingVariables.TrialStateFunction < 0.0) {
	  PlasticityActive = false;
       }
       else {
	  PlasticityActive = true;
       }

       this->CalculateOneExplicitStep(DeltaHenckyStrain, PreviousElasticHenckyStrain, PlasticVariables, NewElasticHenckyStrain, NewStressVector, NewEquivalentPlasticStrain, PlasticityActive, StressErrorMeasure);

       PlasticVariables.EquivalentPlasticStrain = NewEquivalentPlasticStrain;

    }
    else {


    this->CalculateOneExplicitStep(DeltaHenckyStrain, PreviousElasticHenckyStrain, PlasticVariables, NewElasticHenckyStrain, NewStressVector, NewEquivalentPlasticStrain, false, StressErrorMeasure);
    

    //2. Check if the Trial State is elastic
    double ElasticTrialStateFunction = 0.0;
    ElasticTrialStateFunction = mpYieldCriterion->CalculateYieldCondition(ElasticTrialStateFunction, NewStressVector, PlasticVariables.EquivalentPlasticStrain);

    //if ElasticStep
    if (ElasticTrialStateFunction < Tolerance)  {

        if (StressErrorMeasure < Tolerance) {
		// YA LO TENEMOS HECHO
            rReturnMappingVariables.TrialStateFunction = ElasticTrialStateFunction;
        }

        //Maybe elastic but error control did not converge
        else { 
	    //2.b Recalculate with an adaptive timeStep
            this->CalculateExplicitSolution(DeltaHenckyStrain, PreviousElasticHenckyStrain, PlasticVariables, NewElasticHenckyStrain, NewStressVector,  false, Tolerance); 
             
            rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, NewStressVector, PlasticVariables.EquivalentPlasticStrain);
            if (rReturnMappingVariables.TrialStateFunction > Tolerance) {
                  //2b.b From elastic to plastic
                  this->CalculateExplicitSolutionWithChange(DeltaHenckyStrain, PreviousElasticHenckyStrain, PlasticVariables, NewElasticHenckyStrain, NewStressVector, Tolerance);
		  PlasticityActive = true;
            }
        }

    }
    // Elastoplastic Step 
    else {
        if ( (rReturnMappingVariables.TrialStateFunction > Tolerance)  && (ElasticTrialStateFunction > Tolerance )) {
             std::cout << "Consider Drift Correction " << rReturnMappingVariables.TrialStateFunction << std::endl;
                this->CalculateExplicitSolution(DeltaHenckyStrain, PreviousElasticHenckyStrain, PlasticVariables, NewElasticHenckyStrain, NewStressVector, true,Tolerance);
		PlasticityActive = true;
        }	//It is clear that the drift correction must be considered 

        else if ((rReturnMappingVariables.TrialStateFunction < -Tolerance) && (ElasticTrialStateFunction > Tolerance)) {
                this->CalculateExplicitSolutionWithChange(DeltaHenckyStrain, PreviousElasticHenckyStrain, PlasticVariables, NewElasticHenckyStrain, NewStressVector, Tolerance);
		PlasticityActive = true;
        }
	else if (( fabs(rReturnMappingVariables.TrialStateFunction)<Tolerance ) && (ElasticTrialStateFunction > Tolerance)) {

	   double Cosinus = 0.0;
       double AuxNorma1 = 0.0;
	   double AuxNorma2 = 0.0;

       AuxiliarDerivativesStructure AuxiliarDerivatives;
       this->UpdateDerivatives(PreviousElasticHenckyStrain, AuxiliarDerivatives,PlasticVariables.EquivalentPlasticStrain);

	   for (unsigned int i = 0; i<6; i++) {
                 Cosinus   += AuxiliarDerivatives.YieldFunctionD(i)* (NewStressVector(i)-PreviousStressVector(i));
                 AuxNorma1 +=  pow((NewStressVector(i)-PreviousStressVector(i)), 2.0);
                 AuxNorma2 +=  pow(AuxiliarDerivatives.YieldFunctionD(i), 2.0);
           }
	   Cosinus /= sqrt(AuxNorma1*AuxNorma2);
           double AngleTolerance = 0.0;

           if (Cosinus > -AngleTolerance) {
                // Pure Plastic Increment
                std::cout << "INDESEABLE " << std::endl;
                this->CalculateExplicitSolution(DeltaHenckyStrain, PreviousElasticHenckyStrain, PlasticVariables, NewElasticHenckyStrain, NewStressVector, true, Tolerance);
		PlasticityActive = true;
           }
           else {
		//Elasto-Plastic Unloading etc.
                this->CalculateExplicitSolutionWithChange(DeltaHenckyStrain, PreviousElasticHenckyStrain, PlasticVariables, NewElasticHenckyStrain, NewStressVector, Tolerance);
		PlasticityActive = true;
                //Lo tengo solucionado con el tema de que F0 = -1
           }
        }
        else {
	   std::cout << "The Newcastle Group thinks that you can not enter here " << rReturnMappingVariables.TrialStateFunction << " " << ElasticTrialStateFunction << " " << " " << std::endl;
          //No es posible según Newcastle
        }
    }

   }
    rReturnMappingVariables.DeltaGamma = PlasticVariables.EquivalentPlasticStrain;
    rStressMatrix = MathUtils<double>::StressVectorToTensor(NewStressVector);

    rNewElasticLeftCauchyGreen = this->ConvertHenckyStrainToCauchyGreenTensor(NewElasticHenckyStrain);

    rReturnMappingVariables.Options.Set(PLASTIC_REGION,PlasticityActive);
    rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED, true);

    if (PlasticityActive) {
     double Aux;
    Aux = mpYieldCriterion->CalculateYieldCondition(Aux, NewStressVector, PlasticVariables.EquivalentPlasticStrain);
    std::cout << " DRIFT NEEDED TO BE CORRECTED " << Aux << std::endl;
  
    } 
   return PlasticityActive;

}

// ************** CONVERT HENCKY STRAIN TO CAUCHY-GREEN TENSOR  *********
// **********************************************************************

Vector NonAssociativeExplicitPlasticFlowRule::ConvertCauchyGreenTensorToHenckyStrain(const Matrix& rCauchyGreenMatrix)
{
    Matrix EigenVectors;
    Vector EigenValues;
    SolidMechanicsMathUtilities<double>::EigenVectors(rCauchyGreenMatrix, EigenVectors, EigenValues); 

    Matrix Aux = ZeroMatrix(3);
    
    for (unsigned int i = 0; i < 3; ++i)
        Aux(i,i) = (std::log(EigenValues(i)))/2.0;

    Aux = prod(Aux, trans(EigenVectors));
    Aux = prod(trans(EigenVectors), Aux);
   
    Vector Result= ZeroVector(6);

    Result = MathUtils<double>::StrainTensorToVector(Aux);
    return Result;
}
// ********************** AND THE CONTRARY ***************
// **99**88**77**66**
Matrix NonAssociativeExplicitPlasticFlowRule::ConvertHenckyStrainToCauchyGreenTensor(const Vector& rElasticHenckyStrain)
{

     Matrix HenckyStrainMatrix;
     HenckyStrainMatrix = MathUtils<double>::StrainVectorToTensor(rElasticHenckyStrain);

     Matrix EigenVectors;
     Vector EigenValues;
 
     SolidMechanicsMathUtilities<double>::EigenVectors(HenckyStrainMatrix, EigenVectors, EigenValues);

     Matrix Aux = ZeroMatrix(3);

    for (unsigned int i = 0; i < 3; ++i)
        Aux(i,i) = std::exp(2.0*EigenValues(i));

    Aux = prod(Aux, trans(EigenVectors));
    Aux = prod(EigenVectors, Aux);

    return Aux;

}


bool NonAssociativeExplicitPlasticFlowRule::UpdateInternalVariables( RadialReturnVariables & rReturnMappingVariables)
{
   mInternalVariables.EquivalentPlasticStrain = rReturnMappingVariables.DeltaGamma;
   return 0;
}


void NonAssociativeExplicitPlasticFlowRule::UpdateDerivatives(const Vector& rHenckyStrain, AuxiliarDerivativesStructure& rAuxiliarDerivatives, const double& rAlpha)
{
   Vector StressVector = ZeroVector(6);
   this->CalculateKirchhoffStressVector( rHenckyStrain, StressVector);
//   this->CalculatePlasticPotentialDerivatives(StressVector, rAuxiliarDerivatives.PlasticPotentialD, rAuxiliarDerivatives.PlasticPotentialDD);
   mpYieldCriterion->CalculateYieldFunctionDerivative(StressVector, rAuxiliarDerivatives.YieldFunctionD, rAlpha);
   rAuxiliarDerivatives.PlasticPotentialD = rAuxiliarDerivatives.YieldFunctionD;
}



//********************* UPDATES ONE STRAIN INCREMENT THE STRESS ****************
//******************************************************************************
void NonAssociativeExplicitPlasticFlowRule::CalculateOneExplicitStep(const Vector& rHenckyStrainIncrement, const Vector& rPreviousElasticHenckyStrain, InternalVariables& rPlasticVariables, Vector& rNewElasticHenckyStrain, Vector& rNewStressVector, double& rNewEquivalentPlasticStrain, const bool & rElastoPlasticBool, double& rStressErrorMeasure)
{
  
////// ATENCIÓ, AIXÒ ESTÀ FENT ABANS DE QUE RELLEGIS, cOM QUE ES TRACTA DE HYPERELASTIC s = s(epsilon);
// No fa falta calcular matrius tangents
// La pc també té equació explícita 
    
    // Runge Kutta Variables (Bogacki-Shampine)
    unsigned int nRK = 4;
    nRK = 4;
    Vector  CRK = ZeroVector(nRK);
    Matrix ARK = ZeroMatrix(nRK);
    Vector  b1 = ZeroVector(nRK);
    Vector  b2 = ZeroVector(nRK);

    if (nRK == 4) {

    CRK(0) = 0.0;
    CRK(1) = 0.5;
    CRK(2) = 0.75;
    CRK(3) = 1.0;

    ARK(0,0) = 0.0;
    ARK(1,0) = 1.0/2.0;
    ARK(2,0) = 0.0;     ARK(2,1) = 3.0/4.0;
    ARK(3,0) = 2.0/9.0; ARK(3,1) = 1.0/3.0; ARK(3,2) = 4.0/9.0;

    b1(0) = 2.0;
    b1(1) = 3.0;
    b1(2) = 4.0;
    b1(3) = 0.0;
    b1 /= 9.0;
    
    b2(0) = 7.0/24.0; b2(1) = 1.0/4.0; b2(2) = 1.0/3.0; b2(3) = 1.0/8.0;
}
else {
    nRK = 1;
    CRK(0) = 0.0; 
    ARK(0,0) = 0.0;
    b1(0) = 1.0; b2(0) = 1.0;
}
//    Vector  StressVector = ZeroVector(6);
    Vector  ElasticStrainVector = ZeroVector(6);

    //RK increments are saved in a matrix
    Matrix  ElasticStrainMatrix = ZeroMatrix(6, nRK);
    Vector  DeltaPlasticStrainVector = ZeroVector(nRK); 
    Matrix  ElasticMatrix = ZeroMatrix(6);

    Vector StepStrainIncrement = ZeroVector(6);
    Vector ElasticStrainIncrement = ZeroVector(6);
    double StepEquivalentPlasticStrain;
    double DeltaGamma;

    for (unsigned int iRK =0 ; iRK<nRK; ++iRK)
    {
        //Compute the stress and strain (at the Step iRK of RK) 
        ElasticStrainVector = rPreviousElasticHenckyStrain;
        StepEquivalentPlasticStrain = rPlasticVariables.EquivalentPlasticStrain;

        for (unsigned int i = 0; i<iRK; ++i){
	   StepEquivalentPlasticStrain += ARK(iRK,i)*DeltaPlasticStrainVector(i);
	   for (unsigned int j = 0; j < 6; ++j) {
              ElasticStrainVector(j) +=ARK(iRK,i)*ElasticStrainMatrix(j,i);
           }
        }

        //StepStrainIncrement = CRK(iRK)*rHenckyStrainIncrement;
        StepStrainIncrement = rHenckyStrainIncrement; // / 4.0;

        // Compute the ElastoPlastic Matrix, Elastic Incremental deformation and Hardening Parameters
    	this->ComputeElasticMatrix(ElasticStrainVector, ElasticMatrix);
    	if ( rElastoPlasticBool){
       	  AuxiliarDerivativesStructure AuxiliarDerivatives;
       	  this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, StepEquivalentPlasticStrain);

          double H;
          this->ComputePlasticHardeningParameter(ElasticStrainVector, StepEquivalentPlasticStrain, H);

          Vector auxVector;
          auxVector = prod(ElasticMatrix, StepStrainIncrement);
          DeltaGamma = MathUtils<double>::Dot( AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, StepStrainIncrement));
          double auxDenominador = H + MathUtils<double>::Dot( AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD));


	  DeltaGamma /= auxDenominador;
 
          if (DeltaGamma < 0.0)
	      DeltaGamma = 0.0;

          ElasticStrainIncrement = StepStrainIncrement - DeltaGamma*AuxiliarDerivatives.PlasticPotentialD;

        }
	else {
 	   ElasticStrainIncrement = StepStrainIncrement;
        }
	
        //Save the Increments
	for (unsigned int j = 0; j < 6; ++j) {
            ElasticStrainMatrix(j,iRK) = ElasticStrainIncrement(j);
        }
        if (rElastoPlasticBool) {
	       DeltaPlasticStrainVector(iRK) = 0.0;
	       for (unsigned int i = 0; i<3; ++i) {
		  //Update the Volumetric Plastic Strain Increment
		  DeltaPlasticStrainVector(iRK) += (StepStrainIncrement(i)-ElasticStrainIncrement(i) );
                }
        }
     }  //END FOR COMPUTE INCREMENTS

    ElasticStrainVector = ZeroVector(6); 
     
    rNewElasticHenckyStrain = ZeroVector(6); 

    for (unsigned int iRK = 0; iRK<nRK; ++iRK) {
	for (unsigned int j = 0; j < 6; ++j) {
	    ElasticStrainVector(j)     += b1(iRK)*ElasticStrainMatrix(j, iRK);
            rNewElasticHenckyStrain(j) += b2(iRK)*ElasticStrainMatrix(j, iRK);
         }
   }     

   rNewEquivalentPlasticStrain = rPlasticVariables.EquivalentPlasticStrain;

   if (rElastoPlasticBool) {
      for (unsigned int i = 0 ; i < 3; ++i) {
          //StepEquivalentPlasticStrain += (rHenckyStrainIncrement(i)-ElasticStrainIncrement(i));
          rNewEquivalentPlasticStrain += (rHenckyStrainIncrement(i)-rNewElasticHenckyStrain(i));
      }
   }
    ElasticStrainVector     += rPreviousElasticHenckyStrain;
    rNewElasticHenckyStrain += rPreviousElasticHenckyStrain;

   //COMPUTE AN ERROR MEASURE
   rStressErrorMeasure = 0.0;
   double denom = 0.0;
   for (unsigned int j = 0; j<6; ++j) {
       rStressErrorMeasure += pow(ElasticStrainVector(j) - rNewElasticHenckyStrain(j), 2.0 );
       denom += pow(rNewElasticHenckyStrain(j), 2.0);
   }
   rStressErrorMeasure /= (denom + 1e-8);
 
   this->CalculateKirchhoffStressVector(rNewElasticHenckyStrain, rNewStressVector); 

}


void NonAssociativeExplicitPlasticFlowRule::CalculateExplicitSolution(const Vector& rHenckyStrainIncrement, const Vector& rPreviousElasticHenckyStrain, InternalVariables& rPlasticVariables, Vector& rNewElasticHenckyStrain, Vector& rNewStressVector,  const bool& rElastoPlasticBool, const double & rTolerance)
{

   double TimeStep = 0.5;
   double MinTimeStep = 1e-5;
   double DoneTimeStep = 0.0;

   Vector DeltaHenckyStrain = ZeroVector(6);
   Vector PreviousElasticHenckyStrain = rPreviousElasticHenckyStrain;
   double NewEquivalentPlasticStrain;
   

   double StressErrorMeasure = 0.0;

   while (DoneTimeStep < 1.0)
   {
	if (DoneTimeStep + TimeStep >= 1.0) {
	   TimeStep = 1.0 - DoneTimeStep; 
       }
       DeltaHenckyStrain = TimeStep * rHenckyStrainIncrement;

       this->CalculateOneExplicitStep(DeltaHenckyStrain, PreviousElasticHenckyStrain, rPlasticVariables,  rNewElasticHenckyStrain, rNewStressVector, NewEquivalentPlasticStrain, rElastoPlasticBool, StressErrorMeasure);

       // Converge el step, (reasignar)
       if (StressErrorMeasure < rTolerance) {
 	  PreviousElasticHenckyStrain = rNewElasticHenckyStrain;
          DoneTimeStep += TimeStep;
          rPlasticVariables.EquivalentPlasticStrain = NewEquivalentPlasticStrain;

       }
       // El step no converge
       else {
	  if (TimeStep == MinTimeStep) {
	    //No converge el step

	     std::cout << "Explicit Stress Integration did not converge " << std::endl;
 	     PreviousElasticHenckyStrain = rNewElasticHenckyStrain;
             DoneTimeStep += TimeStep;
	     rPlasticVariables.EquivalentPlasticStrain = NewEquivalentPlasticStrain;
          }
       } 
      
       TimeStep *= 0.9* (pow( rTolerance / (StressErrorMeasure+ 1e-6), 1.0/2.0 ));
       TimeStep = std::max(TimeStep, MinTimeStep);

   }

}

// DEBE DE ARMONIARSE CON LOS CAMBIOS HECHOS PREVIAMENTE
void NonAssociativeExplicitPlasticFlowRule::CalculateExplicitSolutionWithChange(const Vector& rHenckyStrainIncrement, const Vector& rPreviousElasticHenckyStrain, InternalVariables& rPlasticVariables, Vector& rNewElasticHenckyStrain, Vector& rNewStressVector, const double& rTolerance)
{
    // To perform the bisection in order to obtain the 
    unsigned int NumberOfSubsteps = 10;
    double IncrementalAlpha = 1/double(NumberOfSubsteps);

    double AlphaEndOfElastic = 0.0;
    Vector StrainAtAlpha = rPreviousElasticHenckyStrain;


    double StateFunctionEndSubstep;     

    for (unsigned int iSubstep = 0; iSubstep<NumberOfSubsteps; ++iSubstep)
    {

	this->CalculateExplicitSolution(IncrementalAlpha*rHenckyStrainIncrement, StrainAtAlpha, rPlasticVariables,  rNewElasticHenckyStrain, rNewStressVector, false, rTolerance); 
        StateFunctionEndSubstep = mpYieldCriterion->CalculateYieldCondition(StateFunctionEndSubstep, rNewStressVector, rPlasticVariables.EquivalentPlasticStrain);

	//Elastic Substep, add and continue
        if (StateFunctionEndSubstep < rTolerance)
        {
	    AlphaEndOfElastic += IncrementalAlpha;
	    StrainAtAlpha = rNewElasticHenckyStrain;

	    if ( iSubstep == (NumberOfSubsteps -1) ) {
		return;
	    }
	}

	//ElastoPlastic Substep
	else {

   	  double StateFunction0 = -1.0; 
	  double StateFunctionHalf; 

	  double HalfSubstep = 0.5*(IncrementalAlpha);
	
	  bool Convergence = false;

	//1. Bisection, calculate Ee just at the yield surface
	  while (! Convergence) {
	
	     this->CalculateExplicitSolution( HalfSubstep*rHenckyStrainIncrement, StrainAtAlpha, rPlasticVariables, rNewElasticHenckyStrain, rNewStressVector, false, rTolerance); 
 
             StateFunctionHalf = mpYieldCriterion->CalculateYieldCondition(StateFunctionHalf, rNewStressVector, rPlasticVariables.EquivalentPlasticStrain);

	     if ( StateFunctionHalf < rTolerance)   {

	       AlphaEndOfElastic += HalfSubstep;
	       StrainAtAlpha = rNewElasticHenckyStrain;
	       StateFunction0 = StateFunctionHalf;
	       HalfSubstep *= 0.5;

             }
	     else {
               StateFunctionEndSubstep = StateFunctionHalf;
               HalfSubstep *= 0.5;
	     }

	     if ((HalfSubstep < 1E-5) || (fabs(StateFunction0) < rTolerance))  {
	        Convergence = true;
	        //We are at the yield surface, let's continue elastoplastically//
	     }
          } //End While Bisection 

StrainAtAlpha = rNewElasticHenckyStrain;
Vector Aux; this->CalculateKirchhoffStressVector(StrainAtAlpha, Aux);

	  this->CalculateExplicitSolution( (1-AlphaEndOfElastic)*rHenckyStrainIncrement, StrainAtAlpha, rPlasticVariables,  rNewElasticHenckyStrain, rNewStressVector,  true, rTolerance);
          return;
       }

    } //End for gran

    std::cout << "NOt supposed to be here " << std::endl;

}

   
void NonAssociativeExplicitPlasticFlowRule::ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rLeftCauchyGreenMatrix, const double& rAlpha, Matrix& rElasticMatrix)
{

   Vector ElasticStrainVector = ZeroVector(6);
   ElasticStrainVector = this->ConvertCauchyGreenTensorToHenckyStrain(rLeftCauchyGreenMatrix);
   this->ComputeElasticMatrix(ElasticStrainVector, rElasticMatrix);


   //if (rReturnMappingVariables.Control.PlasticRegion) {
   if (rReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION) ) {

      AuxiliarDerivativesStructure AuxiliarDerivatives;

      this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, rAlpha);

      double H;
      this->ComputePlasticHardeningParameter(ElasticStrainVector, rAlpha, H);
 
      Vector AuxVectorF;
      Vector AuxVectorG;

      AuxVectorF = prod(trans(AuxiliarDerivatives.YieldFunctionD), rElasticMatrix);
      AuxVectorG = prod( rElasticMatrix, AuxiliarDerivatives.PlasticPotentialD);

      rElasticMatrix -= 1.0*MyCrossProduct(trans(AuxVectorG), AuxVectorF) / ( H + MathUtils<double>::Dot(AuxVectorF, AuxiliarDerivatives.PlasticPotentialD));

  }

}


Matrix NonAssociativeExplicitPlasticFlowRule::MyCrossProduct(const Vector& rA, const Vector& rB)
{
    Matrix Result = ZeroMatrix(6);
    for (unsigned int i = 0; i<6; ++i) {
        for ( unsigned int j = 0; j<6; ++j) {
            Result(i,j) = rA(i)*rB(j);
        }
     }
     return Result;

}


} //End Namepace Kratos

//

#include <iostream>
#include<cmath>

#include "../PfemSolidMechanicsApplication/custom_constitutive/custom_flow_rules/non_associative_explicit_flow_rule.hpp"
#include "utilities/math_utils.h"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"
#include "includes/ublas_interface.h"
#include "solid_mechanics_application.h"
namespace Kratos
{


NonAssociativeExplicitPlasticFlowRule::NonAssociativeExplicitPlasticFlowRule()
  :FlowRule()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

NonAssociativeExplicitPlasticFlowRule::NonAssociativeExplicitPlasticFlowRule(YieldCriterionPointer pYieldCriterion)
	:FlowRule(pYieldCriterion)
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

bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{

   //std::cout << " COMPUTING STRESS, Initial rNewElasticLeftCauchyGreen" << rNewElasticLeftCauchyGreen << std::endl;
   bool PlasticityActive = false;
   bool Explicit  = true;

   if (Explicit) {
       PlasticityActive = CalculateReturnMappingExpl( rReturnMappingVariables, rIncrementalDeformationGradient, rStressMatrix, rNewElasticLeftCauchyGreen);
   }
   else {
       PlasticityActive = CalculateReturnMappingImpl( rReturnMappingVariables, rIncrementalDeformationGradient, rStressMatrix, rNewElasticLeftCauchyGreen);
   }

   return PlasticityActive;


}

/*bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMappingImpl2(RadialReturnVariables& rReturnMappingVariables, const Matrix& rDeltaDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{
    bool PlasticityActive = false;
    double Tolerance = 1e-5;

    InternalVariables PlasticVariables = mInternalVariables;
    rReturnMappingVariables.DeltaGamma = mInternalVariables.EquivalentPlasticStrain;
    mIncrementalPlasticShearStrain = 0.0;   
    rReturnMappingVariables.IncrementalPlasticShearStrain = 0.0;

    rReturnMappingVariables.Options.Set(PLASTIC_REGION, false);

    rNewElasticLeftCauchyGreen = prod( rNewElasticLeftCauchyGreen, trans( rDeltaDeformationGradient) );
    rNewElasticLeftCauchyGreen = prod( rDeltaDeformationGradient, rNewElasticLeftCauchyGreen);

    Vector StressVector;
    this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, StressVector);
  
    double ElasticTrialStateFunction;
    ElasticTrialStateFunction = mpYieldCriterion->CalculateYieldCondition( ElasticTrialStateFunction, StressVector, PlasticVariables.EquivalentPlasticStrain);

    if ( ElasticTrialStateFunction > Tolerance)  {
        PlasticityActive = true;
        rReturnMappingVariables.Options.Set(PLASTIC_REGION, true);
        this->ReturnStressToYieldSurface2( rReturnMappingVariables, rNewElasticLeftCauchyGreen, StressVector, ElasticTrialStateFunction, Tolerance);
    }

    rStressMatrix = MathUtils<double>::StressVectorToTensor(StressVector);

    rReturnMappingVariables.Options.Set(PLASTIC_REGION,PlasticityActive);
    rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED, true);

    return PlasticityActive;

}
*/ // TAKE ME OUT

// THIS FUNCTION IS EQUIVALENT OF THE NEXT ONE  BUT IMPLICIT (RETURN MAPPING)
// CORRESPONTS TO THE CUTTING PLANE ALGORITHM (CPA AND NOT THE CPPM)
// TO BE PUT IN A DERIVED CLASS.
bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMappingImpl(RadialReturnVariables& rReturnMappingVariables, const Matrix& rDeltaDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{
    bool PlasticityActive = false;
    double Tolerance = 1e-5;

    InternalVariables PlasticVariables = mInternalVariables;
    rReturnMappingVariables.DeltaGamma = mInternalVariables.EquivalentPlasticStrain;
    mIncrementalPlasticShearStrain = 0.0;   
    rReturnMappingVariables.IncrementalPlasticShearStrain = 0.0;

    rReturnMappingVariables.Options.Set(PLASTIC_REGION, false);

    rNewElasticLeftCauchyGreen = prod( rNewElasticLeftCauchyGreen, trans( rDeltaDeformationGradient) );
    rNewElasticLeftCauchyGreen = prod( rDeltaDeformationGradient, rNewElasticLeftCauchyGreen);

    Vector StressVector;
    this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, StressVector);
  
    double ElasticTrialStateFunction;
    ElasticTrialStateFunction = mpYieldCriterion->CalculateYieldCondition( ElasticTrialStateFunction, StressVector, PlasticVariables.EquivalentPlasticStrain);

    if ( ElasticTrialStateFunction > Tolerance)  {
        PlasticityActive = true;
        rReturnMappingVariables.Options.Set(PLASTIC_REGION, true);
        this->ReturnStressToYieldSurface( rReturnMappingVariables, rNewElasticLeftCauchyGreen, StressVector, ElasticTrialStateFunction, Tolerance);
    }

    rStressMatrix = MathUtils<double>::StressVectorToTensor(StressVector);

    rReturnMappingVariables.Options.Set(PLASTIC_REGION,PlasticityActive);
    rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED, true);

    return PlasticityActive;

}

//bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables, const Matrix& rDeformationGradientF0, const Matrix& rDeltaDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMappingExpl(RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{

   //1- Initialize some variables
   bool PlasticityActive = false;
   rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);

   InternalVariables PlasticVariables = mInternalVariables;
  
   // IncrementalPlasticStrain and Total VolumetricStrain 
   rReturnMappingVariables.IncrementalPlasticShearStrain = 0.0;
   mIncrementalPlasticShearStrain = 0.0;   
   rReturnMappingVariables.DeltaGamma = mInternalVariables.EquivalentPlasticStrain; 

   Matrix PreviousElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
   Vector PreviousStressVector;

   Vector NewStressVector;
   
   double Tolerance = 1E-5;


   //2- Check for yield Condition (at the beggining)
   this->CalculateKirchhoffStressVector( PreviousElasticLeftCauchyGreen, PreviousStressVector);
   rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, PreviousStressVector, rReturnMappingVariables.DeltaGamma );

    if (rReturnMappingVariables.TrialStateFunction > Tolerance) {
//            this->ReturnStressToYieldSurface( rReturnMappingVariables, PreviousElasticLeftCauchyGreen, PreviousStressVector, rReturnMappingVariables.DeltaGamma, rReturnMappingVariables.TrialStateFunction, Tolerance);
           this->ReturnStressToYieldSurface( rReturnMappingVariables, PreviousElasticLeftCauchyGreen, PreviousStressVector, rReturnMappingVariables.TrialStateFunction, Tolerance);
            //THE INITIAL TENSIONAL STATE MAY BE OUTSIDE THE YIELD SURFACE??
    }

    //3. Compute an Elastic Trial State
    ExplicitStressUpdateInformation  StressUpdateInformation;

    this->CalculateOneExplicitStep(rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, false, StressUpdateInformation);

    double ElasticTrialStateFunction;
    ElasticTrialStateFunction = mpYieldCriterion->CalculateYieldCondition( ElasticTrialStateFunction, NewStressVector, PlasticVariables.EquivalentPlasticStrain);


    //4a. COMPLETELY ELASTIC STEP
    if ( ElasticTrialStateFunction <= Tolerance)   {
        PlasticityActive = false;
    }

    //4b. SOME PART OF THE STEP IS PLASTIC
    else {

         PlasticityActive = true; 

         if ( ( rReturnMappingVariables.TrialStateFunction < -Tolerance) && (ElasticTrialStateFunction > Tolerance) ) {
            //2b.1 Transition from elastic to plastic
             //std::cout << " ELASTIC AND PLASTIC STEP " << std::endl;
            this->CalculateExplicitSolutionWithChange( rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, Tolerance);
          
         }
         else {         

            bool UnloadingCondition;
            UnloadingCondition =  this->EvaluateElastoPlasticUnloadingCondition( UnloadingCondition, PreviousElasticLeftCauchyGreen, rIncrementalDeformationGradient, PlasticVariables, Tolerance); 

            if (UnloadingCondition ) {
                 //2c.1 ElastoPlastic Unloading
	         //std::cout << " ELASTIC UNLOADING STEP " << std::endl;
                 this->CalculateExplicitSolutionWithChange( rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, Tolerance);

            }
            else {
	         //std::cout << " PLASTIC STEP " << std::endl;
                //2c. 2 Completedly ELastoPlastic Step 
                this->CalculateExplicitSolution( rIncrementalDeformationGradient, PreviousElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, true, Tolerance);

            }
         }
 

    }


    // WHAT??    rReturnMappingVariables.NormIsochoricStress = PlasticVariables.EquivalentPlasticStrain;
    if (PlasticityActive) {

       double DriftViolation;
     
       DriftViolation = mpYieldCriterion->CalculateYieldCondition(DriftViolation, NewStressVector, PlasticVariables.EquivalentPlasticStrain);
       //rReturnMappingVariables.EigenVectors = NewElasticLeftCauchyGreen; 

       if ( fabs(DriftViolation) > Tolerance ) {

            this->ReturnStressToYieldSurface( rReturnMappingVariables, rNewElasticLeftCauchyGreen, NewStressVector, DriftViolation, Tolerance);
       }


    }

    rStressMatrix = MathUtils<double>::StressVectorToTensor(NewStressVector);

    rReturnMappingVariables.Options.Set(PLASTIC_REGION,PlasticityActive);
    rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED, true);

   //rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, NewStressVector, PlasticVariables.EquivalentPlasticStrain);
    //rNewElasticLeftCauchyGreen = NewElasticLeftCauchyGreen; //WHY?? 

   return PlasticityActive;

}


bool& NonAssociativeExplicitPlasticFlowRule::EvaluateElastoPlasticUnloadingCondition( bool& rUnloadingCondition, const Matrix& rElasticLeftCauchyGreen, const Matrix& rDeltaDeformationGradient, const InternalVariables& rPlasticVariables, const double& rTolerance)
{

     rUnloadingCondition = false;
    
     Vector StressVector;
     Vector YieldFunctionD;
     this->CalculateKirchhoffStressVector( rElasticLeftCauchyGreen, StressVector);
     mpYieldCriterion->CalculateYieldFunctionDerivative(StressVector, YieldFunctionD, rPlasticVariables.EquivalentPlasticStrain);

     Matrix ElasticMatrix;
     Vector ElasticStrainVector ;
     Vector IncrementalHencky;

     ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rElasticLeftCauchyGreen);
     IncrementalHencky = ConvertCauchyGreenTensorToHenckyStrain( prod(rDeltaDeformationGradient, trans(rDeltaDeformationGradient) ));
     this->ComputeElasticMatrix(ElasticStrainVector, ElasticMatrix);

     Vector IncrementalElasticStress;
     IncrementalElasticStress = prod(ElasticMatrix, IncrementalHencky);

     double Numerador = 0.0;
     double Denom1 = 0.0;
     double Denom2 = 0.0;

     for (unsigned int i = 0; i < 6; ++i) {
         Numerador += YieldFunctionD(i)*IncrementalElasticStress(i);
         Denom1 += pow( YieldFunctionD(i), 2.0);
         Denom2 += pow( IncrementalElasticStress(i), 2.0);
      }

      Denom1 = pow( Denom1*Denom2, 0.5);
    
      Numerador /= Denom1;

      if (Numerador < rTolerance) 
          rUnloadingCondition = true;

       return rUnloadingCondition;

}



// *************** YIELD SURFACE VIOLATION DRIFT CORRECTION ********** TAKE ME OUT
// *******************************************************************
//void NonAssociativeExplicitPlasticFlowRule::ReturnStressToYieldSurface2( RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rStressVector, double& rDrift, const double& rTolerance)
//{
/*     Vector TrialElasticHenckyStrain = this->ConvertCauchyGreenTensorToHenckyStrain( rNewElasticLeftCauchyGreen);
     Vector ActualElasticHenckyStrain = TrialElasticHenckyStrain;
     Vector PlasticDeformation = ZeroVector(6);
     double TrialPlasticHardening = rReturnMappingVariables.EquivalentPlasticStrain;
     double ActualPlasticHardening = TrialPlasticHardening;
     double DeltaGamma = 0.0;

     Vector StessVector;
     Vector Residual = ZeroVector(8);
     AuxiliarDerivativesStructure  AuxiliarDerivatives;

     Matrix = ElasticMatrix;

     for (unsigned int i = 0; i < 150; ++i)
   
        //1. Compute The Residual
        PlasticDeformation = DeltaGamma * AuxiliarDerivatives.PlasticPotentialD ;
        ActualElasticHenckyStrain = TrialElasticHenckyStrain - PlasticDeformation;

        this->UpdateDerivatives(ActualElasticHenckyStrain, AuxiliarDerivatives, ActualPlasticHardening);

        this->CalculateKirchhoffStressVector( CorrectedLeftCauchyGreen, rStressVector);
        rDrift = mpYieldCriterion->CalculateYieldCondition( rDrift, rStressVector, ActualPlasticHardening);
        
        for (unsigned j = 0, j < 6; ++j)
            Residual(j) = -PlasticDeformation(j) + DeltaGamma * AuxiliarDerivatives.PlasticPotentialD(j);  

        Residual(6) = -ActualPlasticHardening + TrialPlasticHardening;
        for (unsigned j = 0 ; j < 3; ++j)
            Residual(6) += DeltaGamma * AuxiliarDerivatives.PlasticPotentialD(j);

        for (unsigned j = 0; j < 6; ++j)
            Residual(j) = -ActualElasticHenckyStrain(j) + TrialElasticHenckyStrain(j) - DeltaGamma * AuxiliarDerivatives.PlasticPotentialD(j);

        Residual(6) = -ActualPlasticHardening + TrialPlasticHardening;
        for (unsigned j = 0; j < 3; ++j)
             Residual(6) += DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(j);

        this->ComputeElasticMatrix(ActualElasticHenckyStrain, ElasticMatrix);

*/
//}



            
void NonAssociativeExplicitPlasticFlowRule::ReturnStressToYieldSurface( RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rStressVector, double& rDrift, const double& rTolerance)
{

/*     Vector ActualElasticHenckyStrain;
     ActualElasticHenckyStrain = this->ConvertCauchyGreenTensorToHenckyStrain( rNewElasticLeftCauchyGreen);

     this->ReturnStressToYieldSurface( ActualElasticHenckyStrain, rStressVector, rAlpha, rDrift, rTolerance);

     rNewElasticLeftCauchyGreen = this->ConvertHenckyStrainToCauchyGreenTensor( ActualElasticHenckyStrain);
 */    // HAY ALGO MAL PORQUE DEJA DE CONVERGER LO BIEN QUE LO HACIA 
 
    AuxiliarDerivativesStructure  AuxiliarDerivatives;
    Matrix ElasticMatrix;
    double H;
    double DeltaGamma;
    Vector ElasticCorrection;
    Matrix  CorrectedLeftCauchyGreen;
    Matrix  UpdateMatrix;
    Vector ActualElasticHenckyStrain;
    Vector StressVector;
    double Alpha = rReturnMappingVariables.DeltaGamma; 
    rReturnMappingVariables.DeltaBeta = 0.0;

    for (unsigned int i = 0; i < 150; ++i) {

        ActualElasticHenckyStrain = this->ConvertCauchyGreenTensorToHenckyStrain( rNewElasticLeftCauchyGreen);

        this->UpdateDerivatives(ActualElasticHenckyStrain, AuxiliarDerivatives, Alpha);
        this->ComputeElasticMatrix(ActualElasticHenckyStrain, ElasticMatrix);
        this->ComputePlasticHardeningParameter(ActualElasticHenckyStrain, Alpha, H);

        DeltaGamma = rDrift;
        DeltaGamma /= ( H + MathUtils<double>::Dot(AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD)));

        Vector AuxF = prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD);

        double denominador = 0.0;
        for (unsigned int k = 0; k < 6; ++k)
             denominador += AuxiliarDerivatives.YieldFunctionD(k) * AuxF(k);

        denominador += H;
        Vector MuyAuxiliar = -DeltaGamma*AuxiliarDerivatives.PlasticPotentialD/ 2.0;
        UpdateMatrix = this->ConvertHenckyStrainToCauchyGreenTensor( MuyAuxiliar);

 
        CorrectedLeftCauchyGreen =  prod((UpdateMatrix), rNewElasticLeftCauchyGreen);
        CorrectedLeftCauchyGreen =  prod( CorrectedLeftCauchyGreen, trans(UpdateMatrix));
    
        double AlphaUpdate = 0.0;
        for (unsigned int j = 0; j < 3; ++j)
            AlphaUpdate += DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(j);

        for (unsigned int j = 0; j < 3; ++j)
            rReturnMappingVariables.IncrementalPlasticShearStrain += pow(DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(j) - AlphaUpdate/3.0, 2.0);

        for (unsigned int j = 3; j < 6; ++j)
            rReturnMappingVariables.IncrementalPlasticShearStrain += pow(DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(j), 2.0);


        Alpha += AlphaUpdate;
        rReturnMappingVariables.DeltaBeta += DeltaGamma ;
        this->CalculateKirchhoffStressVector( CorrectedLeftCauchyGreen, rStressVector);

        rDrift = mpYieldCriterion->CalculateYieldCondition( rDrift, rStressVector, Alpha);

        rNewElasticLeftCauchyGreen = CorrectedLeftCauchyGreen;
        rReturnMappingVariables.DeltaGamma = Alpha;
//        this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, StressVector);

//        std::cout << i << " Drift " <<  rDrift << " rAlpha " << rAlpha <<  "DELTA " << DeltaGamma  << std::endl;
//        std::cout << rNewElasticLeftCauchyGreen << std::endl;
//        std::cout << " stress " << StressVector << std::endl;
//        std::cout << " AlphPrev " << rAlpha-AlphaUpdate << " UPDATE " << AlphaUpdate << std::endl;
//        std::cout << " h  " << H <<  " " << MathUtils<double>::Dot(AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD)) << std::endl;
//        std::cout << " " << std::endl;
        if( fabs(rDrift) < 10.0*rTolerance) {
//            std::cout << " DRIFT LEAVING" << rDrift << " " << rTolerance << std::endl;
            return;
        }

    }

    if ( fabs(rDrift) > 500.0*rTolerance) {
       std::cout<< " " << std::endl;
       std::cout << "Leaving drift correction WITHOUT converging " << rDrift << std::endl;
       std::cout << " StressVector " << rStressVector << std::endl;
}


}


void NonAssociativeExplicitPlasticFlowRule::ReturnStressToYieldSurface(Vector& rElasticHenckyStrainVector, Vector& rStressVector, double& rAlpha, double& rDrift, const double& rTolerance)
{
    AuxiliarDerivativesStructure  AuxiliarDerivatives;
    Matrix ElasticMatrix;
    double H;
    double DeltaGamma;
    Vector ElasticCorrection;
        

    for (unsigned int i = 0; i < 15; ++i)     {


        this->UpdateDerivatives(rElasticHenckyStrainVector, AuxiliarDerivatives, rAlpha);
        this->ComputeElasticMatrix(rElasticHenckyStrainVector, ElasticMatrix);
        this->ComputePlasticHardeningParameter(rElasticHenckyStrainVector, rAlpha, H);

        DeltaGamma = rDrift;
        DeltaGamma /= ( H + MathUtils<double>::Dot(AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD)));

        ElasticCorrection = -DeltaGamma*AuxiliarDerivatives.PlasticPotentialD;

        rElasticHenckyStrainVector += ElasticCorrection;

        for (unsigned int j = 0; j<3; ++j)
           rAlpha -= ElasticCorrection(j);

        this->CalculateKirchhoffStressVector( rElasticHenckyStrainVector, rStressVector);

        rDrift = mpYieldCriterion->CalculateYieldCondition( rDrift, rStressVector, rAlpha);
    
        if (fabs(rDrift) < rTolerance)  {
            return;
        }
     }

     if (fabs(rDrift) > rTolerance)
         std::cout << "Leaving Drift Correction whitout converging! " << rDrift << std::endl;


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

    Aux = prod(Aux, (EigenVectors));
    Aux = prod(trans(EigenVectors), Aux);
    Vector Result; 
    Result = MathUtils<double>::StrainTensorToVector(Aux, 6);
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

    Aux = prod(Aux, (EigenVectors));
    Aux = prod(trans(EigenVectors), Aux);

    // COMPROVACION
/*    Matrix A = ZeroMatrix(3);
    for (unsigned int i = 0; i < 3; ++ i ) 
        A(i,i) = EigenValues(i);

    A = prod(trans(EigenVectors), A);
    A = prod(A, EigenVectors);

    A = A - HenckyStrainMatrix; 
*/
    return Aux;

}


bool NonAssociativeExplicitPlasticFlowRule::UpdateInternalVariables( RadialReturnVariables & rReturnMappingVariables)
{
   mInternalVariables.EquivalentPlasticStrain = rReturnMappingVariables.DeltaGamma;
   //mInternalVariables.EquivalentPlasticStrainOld = rReturnMappingVariables.LameMu ; 
   mIncrementalPlasticShearStrain = rReturnMappingVariables.IncrementalPlasticShearStrain;   
   mInternalVariables.EquivalentPlasticStrainOld = rReturnMappingVariables.IncrementalPlasticShearStrain;
   mInternalVariables.DeltaPlasticStrain += rReturnMappingVariables.IncrementalPlasticShearStrain; 
   return 0;
}


void NonAssociativeExplicitPlasticFlowRule::UpdateDerivatives(const Vector& rHenckyStrain, AuxiliarDerivativesStructure& rAuxiliarDerivatives, const double& rAlpha)
{
   Vector StressVector = ZeroVector(6);
   this->CalculateKirchhoffStressVector( rHenckyStrain, StressVector);
   this->CalculatePlasticPotentialDerivatives(StressVector, rAuxiliarDerivatives.PlasticPotentialD, rAuxiliarDerivatives.PlasticPotentialDD);
   mpYieldCriterion->CalculateYieldFunctionDerivative(StressVector, rAuxiliarDerivatives.YieldFunctionD, rAlpha);
   rAuxiliarDerivatives.PlasticPotentialD = rAuxiliarDerivatives.YieldFunctionD;
}

//****** COMPUTE INCREMENTAL DEFORMATION GRADIENT  ********* 
// FROM BATHE NOTATION ^{rFinalConfiguration} _{rReferenceConfiguration} X  from  ^{t+1}_{t} X
//*********************************************************
void NonAssociativeExplicitPlasticFlowRule::ComputeSubstepIncrementalDeformationGradient( const Matrix& rIncrementalDeformationGradient, const double& rReferenceConfiguration, const double& rFinalConfiguration, Matrix& rSubstepIncrementalDeformationGradient)
{

    
    Matrix DeformationGradientReference;
    Matrix DeformationGradientFinal;

    Matrix IdentityMatrix = ZeroMatrix(3);
    for (unsigned int i = 0; i < 3; ++i)
         IdentityMatrix(i,i) = 1.0;   
 
    DeformationGradientReference = rReferenceConfiguration*rIncrementalDeformationGradient  + (1.0 - rReferenceConfiguration)*IdentityMatrix;
    DeformationGradientFinal     =     rFinalConfiguration*rIncrementalDeformationGradient  + (1.0 -     rFinalConfiguration)*IdentityMatrix;

    double Det; 
    MathUtils<double>::InvertMatrix(DeformationGradientReference, rSubstepIncrementalDeformationGradient, Det);
    rSubstepIncrementalDeformationGradient = prod( DeformationGradientFinal, rSubstepIncrementalDeformationGradient);
   
}




//********************* UPDATES ONE ELASTOPLASTIC STRAIN INCREMENT THE STRESS ****************
//********************************************************************************************

void NonAssociativeExplicitPlasticFlowRule::CalculateOneExplicitPlasticStep(const Matrix& rDeltaDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen, const double& rPreviousEquivalentPlasticStrain, Matrix& rNewElasticLeftCauchyGreen, double& rNewEquivalentPlasticStrain, double& rNewPlasticShearStrain)
{
        rNewEquivalentPlasticStrain = rPreviousEquivalentPlasticStrain;
        Matrix ElasticMatrix;
        Vector ElasticStrainVector ;

        ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rPreviousElasticLeftCauchyGreen);
        this->ComputeElasticMatrix(ElasticStrainVector, ElasticMatrix);

       	AuxiliarDerivativesStructure AuxiliarDerivatives;
       	this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, rNewEquivalentPlasticStrain);

        double H;
        this->ComputePlasticHardeningParameter(ElasticStrainVector, rNewEquivalentPlasticStrain, H);


        ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( prod(rDeltaDeformationGradient, trans(rDeltaDeformationGradient)) );


        Vector auxVector;
        auxVector = prod(ElasticMatrix, ElasticStrainVector) ;
        double DeltaGamma = 0.0;
        for (unsigned int i = 0; i<6; ++i)
             DeltaGamma += auxVector(i)*AuxiliarDerivatives.YieldFunctionD(i);

        double auxDenominador = H + MathUtils<double>::Dot( AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD));

	DeltaGamma /= auxDenominador;

        if (DeltaGamma < 0 )
                DeltaGamma =0;

        Vector MuyAuxiliar;
        MuyAuxiliar = -DeltaGamma * AuxiliarDerivatives.PlasticPotentialD / 2.0; 

        Matrix UpdateMatrix;
        UpdateMatrix = this->ConvertHenckyStrainToCauchyGreenTensor (MuyAuxiliar);
        UpdateMatrix = prod( rDeltaDeformationGradient, UpdateMatrix);

        rNewElasticLeftCauchyGreen = prod(UpdateMatrix, rPreviousElasticLeftCauchyGreen);
        rNewElasticLeftCauchyGreen = prod( rNewElasticLeftCauchyGreen, trans(UpdateMatrix));

        double PlasticPotentialP = 0.0;
        for (unsigned int i = 0; i < 3; ++i)
             PlasticPotentialP += DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(i);

        rNewEquivalentPlasticStrain +=  PlasticPotentialP ;// 3.0; 

        rNewPlasticShearStrain = 0.0;
        for (unsigned int i = 0; i < 3; ++i)
            rNewPlasticShearStrain += pow( DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(i) - rNewEquivalentPlasticStrain, 2.0);

        for (unsigned int i = 3; i < 6; ++i)
            rNewPlasticShearStrain += pow( DeltaGamma*AuxiliarDerivatives.PlasticPotentialD(i), 2.0);

        rNewPlasticShearStrain = pow( rNewPlasticShearStrain, 1.0/2.0) ;
       

}

//********************* UPDATES ONE STRAIN INCREMENT THE STRESS ****************
//******************************************************************************
void NonAssociativeExplicitPlasticFlowRule::CalculateOneExplicitStep(const Matrix& rDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen, const RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector, const bool& rElastoPlasticBool, ExplicitStressUpdateInformation& rStressUpdateInformation)
{
        
    if ( rElastoPlasticBool)  {

        double NewEquivalentPlasticStrain;
        double NewPlasticShearStrain;
        Vector FirstApproxStressVector;
        // ONE STEP UPDATE
        this->CalculateOneExplicitPlasticStep( rDeformationGradient, rPreviousElasticLeftCauchyGreen, rReturnMappingVariables.DeltaGamma, rNewElasticLeftCauchyGreen, NewEquivalentPlasticStrain, NewPlasticShearStrain);

        this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, FirstApproxStressVector); 

        //TWO STEP UPDATE

        Matrix HalfStepIncrementalDeformation;
        ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, 0.0, 0.5, HalfStepIncrementalDeformation);
        
        Matrix SecondApproxLeftCauchyGreen; 

        this->CalculateOneExplicitPlasticStep( HalfStepIncrementalDeformation, rPreviousElasticLeftCauchyGreen, rReturnMappingVariables.DeltaGamma, SecondApproxLeftCauchyGreen, NewEquivalentPlasticStrain, NewPlasticShearStrain);

        ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, 0.5, 1.0, HalfStepIncrementalDeformation);
        rStressUpdateInformation.IncrementPlasticStrain = NewPlasticShearStrain; 

        this->CalculateOneExplicitPlasticStep( HalfStepIncrementalDeformation, SecondApproxLeftCauchyGreen, NewEquivalentPlasticStrain, rNewElasticLeftCauchyGreen, NewEquivalentPlasticStrain, NewPlasticShearStrain);

        rStressUpdateInformation.IncrementPlasticStrain += NewPlasticShearStrain; 
        rStressUpdateInformation.NewEquivalentPlasticStrain = NewEquivalentPlasticStrain;
         this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, rNewStressVector);

         rStressUpdateInformation.StressErrorMeasure = 0.0;
         double Denominador = 0.0;
         for (unsigned int i = 0; i<6; ++i) {
            rStressUpdateInformation.StressErrorMeasure += pow( rNewStressVector(i) - FirstApproxStressVector(i), 2.0);
            Denominador  += pow(rNewStressVector(i), 2.0);
         }
          rStressUpdateInformation.StressErrorMeasure = 1.0* pow( rStressUpdateInformation.StressErrorMeasure/Denominador, 0.5);

         
   } 
   else  {
    
      rNewElasticLeftCauchyGreen = prod( rPreviousElasticLeftCauchyGreen, trans( rDeformationGradient) );
      rNewElasticLeftCauchyGreen = prod( rDeformationGradient, rNewElasticLeftCauchyGreen );
  
      rStressUpdateInformation.StressErrorMeasure = 0.0;
//      rNewPlasticShearStrain = 0.0;
//      rNewEquivalentPlasticStrain = rPlasticVariables.EquivalentPlasticStrain;   
 
      this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, rNewStressVector);
     
   }


}


//***************** CalculateStress From LEFT CAUCHY GREEN MATRIX **************
void NonAssociativeExplicitPlasticFlowRule::CalculateKirchhoffStressVector( const Matrix& rElasticLeftCauchyGreen, Vector& rNewStressVector)
{

    Vector ElasticHenckyStrainVector;
    ElasticHenckyStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rElasticLeftCauchyGreen); 
    this->CalculateKirchhoffStressVector(ElasticHenckyStrainVector, rNewStressVector);

}


           
void NonAssociativeExplicitPlasticFlowRule::UpdateRadialReturnVariables( RadialReturnVariables& rReturnMappingVariables, const ExplicitStressUpdateInformation& rStressUpdateInformation)
{
     rReturnMappingVariables.DeltaGamma = rStressUpdateInformation.NewEquivalentPlasticStrain;
     rReturnMappingVariables.IncrementalPlasticShearStrain += rStressUpdateInformation.IncrementPlasticStrain;
}



void NonAssociativeExplicitPlasticFlowRule::CalculateExplicitSolution( const Matrix& rIncrementalDeformationGradient, const Matrix& rPreviousElasticCauchyGreen, RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector, const bool& rElastoPlasticBool, const double& rTolerance)
{
   double TimeStep = 0.5;
   double MinTimeStep = 1.0e-3;
   double DoneTimeStep = 0.0;
   double MaxTimeStep = 0.5;
   TimeStep = MaxTimeStep;

   //Matrix ActualDeformationGradient = rDeformationGradientF0;
   Matrix ActualElasticLeftCauchyGreen = rPreviousElasticCauchyGreen;
   Matrix SubstepDeformationGradient; 

   //double NewEquivalentPlasticStrain;
   //double NewPlasticShearStrain = 0.0;

   ExplicitStressUpdateInformation  StressUpdateInformation;
   bool MayBeLast = false;
   while (DoneTimeStep < 1.0)  {
      //std::cout << " TIME " << DoneTimeStep << " " << TimeStep  << std::endl;
       if (DoneTimeStep + TimeStep >= 1.0) {
           TimeStep = 1.0 - DoneTimeStep;
           MayBeLast = true;
       }

       this->ComputeSubstepIncrementalDeformationGradient( rIncrementalDeformationGradient, DoneTimeStep, DoneTimeStep + TimeStep, SubstepDeformationGradient);

       this->CalculateOneExplicitStep( SubstepDeformationGradient, ActualElasticLeftCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, rNewStressVector,  rElastoPlasticBool, StressUpdateInformation);


       if ( StressUpdateInformation.StressErrorMeasure < 10.0*rTolerance ) {
           // Se acepta el paso
           //std::cout << "   SUbstepping " << DoneTimeStep << " dt " << TimeStep << "Error " << StressErrorMeasure << std::endl;
           ActualElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;

           this->UpdateRadialReturnVariables( rReturnMappingVariables, StressUpdateInformation);
           DoneTimeStep += TimeStep;
           if (MayBeLast == true)
               return;

       }
       else {
          if (TimeStep <= MinTimeStep) {
              if ( StressUpdateInformation.StressErrorMeasure > rTolerance) {
                 std::cout << "ExplicitStressIntegrationDidNotConverged" << StressUpdateInformation.StressErrorMeasure << " " << TimeStep << std::endl;
                 std::cout << "    stress " << rNewStressVector << std::endl;
               /*  std::cout << "    PreviousEGCG" << rPreviousElasticLeftCauchyGreen << std::endl;
                 std::cout << "    Actual  EGCG" <<    ActualElasticLeftCauchyGreen << std::endl;
                 std::cout << "    New     EGCG" << rNewElasticLeftCauchyGreen << std::endl;
                 std::cout << "    rAlpha      " << NewEquivalentPlasticStrain << std::endl;
                 double Hola;
                 Hola = mpYieldCriterion->CalculateYieldCondition(Hola, rNewStressVector, NewEquivalentPlasticStrain);
                 std::cout << "    rYieldValue " << Hola << std::endl;*/
              }

              DoneTimeStep += TimeStep;

              ActualElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
              this->UpdateRadialReturnVariables( rReturnMappingVariables, StressUpdateInformation);
              //rPlasticVariables.EquivalentPlasticStrain = NewEquivalentPlasticStrain;
              //rPlasticVariables.EquivalentPlasticStrainOld += NewPlasticShearStrain;
              if (MayBeLast == true )
                    return;
          }
 
 
      }

      TimeStep *= 0.9* (pow( 10.0*rTolerance / (StressUpdateInformation.StressErrorMeasure+ 1e-8), 1.0/2.0 ));
      TimeStep = std::max(TimeStep, MinTimeStep);
      TimeStep = std::min(TimeStep, MaxTimeStep);
  }

}





// Search for the change and then continue
void NonAssociativeExplicitPlasticFlowRule::CalculateExplicitSolutionWithChange(const Matrix& rDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen, RadialReturnVariables& rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector,  const double& rTolerance)
{

   //Matrix FinalDeformationGradient = prod( rDeltaDeformationGradient, rDeformationGradientF0);

// bisecction in elastic Regime

   double InitialPosition = 0.0;
   double EndPosition = 1.0;
   double HalfPosition;

   double InitialStateFunction = -1.0;
   double EndStateFunction = 1.0;
   double HalfPositionStateFunction;

   Matrix HalfPositionDeformationGradient;
   Matrix HalfPositionElasticCauchyGreen;
    
   ExplicitStressUpdateInformation StressUpdateInformation;

   //double NewEquivalentPlasticStrain = rReturnMappingVariables.train ;
   //double NewPlasticShearStrain  = 0.0;
   double ErrorMeasure1;
   double ErrorMeasure2;



   for (unsigned int i = 0 ; i< 150; ++i) {

        HalfPosition = 0.5*(InitialPosition + EndPosition);
  
        this->ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, 0.0, HalfPosition, HalfPositionDeformationGradient);

        this->CalculateOneExplicitStep( HalfPositionDeformationGradient, rPreviousElasticLeftCauchyGreen, rReturnMappingVariables, HalfPositionElasticCauchyGreen, rNewStressVector, false, StressUpdateInformation);  

        HalfPositionStateFunction = mpYieldCriterion->CalculateYieldCondition( HalfPositionStateFunction, rNewStressVector, rReturnMappingVariables.DeltaGamma);
        

       if ( HalfPositionStateFunction < 0.0)  {
           InitialStateFunction = HalfPositionStateFunction;
           InitialPosition = HalfPosition;
       }
       else {
           EndPosition = HalfPosition;
           EndStateFunction = HalfPositionStateFunction;
       }
       ErrorMeasure1 = fabs(InitialStateFunction-EndStateFunction);
       ErrorMeasure2 = fabs(InitialPosition - EndPosition);

       if ( (ErrorMeasure1 < rTolerance ) && (ErrorMeasure2 < rTolerance)) {
           break;
       }
 
   }
   //std::cout << " Half Position Found " << HalfPosition << " SF " << HalfStateFunction << std::endl;

   // COMPUTE ELASTIC STEP;

   this->ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, 0.0, HalfPosition, HalfPositionDeformationGradient);
  
   this->CalculateOneExplicitStep( HalfPositionDeformationGradient, rPreviousElasticLeftCauchyGreen, rReturnMappingVariables, HalfPositionElasticCauchyGreen, rNewStressVector, false, StressUpdateInformation);

   // COMPUTE ELASTOPLASTIC STEP

   this->ComputeSubstepIncrementalDeformationGradient( rDeformationGradient, HalfPosition, 1, HalfPositionDeformationGradient); 

   this->CalculateExplicitSolution( HalfPositionDeformationGradient, HalfPositionElasticCauchyGreen, rReturnMappingVariables, rNewElasticLeftCauchyGreen, rNewStressVector,  true, rTolerance);


}

   
void NonAssociativeExplicitPlasticFlowRule::ComputeElastoPlasticTangentMatrix(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rLeftCauchyGreenMatrix, const double& rAlpha, Matrix& rElasticMatrix)
{

   Vector ElasticStrainVector = ZeroVector(6);
       
   ElasticStrainVector = this->ConvertCauchyGreenTensorToHenckyStrain(rLeftCauchyGreenMatrix);

   this->ComputeElasticMatrix(ElasticStrainVector, rElasticMatrix);

   //if (rReturnMappingVariables.Control.PlasticRegion) {
   if (rReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION) ) {


   // MATRIX 1.
      AuxiliarDerivativesStructure AuxiliarDerivatives;
      this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, rAlpha);

      double H;
      this->ComputePlasticHardeningParameter(ElasticStrainVector, rAlpha, H);
 
      Vector AuxVectorF;
      Vector AuxVectorG;

      AuxVectorF = prod(trans(AuxiliarDerivatives.YieldFunctionD), rElasticMatrix);
      AuxVectorG = prod( rElasticMatrix, AuxiliarDerivatives.PlasticPotentialD);

      Matrix PlasticUpdate;
      PlasticUpdate = MyCrossProduct(rElasticMatrix, AuxiliarDerivatives.PlasticPotentialD, AuxiliarDerivatives.YieldFunctionD);

      rElasticMatrix -= 1.0*PlasticUpdate / ( H + MathUtils<double>::Dot(AuxVectorF, AuxiliarDerivatives.PlasticPotentialD));
//      std::cout << " MATRIX 1 " << rElasticMatrix << std::endl;
      return;
      // TRY TO COMPUTE THE CONSISTENT TANGENT MATRIX.
      if ( rReturnMappingVariables.DeltaBeta == 0.0) {
         return;
      }
      // MATRIX 2
/*
      ElasticStrainVector = this->ConvertCauchyGreenTensorToHenckyStrain( rReturnMappingVariables.EigenVectors);

      Matrix ElasticMatrixTrial;
      this->ComputeElasticMatrix(ElasticStrainVector, ElasticMatrixTrial);
      this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, rReturnMappingVariables.NormIsochoricStress);

//      double H;
      this->ComputePlasticHardeningParameter(ElasticStrainVector, rReturnMappingVariables.NormIsochoricStress, H);
 
      Matrix AuxiliarInverse = ElasticMatrixTrial;

      int Singular = 0;
      Singular = SolidMechanicsMathUtilities<double>::InvertMatrix( ElasticMatrixTrial, AuxiliarInverse);
  //    std::cout << Hola << " HOLA 1 : " << rElasticMatrix << std::endl;
      AuxiliarInverse += rReturnMappingVariables.DeltaBeta * AuxiliarDerivatives.PlasticPotentialDD;

//      std::cout<< " AUX INV " << AuxiliarInverse << " UPD " << rReturnMappingVariables.DeltaBeta * AuxiliarDerivatives.PlasticPotentialDD << std::endl;

      Singular = SolidMechanicsMathUtilities<double>::InvertMatrix( AuxiliarInverse,ElasticMatrixTrial); 
  //    std::cout << Hola << " HOLA 2 : " << rElasticMatrix << std::endl;

//      Vector AuxVectorF;
//      Vector AuxVectorG;

      AuxVectorF = prod(trans(AuxiliarDerivatives.YieldFunctionD), ElasticMatrixTrial);
      AuxVectorG = prod( ElasticMatrixTrial, AuxiliarDerivatives.PlasticPotentialD);

//      Matrix PlasticUpdate;
      PlasticUpdate = MyCrossProduct(ElasticMatrixTrial, AuxiliarDerivatives.PlasticPotentialD, AuxiliarDerivatives.YieldFunctionD);

      ElasticMatrixTrial -= 1.0*PlasticUpdate / ( H + MathUtils<double>::Dot(AuxVectorF, AuxiliarDerivatives.PlasticPotentialD));
//      rElasticMatrix = ElasticMatrixTrial; 
//      std::cout << " MATRIX 2 " << ElasticMatrixTrial << std::endl;
//      std::cout << "    DBeta " << rReturnMappingVariables.DeltaBeta << std::endl;
//      std::cout << std::endl;
*/
  }

}


Matrix NonAssociativeExplicitPlasticFlowRule::MyCrossProduct(const Matrix& rM, const Vector& rA, const Vector& rB)
{
 

    Vector A = rA;
    Vector B = rB;

    for (unsigned int i = 3; i<6; ++i) {
//      A(i) /= 2.0;
//      B(i) /= 2.0;
    }

    A = prod(rM, A);
    B = prod(trans(B), rM); 

    Matrix Result = ZeroMatrix(6);
    for (unsigned int i = 0; i<6; ++i) {
        for ( unsigned int j = 0; j<6; ++j) {
            Result(i,j) = A(i)*B(j);
        }
     }
     return Result;

}

void NonAssociativeExplicitPlasticFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
}

void NonAssociativeExplicitPlasticFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )
}



} //End Namepace Kratos

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




//NonAssociativeExplicitPlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables,  Matrix & rStressMatrix)
//NonAssociatievExplicitPlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables, const Matrix& DeltaHenckyStrain, ) 

bool NonAssociativeExplicitPlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables, const Matrix& rDeformationGradientF0, const Matrix& rDeltaDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{


   bool PlasticityActive = false;
   rReturnMappingVariables.Options.Set(PLASTIC_REGION,false);

   //1- Initialize some variables
   Vector StressVector = ZeroVector(6);
   InternalVariables PlasticVariables = mInternalVariables;
   PlasticVariables.DeltaPlasticStrain = 0.0;
   
   Matrix PreviousElasticLeftCauchyGreen = ZeroMatrix(3);
   PreviousElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
   
   Matrix NewElasticLeftCauchyGreen;

   double Tolerance = 1E-5;


   //2- Check for yield Condition (at the beggining)
   Vector NewStressVector;
   Vector PreviousStressVector;
   //TO DO SOMETHING WITH THE STR VEC

   this->CalculateKirchhoffStressVector( PreviousElasticLeftCauchyGreen, PreviousStressVector); 
   rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, PreviousStressVector, PlasticVariables.EquivalentPlasticStrain);
   std::cout << "Beggin TrialStateFunction " << rReturnMappingVariables.TrialStateFunction << std::endl;


    //1. Compute an Elastic Trial State
    double StressErrorMeasure;
    double NewEquivalentPlasticStrain;

    this->CalculateOneExplicitStep(rDeltaDeformationGradient, rDeformationGradientF0, PreviousElasticLeftCauchyGreen, PlasticVariables, NewElasticLeftCauchyGreen, NewStressVector, NewEquivalentPlasticStrain, false, StressErrorMeasure);

    double ElasticTrialStateFunction;
    ElasticTrialStateFunction = mpYieldCriterion->CalculateYieldCondition( ElasticTrialStateFunction, NewStressVector, PlasticVariables.EquivalentPlasticStrain);

    //2a. It is an elastic step
    if ( ElasticTrialStateFunction <= Tolerance)   {
        //Error does not matter because we are analitically in an hyperelastic bla bla bla
        PlasticityActive = false;
    }
    //2b. It seems an elasto-plastic step
    else {
         if ( ( rReturnMappingVariables.TrialStateFunction < -Tolerance) && (ElasticTrialStateFunction > Tolerance) ) {
            //2b.2 Transition from elastic to plastic
            this->CalculateExplicitSolutionWithChange( rDeltaDeformationGradient, rDeformationGradientF0, PreviousElasticLeftCauchyGreen, PlasticVariables, NewElasticLeftCauchyGreen, NewStressVector, NewEquivalentPlasticStrain,  Tolerance);
          
            PlasticityActive = true;
         }

//        else if ( (ElasticTrialStateFunction > Tolerance) && ( (rReturnMappingVariables.TrialStateFunction) > -Tolerance) ) {
            //2b.1 The previous drift correction didn't work out, don't worry. Completeldly PlasticStep.
         // AQUI QUEDAN TODOS LOS CASOS QUE DISCUTEN; PERO DE MOMENTO PASO
         else {         
            this->CalculateExplicitSolution( rDeltaDeformationGradient, rDeformationGradientF0, PreviousElasticLeftCauchyGreen, PlasticVariables, NewElasticLeftCauchyGreen, NewStressVector, true, Tolerance);

            PlasticityActive = true;

         }
 

    }


   PlasticVariables.EquivalentPlasticStrain = NewEquivalentPlasticStrain;
   rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, NewStressVector, PlasticVariables.EquivalentPlasticStrain);
   std::cout << "Exit  TrialStateFunction " << rReturnMappingVariables.TrialStateFunction << std::endl;



//    if (PlasticityActive) {
//
//       double DriftViolation;
//       DriftViolation = mpYieldCriterion->CalculateYieldCondition(DriftViolation, NewStressVector, PlasticVariables.EquivalentPlasticStrain);
//
//       if ( fabs( DriftViolation) > Tolerance ) {
//            this->ReturnStressToYieldSurface( NewElasticHenckyStrain, NewStressVector, PlasticVariables.EquivalentPlasticStrain, DriftViolation, Tolerance);
//       }
//    } 

    // COMO EL TEMA CONMUTA forque f = eye(3), luego se puede aprovechar lo anterior!!! BIEN!!!!
/*    if (PlasticityActive) {
       double DriftViolation;
       DriftViolation = mpYieldCriterion->CalculateYieldCondition(DriftViolation, NewStressVector, PlasticVariables.EquivalentPlasticStrain);

       if ( fabs(DriftViolation) > Tolerance ) {
            this->ReturnStressToYieldSurface( NewElasticLeftCauchyGreen, NewStressVector, PlasticVariables.EquivalentPlasticStrain, DriftViolation, Tolerance);
       }


    }
*/
    rReturnMappingVariables.DeltaGamma = PlasticVariables.EquivalentPlasticStrain;
    rStressMatrix = MathUtils<double>::StressVectorToTensor(NewStressVector);
    rNewElasticLeftCauchyGreen = NewElasticLeftCauchyGreen; 


    rReturnMappingVariables.Options.Set(PLASTIC_REGION,PlasticityActive);
    rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED, true);


   return PlasticityActive;

}
// *************** YIELD SURFACE VIOLATION DRIFT CORRECTION **********
// *******************************************************************
void NonAssociativeExplicitPlasticFlowRule::ReturnStressToYieldSurface(Matrix& rNewElasticLeftCauchyGreen, Vector& rStressVector, double& rAlpha, double& rDrift, const double& rTolerance)
{
     Vector ActualElasticHenckyStrain;
     ActualElasticHenckyStrain = this->ConvertCauchyGreenTensorToHenckyStrain( rNewElasticLeftCauchyGreen);

     this->ReturnStressToYieldSurface( ActualElasticHenckyStrain, rStressVector, rAlpha, rDrift, rTolerance);

     rNewElasticLeftCauchyGreen = this->ConvertHenckyStrainToCauchyGreenTensor( ActualElasticHenckyStrain);
     // HAY ALGO MAL PORQUE DEJA DE CONVERJER LO BIEN QUE LO HACIA 
}


void NonAssociativeExplicitPlasticFlowRule::ReturnStressToYieldSurface(Vector& rElasticHenckyStrainVector, Vector& rStressVector, double& rAlpha, double& rDrift, const double& rTolerance)
{
    AuxiliarDerivativesStructure  AuxiliarDerivatives;
    Matrix ElasticMatrix;
    double H;
    double DeltaGamma;
    Vector ElasticCorrection;
        
    std::cout << "DriftViolation " << rDrift << std::endl; 

    for (unsigned int i = 0; i < 80; ++i)
    {

        this->UpdateDerivatives(rElasticHenckyStrainVector, AuxiliarDerivatives, rAlpha);
        this->ComputeElasticMatrix(rElasticHenckyStrainVector, ElasticMatrix);
        this->ComputePlasticHardeningParameter(rElasticHenckyStrainVector, rAlpha, H);

        DeltaGamma = rDrift;
        DeltaGamma /= ( H + MathUtils<double>::Dot(AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD)));

        ElasticCorrection = -DeltaGamma*AuxiliarDerivatives.PlasticPotentialD;

        rElasticHenckyStrainVector += ElasticCorrection;
        for (unsigned int j = 0; i<3; ++i)
           rAlpha -= ElasticCorrection(j);

        this->CalculateKirchhoffStressVector( rElasticHenckyStrainVector, rStressVector);

        rDrift = mpYieldCriterion->CalculateYieldCondition( rDrift, rStressVector, rAlpha);
        std::cout << "DriftViolation " << rDrift << std::endl; 
    
        if (fabs(rDrift) < rTolerance)
            std::cout << "    " << std::endl;
            return;

     }

     if (fabs(rDrift) > 10.0*rTolerance)
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

//****** COMPUTE INCREMENTAL DEFORMATION GRADIENT  ********* 
// FROM BATHE NOTATION ^{rFinalConfiguration} _{rReferenceConfiguration} X
//*********************************************************
void NonAssociativeExplicitPlasticFlowRule::ComputeIncrementalDeformationGradient(const Matrix& rInitialDeformationGradient, const Matrix& rFinalDeformationGradient, const double & rReferenceConfiguration, const double & rFinalConfiguration, Matrix& rIncrementalDeformationGradient)
{

    std::cout << "COMPUTEINCREMENTAL " << std::endl;
    Matrix DeformationGradientReference;
    Matrix DeformationGradientFinal;

    DeformationGradientReference = rReferenceConfiguration*rFinalDeformationGradient + ( 1.0 - rReferenceConfiguration) * rInitialDeformationGradient;
    DeformationGradientFinal     =     rFinalConfiguration*rFinalDeformationGradient + ( 1.0 -     rFinalConfiguration) * rInitialDeformationGradient;

//    std::cout << "-1 ! " << rInitialDeformationGradient << std::endl;
//    std::cout << "0 ! " << rFinalDeformationGradient << std::endl;
//    std::cout << "1 ! "<<  DeformationGradientReference << std::endl;
//    std::cout << "2 ! "<<  DeformationGradientFinal << std::endl;
   
    double Det; 
    MathUtils<double>::InvertMatrix(DeformationGradientReference, rIncrementalDeformationGradient, Det);
//    std::cout << "3 ! " << rIncrementalDeformationGradient << std::endl;
    rIncrementalDeformationGradient = prod( DeformationGradientFinal, rIncrementalDeformationGradient);
//    std::cout << "4 ! " << rIncrementalDeformationGradient << std::endl;
   
}




//********************* UPDATES ONE ELASTOPLASTIC STRAIN INCREMENT THE STRESS ****************
//******************************************************************************

void NonAssociativeExplicitPlasticFlowRule::CalculateOneExplicitPlasticStep(const Matrix& rDeltaDeformationGradient, const Matrix& rPreviousElasticLeftCauchyGreen, InternalVariables& rPlasticVariables, Matrix& rNewElasticLeftCauchyGreen, double& rNewEquivalentPlasticStrain)
{
        double StepEquivalentPlasticStrain = 0.0; 
        Matrix ElasticMatrix;
        Vector ElasticStrainVector ;

        ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rPreviousElasticLeftCauchyGreen);
        this->ComputeElasticMatrix(ElasticStrainVector, ElasticMatrix);

       	AuxiliarDerivativesStructure AuxiliarDerivatives;
       	this->UpdateDerivatives(ElasticStrainVector, AuxiliarDerivatives, StepEquivalentPlasticStrain);

        ElasticStrainVector = ConvertCauchyGreenTensorToHenckyStrain( prod(rDeltaDeformationGradient, trans(rDeltaDeformationGradient)) );


        double H;
        this->ComputePlasticHardeningParameter(ElasticStrainVector, StepEquivalentPlasticStrain, H);


        Vector auxVector;
        auxVector = prod(ElasticMatrix, ElasticStrainVector) ;
        double DeltaGamma = 0.0;
        for (unsigned int i = 0; i<6; ++i)
             DeltaGamma += auxVector(i)*AuxiliarDerivatives.YieldFunctionD(i);

        std::cout << "DeltaGamma 1 " << DeltaGamma << " cosa " << ElasticStrainVector << " auxV " << auxVector <<  std::endl;
        double auxDenominador = H + MathUtils<double>::Dot( AuxiliarDerivatives.YieldFunctionD, prod(ElasticMatrix, AuxiliarDerivatives.PlasticPotentialD));

        std::cout << "DEN " << auxDenominador << std::endl;
	DeltaGamma /= auxDenominador;

        if (DeltaGamma < 0 )
                DeltaGamma =0;

        Vector MuyAuxiliar;
        MuyAuxiliar = -DeltaGamma * AuxiliarDerivatives.PlasticPotentialD / 2.0; 

        Matrix UpdateMatrix;
        UpdateMatrix = this->ConvertHenckyStrainToCauchyGreenTensor (MuyAuxiliar);
        std::cout << "UP1 1 " << UpdateMatrix << std::endl;
        UpdateMatrix = prod( rDeltaDeformationGradient, UpdateMatrix);
        std::cout << "UP1 2 " << UpdateMatrix << std::endl;

        rNewElasticLeftCauchyGreen = prod(UpdateMatrix, rPreviousElasticLeftCauchyGreen);
        rNewElasticLeftCauchyGreen = prod( rNewElasticLeftCauchyGreen, trans(UpdateMatrix));


}

//********************* UPDATES ONE STRAIN INCREMENT THE STRESS ****************
//******************************************************************************

void NonAssociativeExplicitPlasticFlowRule::CalculateOneExplicitStep(const Matrix& rDeltaDeformationGradient, const Matrix& rDeformationGradientF0, const Matrix& rPreviousElasticLeftCauchyGreen, InternalVariables& rPlasticVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector, double& rNewEquivalentPlasticStrain, const bool & rElastoPlasticBool, double& rStressErrorMeasure)
{
        
    double DeltaGamma = 0.0;

    if ( rElastoPlasticBool)  {

        this->CalculateOneExplicitPlasticStep( rDeltaDeformationGradient, rPreviousElasticLeftCauchyGreen, rPlasticVariables, rNewElasticLeftCauchyGreen, rNewEquivalentPlasticStrain);

        this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, rNewStressVector); 

        //ADAPTIVE SECOND PART
        Matrix FinalDeformationGradient = prod(rDeltaDeformationGradient, rDeformationGradientF0);

        Matrix HalfStepIncrementalDeformation;
        ComputeIncrementalDeformationGradient(rDeformationGradientF0, FinalDeformationGradient, 0.0, 0.5, HalfStepIncrementalDeformation); 
        
        Matrix SecondApproxLeftCauchyGreen; 

        this->CalculateOneExplicitPlasticStep( HalfStepIncrementalDeformation, rPreviousElasticLeftCauchyGreen, rPlasticVariables, SecondApproxLeftCauchyGreen, rNewEquivalentPlasticStrain);

        std::cout << " hola " << HalfStepIncrementalDeformation << std::endl;
        ComputeIncrementalDeformationGradient(rDeformationGradientF0, FinalDeformationGradient, 0.5, 1.0, HalfStepIncrementalDeformation);
        std::cout << " hola " << HalfStepIncrementalDeformation << std::endl;
        this->CalculateOneExplicitPlasticStep( HalfStepIncrementalDeformation, SecondApproxLeftCauchyGreen, rPlasticVariables, SecondApproxLeftCauchyGreen, rNewEquivalentPlasticStrain);
 
         Vector SecondApproxStressVector;
         this->CalculateKirchhoffStressVector( SecondApproxLeftCauchyGreen, SecondApproxStressVector);

         rStressErrorMeasure = 0.0;
         double Denominador = 0.0;
         for (unsigned int i = 0; i<6; ++i) {
            rStressErrorMeasure += pow( rNewStressVector(i) - SecondApproxStressVector(i), 2.0);
            Denominador  += pow(rNewStressVector(i), 2.0);
         }
          rStressErrorMeasure = pow( rStressErrorMeasure/Denominador, 0.5);
         
         std::cout << "ST1 1 " << rNewStressVector << std::endl;
         std::cout << "ST1 2 " << SecondApproxStressVector << std::endl;
   } 
   else  {
    
      rNewElasticLeftCauchyGreen = prod( rPreviousElasticLeftCauchyGreen, trans( rDeltaDeformationGradient) );
      rNewElasticLeftCauchyGreen = prod( rDeltaDeformationGradient, rNewElasticLeftCauchyGreen );
  
      rStressErrorMeasure = 0.0;
       
      this->CalculateKirchhoffStressVector( rNewElasticLeftCauchyGreen, rNewStressVector); 
   }

   std::cout << " Gauss Point abstract " << rElastoPlasticBool << std::endl;
   std::cout << "Previous LCG " << rPreviousElasticLeftCauchyGreen << std::endl;
   std::cout << "INcrement    " << rDeltaDeformationGradient << std::endl;
   std::cout << "ALLNEW   LCG " << rNewElasticLeftCauchyGreen << std::endl;
   std::cout << "Stress       " << rNewStressVector << std::endl;
   std::cout << "Error        " << rStressErrorMeasure << std::endl;
   std::cout << "Landa        " << DeltaGamma << std::endl;
   std::cout << " " << std::endl;

}


//***************** CalculateStress From LEFT CAUCHY GREEN MATRIX **************
void NonAssociativeExplicitPlasticFlowRule::CalculateKirchhoffStressVector( const Matrix& rElasticLeftCauchyGreen, Vector& rNewStressVector)
{

    Vector ElasticHenckyStrainVector;
    ElasticHenckyStrainVector = ConvertCauchyGreenTensorToHenckyStrain( rElasticLeftCauchyGreen); 

    this->CalculateKirchhoffStressVector(ElasticHenckyStrainVector, rNewStressVector);
}




void NonAssociativeExplicitPlasticFlowRule::CalculateExplicitSolution( const Matrix & rDeltaDeformationGradient, const Matrix& rDeformationGradientF0, const Matrix& rPreviousElasticLeftCauchyGreen, InternalVariables& rPlasticVariables, Matrix&  rNewElasticLeftCauchyGreen, Vector& rNewStressVector, const bool& rElastoPlasticBool, const double& rTolerance)
{
   double TimeStep = 0.5;
   double MinTimeStep = 1e-5;
   double DoneTimeStep = 0.0;

   Matrix FinalDeformationGradient = prod( rDeltaDeformationGradient, rDeformationGradientF0);
  
   Matrix ActualDeformationGradient = rDeformationGradientF0;
   Matrix ActualElasticLeftCauchyGreen = rPreviousElasticLeftCauchyGreen;
   Matrix IncrementalDeformationGradientF;

   double StressErrorMeasure ;
   double NewEquivalentPlasticStrain = rPlasticVariables.EquivalentPlasticStrain;

   while (DoneTimeStep < 1.0)  {
 
       if (DoneTimeStep + TimeStep >= 1.0) {
           TimeStep = 1.0 - DoneTimeStep;
       }

       this->ComputeIncrementalDeformationGradient( rDeformationGradientF0, FinalDeformationGradient, DoneTimeStep, DoneTimeStep + TimeStep, IncrementalDeformationGradientF);
 
       this->CalculateOneExplicitStep( IncrementalDeformationGradientF, ActualDeformationGradient, ActualElasticLeftCauchyGreen, rPlasticVariables, rNewElasticLeftCauchyGreen, rNewStressVector, NewEquivalentPlasticStrain, rElastoPlasticBool, StressErrorMeasure);

       std::cout << "EXPLICITCOMPUTE " << TimeStep << " " << DoneTimeStep << " " << StressErrorMeasure << std::endl;

       if ( StressErrorMeasure < 100.0*rTolerance ) {
           // Se acepta el paso
           ActualElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
           ActualDeformationGradient = prod( IncrementalDeformationGradientF, ActualDeformationGradient);
           DoneTimeStep += TimeStep;
           rPlasticVariables.EquivalentPlasticStrain = NewEquivalentPlasticStrain;
       }
       else {
          if (TimeStep == MinTimeStep) {
              std::cout << "ExplicitStressIntegrationDidNotConverged" << StressErrorMeasure << " " << TimeStep << std::endl;
              ActualElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
              ActualDeformationGradient = prod( IncrementalDeformationGradientF, ActualDeformationGradient);
              DoneTimeStep += TimeStep;
              rPlasticVariables.EquivalentPlasticStrain = NewEquivalentPlasticStrain;
          }
 
 
      }

      TimeStep *= 0.9* (pow( rTolerance / (StressErrorMeasure+ 1e-6), 1.0/2.0 ));
      TimeStep = std::max(TimeStep, MinTimeStep);
  }

}





// Search for the change and then continue
void NonAssociativeExplicitPlasticFlowRule::CalculateExplicitSolutionWithChange(const Matrix& rDeltaDeformationGradient, const Matrix& rDeformationGradientF0, const Matrix& rPreviousElasticLeftCauchyGreen, InternalVariables& rPlasticVariables, Matrix& rNewElasticLeftCauchyGreen, Vector& rNewStressVector, double& rNewEquivalentPlasticStrain, const double& rTolerance)
{


   Matrix FinalDeformationGradient = prod( rDeltaDeformationGradient, rDeformationGradientF0);

// bisecction in elastic Regime

   double InitialPosition = 0.0;
   double EndPosition = 1.0;
   double HalfPosition;

   double InitialStateFunction = -1.0;
   double EndStateFunction = 1.0;
   double HalfStateFunction;

   Matrix IncrementalDeformationGradientF;

   double StressErrorMeasure;
   double NewEquivalentPlasticStrain;

   double Measure1;
   double Measure2;

   for (unsigned int i = 0 ; i< 150; ++i) {

        HalfPosition = 0.5*(InitialPosition + EndPosition);
        this->ComputeIncrementalDeformationGradient( rDeformationGradientF0, FinalDeformationGradient, 0, HalfPosition, IncrementalDeformationGradientF);

        this->CalculateOneExplicitStep( IncrementalDeformationGradientF, rDeformationGradientF0, rPreviousElasticLeftCauchyGreen, rPlasticVariables, rNewElasticLeftCauchyGreen, rNewStressVector, NewEquivalentPlasticStrain, false, StressErrorMeasure);
        
        HalfStateFunction = mpYieldCriterion->CalculateYieldCondition( HalfStateFunction, rNewStressVector, rPlasticVariables.EquivalentPlasticStrain);

       if ( HalfStateFunction < 0)  {
           InitialStateFunction = HalfStateFunction;
           InitialPosition = HalfPosition;
       }
       else {
           EndPosition = HalfPosition;
           EndStateFunction = HalfStateFunction;
       }
       Measure1 = fabs(InitialStateFunction-EndStateFunction);
       Measure2 = fabs(InitialPosition - EndPosition);
       std::cout << "BISECTION " << i << " M1 " << Measure1 << " M2 " << Measure2 << std::endl;

       if ( (Measure1 < rTolerance ) && (Measure2 < rTolerance)) {
           std::cout << "BISECTION "  <<i  << "Converged " << std::endl;
           break;
       }
 
   }

   // COMPUTE ELASTIC STEP;

   this->ComputeIncrementalDeformationGradient( rDeformationGradientF0, FinalDeformationGradient, 0, HalfPosition, IncrementalDeformationGradientF);

   this->CalculateOneExplicitStep( IncrementalDeformationGradientF, rDeformationGradientF0, rPreviousElasticLeftCauchyGreen, rPlasticVariables, rNewElasticLeftCauchyGreen, rNewStressVector, NewEquivalentPlasticStrain, false, StressErrorMeasure);

   // COMPUTE ELASTOPLASTIC STEP
   Matrix HalfStepDeformationGradient = prod(IncrementalDeformationGradientF, rDeformationGradientF0);

   this->ComputeIncrementalDeformationGradient( rDeformationGradientF0, FinalDeformationGradient, HalfPosition, 1, IncrementalDeformationGradientF); 

   this->CalculateExplicitSolution( IncrementalDeformationGradientF, HalfStepDeformationGradient, rNewElasticLeftCauchyGreen, rPlasticVariables, rNewElasticLeftCauchyGreen, rNewStressVector,  true, rTolerance);

// AQUI HAY UN PROBLEMA; CREO

}
/*{
   unsigned int NumberOfSubsteps = 10;
   double IncrementalAlpha = 1/double(NumberOfSubsteps);

   Matrix FinalDeformationGradient = prod( rDeltaDeformationGradient, rDeformationGradientF0);
  
   Matrix ActualDeformationGradient = rDeformationGradientF0;
   Matrix ActualElasticLeftCauchyGreen = rPreviousElasticLeftCauchyGreen;
   Matrix IncrementalDeformationGradient;
   
   double StateFunctionEndSubStep;
   double StressErrorMeasure ;

   for (unsigned int iSubstep = 0; iSubstep < NumberOfSubsteps, ++iSubstep)
   {
        //1. Compute Small ElasticIncrements and check for YieldSurface
        this->ComputeIncrementalDeformationGradient( rDeformationGradientF0, FinalDeformationGradient, DoneTimeStep, DoneTimeStep + IncrementalAlpha, IncrementalDeformationGradientF);
        this->CalculateOneExplicitStep( IncrementalDeformationGradientF, ActualDeformationGradient, ActualElasticLeftCauchyGreen, rPlasticVariables, rNewElasticLeftCauchyGreen, rNewStressVector, NewEquivalentPlasticStrain, true, StressErrorMeasure);
        StateFunctionEndSubStep = mpYieldCriterion->CalcualteYieldCondition(StateFunctionEndSubStep, rNewStressVector, rPlasticVariables.EquivalentPlasticStrain);
        //Elastic Substep; add and continue;
        if ( StateFunctionEndSubStep < rTolerance )  {
           ActualElasticLeftCauchyGreen = rNewElasticLeftCauchyGreen;
           ActualDeformationGradient = prod( IncrementalDeformationGradientF, AcualDeformationGradient);
           DoneTimeStep += TimeStep;
           rPlasticVariables.EquivalentPlasticStrain = NewEquivalentPlasticStrain;
        }
        // ElastoPlastic Substep. 
        else {
        }
   }



}

*/


   
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

      Matrix PlasticUpdate;
      PlasticUpdate = MyCrossProduct(rElasticMatrix, AuxiliarDerivatives.PlasticPotentialD, AuxiliarDerivatives.YieldFunctionD);

/*      std::cout <<"PlasMat " << PlasticUpdate << std::endl;
      std::cout <<"DERIV   " << AuxiliarDerivatives.YieldFunctionD << std::endl;
      std::cout <<"Denom   " << H << " " << MathUtils<double>::Dot(AuxVectorF, AuxiliarDerivatives.PlasticPotentialD) << std::endl;
      std::cout <<"ElasticM" << rElasticMatrix << std::endl;
*/
      rElasticMatrix -= 1.0*PlasticUpdate / ( H + MathUtils<double>::Dot(AuxVectorF, AuxiliarDerivatives.PlasticPotentialD));

//      std::cout <<"TOTAL   " << rElasticMatrix << std::endl;
  }

}


Matrix NonAssociativeExplicitPlasticFlowRule::MyCrossProduct(const Matrix& rM, const Vector& rA, const Vector& rB)
{
 
    for (unsigned int i = 3; i<6; ++i) {
    //   rA(i) /= 2.0;
    //   rB(i) /= 2.0;
    }

    Vector A = rA;
    Vector B = rB;

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


} //End Namepace Kratos

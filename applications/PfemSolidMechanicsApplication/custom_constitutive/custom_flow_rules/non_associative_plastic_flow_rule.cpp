//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes

#include "custom_constitutive/custom_flow_rules/non_associative_plastic_flow_rule.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{


NonAssociativePlasticFlowRule::NonAssociativePlasticFlowRule()
   :FlowRule()
{
}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

NonAssociativePlasticFlowRule::NonAssociativePlasticFlowRule(YieldCriterionPointer pYieldCriterion)
	:FlowRule(pYieldCriterion)
{
   
}

NonAssociativePlasticFlowRule::NonAssociativePlasticFlowRule( const NonAssociativePlasticFlowRule & rOther)
  :FlowRule(rOther)
{

} 


NonAssociativePlasticFlowRule::~NonAssociativePlasticFlowRule()
{
}

void NonAssociativePlasticFlowRule::InitializeMaterial(YieldCriterionPointer pYieldCriterionPointer, HardeningLawPointer pHardeningPointer, const Properties& rProp)
{
    FlowRule::InitializeMaterial(pYieldCriterionPointer, pHardeningPointer, rProp);
    //mTrialEigenValues = ZeroVector(3);
    //mEigenVectors = ZeroMatrix(3);
    mElasticPrincipalStrain = ZeroVector(3);

    mLargeStrainBool = true; 

}



// rStressMatrix es la matrix de tensiones que devuelve
// rStrainMatrix son las deformaciones que se imponen
// mElasticPrincipalStrain almacena las deformaciones elasticas en los ejes principales (para luego calcular be)
bool NonAssociativePlasticFlowRule::CalculateReturnMapping( RadialReturnVariables& rReturnMappingVariables, const Matrix& rIncrementalDeformationGradient, Matrix& rStressMatrix, Matrix& rNewElasticLeftCauchyGreen)
{

   bool PlasticityActive = false; 
   rReturnMappingVariables.Options.Set(PLASTIC_REGION,false); 
   //1.-Compute Principal Axis

   //Computes the polar decomposition and Gives the PrincipalStrain
   // mTrialEigenValues and mEigenValues are computed and stored (crec que és una guarrada, pero ...)
   Vector MainStrain      = ZeroVector(3); 
   // Matrix EigenVectors    = ZeroMatrix(3);
   // this-> ComputePrincipalAxisStrain(rReturnMappingVariables, rNewElasticLeftCauchyGreen, MainStrain, EigenVectors);

   for (unsigned int i = 0; i<3; ++i)
     MainStrain[i] = rNewElasticLeftCauchyGreen(i,i);

   double YoungModulus        = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
   double PoissonCoefficient  = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];

   double LameLambda         = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
   double LameMu             =  YoungModulus/(2*(1+PoissonCoefficient));

   //2.-Compute ElasticMatrix &  Trial State
   Matrix ElasticMatrix = ZeroMatrix(3);
   for (unsigned int i = 0; i<3; ++i)
   {  
     for (unsigned int j = 0; j<3; ++j)
       ElasticMatrix(i,j) = LameLambda + 2.0/3.0*LameMu;

     ElasticMatrix(i,i) += 2.0*LameMu;
   }
 
   Vector PrincipalStress = ZeroVector(3);
   PrincipalStress = prod(ElasticMatrix, MainStrain);

   //3.- Check for the yield Condition
   InternalVariables PlasticVariables = mInternalVariables;
   rReturnMappingVariables.TrialStateFunction = 0.0;
   rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, PrincipalStress, PlasticVariables.EquivalentPlasticStrain );

   if (rReturnMappingVariables.TrialStateFunction < 0.0)
   {
      PlasticVariables.DeltaPlasticStrain = 0;
      mElasticPrincipalStrain = MainStrain;
      PlasticityActive = false; 
      rReturnMappingVariables.Options.Set(PLASTIC_REGION,false); 
   }
   else
   {
      // First of all compute the inverse elastic matrix, since it's a tedious thing to do.
      Matrix InverseElasticMatrix = ZeroMatrix(3);
      CalculateInverseElasticMatrix(rReturnMappingVariables, InverseElasticMatrix);

      mElasticPrincipalStrain = MainStrain;
      bool converged = this->CalculateConsistencyCondition(rReturnMappingVariables, PlasticVariables, InverseElasticMatrix, PrincipalStress, mElasticPrincipalStrain);
      if (!converged)
         std::cout<<" Constit Law Did Not Converge "<<std::endl;

      PlasticityActive = true; 
      rReturnMappingVariables.Options.Set(PLASTIC_REGION,true); 
   }

   this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, PrincipalStress, rStressMatrix);

   rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED,true);
      
   return PlasticityActive; 


}
      

void NonAssociativePlasticFlowRule::CalculateInverseElasticMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rInverseElasticMatrix)
{

        double YoungModulus        = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
	double PoissonCoefficient  = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];

	double LameLambda         = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
	double LameMu             =  YoungModulus/(2*(1+PoissonCoefficient));

	double Diagonal    = (LameLambda + LameMu)/(LameMu*(3.0*LameLambda+2.0*LameMu));
	double NonDiagonal = (-LameLambda)/( 2.0*LameMu*(3.0*LameLambda + 2.0*LameMu));

	for (unsigned int i = 0; i<3; ++i)
	{
	  for (unsigned int j = 0; j<3; ++j)
  	  {
	     if (i == j) {
		rInverseElasticMatrix(i,i) = Diagonal;
	     }
	     else {
		rInverseElasticMatrix(i,j) = NonDiagonal;
	     }
	   }
 	}
}

bool NonAssociativePlasticFlowRule::CalculateConsistencyCondition(RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables,const Matrix& rInverseElasticMatrix, Vector& rPrincipalStress, Vector& rPrincipalStrain)
{

	//Set convergence parameters
	unsigned int iter   = 0;
	double Tolerance    = 1e-3;
        double MaxIterations= 80;

	// start
	double DeltaGamma = 0.0;
        double DeltaDeltaGamma    = 0.0;
	//double DeltaStateFunction = 0;
        rReturnMappingVariables.DeltaGamma  = 0;

	double StateFunction     = rReturnMappingVariables.TrialStateFunction;
	//double InitialEquivalentPlasticStrain = rPlasticVariables.EquivalentPlasticStrain;
	
	Vector PlasticStrain = ZeroVector(3);
	Vector DeltaPlasticStrain = ZeroVector(3);
	Vector DeltaStress = ZeroVector(3);
	Vector Residual = ZeroVector(3);
        AuxiliarDerivativesStructure AuxiliarDerivatives;

	while (fabs(StateFunction) >= Tolerance && iter<=MaxIterations)
	{
		// 1. Compute DeltaDeltaGamma
		this->UpdateDerivatives(rPrincipalStress, AuxiliarDerivatives);
		
		Residual = - PlasticStrain + DeltaGamma*AuxiliarDerivatives.PlasticPotentialD;
		Matrix A = ZeroMatrix(3);
		A = rInverseElasticMatrix + DeltaGamma*AuxiliarDerivatives.PlasticPotentialDD;
		Matrix Ainverse = ZeroMatrix(3);
		double detA;
		MathUtils<double>::InvertMatrix(A, Ainverse, detA);
		Vector AuxiliarVector = ZeroVector(3);
		
		AuxiliarVector  = prod(AuxiliarDerivatives.YieldFunctionD, Ainverse);
		
		
		DeltaDeltaGamma = StateFunction -  MathUtils<double>::Dot3(AuxiliarVector, Residual);
		DeltaDeltaGamma /= MathUtils<double>::Dot3(AuxiliarVector, AuxiliarDerivatives.PlasticPotentialD);
		// 2. Compute New Stress and Strain
		AuxiliarVector = Residual + DeltaDeltaGamma*AuxiliarDerivatives.PlasticPotentialD;
		DeltaStress = prod(Ainverse, AuxiliarVector);
		DeltaPlasticStrain = prod( rInverseElasticMatrix, DeltaStress);

		// 3.  ACUMULAR Y CERRAR.
		DeltaGamma   += DeltaDeltaGamma; 
		PlasticStrain += DeltaPlasticStrain;
		rPrincipalStress -= DeltaStress;
		iter += 1;
		StateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, rPrincipalStress, rPlasticVariables.EquivalentPlasticStrain); 
	}
	rReturnMappingVariables.DeltaGamma = DeltaGamma;
        rPrincipalStrain -= PlasticStrain;
	if (iter> MaxIterations)
		return false;
	
	return true;
 
}



  void NonAssociativePlasticFlowRule::ComputePrincipalAxisStrain(RadialReturnVariables& rReturnMappingVariables, const Matrix& rStrainMatrix, Vector& rPrincipalStrain, Matrix& rEigenVectors)
{

   Matrix StrainMatrix = ZeroMatrix(3);
   if (mLargeStrainBool)
   {
      StrainMatrix = rStrainMatrix;
   } else
   {
   // Aquí sha de posar StrainMatrix = epsilon;
   }
   
   rEigenVectors    =  ZeroMatrix(3);
   Vector EigenValues     =  ZeroVector(3);

   SolidMechanicsMathUtilities<double>::EigenVectors(StrainMatrix, rEigenVectors, EigenValues);

   if (mLargeStrainBool)
   {
      for (unsigned int i = 0; i<3; ++i)
         rPrincipalStrain[i] = 0.50*std::log(EigenValues[i]);
   
   } else 
   {
       for (unsigned int i = 0; i<3; ++i)
          rPrincipalStrain[i] = EigenValues[i];
   }


}


void NonAssociativePlasticFlowRule::ReturnStressFromPrincipalAxis(const Matrix& rEigenVectors, const Vector& rPrincipalStress, Matrix& rStressMatrix)
{
   rStressMatrix = ZeroMatrix(3); //Esto seguro que está prohibido
   Vector auxN = ZeroVector(3);
   Matrix auxM = ZeroMatrix(3);
   for (unsigned int i = 0; i<3; ++i)
   {
      for (unsigned int j = 0; j<3; ++j)
         auxN(j) = rEigenVectors(j,i);
      auxM = MathUtils<double>::TensorProduct3(auxN, auxN);
      rStressMatrix += rPrincipalStress(i)*auxM;
   }
}

void NonAssociativePlasticFlowRule::CalculatePrincipalAxisAlgorithmicTangent(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rStressMatrix, Matrix& rConsistMatrix)
{
  Matrix InverseElasticMatrix = ZeroMatrix(3);
  this->CalculateInverseElasticMatrix(rReturnMappingVariables, InverseElasticMatrix);
  AuxiliarDerivativesStructure AuxiliarDerivatives;
  Vector PrincipalStress = ZeroVector(3);
  {
	Matrix aux=ZeroMatrix(3);
	aux = prod(trans(rReturnMappingVariables.MainDirections), rStressMatrix);
	aux = prod(aux, rReturnMappingVariables.MainDirections);
	for (unsigned int i = 0; i<3; ++i)
           PrincipalStress(i) = aux(i,i);
  }
  this->UpdateDerivatives(PrincipalStress, AuxiliarDerivatives);

  InverseElasticMatrix += rReturnMappingVariables.DeltaGamma * AuxiliarDerivatives.PlasticPotentialDD;
  Matrix ConsistElastic = ZeroMatrix(3);
  double detA;
  MathUtils<double>::InvertMatrix(InverseElasticMatrix, rConsistMatrix, detA);

  Vector auxF = ZeroVector(3);
  Vector auxG = ZeroVector(3);
 
  auxG = prod(rConsistMatrix, AuxiliarDerivatives.PlasticPotentialD);
  auxF = prod(AuxiliarDerivatives.YieldFunctionD, rConsistMatrix);
  double denominador;
  denominador = MathUtils<double>::Dot3(auxF, AuxiliarDerivatives.PlasticPotentialD);
  if (fabs(denominador) > 1){
	  rConsistMatrix -= (1/denominador)*MathUtils<double>::TensorProduct3(auxG, auxF);
  }
 
}


//Vector& NonAssociativePlasticFlowRule::GetStressVectorFromMatrix(const Matrix& rStressMatrix)
//{
//   Vector OutPut = ZeroVector(3);
//   Matrix auxMatrix = ZeroMatrix(3);
//   auxMatrix = prod( rStressMatrix, mEigenVectors);
//   auxMatrix = prod( trans(mEigenVectors), auxMatrix);
//   for (unsigned int i = 0; i<3; ++i)
//      OutPut(i) = auxMatrix(i,i);
//   return OutPut;
//
//}


Matrix NonAssociativePlasticFlowRule::GetElasticLeftCauchyGreen(const RadialReturnVariables& rReturnMappingVariables)
{

   Vector Landa2 = ZeroVector(3);
   if (mLargeStrainBool)
   {   
       for (unsigned int i = 0; i<3; ++i)
            Landa2(i) = std::exp(2.0*mElasticPrincipalStrain(i));
   } 
   else {
       Landa2 = mElasticPrincipalStrain;
   }
   Matrix OutPut = ZeroMatrix(3);
   this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.MainDirections, Landa2, OutPut);
   return OutPut;

}

void NonAssociativePlasticFlowRule::UpdateDerivatives(const Vector& rPrincipalStress, AuxiliarDerivativesStructure& rAuxiliarDerivatives)
{
   this->CalculatePlasticPotentialDerivatives(rPrincipalStress, rAuxiliarDerivatives.PlasticPotentialD, rAuxiliarDerivatives.PlasticPotentialDD);
   mpYieldCriterion->CalculateYieldFunctionDerivative(rPrincipalStress, rAuxiliarDerivatives.YieldFunctionD);
}


void NonAssociativePlasticFlowRule::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FlowRule )
}

void NonAssociativePlasticFlowRule::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FlowRule )
}


bool NonAssociativePlasticFlowRule::UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables )
{
	
	//mInternalVariables.EquivalentPlasticStrainOld  = PlasticVariables.EquivalentPlasticStrain;
        mInternalVariables.EquivalentPlasticStrainOld = 0.0;

	mInternalVariables.DeltaPlasticStrain          = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

	mInternalVariables.EquivalentPlasticStrain    += mInternalVariables.DeltaPlasticStrain;

	mInternalVariables.DeltaPlasticStrain         *= ( 1.0/rReturnMappingVariables.DeltaTime );
 	
	return true;
}

} //nameespace Kratos

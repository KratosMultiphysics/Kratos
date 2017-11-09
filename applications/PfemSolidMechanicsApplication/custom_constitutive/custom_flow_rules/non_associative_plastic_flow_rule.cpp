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
bool NonAssociativePlasticFlowRule::CalculateReturnMapping(RadialReturnVariables& rReturnMappingVariables, Matrix& rStressMatrix, const Vector& rPrincipalStrainTrial)
{

   bool PlasticityActive = false; 
   rReturnMappingVariables.Options.Set(PLASTIC_REGION,false); 
   //1.-Compute Principal Axis
//   Vector PrincipalStrain = ZeroVector(3); 

   //Computes the polar decomposition and Gives the PrincipalStrain
   // mTrialEigenValues and mEigenValues are computed and stored (crec que és una guarrada, pero ...)
//   this->ComputePrincipalAxisStrain(rReturnMappingVariables, rStrainMatrix, PrincipalStrain);

   //2.-Compute ElasticMatrix &  Trial State
   Matrix ElasticMatrix = ZeroMatrix(3,3);
   const double& Young = rReturnMappingVariables.YoungModulus;
   const double& Nu    = rReturnMappingVariables.PoissonCoefficient;
   for (unsigned int i = 0; i<3; ++i)   { 
        for (unsigned int j = 0; j<3; ++j)  {
           if (i == j)   {
               ElasticMatrix(i,j) = Young / (1.0+Nu) / (1.0-2.0*Nu)*(1.0-Nu)   ;
           }
           else {
               ElasticMatrix(i,j) = Young / (1.0+Nu) / (1.0-2.0*Nu)*(Nu);
           }
        }
   }
 
   Vector PrincipalStress = ZeroVector(3);
   PrincipalStress = prod(ElasticMatrix, rPrincipalStrainTrial);

   //3.- Check for the yield Condition
   InternalVariables PlasticVariables = mInternalVariables;
   rReturnMappingVariables.TrialStateFunction = 0.0;
   rReturnMappingVariables.TrialStateFunction = mpYieldCriterion->CalculateYieldCondition(rReturnMappingVariables.TrialStateFunction, PrincipalStress, PlasticVariables.EquivalentPlasticStrain );

   if (rReturnMappingVariables.TrialStateFunction < 0.0)
   {
      PlasticVariables.DeltaPlasticStrain = 0;
      mElasticPrincipalStrain = rPrincipalStrainTrial;
      PlasticityActive = false; 
      rReturnMappingVariables.Options.Set(PLASTIC_REGION,false); 
   }
   else
   {
      // First of all compute the inverse elastic matrix, since it's a tedious thing to do.
      Matrix InverseElasticMatrix = ZeroMatrix(3,3);
      CalculateInverseElasticMatrix(rReturnMappingVariables, InverseElasticMatrix);

      mElasticPrincipalStrain = rPrincipalStrainTrial;
      bool converged = this->CalculateConsistencyCondition(rReturnMappingVariables, PlasticVariables, InverseElasticMatrix, PrincipalStress, mElasticPrincipalStrain);
      if (!converged)
         std::cout<<" Constit Law Did Not Converge "<<std::endl;

      PlasticityActive = true; 
      rReturnMappingVariables.Options.Set(PLASTIC_REGION,true); 
   }

   this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.EigenVectors, PrincipalStress, rStressMatrix);

   rReturnMappingVariables.Options.Set(RETURN_MAPPING_COMPUTED,true);

   return PlasticityActive; 


}
      

void NonAssociativePlasticFlowRule::CalculateInverseElasticMatrix(const RadialReturnVariables& rReturnMappingVariables, Matrix& rInverseElasticMatrix)
{
	const double& Young = rReturnMappingVariables.YoungModulus;
	const double& Nu = rReturnMappingVariables.PoissonCoefficient;
	double Diagonal = 1.0/Young;
	double NonDiagonal = -Nu/Young;

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
	double Tolerance    = 1e-6;
        double MaxIterations= 80;

	// start
	double DeltaGamma = 0.0;
        double DeltaDeltaGamma    = 0.0;
        double DeltaStateFunction = 0;
        rReturnMappingVariables.DeltaGamma  = 0;

	double StateFunction     = rReturnMappingVariables.TrialStateFunction;
	double InitialEquivalentPlasticStrain = rPlasticVariables.EquivalentPlasticStrain;
	
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
		Matrix A = ZeroMatrix(3,3);
		A = rInverseElasticMatrix + DeltaGamma*AuxiliarDerivatives.PlasticPotentialDD;
		Matrix Ainverse = ZeroMatrix(3,3);
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



void NonAssociativePlasticFlowRule::ComputePrincipalAxisStrain(RadialReturnVariables& rReturnMappingVariables, const Matrix& rStrainMatrix, Vector& rPrincipalStrain)
{

   Matrix StrainMatrix = ZeroMatrix(3,3);
   mLargeStrainBool = true;
   if (mLargeStrainBool)
   {
      StrainMatrix = rStrainMatrix;
   } else
   {
   // Aquí sha de posar StrainMatrix = epsilon;
   }
   
   SolidMechanicsMathUtilities<double>::EigenVectors(StrainMatrix, rReturnMappingVariables.EigenVectors, rReturnMappingVariables.TrialEigenValues);
   if (mLargeStrainBool)
   {
      for (unsigned int i = 0; i<3; ++i)
         rPrincipalStrain(i) = 0.50*std::log(rReturnMappingVariables.TrialEigenValues(i));
   
   } else 
   {
       for (unsigned int i = 0; i<3; ++i)
          rPrincipalStrain(i) = rReturnMappingVariables.TrialEigenValues(i);
   }


}


void NonAssociativePlasticFlowRule::ReturnStressFromPrincipalAxis(const Matrix& rEigenVectors, const Vector& rPrincipalStress, Matrix& rStressMatrix)
{
    rStressMatrix = ZeroMatrix(3,3); //Esto seguro que está prohibido
    Matrix Aux = ZeroMatrix(3,3);
 
    for (unsigned int j = 0; j<3; ++j)
        Aux(j,j) = rPrincipalStress(j);

    rStressMatrix = prod(Aux, (rEigenVectors));
    rStressMatrix = prod(trans(rEigenVectors), rStressMatrix);

}


void NonAssociativePlasticFlowRule::UpdateDerivatives(const Vector& rPrincipalStress, AuxiliarDerivativesStructure& rAuxiliarDerivatives)
{
   this->CalculatePlasticPotentialDerivatives(rPrincipalStress, rAuxiliarDerivatives.PlasticPotentialD, rAuxiliarDerivatives.PlasticPotentialDD);
   mpYieldCriterion->CalculateYieldFunctionDerivative(rPrincipalStress, rAuxiliarDerivatives.YieldFunctionD);
   rAuxiliarDerivatives.YieldFunctionD = rAuxiliarDerivatives.PlasticPotentialD;
}

void NonAssociativePlasticFlowRule::UpdateDerivatives(const Matrix& rEigenVectors, const Matrix& rStressMatrix, AuxiliarDerivativesStructure& rAuxiliarDerivatives)
{

    // STEP 1. Compute the Stress In Principal Axis
    Matrix AuxiliarStressMatrix ;
    AuxiliarStressMatrix = prod( rStressMatrix, trans(rEigenVectors) );
    AuxiliarStressMatrix = prod( rEigenVectors, AuxiliarStressMatrix);

    Vector StressVector = ZeroVector(3);

    for (unsigned int i = 0; i < 3; ++i)
        StressVector(i) = AuxiliarStressMatrix(i,i);
 
    // STEP 2. Compute Principal Axis Derivatives
    AuxiliarDerivativesStructure PrincipalAxisDerivatives;
    this->UpdateDerivatives(StressVector, PrincipalAxisDerivatives);

    // STEP3. Reconvert to the Cartesian Axis

    Matrix YieldDerivative = ZeroMatrix(3,3);
    for (unsigned int i = 0; i <3; ++i)
        YieldDerivative(i,i) = PrincipalAxisDerivatives.YieldFunctionD(i);

    YieldDerivative = prod( YieldDerivative, rEigenVectors);
    YieldDerivative = prod( trans(rEigenVectors), YieldDerivative);

    rAuxiliarDerivatives.YieldFunctionD = MathUtils<double>::StrainTensorToVector( YieldDerivative, 6);
    rAuxiliarDerivatives.PlasticPotentialD = rAuxiliarDerivatives.YieldFunctionD;

}


void NonAssociativePlasticFlowRule::CalculateElasticMatrix(const RadialReturnVariables& rReturnMappingVariables,  Matrix& rElasticMatrix)
{
  
   rElasticMatrix = ZeroMatrix(6,6);

   const double& Young = rReturnMappingVariables.YoungModulus;
   const double& Nu    = rReturnMappingVariables.PoissonCoefficient;

   for (unsigned int i = 0; i<3; ++i)   { 
        for (unsigned int j = 0; j<3; ++j)  {
           if (i == j)   {
               rElasticMatrix(i,j) = Young / (1.0+Nu) / (1.0-2.0*Nu)*(1.0-Nu)   ;
           }
           else {
               rElasticMatrix(i,j) = Young / (1.0+Nu) / (1.0-2.0*Nu)*(Nu);
           }
        }
   }

   for (unsigned int i = 3; i <6; ++i)
          rElasticMatrix(i,i) = Young / ( 1.0 + Nu) / 2.0;

}

void NonAssociativePlasticFlowRule::CalculatePrincipalAxisAlgorithmicTangent(const RadialReturnVariables& rReturnMappingVariables, const Matrix& rStressMatrix, Matrix& rConsistMatrix)
{


  Matrix ElasticMatrix = ZeroMatrix(3,3);
  //this->CalculateInverseElasticMatrix(rReturnMappingVariables, InverseElasticMatrix);
  this->CalculateElasticMatrix(rReturnMappingVariables, ElasticMatrix) ; 
  if ( ! rReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION))  {

    rConsistMatrix = ElasticMatrix;

    return;    

  }

  // TO DO: SECOND DERIVATIVE INVERSE ....

  AuxiliarDerivativesStructure AuxiliarDerivatives; 
  this->UpdateDerivatives( rReturnMappingVariables.EigenVectors, rStressMatrix, AuxiliarDerivatives);

  rConsistMatrix = ElasticMatrix;

  Vector auxF; 
  Vector auxG; 
 
  auxG = prod(rConsistMatrix, AuxiliarDerivatives.PlasticPotentialD);
  auxF = prod(AuxiliarDerivatives.YieldFunctionD, rConsistMatrix);

  double denominador = 0.0;

  for (unsigned int i = 0; i < 6; ++i )
      denominador += auxF(i)*AuxiliarDerivatives.PlasticPotentialD(i);

  //denominador = MathUtils<double>::Dot3(auxF, AuxiliarDerivatives.PlasticPotentialD);

  if (fabs(denominador) > 1e-8){
        for (unsigned int i = 0; i < 6; ++i)  {
           for (unsigned j = 0; j <6; ++j)   {
               rConsistMatrix(i,j) -= (1/denominador)* auxG(i)*auxF(j) ;
           }
       }
	 // rConsistMatrix -= (1/denominador)*MathUtils<double>::TensorProduct3(auxG, auxF);
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
   mLargeStrainBool = true;
   if (mLargeStrainBool)
   {   
       for (unsigned int i = 0; i<3; ++i)
            Landa2(i) = std::exp(2.0*mElasticPrincipalStrain(i));
   } 
   else {
       for (unsigned int i = 0; i<3; ++i)
       Landa2(i) = std::exp(2.0*mElasticPrincipalStrain(i)); //mElasticPrincipalStrain;
   }
   Matrix OutPut = ZeroMatrix(3,3);

   this->ReturnStressFromPrincipalAxis(rReturnMappingVariables.EigenVectors, Landa2, OutPut);

   return OutPut;


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
   
        mInternalVariables.EquivalentPlasticStrain *= 0.9;

	mInternalVariables.DeltaPlasticStrain          = sqrt(2.0/3.0) * rReturnMappingVariables.DeltaGamma;

	mInternalVariables.EquivalentPlasticStrain    += mInternalVariables.DeltaPlasticStrain;

	mInternalVariables.DeltaPlasticStrain         *= ( 1.0/rReturnMappingVariables.DeltaTime );

	return true;
}

} //nameespace Kratos

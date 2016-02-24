//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/non_linear_hencky_plastic_3D_law.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

NonLinearHenckyElasticPlastic3DLaw::NonLinearHenckyElasticPlastic3DLaw()
   : HyperElasticPlastic3DLaw()
{
}


NonLinearHenckyElasticPlastic3DLaw::NonLinearHenckyElasticPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
   : HyperElasticPlastic3DLaw( pFlowRule, pYieldCriterion, pHardeningLaw)
{
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

NonLinearHenckyElasticPlastic3DLaw::NonLinearHenckyElasticPlastic3DLaw(const NonLinearHenckyElasticPlastic3DLaw&  rOther)
  : HyperElasticPlastic3DLaw(rOther)
{
}

 
//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer NonLinearHenckyElasticPlastic3DLaw::Clone() const
{
    NonLinearHenckyElasticPlastic3DLaw::Pointer p_clone(new NonLinearHenckyElasticPlastic3DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

NonLinearHenckyElasticPlastic3DLaw::~NonLinearHenckyElasticPlastic3DLaw()
{
}

//************************************************************************************
//************************************************************************************

void NonLinearHenckyElasticPlastic3DLaw::InitializeMaterial(const Properties& rProps, 
							    const GeometryType& rGeom,
							    const Vector& rShapeFunctionsValues)
{

   mDeterminantF0                = 1;
   mInverseDeformationGradientF0 = identity_matrix<double> (3);
   mElasticLeftCauchyGreen       = identity_matrix<double> (3);

   mpFlowRule->InitializeMaterial(mpYieldCriterion, mpHardeningLaw, rProps);
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//*****************************MATERIAL RESPONSES*************************************
//************************************************************************************

void NonLinearHenckyElasticPlastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const ProcessInfo&  CurrentProcessInfo = rValues.GetProcessInfo();

    const Matrix&   DeformationGradientF   = rValues.GetDeformationGradientF();
    const double&   DeterminantF           = rValues.GetDeterminantF();

    const GeometryType&  DomainGeometry    = rValues.GetElementGeometry ();
    const Vector&        ShapeFunctions    = rValues.GetShapeFunctionsValues ();

    Vector& StrainVector                   = rValues.GetStrainVector();
    Vector& StressVector                   = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix             = rValues.GetConstitutiveMatrix();


    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    ElasticVariables.IdentityMatrix = identity_matrix<double> ( 3 );

    ElasticVariables.SetElementGeometry(DomainGeometry);
    ElasticVariables.SetShapeFunctionsValues(ShapeFunctions);

    FlowRule::RadialReturnVariables ReturnMappingVariables;
    //ReturnMappingVariables.initialize(); //it has to be called at the start
    ReturnMappingVariables.clear();

    // Initialize variables from the process information
    ReturnMappingVariables.DeltaTime = CurrentProcessInfo[DELTA_TIME];

    if(CurrentProcessInfo[IMPLEX] == 1)	
      ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,true);
    else
      ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,false);

    //1.- Lame constants
    // const double& YoungModulus       = MaterialProperties[YOUNG_MODULUS];
    // const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];

    //ElasticVariables.LameLanda      = (YoungModulus*PoissonCoefficient)/((1.0+PoissonCoefficient)*(1.0-2.0*PoissonCoefficient));
    //ElasticVariables.LameMu         =  YoungModulus/(2.0*(1.0+PoissonCoefficient));

    //2.-Determinant of the Total Deformation Gradient
    ElasticVariables.DeterminantF = DeterminantF;

    //3.-Compute Incremental DeformationGradient (in 3D)
    ElasticVariables.DeformationGradientF = DeformationGradientF;
    ElasticVariables.DeformationGradientF = Transform2DTo3D(ElasticVariables.DeformationGradientF);
    ElasticVariables.DeformationGradientF = prod(ElasticVariables.DeformationGradientF,mInverseDeformationGradientF0); 

    //4.-Left Cauchy Green tensor b: (stored in the CauchyGreenMatrix)
    ElasticVariables.CauchyGreenMatrix = mElasticLeftCauchyGreen;

    //4.-Almansi Strain:
    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
    {
        // e= 0.5*(1-invbT*invb)
        this->CalculateAlmansiStrain(ElasticVariables.CauchyGreenMatrix,StrainVector);
    }

 
    //5.-Calculate Total Kirchhoff stress
    Matrix StressMatrix = ZeroMatrix(3);    
    Matrix NewElasticLeftCauchyGreen = mElasticLeftCauchyGreen; 

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
         this->CalculateOnlyDeviatoricPart( ElasticVariables.DeformationGradientF );
         mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, ElasticVariables.DeformationGradientF, StressMatrix, NewElasticLeftCauchyGreen);
         this->CorrectDomainPressure( StressMatrix, ElasticVariables);
    }
    //OPTION 1:
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {
        //Kirchhoff Stress:
        StressVector = MathUtils<double>::StressTensorToVector(StressMatrix, StressVector.size());

    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {

        //initialize constitutive tensors
        ConstitutiveMatrix.clear();

        Matrix ElastoPlasticTangentMatrix;
        double rAlpha = ReturnMappingVariables.DeltaGamma;

        this->CalculateElastoPlasticTangentMatrix( ReturnMappingVariables, NewElasticLeftCauchyGreen, rAlpha, ElastoPlasticTangentMatrix, ElasticVariables);

        //if( Options.Is(ConstitutiveLaw::TOTAL_TENSOR ) )
        ConstitutiveMatrix = this->SetConstitutiveMatrixToAppropiateDimension(ElastoPlasticTangentMatrix);

    }

    if( Options.Is( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ) )
    {
      mpFlowRule->UpdateInternalVariables ( ReturnMappingVariables );

      mElasticLeftCauchyGreen = NewElasticLeftCauchyGreen;

      ElasticVariables.DeformationGradientF = DeformationGradientF;
      ElasticVariables.DeformationGradientF = Transform2DTo3D(ElasticVariables.DeformationGradientF);
      MathUtils<double>::InvertMatrix( ElasticVariables.DeformationGradientF, mInverseDeformationGradientF0, mDeterminantF0);
      mDeterminantF0 = DeterminantF; //special treatment of the determinant     
   }

}

Matrix& NonLinearHenckyElasticPlastic3DLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
{
   if ( rThisVariable == KIRCHHOFF_STRESS_TENSOR )
   {

      Matrix StressMatrix;
      Matrix ElasticLeftCauchyGreen = mElasticLeftCauchyGreen;
      Matrix DeformationGradientF   = ZeroMatrix(3);

      for (unsigned int i = 0; i < 3; ++i)
         DeformationGradientF(i,i) = 1.0;


      FlowRule::RadialReturnVariables ReturnMappingVariables;

      mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, DeformationGradientF, StressMatrix, ElasticLeftCauchyGreen);

      rValue = StressMatrix;

   }
   else if ( rThisVariable == ELASTIC_LEFT_CAUCHY_GREEN_TENSOR)
   {
      rValue = mElasticLeftCauchyGreen;
   }
   else {
      rValue = HyperElasticPlastic3DLaw::GetValue( rThisVariable, rValue);
   }
   return rValue;
}

double& NonLinearHenckyElasticPlastic3DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
{

   // POST-PROCESS LOTS OF VARIABLES. ONLY IF DEFINED IN THE APPROPIATE FILE
   if ( (rThisVariable==VOLUMETRIC_PLASTIC) || (rThisVariable==INCR_SHEAR_PLASTIC) || (rThisVariable==PLASTIC_STRAIN) )
   {

       const FlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
       if ( rThisVariable==VOLUMETRIC_PLASTIC) {
             rValue = InternalVariables.EquivalentPlasticStrain; 
       }
       else if (rThisVariable==PLASTIC_STRAIN) {
            rValue = InternalVariables.DeltaPlasticStrain;
       }
       else {
            rValue = InternalVariables.EquivalentPlasticStrainOld;
       }

       return rValue;
   }
   else if ( (rThisVariable==STRESS_INV_P) || (rThisVariable==STRESS_INV_J2) || (rThisVariable==STRESS_INV_THETA)  ) 
   {
         Matrix StressMatrix;
         Matrix NewElasticLeftCauchyGreen = mElasticLeftCauchyGreen;
         Matrix DeformationGradientF = ZeroMatrix(3);
         for (unsigned int i = 0; i < 3; ++i)
	   DeformationGradientF(i,i) = 1.0;
 
         FlowRule::RadialReturnVariables ReturnMappingVariables;
       
         mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, DeformationGradientF, StressMatrix, NewElasticLeftCauchyGreen);

         double MeanStress = 0.0;
         for (unsigned int i = 0; i < 3; ++i)
              MeanStress += StressMatrix(i,i)/3.0;

         if (rThisVariable==STRESS_INV_P) {
              rValue = -MeanStress;
              return rValue;
         }

         double StressQ = 0.0;
         for (unsigned int i = 0; i <3; ++i) 
              StressQ += pow( StressMatrix(i,i) - MeanStress, 2.0);

          StressQ += 2.0*pow( StressMatrix(0,1) , 2.0);
          StressQ += 2.0*pow( StressMatrix(0,2) , 2.0);
          StressQ += 2.0*pow( StressMatrix(1,2) , 2.0);
    
          if (rThisVariable== STRESS_INV_J2)  {

              rValue = pow( 3.0/2.0*StressQ, 1.0/2.0);
              return rValue; 
          }

          StressQ /= 2.0;

          if (StressQ < 1e-5)
               rValue = 0.0;
         else { 

             for (unsigned int i = 0; i < 3 ; ++i ) 
                 StressMatrix(i,i) -= MeanStress;
 
             double ThirdInvariant = 0.0;
     
             ThirdInvariant = MathUtils<double>::Det( StressMatrix );
 
             ThirdInvariant *= 3.0*pow( 3.0, 1.0/2.0) / 2.0;
             ThirdInvariant /= pow( StressQ, 3.0/2.0);
             double Epsi = 1e-5;
             if (ThirdInvariant > 1.0-Epsi) {
                 rValue = 30.0;
             }
             else if (ThirdInvariant < -1.0+Epsi) {
                 rValue = -30.0;
             }
             else {
                 ThirdInvariant = std::asin( ThirdInvariant);
                 rValue = ThirdInvariant / 3.0 * 180.0 / 3.14159265359; 
             }
          }
   }

   else if (rThisVariable==PRECONSOLIDATION) 
   {
       rValue = 0.0;
       const FlowRule::InternalVariables& InternalVariables=mpFlowRule->GetInternalVariables();
       double Alpha = InternalVariables.EquivalentPlasticStrain;
       rValue = mpHardeningLaw->CalculateHardening(rValue, Alpha);
   }
   else if ( rThisVariable==M_MODULUS)
   {
      double YoungModulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
      double PoissonCoef  = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];
      double K = YoungModulus/ 3.0 / (1.0 - 2.0*PoissonCoef);
      double G = YoungModulus/ 2.0 / ( 1.0 + PoissonCoef);
      rValue = K + 4.0*G / 3.0;
   }
   else if ( rThisVariable==SIMILAR_YOUNG_MODULUS)
   {
      rValue = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
   }
   else {
       rValue = HyperElasticPlastic3DLaw::GetValue( rThisVariable, rValue);
    }
   // THING ABOUT ALL THE VARIABLES.
   return (rValue);

}


void NonLinearHenckyElasticPlastic3DLaw::SetValue(const Variable<double>& rThisVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo)
{
   HyperElasticPlastic3DLaw::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

void NonLinearHenckyElasticPlastic3DLaw::SetValue(const Variable<Vector>& rThisVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
{

   if ( rThisVariable == ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS) {
      // SETS THE VALUE OF THE ELASTIC LEFT CAUCHY GREEN FROM A KIRCHHOFF STRESS VECTOR

      Matrix StressMat = ZeroMatrix(3);
      for (int i = 0; i < 3; ++i)
         StressMat(i,i) = rValue(i);

      StressMat(0, 1) = rValue(3);
      StressMat(1, 0) = rValue(3);
      StressMat(2, 0) = rValue(4);
      StressMat(0, 2) = rValue(4);
      StressMat(1, 2) = rValue(5);  // VALE, està malament
      StressMat(2, 1) = rValue(5);  // VALE, està malament

      Vector EigenStress;
      Matrix  EigenV;
      SolidMechanicsMathUtilities<double>::EigenVectors( StressMat, EigenV, EigenStress);


      Vector ElasticHenckyStrain(3);
      Matrix InverseElastic(3,3);
      double YoungModulus = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
      double PoissonCoef  = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];


      for (unsigned int i = 0; i < 3; ++i) {
         for (unsigned int j = 0; j < 3; ++j) {
            if (i == j ) {
               InverseElastic(i,i) = 1.0/YoungModulus;
            }
            else {
               InverseElastic(i,j) = - PoissonCoef / YoungModulus;
            }
         }
      }

      ElasticHenckyStrain = prod( InverseElastic, EigenStress);

      
      mElasticLeftCauchyGreen = ZeroMatrix(3);
      for (unsigned int i = 0; i < 3; ++i) {
         mElasticLeftCauchyGreen(i,i) = std::exp(2.0*ElasticHenckyStrain(i));
      }

      mElasticLeftCauchyGreen = prod( trans( EigenV), mElasticLeftCauchyGreen);
      mElasticLeftCauchyGreen = prod( mElasticLeftCauchyGreen, (EigenV) );

   }
   else if ( rThisVariable == ELASTIC_LEFT_CAUCHY_GREEN_VECTOR)
   {
      Matrix ElasticMatrix = ZeroMatrix(3);
      for (int i = 0; i < 3; i++)
         ElasticMatrix(i,i) = rValue(i);
      ElasticMatrix(0,1) = rValue(3);
      ElasticMatrix(1,0) = rValue(3);
      ElasticMatrix(0,2) = rValue(4);
      ElasticMatrix(2,0) = rValue(4);
      ElasticMatrix(2,1) = rValue(5);
      ElasticMatrix(1,2) = rValue(5);

      mElasticLeftCauchyGreen = ElasticMatrix;

   }
   else {
      HyperElasticPlastic3DLaw::SetValue( rThisVariable, rValue, rCurrentProcessInfo);   
   }
}

Matrix NonLinearHenckyElasticPlastic3DLaw::SetConstitutiveMatrixToAppropiateDimension(const Matrix& rElastoPlasticTangentMatrix)
{
     return rElastoPlasticTangentMatrix;

}

void NonLinearHenckyElasticPlastic3DLaw::CorrectDomainPressure( Matrix& rStressMatrix, const MaterialResponseVariables & rElasticVariables)
{

}


void NonLinearHenckyElasticPlastic3DLaw::CalculateElastoPlasticTangentMatrix( const FlowRule::RadialReturnVariables & rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& rAlpha, Matrix& rElastoPlasticTangentMatrix, const MaterialResponseVariables& rElasticVariables )
{

     mpFlowRule->ComputeElastoPlasticTangentMatrix( rReturnMappingVariables,  rNewElasticLeftCauchyGreen, rAlpha, rElastoPlasticTangentMatrix);

}
    
void NonLinearHenckyElasticPlastic3DLaw::CalculateOnlyDeviatoricPart( Matrix& rIncrementalDeformationGradient)
{

}
//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void NonLinearHenckyElasticPlastic3DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( FINITE_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);
	
	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

}


} // namespace Kratos

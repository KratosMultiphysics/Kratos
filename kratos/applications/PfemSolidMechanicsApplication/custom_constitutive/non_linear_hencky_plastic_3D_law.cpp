// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hyperelastic_plastic_3D_law.hpp"

#include "../PfemSolidMechanicsApplication/custom_constitutive/non_linear_hencky_plastic_3D_law.hpp"

#include "solid_mechanics_application.h"
#include "pfem_solid_mechanics_application.h"

//#include <iostream>


namespace Kratos
{

NonLinearHenckyElasticPlastic3DLaw::NonLinearHenckyElasticPlastic3DLaw()
   : HyperElasticPlastic3DLaw()
{

}

NonLinearHenckyElasticPlastic3DLaw::NonLinearHenckyElasticPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
   : HyperElasticPlastic3DLaw( pFlowRule, pYieldCriterion, pHardeningLaw)
{

}


NonLinearHenckyElasticPlastic3DLaw::NonLinearHenckyElasticPlastic3DLaw(const NonLinearHenckyElasticPlastic3DLaw&  rOther)
  : HyperElasticPlastic3DLaw(rOther)
{

}

 
NonLinearHenckyElasticPlastic3DLaw::~NonLinearHenckyElasticPlastic3DLaw()
{
}



void NonLinearHenckyElasticPlastic3DLaw::InitializeMaterial(const Properties& rProps, 
	const GeometryType& rGeom, const Vector& rShapeFunctionsValues)
{
   mElasticLeftCauchyGreen = identity_matrix<double> (3);
   mpFlowRule->InitializeMaterial(mpYieldCriterion, mpHardeningLaw, rProps);
}


void NonLinearHenckyElasticPlastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const ProcessInfo&  CurProcessInfo    = rValues.GetProcessInfo();
    // const Properties& MaterialProperties  = rValues.GetMaterialProperties(); not used JM
    const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();
    const double&   DeterminantF          = rValues.GetDeterminantF();

    const GeometryType&  DomainGeometry   = rValues.GetElementGeometry ();
    const Vector&        ShapeFunctions   = rValues.GetShapeFunctionsValues ();

    Vector& StrainVector                  = rValues.GetStrainVector();

    double& DeterminantF0                 = rValues.GetDeterminantF0();

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

    const Matrix&  DeformationGradientF0P= rValues.GetDeformationGradientF0();

    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    ElasticVariables.IdentityMatrix = identity_matrix<double> ( 3 );

    ElasticVariables.SetElementGeometry(DomainGeometry);
    ElasticVariables.SetShapeFunctionsValues(ShapeFunctions);

    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.clear();
    ReturnMappingVariables.DeltaTime = CurProcessInfo[DELTA_TIME];

    MatrixSplit SplitConstitutiveMatrix;

    //1.- Lame constants
    //const double& YoungModulus       = MaterialProperties[YOUNG_MODULUS];
    //const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    // const double YoungModulus = 2.069e5; not used JM
    // const double PoissonCoefficient = 0.3; not used JM

    //ReturnMappingVariables.LameLanda      = (YoungModulus*PoissonCoefficient)/((1.0+PoissonCoefficient)*(1.0-2.0*PoissonCoefficient));
    //ReturnMappingVariables.LameMu          =  YoungModulus/(2.0*(1.0+PoissonCoefficient));
    
    //2.-Determinant of the Total Deformation Gradient
    ElasticVariables.DeterminantF0 = DeterminantF0 * DeterminantF;

    //3.-Compute DeformationGradient (in 3D)
    
    Matrix DeformationGradientFbar = DeformationGradientF;
    DeformationGradientFbar = DeformationGradient3D(DeformationGradientFbar);

    Matrix DeformationGradientF0 = DeformationGradientF0P;
    DeformationGradientF0 = DeformationGradient3D(DeformationGradientF0);

    //4.-Left Cauchy-Green tensor b (without bar) to the new configuration
    Matrix IncrementalDeformationGradient;
//    IncrementalDeformationGradient = prod(DeformationGradientFbar, trans(DeformationGradientFbar));
    IncrementalDeformationGradient = DeformationGradientFbar;
    

/*    //4.-Left Cauchy-Green tensor b_bar to the new configuration
    ElasticVariables.CauchyGreenMatrix = prod(mElasticLeftCauchyGreen,trans(DeformationGradientFbar));
    ElasticVariables.CauchyGreenMatrix = prod(DeformationGradientFbar,ElasticVariables.CauchyGreenMatrix);
    //5.-Calculate trace of Left Cauchy-Green tensor b_bar
    ElasticVariables.traceCG = 0;
    for( unsigned int i=0; i<3; i++)
    {
       ElasticVariables.traceCG += ElasticVariables.CauchyGreenMatrix( i , i );
    }

    ReturnMappingVariables.LameMu_bar = ElasticVariables.LameMu * ( ElasticVariables.traceCG / 3.0  );
*/
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


    //Matrix IsochoricStressMatrix = ZeroMatrix(3);

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
         this->CalculateOnlyDeviatoricPart( IncrementalDeformationGradient );
         mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, IncrementalDeformationGradient, StressMatrix, NewElasticLeftCauchyGreen);
         this->CorrectDomainPressure( StressMatrix, ElasticVariables);

    }
    //OPTION 1:
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {

        //Kirchhoff Stress:
        StressVector = MathUtils<double>::StressTensorToVector(StressMatrix, StressVector.size());

        if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) )
        {
//            StressVector = SplitStressVector.Isochoric;
        }
        else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) )
        {
//            StressVector = SplitStressVector.Volumetric;
        }

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

        if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) )
        {
//            ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric + SplitConstitutiveMatrix.Plastic;
        }
        else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) )
        {
//            ConstitutiveMatrix = SplitConstitutiveMatrix.Volumetric;
        }

    }



    if( Options.Is( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ) )
    {
      mpFlowRule->UpdateInternalVariables ( ReturnMappingVariables );

      mElasticLeftCauchyGreen = NewElasticLeftCauchyGreen;

      //IncrementalDeformationGradient = prod(DeformationGradientFbar, DeformationGradientF0);
      //std::cout << " DEBUG " << IncrementalDeformationGradient(1,1)-1.0 <<" " <<  NewElasticLeftCauchyGreen(1,1)-1.0 << " " << prod(DeformationGradientFbar, DeformationGradientF0) << " STRESS " << StressVector << " rAlpha " << ReturnMappingVariables.DeltaGamma << std::endl;
      //std::cout << " FINALIZE: Stress " << StressVector << " rAlphA " << ReturnMappingVariables.DeltaGamma << std::endl; 
    }

}


double& NonLinearHenckyElasticPlastic3DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
{

   /* POST-PROCESS LOTS OF VARIABLES. ONLY IF DEFINED IN THE APPROPIATE FILE
   if ( (rThisVariable==PLASTIC_VOL) || (rThisVariable==PLASTIC_SHEAR) || (rThisVariable==PLASTIC_STRAIN) )
   {

       const FlowRule::InternalVariables& InternalVariables = mpFlowRule->GetInternalVariables();
       if ( rThisVariable==PLASTIC_VOL) {
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

   else if ( (rThisVariable==STRESS_P) || (rThisVariable==STRESS_Q) || (rThisVariable==STRESS_THETA) || (rThisVariable==STRESS_RATIO ) ) 
   {
         Matrix StressMatrix;
         Matrix NewElasticLeftCauchyGreen = mElasticLeftCauchyGreen;
         Matrix DeformationGradientF0 = ZeroMatrix(3);
         for (unsigned int i = 0; i < 3; ++i)
              DeformationGradientF0(i,i) = 1.0;
         Matrix IncrementalDeformationGradient = DeformationGradientF0;
 
         FlowRule::RadialReturnVariables ReturnMappingVariables;
       
         mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, IncrementalDeformationGradient, StressMatrix, NewElasticLeftCauchyGreen);

         double MeanStress = 0.0;
         for (unsigned int i = 0; i < 3; ++i)
              MeanStress += StressMatrix(i,i)/3.0;

         if (rThisVariable==STRESS_P) {
              rValue = -MeanStress;
              return rValue;
         }

         double StressQ = 0.0;
         for (unsigned int i = 0; i <3; ++i) 
              StressQ += pow( StressMatrix(i,i) - MeanStress, 2.0);

          StressQ += 2.0*pow( StressMatrix(0,1) , 2.0);
          StressQ += 2.0*pow( StressMatrix(0,2) , 2.0);
          StressQ += 2.0*pow( StressMatrix(1,2) , 2.0);
    
          if (rThisVariable== STRESS_Q)  {

              rValue = pow( 3.0/2.0*StressQ, 1.0/2.0);
              return rValue; 
          }
          if ( rThisVariable == STRESS_RATIO) {
              //if ( ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION) ) {
                  rValue = pow( 3.0 / 2.0 * StressQ, 1.0/2.0) / (-MeanStress);
                  return rValue;
             // / * }
             // else {
             //     rValue = 0.0;
             //     return rValue;
             // } * /
          } 
          StressQ /= 2.0;

          for (unsigned int i = 0; i < 3 ; ++i ) 
              StressMatrix(i,i) -= MeanStress;

          double ThirdInvariant = 0.0;
     
          ThirdInvariant = MathUtils<double>::Det( StressMatrix );
 
          ThirdInvariant *= 3.0*pow( 3.0, 1.0/2.0) / 2.0;
          ThirdInvariant /= pow( StressQ, 3.0/2.0);

          ThirdInvariant = std::asin( ThirdInvariant);
          rValue = ThirdInvariant / 3.0 * 180 / 3.14159265359; 

   }

   else if (rThisVariable==PRECONSOLIDATION) 
   {
       rValue = 0.0;
       const FlowRule::InternalVariables& InternalVariables=mpFlowRule->GetInternalVariables();
       double Alpha = InternalVariables.EquivalentPlasticStrain;
       rValue = mpHardeningLaw->CalculateHardening(rValue, Alpha);
   } 
   else {
       rValue = HyperElasticPlastic3DLaw::GetValue( rThisVariable, rValue);
    }
   */ // THING ABOUT ALL THE VARIABLES.
   rValue = HyperElasticPlastic3DLaw::GetValue( rThisVariable, rValue);
   return (rValue);

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


} // namespace Kratos

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hyperelastic_plastic_3D_law.hpp"

#include "../PfemSolidMechanicsApplication/custom_constitutive/hencky_plastic_3d_law.hpp"

#include "solid_mechanics_application.h"
//Molt important, el tema de constructors... etc
namespace Kratos
{

HenckyElasticPlastic3DLaw::HenckyElasticPlastic3DLaw()
   : HyperElasticPlastic3DLaw()
{

}

HenckyElasticPlastic3DLaw::HenckyElasticPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
   : HyperElasticPlastic3DLaw( pFlowRule, pYieldCriterion, pHardeningLaw)
{

}


HenckyElasticPlastic3DLaw::HenckyElasticPlastic3DLaw(const HenckyElasticPlastic3DLaw&  rOther)
  : HyperElasticPlastic3DLaw(rOther)
{

}

 
HenckyElasticPlastic3DLaw::~HenckyElasticPlastic3DLaw()
{
}



void HenckyElasticPlastic3DLaw::InitializeMaterial(const Properties& rProps, 
	const GeometryType& rGeom, const Vector& rShapeFunctionsValues)
{
   mElasticLeftCauchyGreen = identity_matrix<double> (3);
   mpFlowRule->InitializeMaterial(mpYieldCriterion, mpHardeningLaw, rProps);
}


void HenckyElasticPlastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const ProcessInfo&  CurProcessInfo    = rValues.GetProcessInfo();
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();
    const double&   DeterminantF          = rValues.GetDeterminantF();

    const GeometryType&  DomainGeometry   = rValues.GetElementGeometry ();
    const Vector&        ShapeFunctions   = rValues.GetShapeFunctionsValues ();

    Vector& StrainVector                  = rValues.GetStrainVector();

    double& DeterminantF0                 = rValues.GetDeterminantF0();

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    ElasticVariables.IdentityMatrix = identity_matrix<double> ( 3 );

    ElasticVariables.SetElementGeometry(DomainGeometry);
    ElasticVariables.SetShapeFunctionsValues(ShapeFunctions);

    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.clear();
    ReturnMappingVariables.DeltaTime = CurProcessInfo[DELTA_TIME];

    // Initialize Splited Parts: Isochoric and Volumetric stresses and constitutive tensors
    double voigtsize = StressVector.size();
 //   VectorSplit SplitStressVector;
    MatrixSplit SplitConstitutiveMatrix;

    //1.- Lame constants
    const double& YoungModulus       = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    //const double YoungModulus = 2.0e5;
    //const double PoissonCoefficient = 0.3;

//    ReturnMappingVariables.LameLanda      = (YoungModulus*PoissonCoefficient)/((1.0+PoissonCoefficient)*(1.0-2.0*PoissonCoefficient));
//    ReturnMappingVariables.LameMu          =  YoungModulus/(2.0*(1.0+PoissonCoefficient));
   
    ReturnMappingVariables.YoungModulus = YoungModulus;
    ReturnMappingVariables.PoissonCoefficient    = PoissonCoefficient; 
    //2.-Determinant of the Total Deformation Gradient
    ElasticVariables.DeterminantF0 = DeterminantF0 * DeterminantF;

    //3.-Compute DeformationGradient (in 3D)
    
    Matrix DeformationGradientFbar = DeformationGradientF;
    DeformationGradientFbar = DeformationGradient3D(DeformationGradientFbar);


    //4.-Left Cauchy-Green tensor b (without bar) to the new configuration
    ElasticVariables.CauchyGreenMatrix = prod(mElasticLeftCauchyGreen,trans(DeformationGradientFbar));
    ElasticVariables.CauchyGreenMatrix = prod(DeformationGradientFbar,ElasticVariables.CauchyGreenMatrix);

    if (false) {

       Matrix f = ZeroMatrix(3);

       for (unsigned int i = 0; i <3; ++i)
           f(i,i) = 1.0; 

        int nPassos = 500;
        Matrix DeltaDef = ZeroMatrix(3);
        Matrix StressMatrix;
        Vector PrincipalStrainTrial = ZeroVector(3);
        for (int i = 0; i < nPassos; ++ i ) {
            std::cout << " i " << i << std::endl;
            if (i < 50) {
               for (unsigned int j = 0; j < 3; ++ j) {
                  DeltaDef(j,j) = 0.999;
                }
            }
            else {
               for (unsigned int j = 0; j < 3; ++j) {
                  DeltaDef(j,j) = 1.0;
               }
               DeltaDef(0,1) = 0.001;
            }
            std::cout << " " << DeltaDef << std::endl;

            ElasticVariables.CauchyGreenMatrix = prod( mElasticLeftCauchyGreen, trans(DeltaDef));
            ElasticVariables.CauchyGreenMatrix = prod( DeltaDef, ElasticVariables.CauchyGreenMatrix);

            std::cout << " " << ElasticVariables.CauchyGreenMatrix << std::endl;
            this->CalculatePrincipalAxisHenckyTrial(ElasticVariables.CauchyGreenMatrix, ReturnMappingVariables, PrincipalStrainTrial);
      
            //2. compute trial stress and perform return mapping
            //returnmappingvariables.control.plasticregion = mpflowrule->calculatereturnmapping( returnmappingvariables, stressmatrix, principalstraintrial);
            mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, StressMatrix, PrincipalStrainTrial);

            mpFlowRule->UpdateInternalVariables ( ReturnMappingVariables );

            mElasticLeftCauchyGreen = mpFlowRule->GetElasticLeftCauchyGreen(ReturnMappingVariables);
            f = prod( DeltaDef, f);
            std::cout << "DEBUG stre " << StressMatrix << std::endl;
            std::cout << "DEBUG elcg " << mElasticLeftCauchyGreen << std::endl;
            std::cout << "DEBUG fff  " << f << std::endl;
 
        }

        KRATOS_ERROR( std::logic_error, "finishing the constitutive test debug zone blah blah blah", "" );
    }






    //4.-Almansi Strain:
    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
    {
        // e= 0.5*(1-invbT*invb)
        this->CalculateAlmansiStrain(ElasticVariables.CauchyGreenMatrix,StrainVector);
    }

 
    //5.-Calculate Total Kirchhoff stress
    Matrix StressMatrix = ZeroMatrix(3);    

    //Matrix IsochoricStressMatrix = ZeroMatrix(3);

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
      Vector PrincipalStrainTrial= ZeroVector(3);
      
      //2. Compute Trial Stress and perform return mapping

      this->CalculatePrincipalAxisHenckyTrial(ElasticVariables.CauchyGreenMatrix, ReturnMappingVariables, PrincipalStrainTrial);
      
      mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, StressMatrix, PrincipalStrainTrial);
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

        //ElasticVariables.CauchyGreenMatrix = ElasticVariables.IdentityMatrix;

	Matrix PrincipalTangentMatrix;
	if (ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION)  ) {
	   mpFlowRule->CalculatePrincipalAxisAlgorithmicTangent(ReturnMappingVariables, StressMatrix, PrincipalTangentMatrix);
        }
        else  {
           mpFlowRule->CalculatePrincipalAxisAlgorithmicTangent(ReturnMappingVariables, StressMatrix, PrincipalTangentMatrix);
        }

        this->ConvertConstitutiveMatrixToAppropiateDimension( PrincipalTangentMatrix);
        ConstitutiveMatrix = PrincipalTangentMatrix; 

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

      mElasticLeftCauchyGreen = mpFlowRule->GetElasticLeftCauchyGreen(ReturnMappingVariables);
    }

}


void HenckyElasticPlastic3DLaw::ConvertConstitutiveMatrixToAppropiateDimension(Matrix&  PrincipalTangentMatrix)
{

}



Vector& HenckyElasticPlastic3DLaw::GetStressVectorFromMatrix(const Matrix& rStressMatrix, Vector& rPrincipalStress, const Matrix& rEigenVectors)
{
//   EigenVectors = mpFlowRule->GetEigenVectors();
//Ll:
   Matrix auxMatrix = ZeroMatrix(3);
   auxMatrix = prod( rStressMatrix, trans(rEigenVectors));
   auxMatrix = prod( (rEigenVectors), auxMatrix);
   rPrincipalStress = ZeroVector(3);
   for (unsigned int i = 0; i<3; i++)
      rPrincipalStress(i) = auxMatrix(i,i);

   return rPrincipalStress;
}

void HenckyElasticPlastic3DLaw::FinalizeMaterialResponseKirchhoff (Parameters & rValues)
{
  HyperElasticPlastic3DLaw::FinalizeMaterialResponseKirchhoff(rValues);
}



void HenckyElasticPlastic3DLaw::CalculatePrincipalAxisHenckyTrial(const Matrix& rCauchyGreenMatrix, FlowRule::RadialReturnVariables& rReturnMappingVariables, Vector& rPrincipalStrain)
{
    //Vector rPrincipalStrain = ZeroVector(3);
    rReturnMappingVariables.EigenVectors = ZeroMatrix(3);
    rReturnMappingVariables.TrialEigenValues = ZeroVector(3);
    SolidMechanicsMathUtilities<double>::EigenVectors(rCauchyGreenMatrix, rReturnMappingVariables.EigenVectors, rReturnMappingVariables.TrialEigenValues);
    for (unsigned int i = 0; i<3; ++i)
           rPrincipalStrain(i) = 0.50*std::log(rReturnMappingVariables.TrialEigenValues(i));
}


} // namespace Kratos

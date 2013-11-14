// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hyperelastic_plastic_3D_law.hpp"

#include "../PfemSolidMechanicsApplication/custom_constitutive/non_linear_hencky_plastic_3d_law.hpp"

#include "solid_mechanics_application.h"
//Molt important, el tema de constructors... etc
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
    //const double& YoungModulus       = MaterialProperties[YOUNG_MODULUS];
    //const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    const double YoungModulus = 2.069e5;
    const double PoissonCoefficient = 0.3;

    ReturnMappingVariables.LameLanda      = (YoungModulus*PoissonCoefficient)/((1.0+PoissonCoefficient)*(1.0-2.0*PoissonCoefficient));
    ReturnMappingVariables.LameMu          =  YoungModulus/(2.0*(1.0+PoissonCoefficient));
    
    //2.-Determinant of the Total Deformation Gradient
    ElasticVariables.DeterminantF0 = DeterminantF0 * DeterminantF;

    //3.-Compute DeformationGradient (in 3D)
    
    Matrix DeformationGradientFbar = DeformationGradientF;
    DeformationGradientFbar = DeformationGradient3D(DeformationGradientFbar);


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
        //ReturnMappingVariables.Control.PlasticRegion = mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, IncrementalDeformationGradient, StressMatrix, NewElasticLeftCauchyGreen);
         mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, IncrementalDeformationGradient, StressMatrix, NewElasticLeftCauchyGreen);

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
        double rAlpha = 0.0;
        
        mpFlowRule->ComputeElastoPlasticTangentMatrix( ReturnMappingVariables,  NewElasticLeftCauchyGreen, rAlpha, ElastoPlasticTangentMatrix);



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
       
 
    }

}


Matrix NonLinearHenckyElasticPlastic3DLaw::SetConstitutiveMatrixToAppropiateDimension(const Matrix& rElastoPlasticTangentMatrix)
{
     return rElastoPlasticTangentMatrix;

}
} // namespace Kratos

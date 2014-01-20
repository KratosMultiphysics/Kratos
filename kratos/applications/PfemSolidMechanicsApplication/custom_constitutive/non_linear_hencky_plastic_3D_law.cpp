// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hyperelastic_plastic_3D_law.hpp"

#include "../PfemSolidMechanicsApplication/custom_constitutive/non_linear_hencky_plastic_3D_law.hpp"
#include "custom_constitutive/custom_flow_rules/non_associative_explicit_flow_rule.hpp"

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

    // Initialize Splited Parts: Isochoric and Volumetric stresses and constitutive tensors
    // double voigtsize = StressVector.size(); not used JM
 //   VectorSplit SplitStressVector;
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


// ONE POINT GAUSS INTEGRATION TO CHECK THE CONSTITUTIVE EQUATION
// DEBUGG ZONE FOR THE CONSTITUTIVE EQUATIONS
if (false) {

   unsigned int nPassos = 500;
   Vector ThisStressVector;
   Matrix IncrementalDeformationGradient;
   DeformationGradientF0 = ZeroMatrix(3);
   for (unsigned int j = 0 ; j < 3; ++j)
      DeformationGradientF0(j,j) = 1.0;
   DeformationGradientF0(0,0) -= 1e-12;
   NewElasticLeftCauchyGreen = prod(DeformationGradientF0, trans(DeformationGradientF0));

   Matrix RotationMatrix = ZeroMatrix(3);
   RotationMatrix(2,2) = 1.0;
   for (unsigned j = 0; j < 3; ++j)
      RotationMatrix(j,j) = 1.0;
   /*RotationMatrix(0,0) = cos(1.0/ double(nPassos) * 3.14159);
   RotationMatrix(1,1) = RotationMatrix(0,0);
   RotationMatrix(0,1) = - pow(  1- RotationMatrix(0,0)*RotationMatrix(0,0), 0.5);
   RotationMatrix(1,0) = - RotationMatrix(0,1);
*/
   for (unsigned int i = 0; i < nPassos; ++i) {

       // SETTING THE NEW DEFORMATION
       IncrementalDeformationGradient = ZeroMatrix(3);
       for (unsigned int j = 0; j < 3; ++j) 
           IncrementalDeformationGradient(j,j) = 1.0;

       double Inc = (double) rand()/ (double) RAND_MAX;
       IncrementalDeformationGradient(1,1) = 1.0 - 0.01*(Inc - 0.3);
       std::cout << "Hello"  <<  IncrementalDeformationGradient(1,1) << std::endl;;
//       IncrementalDeformationGradient(1,1) = 0.999;
       IncrementalDeformationGradient = prod(RotationMatrix, IncrementalDeformationGradient);
       //COMPUTING
       mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, DeformationGradientF0, IncrementalDeformationGradient, StressMatrix, NewElasticLeftCauchyGreen);

       // FINALIZING STEP
       mpFlowRule->UpdateInternalVariables( ReturnMappingVariables);
       mElasticLeftCauchyGreen = NewElasticLeftCauchyGreen;
       ThisStressVector = MathUtils<double>::StressTensorToVector(StressMatrix, 6);
       DeformationGradientF0 = prod( IncrementalDeformationGradient, DeformationGradientF0);

        // WRITTING
        std::cout << "DebugConst Results " << i << std::endl;
        std::cout << "DebugConst F0   " << DeformationGradientF0 << std::endl;
        std::cout << "DebugConst ELCG " << mElasticLeftCauchyGreen << std::endl;
        std::cout << "DebugConst INCR " << IncrementalDeformationGradient << std::endl;
        std::cout << "DebugConst STRE "<< ThisStressVector << std::endl;
        std::cout << " " << std::endl;
    }  
    KRATOS_ERROR( std::logic_error, "FINISHING THE CONSTITUTIVE TEST DEBUG ZONE BLAH BLAH BLAH", "" );

}


    //Matrix IsochoricStressMatrix = ZeroMatrix(3);

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        //ReturnMappingVariables.Control.PlasticRegion = mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, IncrementalDeformationGradient, StressMatrix, NewElasticLeftCauchyGreen);
         mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, DeformationGradientF0, IncrementalDeformationGradient, StressMatrix, NewElasticLeftCauchyGreen);

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

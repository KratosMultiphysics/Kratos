//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
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
#include "custom_constitutive/hencky_plastic_3d_law.hpp"

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************
HenckyElasticPlastic3DLaw::HenckyElasticPlastic3DLaw()
   : HyperElasticPlastic3DLaw()
{

}



HenckyElasticPlastic3DLaw::HenckyElasticPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
   : HyperElasticPlastic3DLaw( pFlowRule, pYieldCriterion, pHardeningLaw)
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyElasticPlastic3DLaw::HenckyElasticPlastic3DLaw(const HenckyElasticPlastic3DLaw&  rOther)
  : HyperElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyElasticPlastic3DLaw::Clone() const
{
    HenckyElasticPlastic3DLaw::Pointer p_clone(new HenckyElasticPlastic3DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlastic3DLaw::~HenckyElasticPlastic3DLaw()
{
}


//************************************************************************************
//************************************************************************************


void HenckyElasticPlastic3DLaw::InitializeMaterial(const Properties& rProps, 
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


void HenckyElasticPlastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const ProcessInfo& CurrrentProcessInfo = rValues.GetProcessInfo();
    const Properties& MaterialProperties   = rValues.GetMaterialProperties();

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
    ReturnMappingVariables.DeltaTime = CurProcessInfo[DELTA_TIME];

    if(CurrentProcessInfo[IMPLEX] == 1)	
      ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,true);
    else
      ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,false);


    //1.- Lame constants
    // const double& YoungModulus       = MaterialProperties[YOUNG_MODULUS];
    // const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];

    //ElasticVariables.LameLanda      = (YoungModulus*PoissonCoefficient)/((1.0+PoissonCoefficient)*(1.0-2.0*PoissonCoefficient));
    //ElasticVariables.LameMu         =  YoungModulus/(2.0*(1.0+PoissonCoefficient));
   
    // lluis
    // ReturnMappingVariables.YoungModulus          =  MaterialProperties[YOUNG_MODULUS];
    // ReturnMappingVariables.PoissonCoefficient    =  MaterialProperties[POISSON_RATIO]; 


    //2.-Determinant of the Total Deformation Gradient
    ElasticVariables.DeterminantF = DeterminantF;

    //3.-Compute Incremental DeformationGradient (in 3D)
    ElasticVariables.DeformationGradientF = DeformationGradientF;
    ElasticVariables.DeformationGradientF = Transform2DTo3D(ElasticVariables.DeformationGradientF);
    ElasticVariables.DeformationGradientF = prod(ElasticVariables.DeformationGradientF,mInverseDeformationGradientF0); 

    //4.-Left Cauchy Green tensor b: (stored in the CauchyGreenMatrix)
    ElasticVariables.CauchyGreenMatrix = prod(mElasticLeftCauchyGreen,trans(ElasticVariables.DeformationGradientF));
    ElasticVariables.CauchyGreenMatrix = prod(ElasticVariables.DeformationGradientF,ElasticVariables.CauchyGreenMatrix);


    //TODO: supress this testing zone -------------------
    if (false) {

       Matrix f = ZeroMatrix(3);

       for (unsigned int i = 0; i <3; ++i)
           f(i,i) = 1.0; 

        int nPassos = 500;
        Matrix DeltaDef = ZeroMatrix(3);
        Matrix StressMatrix;
        Vector MainStrainTrial = ZeroVector(3);
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
            this->CalculateHenckyMainStrain(ElasticVariables.CauchyGreenMatrix, ReturnMappingVariables, MainStrainTrial);
      
            //2. compute trial stress and perform return mapping
            //returnmappingvariables.control.plasticregion = mpflowrule->calculatereturnmapping( returnmappingvariables, stressmatrix, principalstraintrial);`
            //mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, StressMatrix, MainStrainTrial);
	    Matrix NewElasticLeftCauchyGreen = ElasticVariables.CauchyGreenMatrix;
	    mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, ElasticVariables.DeformationGradientF, StressMatrix, NewElasticLeftCauchyGreen);

            mpFlowRule->UpdateInternalVariables ( ReturnMappingVariables );

            mElasticLeftCauchyGreen = mpFlowRule->GetElasticLeftCauchyGreen(ReturnMappingVariables);
            f = prod( DeltaDef, f);
            std::cout << "DEBUG stre " << StressMatrix << std::endl;
            std::cout << "DEBUG elcg " << mElasticLeftCauchyGreen << std::endl;
            std::cout << "DEBUG fff  " << f << std::endl;
 
        }

        KRATOS_THROW_ERROR( std::logic_error, "finishing the constitutive test debug zone blah blah blah", "" );
    }
    //TODO: supress this testing zone -------------------


    //5.-Almansi Strain:
    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
    {
        // e= 0.5*(1-invbT*invb)
        this->CalculateAlmansiStrain(ElasticVariables.CauchyGreenMatrix, StrainVector);
    }

 
    //6.-Calculate Total Kirchhoff stress

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {

      Matrix StressMatrix     = ZeroMatrix(3);     
      Vector HenckyMainStrainVector = ZeroVector(3);

      this->CalculateHenckyMainStrain(ElasticVariables.CauchyGreenMatrix, ReturnMappingVariables, HenckyMainStrainVector);
      // mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, StressMatrix, HenckyMainStrain );

      Matrix HenckyMainStrainMatrix = ZeroMatrix(3);
      for (unsigned int i = 0; i<3; ++i)
	HenckyMainStrainMatrix(i,i) = HenckyMainStrainVector[i];

      mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, ElasticVariables.DeformationGradientF, StressMatrix, );

      StressVector = MathUtils<double>::StressTensorToVector(StressMatrix, StressVector.size());

    }


    //8.-Calculate Constitutive Matrix related to Total Kirchhoff stress
    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        //initialize constitutive tensors
        ConstitutiveMatrix.clear();

        //ElasticVariables.CauchyGreenMatrix = ElasticVariables.IdentityMatrix;

	if (ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION)  ) {
	   mpFlowRule->CalculatePrincipalAxisAlgorithmicTangent(ReturnMappingVariables, StressMatrix, ConstitutiveMatrix);
        }
        else  {
           mpFlowRule->CalculatePrincipalAxisAlgorithmicTangent(ReturnMappingVariables, StressMatrix, ConstitutiveMatrix);
        }

        this->SetConstitutiveMatrixToAppropiateDimension( ConstitutiveMatrix );

    }



    if( Options.Is( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ) )
    {
      mpFlowRule->UpdateInternalVariables ( ReturnMappingVariables );

      mElasticLeftCauchyGreen = mpFlowRule->GetElasticLeftCauchyGreen(ReturnMappingVariables);

      ElasticVariables.DeformationGradientF = DeformationGradientF;
      ElasticVariables.DeformationGradientF = Transform2DTo3D(ElasticVariables.DeformationGradientF);
      MathUtils<double>::InvertMatrix( ElasticVariables.DeformationGradientF, mInverseDeformationGradientF0, mDeterminantF0);
      mDeterminantF0 = DeterminantF; //special treatment of the determinant 
    }

}


//************************************************************************************
//************************************************************************************


void HenckyElasticPlastic3DLaw::SetConstitutiveMatrixToAppropiateDimension(Matrix&  rConstitutiveMatrix)
{


}


//************************************************************************************
//************************************************************************************


Vector& HenckyElasticPlastic3DLaw::GetStressVectorFromMatrix(const Matrix& rStressMatrix, 
							     Vector& rMainStress, 
							     const Matrix& rEigenVectors)
{
//   EigenVectors = mpFlowRule->GetEigenVectors();
//Ll:
   Matrix auxMatrix = ZeroMatrix(3);
   auxMatrix = prod( rStressMatrix, trans(rEigenVectors));
   auxMatrix = prod( (rEigenVectors), auxMatrix);

   rMainStress = ZeroVector(3);
   for (unsigned int i = 0; i<3; i++)
      rMainStress[i] = auxMatrix(i,i);

   return rMainStress;
}



//************************************************************************************
//************************************************************************************


  void HenckyElasticPlastic3DLaw::CalculateHenckyMainStrain(const Matrix& rCauchyGreenMatrix, 
							    FlowRule::RadialReturnVariables& rReturnMappingVariables, 
							    Vector& rMainStrain)
{
    Matrix EigenVectors  = ZeroMatrix(3);
    Vector EigenValues   = ZeroVector(3);
    SolidMechanicsMathUtilities<double>::EigenVectors(rCauchyGreenMatrix, EigenVectors, EigenValues);

    // lluis
    rReturnMappingVariables.MainDirections     = EigenVectors;
    // rReturnMappingVariables.TrialEigenValues = EigenValues;

    for (unsigned int i = 0; i<3; ++i)
      rMainStrain[i] = 0.50 * std::log(EigenValues[i]);
}


} // namespace Kratos

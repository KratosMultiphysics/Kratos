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
    ElasticVariables.CauchyGreenMatrix = prod(mElasticLeftCauchyGreen,trans(DeformationGradientFbar));
    ElasticVariables.CauchyGreenMatrix = prod(DeformationGradientFbar,ElasticVariables.CauchyGreenMatrix);
    

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
    // No només fa el Return Mapping sinó que també calcula el trial (i la descomposició). D'aquesta manera netegem coses.
    //1. Compute Trial Hencky Strian in the principal Axis (and polar decomposition) 
    //      Vector PrincipalStrainTrial = ZeroVector(3);
    //      this->CalculatePrincipalAxisHenckyTrial(ElasticVariables.CauchyGreenMatrix, ReturnMappingVariables, PrincipalStrainTrial); 

      Vector PrincipalStrainTrial= ZeroVector(3);
      this->CalculatePrincipalAxisHenckyTrial(ElasticVariables.CauchyGreenMatrix, ReturnMappingVariables, PrincipalStrainTrial);
      
      //2. Compute Trial Stress and perform return mapping
      //ReturnMappingVariables.Control.PlasticRegion = mpFlowRule->CalculateReturnMapping( ReturnMappingVariables, StressMatrix, PrincipalStrainTrial);
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
        SplitConstitutiveMatrix.EigenValues  = ConstitutiveMatrix;
        SplitConstitutiveMatrix.EigenVectors = ConstitutiveMatrix;

        ElasticVariables.CauchyGreenMatrix = ElasticVariables.IdentityMatrix;

	Matrix PrincipalTangentMatrix;
	if (ReturnMappingVariables.Options.Is(FlowRule::PLASTIC_REGION)  ) {
	   mpFlowRule->CalculatePrincipalAxisAlgorithmicTangent(ReturnMappingVariables, StressMatrix, PrincipalTangentMatrix);
        }
        else  {
           mpFlowRule->CalculatePrincipalAxisAlgorithmicTangent(ReturnMappingVariables, StressMatrix, PrincipalTangentMatrix);
	   //Ll: Poner Aquí la matrix elastica normal en ejes principales.
        }
	//Matrix& EigenVectors = mpFlowRule->GetEigenVectors();
        //Vector& EigenValues  = mpFlowRule->GetEigenValues();
        //Ll:
       Matrix& EigenVectors = ReturnMappingVariables.EigenVectors;
       Vector& EigenValues  = ReturnMappingVariables.TrialEigenValues;

	if ( Options.Is( ConstitutiveLaw::LAST_KNOWN_CONFIGURATION ) || Options.Is( ConstitutiveLaw::INITIAL_CONFIGURATION ) ){
	  
	  Matrix InverseDeformationGradientF ( 3, 3 );
	  double DetInvF=0;
	  MathUtils<double>::InvertMatrix( DeformationGradientFbar, InverseDeformationGradientF, DetInvF);

	  this->CalculateEigenValuesConstitutiveMatrix(PrincipalTangentMatrix, EigenVectors, InverseDeformationGradientF, SplitConstitutiveMatrix.EigenValues);

	  this->CalculateEigenVectorsConstitutiveMatrix(EigenVectors, EigenValues, StressMatrix, InverseDeformationGradientF, SplitConstitutiveMatrix.EigenVectors);

	}
	else{

	  this->CalculateEigenValuesConstitutiveMatrix(PrincipalTangentMatrix, EigenVectors, SplitConstitutiveMatrix.EigenValues);

	  this->CalculateEigenVectorsConstitutiveMatrix(EigenVectors, EigenValues, StressMatrix, SplitConstitutiveMatrix.EigenVectors);
	  
	  
	}


        //if( Options.Is(ConstitutiveLaw::TOTAL_TENSOR ) )
        ConstitutiveMatrix = SplitConstitutiveMatrix.EigenValues+ SplitConstitutiveMatrix.EigenVectors;
        ConstitutiveMatrix(2,2) = 100.0; 

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

      //mElasticLeftCauchyGreen  = ( IsochoricStressMatrix * ( 1.0 / ElasticVariables.LameMu ) );
      //mElasticLeftCauchyGreen += ( ElasticVariables.traceCG/3.0) * ElasticVariables.IdentityMatrix;
      mElasticLeftCauchyGreen = mpFlowRule->GetElasticLeftCauchyGreen(ReturnMappingVariables);
 
    }

}


void HenckyElasticPlastic3DLaw::CalculateEigenValuesConstitutiveMatrix(const Matrix & rPrincipalTangent, const Matrix & rEigenVectors, const Matrix& rInverseDeformationGradientF, Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    static const unsigned int msIndexVoigt3D [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
	  rConstitutiveMatrix( i, j ) = EigenValuesConstitutiveComponent(rConstitutiveMatrix( i, j ),
					  rPrincipalTangent, rEigenVectors, rInverseDeformationGradientF, 
                                          msIndexVoigt3D[i][0], msIndexVoigt3D[i][1], msIndexVoigt3D[j][0], msIndexVoigt3D[j][1]);
        }

    }


}


void HenckyElasticPlastic3DLaw::CalculateEigenVectorsConstitutiveMatrix(const Matrix& rEigenVectors, const Vector& rEigenValues, const Matrix& rStressMatrix, const Matrix& rInverseDeformationGradientF, Matrix& rConstitutiveMatrix)
{
    Vector PrincipalStress;
    PrincipalStress = this->GetStressVectorFromMatrix(rStressMatrix, PrincipalStress);
    
    rConstitutiveMatrix.clear();


    static const unsigned int msIndexVoigt3D [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = EigenVectorsConstitutiveComponent(rConstitutiveMatrix( i, j ), 
					  rEigenVectors, rEigenValues, PrincipalStress, rInverseDeformationGradientF, 
                                          msIndexVoigt3D[i][0], msIndexVoigt3D[i][1], msIndexVoigt3D[j][0], msIndexVoigt3D[j][1]);
        }

    }

}


double& HenckyElasticPlastic3DLaw::EigenValuesConstitutiveComponent(double & rCabcd,
	const Matrix & rPrincipalTangent, const Matrix& rEigenVectors, const Matrix&  rInverseDeformationGradientF, 
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{
    rCabcd = 0;
    double Cijkl=0;

    unsigned int dimension = rInverseDeformationGradientF.size1();

    //Cabcd
    for(unsigned int j=0; j<dimension; j++)
    {
        for(unsigned int l=0; l<dimension; l++)
        {
            for(unsigned int k=0; k<dimension; k++)
            {
                for(unsigned int i=0; i<dimension; i++)
                {
                    //Cijkl
		  rCabcd +=rInverseDeformationGradientF(a,i)*rInverseDeformationGradientF(b,j)*rInverseDeformationGradientF(c,k)*rInverseDeformationGradientF(d,l)*EigenValuesConstitutiveComponent(Cijkl,rPrincipalTangent,rEigenVectors,i,j,k,l);
                }
            }
        }
    }

    return rCabcd;
}


double& HenckyElasticPlastic3DLaw::EigenVectorsConstitutiveComponent(double & rCabcd,
	const Matrix & rEigenVectors, const Vector& rEigenValues, const Vector&  rPrincipalStress,
	const Matrix & rInverseDeformationGradientF,  
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{
    rCabcd = 0;
    double Cijkl=0;

    unsigned int dimension = rInverseDeformationGradientF.size1();

    //Cabcd
    for(unsigned int j=0; j<dimension; j++)
    {
        for(unsigned int l=0; l<dimension; l++)
        {
            for(unsigned int k=0; k<dimension; k++)
            {
                for(unsigned int i=0; i<dimension; i++)
                {
                    //Cijkl
		  rCabcd +=rInverseDeformationGradientF(a,i)*rInverseDeformationGradientF(b,j)*rInverseDeformationGradientF(c,k)*rInverseDeformationGradientF(d,l)*EigenVectorsConstitutiveComponent(Cijkl, rEigenVectors, rEigenValues, rPrincipalStress, i,j,k,l);

                }
            }
        }
    }

    return rCabcd;
}


void HenckyElasticPlastic3DLaw::CalculateEigenValuesConstitutiveMatrix(const Matrix & rPrincipalTangent, const Matrix & rEigenVectors, Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();


    static const unsigned int msIndexVoigt3D [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = EigenValuesConstitutiveComponent(rConstitutiveMatrix( i, j ), 
					  rPrincipalTangent, rEigenVectors, 
                                          msIndexVoigt3D[i][0], msIndexVoigt3D[i][1], msIndexVoigt3D[j][0], msIndexVoigt3D[j][1]);
        }

    }

}


void HenckyElasticPlastic3DLaw::CalculateEigenVectorsConstitutiveMatrix(const Matrix& rEigenVectors, const Vector& rEigenValues, const Matrix& rStressMatrix, Matrix& rConstitutiveMatrix)
{
    Vector PrincipalStress;
    PrincipalStress = this->GetStressVectorFromMatrix(rStressMatrix, PrincipalStress);
    
    rConstitutiveMatrix.clear();


    static const unsigned int msIndexVoigt3D [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = EigenVectorsConstitutiveComponent(rConstitutiveMatrix( i, j ), 
					  rEigenVectors, rEigenValues, PrincipalStress,
                                          msIndexVoigt3D[i][0], msIndexVoigt3D[i][1], msIndexVoigt3D[j][0], msIndexVoigt3D[j][1]);
        }

    }

}



double& HenckyElasticPlastic3DLaw::EigenVectorsConstitutiveComponent(double & rCabcd,
	const Matrix & rEigenVectors, const Vector& rEigenValues,const Vector&  rPrincipalStress, 
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{
	const Matrix& rEV = rEigenVectors;
        rCabcd = 0.0;
	for (unsigned int m = 0; m<3; ++m)
	{
	   for (unsigned int n = 0; n<3; ++n)
           {
	      if ( m == n) {
	      }
	      else {
                 if ( abs(rEigenValues(m)-rEigenValues(n))>10E-5)
                 {
		     rCabcd += (rPrincipalStress(m)-rPrincipalStress(n)) / ( rEigenValues(m)-rEigenValues(n)) * (rEigenValues(m)*rEV(a,m)*rEV(b,n)*rEV(c,m)*rEV(d,n) + rEigenValues(n)*rEV(a,m)*rEV(b,n)*rEV(c,n)*rEV(d,n));
	
	         }
              }
            }
	}
        return rCabcd;
                  
}

double& HenckyElasticPlastic3DLaw::EigenValuesConstitutiveComponent(double & rCabcd,
        const Matrix & rPrincipalTangent, const Matrix & rEigenVectors,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{
    rCabcd = 0.0;
    //Isochoric part of the hyperelastic constitutive tensor component
    for (unsigned int i =0; i<3; ++i) {
       for (unsigned int j = 0; j<3; ++j) {
          rCabcd += rPrincipalTangent(i,j)*rEigenVectors(a,i)*rEigenVectors(b,i)*rEigenVectors(c,j)*rEigenVectors(d,j);
       }
    }

    return rCabcd;
}





Vector& HenckyElasticPlastic3DLaw::GetStressVectorFromMatrix(const Matrix& rStressMatrix, Vector& rPrincipalStress)
{
   Matrix EigenVectors= ZeroMatrix(3);
//   EigenVectors = mpFlowRule->GetEigenVectors();
//Ll:
   Matrix auxMatrix = ZeroMatrix(3);
   auxMatrix = prod( rStressMatrix, EigenVectors);
   auxMatrix = prod( trans(EigenVectors), auxMatrix);
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
    SD_MathUtils<double>::EigenVectors(rCauchyGreenMatrix, rReturnMappingVariables.EigenVectors, rReturnMappingVariables.TrialEigenValues);
    for (unsigned int i = 0; i<3; ++i)
           rPrincipalStrain(i) = 0.50*std::log(rReturnMappingVariables.TrialEigenValues(i));
}


} // namespace Kratos

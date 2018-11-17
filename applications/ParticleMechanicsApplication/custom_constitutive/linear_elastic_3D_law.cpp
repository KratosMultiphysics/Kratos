//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//
//  References:      This class is adapted from applications/SolidMechanicsApplication/custom_constitutive/linear_elastic_3D_law.cpp


// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_elastic_3D_law.hpp"

#include "particle_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearElastic3DLaw::LinearElastic3DLaw()
    : HyperElastic3DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LinearElastic3DLaw::LinearElastic3DLaw(const LinearElastic3DLaw& rOther)
    : HyperElastic3DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearElastic3DLaw::Clone() const
{
    return Kratos::make_shared<LinearElastic3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LinearElastic3DLaw::~LinearElastic3DLaw()
{
}


//************************************************************************************
//************************************************************************************


//******************CALCULATE VALUE: DOUBLE - VECTOR - MATRIX*************************
//************************************************************************************

double& LinearElastic3DLaw::CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue )
{
    return (this->GetValue(rThisVariable,rValue ));
}


double& LinearElastic3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if (rThisVariable == STRAIN_ENERGY)
    {
        rValue = mStrainEnergy;
    }

    return( rValue );
}

//*****************************MATERIAL RESPONSES*************************************
//************************************************************************************

void  LinearElastic3DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
{

    //-----------------------------//

    // a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    // b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();
    mStrainEnergy = 0.0; //When it is not calculated, a zero will be returned

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();

    Vector& StrainVector                  = rValues.GetStrainVector();
    Vector& StressVector                  = rValues.GetStressVector();

    //-----------------------------//

    // 0.- Lame constants
    const double& YoungModulus       = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];

    if(Options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) //large strains
    {
        // 1.-Compute total deformation gradient
        const Matrix& DeformationGradientF = rValues.GetDeformationGradientF();

        // 2.-Right Cauchy-Green tensor C
        Matrix RightCauchyGreen = prod(trans(DeformationGradientF),DeformationGradientF);

        // 3.-Green-Lagrange Strain:
        // E= 0.5*(FT*F-1)
        this->CalculateGreenLagrangeStrain(RightCauchyGreen,StrainVector);

    }

    //4.-Calculate Total PK2 stress
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {
        if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
        {
            Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
            this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );
        }
        else
        {
            Matrix ConstitutiveMatrix(StrainVector.size(),StrainVector.size());
            noalias(ConstitutiveMatrix) = ZeroMatrix(StrainVector.size());
            this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
            this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );
        }

    }
    else if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
    }

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRAIN_ENERGY ) )
    {
        if( Options.IsNot( ConstitutiveLaw::COMPUTE_STRESS ) )
        {
            if(Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ))
            {
                Matrix ConstitutiveMatrix( StrainVector.size(), StrainVector.size() );
                noalias(ConstitutiveMatrix) = ZeroMatrix(StrainVector.size());

                this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
                this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );
            }
            else
            {
                Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
                this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );
            }
        }

        mStrainEnergy = 0.5 * inner_prod(StrainVector,StressVector); //Belytschko Nonlinear Finite Elements pag 226 (5.4.3) : w = 0.5*E:C:E
    }

}


//************************************************************************************
//************************************************************************************


void LinearElastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    //-----------------------------//

    // a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    // b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();

    Vector& StrainVector                  = rValues.GetStrainVector();
    Vector& StressVector                  = rValues.GetStressVector();

    //-----------------------------//

    // 0. Lame constants
    const double& YoungModulus          = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient    = MaterialProperties[POISSON_RATIO];

    if(Options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) //large strains
    {
	    // 1.-Compute total deformation gradient
        const Matrix& DeformationGradientF  = rValues.GetDeformationGradientF();
        const double& DeterminantF = rValues.GetDeterminantF();

        // 2.-Left Cauchy-Green tensor b
        Matrix LeftCauchyGreenMatrix = prod(DeformationGradientF,trans(DeformationGradientF));

        // 3.-Almansi Strain:
        // e= 0.5*(1-invFT*invF)
        this->CalculateAlmansiStrain(LeftCauchyGreenMatrix,StrainVector);

        // 4. Kirchhoff model is set with S = CE
        this->CalculateMaterialResponsePK2 (rValues);

        // 5.-Calculate Total Kirchhoff stress
        if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
        {
            TransformStresses(StressVector, DeformationGradientF, DeterminantF, StressMeasure_PK2, StressMeasure_Kirchhoff);
        }

        if( Options.Is( ConstitutiveLaw::COMPUTE_STRAIN_ENERGY ) )
        {
            mStrainEnergy *= DeterminantF;
        }

    }
    else //small strains
    {

        // 1.-Calculate total Kirchhoff stress

        if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){

            if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
            {
	            Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();

                this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
                this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );
            }
	        else
            {
                Matrix ConstitutiveMatrix( StrainVector.size() ,StrainVector.size());
                noalias(ConstitutiveMatrix) = ZeroMatrix(StrainVector.size());

                this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
                this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );

            }

        }
        else if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
        {
            Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
            this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
        }

        if( Options.Is( ConstitutiveLaw::COMPUTE_STRAIN_ENERGY ) )
        {
        	    if( Options.IsNot( ConstitutiveLaw::COMPUTE_STRESS ) )
                {
                    if(Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ))
                    {
                        Matrix ConstitutiveMatrix( StrainVector.size(), StrainVector.size());
                        noalias(ConstitutiveMatrix) = ZeroMatrix(StrainVector.size());

                        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
                        this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );
                    }
                    else
                    {
                        Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
                        this->CalculateStress( StrainVector, ConstitutiveMatrix, StressVector );
                    }
                }

            mStrainEnergy = 0.5 * inner_prod(StrainVector,StressVector); //Belytschko Nonlinear Finite Elements pag 226 (5.4.3) : w = 0.5*E:C:E
        }

    }

}


//***********************COMPUTE TOTAL STRESS PK2*************************************
//************************************************************************************


void LinearElastic3DLaw::CalculateStress( const Vector & rStrainVector,
					  const Matrix & rConstitutiveMatrix,
					  Vector& rStressVector )
{
    // 2nd Piola Kirchhoff StressVector or Cauchy StressVector
    if( rStressVector.size() != rStrainVector.size() )
        rStressVector.resize(rStrainVector.size(),false);

    noalias(rStressVector) = prod(rConstitutiveMatrix,rStrainVector);

}



//***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
//************************************************************************************


void LinearElastic3DLaw::CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,
        const double &rYoungModulus,
        const double &rPoissonCoefficient )
{
    rConstitutiveMatrix.clear();

    // 3D linear elastic constitutive matrix
    rConstitutiveMatrix ( 0 , 0 ) = (rYoungModulus*(1.0-rPoissonCoefficient)/((1.0+rPoissonCoefficient)*(1.0-2.0*rPoissonCoefficient)));
    rConstitutiveMatrix ( 1 , 1 ) = rConstitutiveMatrix ( 0 , 0 );
    rConstitutiveMatrix ( 2 , 2 ) = rConstitutiveMatrix ( 0 , 0 );

    rConstitutiveMatrix ( 3 , 3 ) = rConstitutiveMatrix ( 0 , 0 )*(1.0-2.0*rPoissonCoefficient)/(2.0*(1.0-rPoissonCoefficient));
    rConstitutiveMatrix ( 4 , 4 ) = rConstitutiveMatrix ( 3 , 3 );
    rConstitutiveMatrix ( 5 , 5 ) = rConstitutiveMatrix ( 3 , 3 );

    rConstitutiveMatrix ( 0 , 1 ) = rConstitutiveMatrix ( 0 , 0 )*rPoissonCoefficient/(1.0-rPoissonCoefficient);
    rConstitutiveMatrix ( 1 , 0 ) = rConstitutiveMatrix ( 0 , 1 );

    rConstitutiveMatrix ( 0 , 2 ) = rConstitutiveMatrix ( 0 , 1 );
    rConstitutiveMatrix ( 2 , 0 ) = rConstitutiveMatrix ( 0 , 1 );

    rConstitutiveMatrix ( 1 , 2 ) = rConstitutiveMatrix ( 0 , 1 );
    rConstitutiveMatrix ( 2 , 1 ) = rConstitutiveMatrix ( 0 , 1 );

}


//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void LinearElastic3DLaw::GetLawFeatures(Features& rFeatures)
{
    // Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	// Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	// Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	// Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

}

//******************CHECK CONSISTENCY IN THE CONSTITUTIVE LAW*************************
//************************************************************************************


int LinearElastic3DLaw::Check(const Properties& rMaterialProperties,
                              const GeometryType& rElementGeometry,
                              const ProcessInfo& rCurrentProcessInfo)
{
    HyperElastic3DLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    return 0;
}


} // Namespace Kratos

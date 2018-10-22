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
//  References:      This class is adapted from applications/SolidMechanicsApplication/custom_constitutive/hyperelastic_U_P_3D_law.cpp


// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_UP_3D_law.hpp"

#include "particle_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

HyperElasticUP3DLaw::HyperElasticUP3DLaw()
    : HyperElastic3DLaw()
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HyperElasticUP3DLaw::HyperElasticUP3DLaw(const HyperElasticUP3DLaw& rOther)
    : HyperElastic3DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HyperElasticUP3DLaw::Clone() const
{
    return Kratos::make_shared<HyperElasticUP3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HyperElasticUP3DLaw::~HyperElasticUP3DLaw()
{
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//*****************************MATERIAL RESPONSES*************************************
//************************************************************************************


void HyperElasticUP3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();
    const double&   DeterminantF          = rValues.GetDeterminantF();

    const GeometryType&  DomainGeometry   = rValues.GetElementGeometry ();
    const Vector&        ShapeFunctions   = rValues.GetShapeFunctionsValues ();

    Vector& StrainVector                  = rValues.GetStrainVector();
    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    ElasticVariables.Identity = identity_matrix<double> ( 3 );

    ElasticVariables.SetElementGeometry(DomainGeometry);
    ElasticVariables.SetShapeFunctionsValues(ShapeFunctions);

    // Initialize Splited Parts: Isochoric and Volumetric stresses and constitutive tensors
    double voigtsize = StressVector.size();
    VectorSplit SplitStressVector;
    MatrixSplit SplitConstitutiveMatrix;

    //1.- Lame constants
    const double& YoungModulus        = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient  = MaterialProperties[POISSON_RATIO];

    ElasticVariables.LameLambda      = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    ElasticVariables.LameMu          =  YoungModulus/(2*(1+PoissonCoefficient));

    //2.- Thermal constants
    if( MaterialProperties.Has(THERMAL_EXPANSION_COEFFICIENT) )
      ElasticVariables.ThermalExpansionCoefficient = MaterialProperties[THERMAL_EXPANSION_COEFFICIENT];
    else
      ElasticVariables.ThermalExpansionCoefficient = 0;

    if( MaterialProperties.Has(REFERENCE_TEMPERATURE) )
      ElasticVariables.ReferenceTemperature = MaterialProperties[REFERENCE_TEMPERATURE];
    else
      ElasticVariables.ReferenceTemperature = 0;

    //3.-Total DeformationGradientF Tensor 3D
    ElasticVariables.DeformationGradientF = DeformationGradientF;
    ElasticVariables.DeformationGradientF = Transform2DTo3D( ElasticVariables.DeformationGradientF );

    //4.-Determinant of the Total DeformationGradientF
    ElasticVariables.DeterminantF         = DeterminantF;

    //5.-Left Cauchy Green tensor b: (stored in the CauchyGreenMatrix)
    ElasticVariables.CauchyGreenMatrix.resize(3,3,false);
    noalias(ElasticVariables.CauchyGreenMatrix) = prod(ElasticVariables.DeformationGradientF,trans(ElasticVariables.DeformationGradientF));

    //6.-Calculate trace of Left Cauchy-Green tensor b
    ElasticVariables.traceCG = 0;
    for( unsigned int i=0; i<3; i++)
    {
        ElasticVariables.traceCG += ElasticVariables.CauchyGreenMatrix( i , i );
    }


    //8.-Almansi Strain:
    if(Options.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ))
    {
        // e= 0.5*(1-invbT*invb)
        this->CalculateAlmansiStrain(ElasticVariables.CauchyGreenMatrix,StrainVector);
    }

    //9.-Calculate Total kirchhoff stress
    SplitStressVector.Isochoric.resize(voigtsize,false);
    noalias(SplitStressVector.Isochoric) = ZeroVector(voigtsize);

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
      this->CalculateIsochoricStress( ElasticVariables, StressMeasure_Kirchhoff, SplitStressVector.Isochoric );

    Vector IsochoricStressVector = SplitStressVector.Isochoric;

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {

        SplitStressVector.Volumetric.resize(voigtsize,false);
        noalias(SplitStressVector.Volumetric) = ZeroVector(voigtsize);

        ElasticVariables.CauchyGreenMatrix = ElasticVariables.Identity;

        this->CalculateVolumetricStress ( ElasticVariables, SplitStressVector.Volumetric );

        //Kirchhoff Stress:
        StressVector = SplitStressVector.Isochoric + SplitStressVector.Volumetric;

        //std::cout<<" StressVector.Isochoric "<<SplitStressVector.Isochoric<<std::endl;
        //std::cout<<" StressVector.Volumetric "<<SplitStressVector.Volumetric<<std::endl;

        if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) )
        {
            StressVector = SplitStressVector.Isochoric;
        }
        else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) )
        {
            StressVector = SplitStressVector.Volumetric;
        }

    }

    //10.-Calculate Constitutive Matrix related to Total kirchhoff stress
    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {

        //initialize constitutive tensors
        ConstitutiveMatrix.clear();
        SplitConstitutiveMatrix.Isochoric  = ConstitutiveMatrix;
        SplitConstitutiveMatrix.Volumetric = ConstitutiveMatrix;

        ElasticVariables.CauchyGreenMatrix = ElasticVariables.Identity;

	Matrix IsoStressMatrix = MathUtils<double>::StressVectorToTensor( IsochoricStressVector );

        this->CalculateIsochoricConstitutiveMatrix ( ElasticVariables, IsoStressMatrix, SplitConstitutiveMatrix.Isochoric );

	this->CalculateVolumetricConstitutiveMatrix ( ElasticVariables, SplitConstitutiveMatrix.Volumetric );

        ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric + SplitConstitutiveMatrix.Volumetric;

        if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) )
        {
            ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric;
        }
        else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) )
        {
            ConstitutiveMatrix = SplitConstitutiveMatrix.Volumetric;
        }
    }

    // std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
    // std::cout<<" Stress "<<StressVector<<std::endl;

    //-----------------------------//


}


//******************************* COMPUTE DOMAIN PRESSURE  ***************************
//************************************************************************************


double &  HyperElasticUP3DLaw::CalculateVolumetricPressure (const MaterialResponseVariables & rElasticVariables,
							    double & rPressure)
{

    const GeometryType&  DomainGeometry =  rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues  =  rElasticVariables.GetShapeFunctionsValues();

    const unsigned int number_of_nodes  =  DomainGeometry.size();

    rPressure = 0;
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        rPressure += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(PRESSURE);
    }

    return rPressure;

}

//************************* COMPUTE DOMAIN PRESSURE FACTORS***************************
//************************************************************************************

Vector&  HyperElasticUP3DLaw::CalculateVolumetricPressureFactors (const MaterialResponseVariables & rElasticVariables,
							      Vector & rFactors)

{
    double Pressure = 0;
    Pressure = this->CalculateVolumetricPressure( rElasticVariables, Pressure );

    if(rFactors.size()!=3) rFactors.resize(3);

    rFactors[0] =  1.0;
    rFactors[1] =  2.0;
    rFactors[2] =  Pressure*rElasticVariables.DeterminantF;

    return rFactors;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void HyperElasticUP3DLaw::GetLawFeatures(Features& rFeatures)
{
    	//Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( FINITE_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );
	rFeatures.mOptions.Set( U_P_LAW );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

}




} // Namespace Kratos

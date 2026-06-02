//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Contri Alessandro
//
//  References:      This class is adapted from applications/ParticleMechanicsApplication/custom_constitutive/hyperelastic_U_P_3D_law.cpp


// System includes

// External includes

// Project includes
#include "custom_constitutive/displacement_newtonian_fluid_UP_3D_law.hpp"
#include "includes/cfd_variables.h"

#include "mpm_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

DispNewtonianFluidUP3DLaw::DispNewtonianFluidUP3DLaw()
    : DispNewtonianFluid3DLaw()
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

DispNewtonianFluidUP3DLaw::DispNewtonianFluidUP3DLaw(const DispNewtonianFluidUP3DLaw& rOther)
    : DispNewtonianFluid3DLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer DispNewtonianFluidUP3DLaw::Clone() const
{
    return Kratos::make_shared<DispNewtonianFluidUP3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

DispNewtonianFluidUP3DLaw::~DispNewtonianFluidUP3DLaw()
{
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//*****************************MATERIAL RESPONSES*************************************
//************************************************************************************


void DispNewtonianFluidUP3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{

    //-----------------------------//

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags &Options=rValues.GetOptions();

    const Properties& MaterialProperties=rValues.GetMaterialProperties();
    const Matrix& DeformationGradientF=rValues.GetDeformationGradientF();
    const double& DeterminantF=rValues.GetDeterminantF();

    Vector& StrainVector=rValues.GetStrainVector();
    Vector& StressVector=rValues.GetStressVector();
    Matrix& ConstitutiveMatrix=rValues.GetConstitutiveMatrix();

    const GeometryType& domain_geometry = rValues.GetElementGeometry();
    const Vector& shape_functions = rValues.GetShapeFunctionsValues();
    const ProcessInfo& current_process_info = rValues.GetProcessInfo();

    // Matrix& VolumetricConstitutiveMatrix  = rValues.GetVolumetricConstitutiveMatrix();


    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ViscousVariables;
    ViscousVariables.Identity=IdentityMatrix(3);

    ViscousVariables.SetElementGeometry(domain_geometry);
    ViscousVariables.SetShapeFunctionsValues(shape_functions);

    // Initialize Splited Parts: Isochoric and Volumetric stresses and constitutive tensors
    double voigtsize = StressVector.size();
    VectorSplit SplitStressVector;
    MatrixSplit SplitConstitutiveMatrix;

    //1.- Material constants
    ViscousVariables.Mu = MaterialProperties[DYNAMIC_VISCOSITY];
    ViscousVariables.BulkModulus = MaterialProperties[BULK_MODULUS];

    ViscousVariables.DeltaTime = current_process_info[DELTA_TIME];

    //3.-Total DeformationGradientF Tensor 3D
    ViscousVariables.DeformationGradientF = DeformationGradientF;
    ViscousVariables.DeformationGradientF = Transform2DTo3D(ViscousVariables.DeformationGradientF);

    //4.-Determinant of the Total DeformationGradientF
    ViscousVariables.DeterminantF = DeterminantF;

    //5.-Left Cauchy Green tensor b: (stored in the CauchyGreenMatrix)
    ViscousVariables.CauchyGreenMatrix.resize(3, 3, false);
    noalias(ViscousVariables.CauchyGreenMatrix) = prod(ViscousVariables.DeformationGradientF, trans(ViscousVariables.DeformationGradientF));

    //6.-Calculate deformation rate
    this->CalculateDeformationRate(ViscousVariables);

    //7.-Almansi Strain:
    if (Options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
    {
        // e= 0.5*(1-invbT*invb)
        this->CalculateAlmansiStrain(ViscousVariables.CauchyGreenMatrix,StrainVector);
    }

    //8.-Calculate Total Cauchy stress
    SplitStressVector.Isochoric.resize(voigtsize,false);
    noalias(SplitStressVector.Isochoric) = ZeroVector(voigtsize);

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
      this->CalculateIsochoricStress( ViscousVariables, StressMeasure_Cauchy, SplitStressVector.Isochoric );

    Vector IsochoricStressVector = SplitStressVector.Isochoric;

    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {

        SplitStressVector.Volumetric.resize(voigtsize,false);
        noalias(SplitStressVector.Volumetric) = ZeroVector(voigtsize);

        this->CalculateVolumetricStress ( ViscousVariables, SplitStressVector.Volumetric );

        //Cauchy Stress:
        StressVector = SplitStressVector.Isochoric + SplitStressVector.Volumetric;

        if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) )
        {
            StressVector = SplitStressVector.Isochoric;
        }
        else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) )
        {
            StressVector = SplitStressVector.Volumetric;
        }

    }

    //9.-Calculate Constitutive Matrix related to Total Cauchy stress
    if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
    {

        //initialize constitutive tensors
        ConstitutiveMatrix.clear();
        SplitConstitutiveMatrix.Isochoric  = ConstitutiveMatrix;
        SplitConstitutiveMatrix.Volumetric = ConstitutiveMatrix;

	    Matrix IsoStressMatrix = MathUtils<double>::StressVectorToTensor( IsochoricStressVector );

        this->CalculateIsochoricConstitutiveMatrix ( ViscousVariables, IsoStressMatrix, SplitConstitutiveMatrix.Isochoric );

	    this->CalculateVolumetricConstitutiveMatrix ( ViscousVariables, SplitConstitutiveMatrix.Volumetric );

        ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric + SplitConstitutiveMatrix.Volumetric;

        // VolumetricConstitutiveMatrix = SplitConstitutiveMatrix.Volumetric

        if( Options.Is(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) )
        {
            ConstitutiveMatrix = SplitConstitutiveMatrix.Isochoric;
        }
        else if( Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) )
        {
            ConstitutiveMatrix = SplitConstitutiveMatrix.Volumetric;
        }
    }
}


//******************************* COMPUTE DOMAIN PRESSURE  ***************************
//************************************************************************************


double &  DispNewtonianFluidUP3DLaw::CalculateVolumetricPressure (const MaterialResponseVariables & rViscousVariables,
							    double & rPressure)
{

    const GeometryType&  DomainGeometry =  rViscousVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues  =  rViscousVariables.GetShapeFunctionsValues();

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

Vector&  DispNewtonianFluidUP3DLaw::CalculateVolumetricPressureFactors (const MaterialResponseVariables & rViscousVariables,
							      Vector & rFactors)

{
    double Pressure = 0;
    Pressure = this->CalculateVolumetricPressure( rViscousVariables, Pressure );

    if(rFactors.size()!=3) rFactors.resize(3,false);

    rFactors[0] =  1.0;
    rFactors[1] =  2.0;
    rFactors[2] =  Pressure;

    return rFactors;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void DispNewtonianFluidUP3DLaw::GetLawFeatures(Features& rFeatures)
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

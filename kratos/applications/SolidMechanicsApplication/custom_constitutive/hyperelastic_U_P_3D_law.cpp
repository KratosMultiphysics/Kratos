//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hyperelastic_U_P_3D_law.hpp"

#include "solid_mechanics_application.h"

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
    HyperElasticUP3DLaw::Pointer p_clone(new HyperElasticUP3DLaw(*this));
    return p_clone;
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


void  HyperElasticUP3DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
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
    Matrix& DeformationGradientF0         = rValues.GetDeformationGradientF0();
    double& DeterminantF0                 = rValues.GetDeterminantF0();

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    ElasticVariables.IdentityMatrix = identity_matrix<double> ( 3 );

    // Initialize Splited Parts: Isochoric and Volumetric stresses and constitutive tensors
    double voigtsize = StressVector.size();
    VectorSplit SplitStressVector;
    MatrixSplit SplitConstitutiveMatrix;

    //1.- Lame constants
    const double& YoungModulus          = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient    = MaterialProperties[POISSON_RATIO];

    ElasticVariables.LameLambda      = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    ElasticVariables.LameMu          =  YoungModulus/(2*(1+PoissonCoefficient));


    //-----------------------------//
    //OPTION 1: ( initial configuration )
    if( Options.Is( ConstitutiveLaw::INITIAL_CONFIGURATION ) )
    {

        //2.-Total Deformation Gradient
        Matrix TotalDeformationGradientF0  = DeformationGradientF;
        TotalDeformationGradientF0         = DeformationGradient3D( TotalDeformationGradientF0 );

        //3.-Determinant of the Total Deformation Gradient
        ElasticVariables.DeterminantF0 = DeterminantF0 * DeterminantF;

        //4.-Right Cauchy Green
        Matrix RightCauchyGreen = prod(trans(TotalDeformationGradientF0),TotalDeformationGradientF0);

        //5.-Inverse of the Right Cauchy-Green tensor C: (stored in the CauchyGreenMatrix)
        ElasticVariables.traceCG = 0;
        ElasticVariables.CauchyGreenMatrix( 3, 3 );
        MathUtils<double>::InvertMatrix( RightCauchyGreen, ElasticVariables.CauchyGreenMatrix, ElasticVariables.traceCG);

        //6.-Green-Lagrange Strain:
        if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
        {
            this->CalculateGreenLagrangeStrain(RightCauchyGreen, StrainVector);
        }

        //7.-Calculate Total PK2 stress
        SplitStressVector.Isochoric = ZeroVector(voigtsize);

        if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
            this->CalculateIsochoricStress( ElasticVariables, StressMeasure_PK2, SplitStressVector.Isochoric );

        Vector IsochoricStressVector = SplitStressVector.Isochoric;

        if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
        {

            SplitStressVector.Volumetric = ZeroVector(voigtsize);

            this->CalculateVolumetricStress ( ElasticVariables, DomainGeometry, ShapeFunctions, SplitStressVector.Volumetric );

            //PK2 Stress:
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

        if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
        {

            //initialize constitutive tensors
            ConstitutiveMatrix.clear();
            SplitConstitutiveMatrix.Isochoric  = ConstitutiveMatrix;
            SplitConstitutiveMatrix.Volumetric = ConstitutiveMatrix;

            this->CalculateIsochoricConstitutiveMatrix ( ElasticVariables, IsochoricStressVector, SplitConstitutiveMatrix.Isochoric );

            this->CalculateVolumetricConstitutiveMatrix ( ElasticVariables, DomainGeometry, ShapeFunctions, SplitConstitutiveMatrix.Volumetric );

            //if( Options.Is(ConstitutiveLaw::TOTAL_TENSOR ) )
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
    }


    //-----------------------------//
    //OPTION 2: ( last known configuration : updated lagrangian approach only )
    if( Options.Is( ConstitutiveLaw::LAST_KNOWN_CONFIGURATION ) || Options.Is( ConstitutiveLaw::FINAL_CONFIGURATION ) )
    {

        //Determinant of the Total Deformation Gradient
        ElasticVariables.DeterminantF0 = DeterminantF0 * DeterminantF;

        //Left Cauchy-Green tensor b
        Matrix TotalDeformationGradientF0  = prod(DeformationGradientF, DeformationGradientF0);
        TotalDeformationGradientF0         = DeformationGradient3D( TotalDeformationGradientF0 );
        ElasticVariables.CauchyGreenMatrix = prod(TotalDeformationGradientF0,trans(TotalDeformationGradientF0));

        //Calculate trace of Left Cauchy-Green tensor b
        ElasticVariables.traceCG = 0;
        for( unsigned int i=0; i<3; i++)
        {
            ElasticVariables.traceCG += ElasticVariables.CauchyGreenMatrix( i , i );
        }

        SplitStressVector.Isochoric = ZeroVector(voigtsize);

        if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
            this->CalculateIsochoricStress( ElasticVariables, StressMeasure_Kirchhoff, SplitStressVector.Isochoric );

        Vector IsochoricStressVector = SplitStressVector.Isochoric;

        TransformStresses(SplitStressVector.Isochoric, DeformationGradientF, DeterminantF, StressMeasure_Kirchhoff, StressMeasure_PK2); //2nd PK Stress in the last known configuration


        if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
        {

            SplitStressVector.Volumetric = ZeroVector(voigtsize);

            ElasticVariables.CauchyGreenMatrix = ElasticVariables.IdentityMatrix;

            this->CalculateVolumetricStress ( ElasticVariables, DomainGeometry, ShapeFunctions, SplitStressVector.Volumetric );

            TransformStresses(SplitStressVector.Volumetric, DeformationGradientF, DeterminantF, StressMeasure_Kirchhoff, StressMeasure_PK2); //2nd PK Stress in the last known configuration

            //PK2 Stress in the last known configuration:
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

        if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
        {


            Matrix InverseDeformationGradientF ( 3, 3 );
            double DetInvF=0;
            MathUtils<double>::InvertMatrix( DeformationGradientF, InverseDeformationGradientF, DetInvF);

            ElasticVariables.CauchyGreenMatrix = ElasticVariables.IdentityMatrix;


            //initialize constitutive tensors
            ConstitutiveMatrix.clear();
            SplitConstitutiveMatrix.Isochoric  = ConstitutiveMatrix;
            SplitConstitutiveMatrix.Volumetric = ConstitutiveMatrix;


            this->CalculateIsochoricConstitutiveMatrix ( ElasticVariables, IsochoricStressVector, InverseDeformationGradientF, SplitConstitutiveMatrix.Isochoric );

            this->CalculateVolumetricConstitutiveMatrix ( ElasticVariables, InverseDeformationGradientF, DomainGeometry, ShapeFunctions, SplitConstitutiveMatrix.Volumetric );

            SplitConstitutiveMatrix.Isochoric  *= DeterminantF;
            SplitConstitutiveMatrix.Volumetric *= DeterminantF;

            //if( Options.Is(ConstitutiveLaw::TOTAL_TENSOR ) )
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

    }

    // std::cout<<" Constitutive "<<ConstitutiveMatrix<<std::endl;
    // std::cout<<" Stress "<<StressVector<<std::endl;

}

//************************************************************************************
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
    Matrix& DeformationGradientF0         = rValues.GetDeformationGradientF0();
    double& DeterminantF0                 = rValues.GetDeterminantF0();

    Vector& StressVector                  = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix            = rValues.GetConstitutiveMatrix();

    //-----------------------------//

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;
    ElasticVariables.IdentityMatrix = identity_matrix<double> ( 3 );


    // Initialize Splited Parts: Isochoric and Volumetric stresses and constitutive tensors
    double voigtsize = StressVector.size();
    VectorSplit SplitStressVector;
    MatrixSplit SplitConstitutiveMatrix;

    //1.- Lame constants
    const double& YoungModulus        = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient  = MaterialProperties[POISSON_RATIO];

    ElasticVariables.LameLambda      = (YoungModulus*PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    ElasticVariables.LameMu          =  YoungModulus/(2*(1+PoissonCoefficient));


    //2.-Determinant of the Total Deformation Gradient
    ElasticVariables.DeterminantF0 = DeterminantF0 * DeterminantF;

    //3.-Push-Forward Left Cauchy-Green tensor b to the new configuration
    Matrix TotalDeformationGradientF0  = prod(DeformationGradientF, DeformationGradientF0);
    TotalDeformationGradientF0         = DeformationGradient3D( TotalDeformationGradientF0 );
    ElasticVariables.CauchyGreenMatrix = prod(TotalDeformationGradientF0,trans(TotalDeformationGradientF0));

    //std::cout<<" TotalDeformationGradientF0 "<<TotalDeformationGradientF0<<std::endl;

    //4.-Calculate trace of Left Cauchy-Green tensor b
    ElasticVariables.traceCG = 0;
    for( unsigned int i=0; i<3; i++)
    {
        ElasticVariables.traceCG += ElasticVariables.CauchyGreenMatrix( i , i );
    }

    //4.-Almansi Strain:
    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
    {
        // e= 0.5*(1-invbT*invb)
        this->CalculateAlmansiStrain(ElasticVariables.CauchyGreenMatrix,StrainVector);
    }

    //5.-Calculate Total PK2 stress
    SplitStressVector.Isochoric = ZeroVector(voigtsize);

    if( Options.Is(ConstitutiveLaw::COMPUTE_STRESS ) || Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
        this->CalculateIsochoricStress( ElasticVariables, StressMeasure_Kirchhoff, SplitStressVector.Isochoric );

    Vector IsochoricStressVector = SplitStressVector.Isochoric;

    //OPTION 1:
    if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) )
    {

        SplitStressVector.Volumetric = ZeroVector(voigtsize);

        ElasticVariables.CauchyGreenMatrix = ElasticVariables.IdentityMatrix;

        this->CalculateVolumetricStress ( ElasticVariables, DomainGeometry, ShapeFunctions, SplitStressVector.Volumetric );

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


    if( Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {

        //initialize constitutive tensors
        ConstitutiveMatrix.clear();
        SplitConstitutiveMatrix.Isochoric  = ConstitutiveMatrix;
        SplitConstitutiveMatrix.Volumetric = ConstitutiveMatrix;

        ElasticVariables.CauchyGreenMatrix = ElasticVariables.IdentityMatrix;

        this->CalculateIsochoricConstitutiveMatrix ( ElasticVariables, IsochoricStressVector, SplitConstitutiveMatrix.Isochoric );

        this->CalculateVolumetricConstitutiveMatrix ( ElasticVariables, DomainGeometry, ShapeFunctions, SplitConstitutiveMatrix.Volumetric );

        //if( Options.Is(ConstitutiveLaw::TOTAL_TENSOR ) )
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



}


//******************************* COMPUTE ISOCHORIC STRESS  **************************
//************************************************************************************



void HyperElasticUP3DLaw::CalculateIsochoricStress( const MaterialResponseVariables & rElasticVariables,
        StressMeasure rStressMeasure,
        Vector& rIsoStressVector )
{

    //1.-Identity build
    Matrix IsoStressMatrix ( 3, 3 );


    if(rStressMeasure == StressMeasure_PK2)
    {

        //rElasticVariables.CauchyGreenMatrix is InverseRightCauchyGreen

        //2.-Incompressible part of the 2nd Piola Kirchhoff Stress Matrix
        IsoStressMatrix  = (rElasticVariables.IdentityMatrix - (rElasticVariables.traceCG/3.0)*rElasticVariables.CauchyGreenMatrix );
        IsoStressMatrix *= rElasticVariables.LameMu*pow(rElasticVariables.DeterminantF0,(-2.0/3.0));

        //std::cout<<" PK2 "<<std::endl;
    }

    if(rStressMeasure == StressMeasure_Kirchhoff)
    {

        //rElasticVariables.CauchyGreenMatrix is LeftCauchyGreen

        //2.-Incompressible part of the Kirchhoff Stress Matrix
        IsoStressMatrix  = (rElasticVariables.CauchyGreenMatrix - (rElasticVariables.traceCG/3.0)*rElasticVariables.IdentityMatrix );
        IsoStressMatrix *= rElasticVariables.LameMu*pow(rElasticVariables.DeterminantF0,(-2.0/3.0));

        //std::cout<<" Kirchooff "<<std::endl;

    }


    rIsoStressVector = MathUtils<double>::StressTensorToVector(IsoStressMatrix,rIsoStressVector.size());

}

//******************************* COMPUTE DOMAIN PRESSURE  ***************************
//************************************************************************************


double &  HyperElasticUP3DLaw::CalculateDomainPressure (const GeometryType& rElementGeometry,
        const Vector & rShapeFunctions,
        double & rPressure)
{

    const unsigned int number_of_nodes = rElementGeometry.size();

    rPressure = 0;
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        rPressure += rShapeFunctions[j] * rElementGeometry[j].GetSolutionStepValue(PRESSURE);
    }

    return rPressure;
}


//******************************* COMPUTE VOLUMETRIC STRESS  *************************
//************************************************************************************

void HyperElasticUP3DLaw::CalculateVolumetricStress(const MaterialResponseVariables & rElasticVariables,
        const GeometryType& rElementGeometry,
        const Vector & rShapeFunctions,
        Vector& rVolStressVector )
{

    //1.- Declaration
    Matrix VolStressMatrix ( 3 , 3 );

    double Pressure = 0;

    Pressure = CalculateDomainPressure (rElementGeometry, rShapeFunctions, Pressure);

    //2.- Volumetric part of the Kirchhoff StressMatrix from nodal pressures
    VolStressMatrix = rElasticVariables.DeterminantF0 * Pressure * rElasticVariables.CauchyGreenMatrix;


    //3.- Volumetric part of the  Kirchhoff StressMatrix

    /*
    //(J²-1)/2
    //double auxiliar1 =  0.5*(rdetF0*rdetF0-1);

    //(ln(J))
    double auxiliar1 =  std::log(rdetF0);

    double BulkModulus= rLameLambda + (2.0/3.0) * rLameMu;

    //2.-Volumetric part of the Kirchhoff Stress Matrix

     VolStressMatrix  = BulkModulus * auxiliar1 * rMatrixIC;
    */

    rVolStressVector = MathUtils<double>::StressTensorToVector(VolStressMatrix,rVolStressVector.size());

}

//***********************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX************************
//************************************************************************************

void HyperElasticUP3DLaw::CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
        const Vector & rIsoStressVector,
        Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    Matrix IsoStressMatrix = MathUtils<double>::StressVectorToTensor( rIsoStressVector );

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = IsochoricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, IsoStressMatrix,
                                          this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}

//***********************COMPUTE VOLUMETRIC CONSTITUTIVE MATRIX***********************
//************************************************************************************


void HyperElasticUP3DLaw::CalculateVolumetricConstitutiveMatrix ( const MaterialResponseVariables & rElasticVariables,
        const GeometryType& rElementGeometry,
        const Vector & rShapeFunctions,
        Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    double Pressure = 0;
    Pressure = CalculateDomainPressure ( rElementGeometry, rShapeFunctions, Pressure);


    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, Pressure,
                                          this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}


//******************COMPUTE ISOCHORIC CONSTITUTIVE MATRIX PULL-BACK*******************
//************************************************************************************

void HyperElasticUP3DLaw::CalculateIsochoricConstitutiveMatrix ( const MaterialResponseVariables & rElasticVariables,
        const Vector & rIsoStressVector,
        const Matrix & rInverseDeformationGradientF,
        Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    Matrix IsoStressMatrix = MathUtils<double>::StressVectorToTensor( rIsoStressVector );

    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = IsochoricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, IsoStressMatrix, rInverseDeformationGradientF,
                                          this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}

//***************COMPUTE VOLUMETRIC CONSTITUTIVE MATRIX PUSH-FORWARD******************
//************************************************************************************

void HyperElasticUP3DLaw::CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables & rElasticVariables,
        const Matrix & rInverseDeformationGradientF,
        const GeometryType& rElementGeometry,
        const Vector & rShapeFunctions,
        Matrix& rConstitutiveMatrix)
{

    rConstitutiveMatrix.clear();

    double Pressure = 0;
    Pressure = CalculateDomainPressure ( rElementGeometry, rShapeFunctions, Pressure);


    for(unsigned int i=0; i<6; i++)
    {
        for(unsigned int j=0; j<6; j++)
        {
            rConstitutiveMatrix( i, j ) = VolumetricConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, rInverseDeformationGradientF, Pressure,
                                          this->msIndexVoigt3D6C[i][0], this->msIndexVoigt3D6C[i][1], this->msIndexVoigt3D6C[j][0], this->msIndexVoigt3D6C[j][1]);
        }

    }


}


//********************CONSTITUTIVE MATRIX ISOCHORIC COMPONENT*************************
//************************************************************************************


double& HyperElasticUP3DLaw::IsochoricConstitutiveComponent(double & rCabcd,
        const MaterialResponseVariables & rElasticVariables,
        const Matrix & rIsoStressMatrix,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{

    //Isochoric part of the hyperelastic constitutive tensor component: (J²-1)/2  -  (ln(J)/J)

    rCabcd  = (1.0/3.0)*(rElasticVariables.CauchyGreenMatrix(a,b)*rElasticVariables.CauchyGreenMatrix(c,d));
    rCabcd -= (0.5*(rElasticVariables.CauchyGreenMatrix(a,c)*rElasticVariables.CauchyGreenMatrix(b,d)+rElasticVariables.CauchyGreenMatrix(a,d)*rElasticVariables.CauchyGreenMatrix(b,c)));
    rCabcd *= pow(rElasticVariables.DeterminantF0,(-2.0/3.0))*rElasticVariables.traceCG*rElasticVariables.LameMu;

    rCabcd += (rElasticVariables.CauchyGreenMatrix(c,d)*rIsoStressMatrix(a,b) + rIsoStressMatrix(c,d)*rElasticVariables.CauchyGreenMatrix(a,b));
    rCabcd *= (-2.0/3.0);

    return rCabcd;
}


//********************CONSTITUTIVE MATRIX VOLUMETRIC COMPONENT************************
//************************************************************************************


double& HyperElasticUP3DLaw::VolumetricConstitutiveComponent(double & rCabcd,
        const MaterialResponseVariables & rElasticVariables,
        const double & rPressure,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{

    //Volumetric part of the hyperelastic constitutive tensor component: (J²-1)/2  -  (ln(J)/J)

    rCabcd = (rElasticVariables.CauchyGreenMatrix(a,b)*rElasticVariables.CauchyGreenMatrix(c,d));
    rCabcd -= 2*(0.5*(rElasticVariables.CauchyGreenMatrix(a,c)*rElasticVariables.CauchyGreenMatrix(b,d)+rElasticVariables.CauchyGreenMatrix(a,d)*rElasticVariables.CauchyGreenMatrix(b,c)));

    rCabcd *= rPressure*rElasticVariables.DeterminantF0;


    /*
    double BulkModulus= rLameLambda + (2.0/3.0) * rLameMu;

    //(J²-1)/2
    //double auxiliar1 =  rElasticVariables.DeterminantF0*rElasticVariables.DeterminantF0;
    //double auxiliar2 =  (rElasticVariables.DeterminantF0*rElasticVariables.DeterminantF0-1);
    //double auxiliar3 =  BulkModulus;

    //(ln(J))
    double auxiliar1 =  1.0;
    double auxiliar2 =  (2.0*std::log(rVariables.DeterminantF0));
    double auxiliar3 =  BulkModulus/rElasticVariables.DeterminantF0;

    //1.Volumetric Elastic constitutive tensor component:
    rCabcd = auxiliar1*(rElasticVariables.CauchyGreenMatrix(a,b)*rElasticVariables.CauchyGreenMatrix(c,d));
    rCabcd -= auxiliar2*(0.5*(rElasticVariables.CauchyGreenMatrix(a,c)*rElasticVariables.CauchyGreenMatrix(b,d)+rElasticVariables.CauchyGreenMatrix(a,d)*rElasticVariables.CauchyGreenMatrix(b,c)));
    rCabcd *= auxiliar3;
    */



    return rCabcd;
}


//**************CONSTITUTIVE MATRIX ISOCHORIC COMPONENT PUSH-FORWARD******************
//************************************************************************************


double& HyperElasticUP3DLaw::IsochoricConstitutiveComponent(double & rCabcd,
        const MaterialResponseVariables & rElasticVariables,
        const Matrix & rIsoStressMatrix,
        const Matrix & rInverseDeformationGradientF,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{

    rCabcd = 0;
    double Cijkl=0;

    unsigned int dimension =  rInverseDeformationGradientF.size1();

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
                    rCabcd +=rInverseDeformationGradientF(a,i)*rInverseDeformationGradientF(b,j)*rInverseDeformationGradientF(c,k)*rInverseDeformationGradientF(d,l)*IsochoricConstitutiveComponent(Cijkl, rElasticVariables, rIsoStressMatrix, i, j, k, l);
                }
            }
        }
    }

    return rCabcd;


}


//**************CONSTITUTIVE MATRIX VOLUMETRIC COMPONENT PUSH-FORWARD*****************
//************************************************************************************

double& HyperElasticUP3DLaw::VolumetricConstitutiveComponent(double & rCabcd,
        const MaterialResponseVariables & rElasticVariables,
        const Matrix & rInverseDeformationGradientF,
        const double & rPressure,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d)
{



    rCabcd = 0;
    double Cijkl = 0;
    unsigned int dimension =  rInverseDeformationGradientF.size1();

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
                    rCabcd +=rInverseDeformationGradientF(a,i)*rInverseDeformationGradientF(b,j)*rInverseDeformationGradientF(c,k)*rInverseDeformationGradientF(d,l)*VolumetricConstitutiveComponent(Cijkl, rElasticVariables, rPressure, i, j, k, l);
                }
            }
        }
    }

    return rCabcd;

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

	//Set strain measure requires by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);
	
	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

}




} // Namespace Kratos

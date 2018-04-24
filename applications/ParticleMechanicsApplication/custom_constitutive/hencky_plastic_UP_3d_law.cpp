//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hencky_plastic_UP_3d_law.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"
#include "particle_mechanics_application.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************
HenckyElasticPlasticUP3DLaw::HenckyElasticPlasticUP3DLaw()
    : HenckyElasticPlastic3DLaw()
{

}



HenckyElasticPlasticUP3DLaw::HenckyElasticPlasticUP3DLaw(MPMFlowRulePointer pMPMFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : HenckyElasticPlastic3DLaw( pMPMFlowRule, pYieldCriterion, pHardeningLaw)
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

HenckyElasticPlasticUP3DLaw::HenckyElasticPlasticUP3DLaw(const HenckyElasticPlasticUP3DLaw&  rOther)
    : HenckyElasticPlastic3DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer HenckyElasticPlasticUP3DLaw::Clone() const
{
    HenckyElasticPlasticUP3DLaw::Pointer p_clone(new HenckyElasticPlasticUP3DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

HenckyElasticPlasticUP3DLaw::~HenckyElasticPlasticUP3DLaw()
{
}


//************************************************************************************
//************************************************************************************




//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//*****************************MATERIAL RESPONSES*************************************
//************************************************************************************
void HenckyElasticPlasticUP3DLaw::CalculatePrincipalStressTrial(const MaterialResponseVariables & rElasticVariables, Parameters& rValues, const MPMFlowRule::RadialReturnVariables & rReturnMappingVariables, Matrix& rNewElasticLeftCauchyGreen, Matrix& rStressMatrix)
{

    const Properties& MaterialProperties   = rValues.GetMaterialProperties();
    const double& Young       = MaterialProperties[YOUNG_MODULUS];
    const double& Nu = MaterialProperties[POISSON_RATIO];

    double ShearModulus = Young/(2*(1 + Nu));

    //1- calculate the deviatoric elastic streches eigenvalues
    Vector MainStrain      = ZeroVector(3);
    //const GeometryType&  DomainGeometry =  rElasticVariables.GetElementGeometry();
    for (unsigned int i = 0; i<3; ++i)
    {
        MainStrain[i] = rNewElasticLeftCauchyGreen(i,i);
    }
    //std::cout<<" rMainStrain in up cl "<<MainStrain<<std::endl;
    Vector DeviatoricMainStrain = ZeroVector(3);

    double TracePrincipalStrain = MainStrain[0] + MainStrain[1] + MainStrain[2];
    //std::cout<<" TracePrincipalStrain in up cl "<<TracePrincipalStrain<<std::endl;
    Vector DeviatoricPrincipalStress = ZeroVector(3);


    for (unsigned int i=0; i<3 ; i++)
    {
        DeviatoricMainStrain[i] = MainStrain[i] -TracePrincipalStrain/3;
        DeviatoricPrincipalStress[i] = 2 * ShearModulus * DeviatoricMainStrain[i];
    }
    //std::cout<<" DeviatoricMainStrain in up cl "<<DeviatoricMainStrain<<std::endl;
    //std::cout<<" DeviatoricPrincipalStress in up cl "<<DeviatoricPrincipalStress<<std::endl;
    //now I have to transform the principal deviaric stress in cartesian stress
    //such that I can add them to the volumetric stress


    Vector auxN = ZeroVector(3);
    Matrix auxM = ZeroMatrix(3,3);
    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            auxN(j) = rReturnMappingVariables.MainDirections(i,j);
        }
        //std::cout<<" auxN in up cl "<<auxN<<std::endl;
        auxM = MathUtils<double>::TensorProduct3(auxN, auxN);
        rStressMatrix += DeviatoricPrincipalStress(i)*auxM;
        //std::cout<<" rStressMatrix in up cl "<<rStressMatrix<<std::endl;
    }

    //if(DomainGeometry[0].Id() == 2295 && DomainGeometry[1].Id() == 2315 && DomainGeometry[2].Id() == 2313)
    //{

    ////std::cout<<" MainStrain "<<MainStrain<<std::endl;
    //std::cout<<" rStressMatrix dev"<<rStressMatrix<<std::endl;

    //}
    //std::cout<<" rStressMatrix dev"<<rStressMatrix<<std::endl;
    double Pressure = 0;
    GetDomainPressure( Pressure, rElasticVariables);
    //std::cout<<" Pressure in up cl "<<Pressure<<std::endl;

    for (unsigned int i = 0; i < 3; ++i)
        rStressMatrix(i,i) += Pressure * rElasticVariables.DeterminantF;
    //std::cout<<" total rStressMatrix in up cl "<<rStressMatrix<<std::endl;
    //if(DomainGeometry[0].Id() == 2295 && DomainGeometry[1].Id() == 2315 && DomainGeometry[2].Id() == 2313)
    //{
    //std::cout<<" Pressure "<<Pressure<<std::endl;
    //std::cout<<" rStressMatrix tot"<<rStressMatrix<<std::endl;

    //}

    //Now I have to apply the spectral theorem
    Matrix EigenVectors  = ZeroMatrix(3,3);
    Vector EigenValues   = ZeroVector(3);

    double tol = 1e-9;
    int iter = 100;
    //std::cout<<" rCauchyGreenMatrix "<<rCauchyGreenMatrix<<std::endl;
    SolidMechanicsMathUtilities<double>::EigenVectors(rStressMatrix, EigenVectors, EigenValues, tol, iter);
    //std::cout<<" EigenValues in up cl "<<EigenValues<<std::endl;


    rStressMatrix.clear();
    for(unsigned int i=0; i<3; i++)
    {
        rStressMatrix(i,i) = EigenValues(i);
    }
    //if(DomainGeometry[0].Id() == 2304 && DomainGeometry[1].Id() == 2326 && DomainGeometry[2].Id() == 2311)
    //{

    //std::cout<<" MainStrain "<<MainStrain<<std::endl;
    //std::cout<<" rStressMatrix"<<rStressMatrix<<std::endl;

    //}
    //std::cout<<" principal rStressMatrix in up cl "<<rStressMatrix<<std::endl;
}

void HenckyElasticPlasticUP3DLaw::CorrectDomainPressure( Matrix& rStressMatrix, const MaterialResponseVariables & rElasticVariables)
{
    //const GeometryType&  DomainGeometry =  rElasticVariables.GetElementGeometry();
    double MeanPressure = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
        MeanPressure += rStressMatrix(i,i);

    MeanPressure /=3.0;
    //if ( fabs(MeanPressure) > 1.0E-4)
    //   std::cout << " UNCORRECTED PRESSURE " << MeanPressure << std::endl;

    for (unsigned int i = 0; i < 3; ++i)
        rStressMatrix(i,i) -= MeanPressure;


    double Pressure = 0;
    GetDomainPressure( Pressure, rElasticVariables);

    for (unsigned int i = 0; i < 3; ++i)
        rStressMatrix(i,i) += Pressure * rElasticVariables.DeterminantF;



    //std::cout << " THIS DET " << rElasticVariables.DeterminantF << std::endl;

}

void HenckyElasticPlasticUP3DLaw::GetDomainPressure( double& rPressure, const MaterialResponseVariables& rElasticVariables)
{

    rPressure = 0.0;
    const GeometryType&  DomainGeometry =  rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues  =  rElasticVariables.GetShapeFunctionsValues();

    const unsigned int number_of_nodes  =  DomainGeometry.size();

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        rPressure += ShapeFunctionsValues[j] * DomainGeometry[j].FastGetSolutionStepValue(PRESSURE); //NOOOOO
    }

}

void HenckyElasticPlasticUP3DLaw::CalculateElastoPlasticTangentMatrix( const MPMFlowRule::RadialReturnVariables & rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& rAlpha, Matrix& rElastoPlasticTangentMatrix, const MaterialResponseVariables& rElasticVariables )
{

    mpMPMFlowRule->ComputeElastoPlasticTangentMatrix( rReturnMappingVariables,  rNewElasticLeftCauchyGreen, rAlpha, rElastoPlasticTangentMatrix);
    // ADDING THE K TERMS
    double Pressure;
    //GetDomainPressure( Pressure, rElasticVariables);
    GetDomainPressure( Pressure, rElasticVariables);

    Pressure *= rElasticVariables.DeterminantF;

    double Young = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    double Nu = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];

    double K = Young / (3.0 * (1.0 - 2.0*Nu));

    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = 0; j < 3 ; ++j)
        {
            rElastoPlasticTangentMatrix(i,j)  -= K;
        }
    }

    Matrix FourthOrderIdentity = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i)
        FourthOrderIdentity(i,i) = 1.0;

    for (unsigned int i = 3; i<6; ++i)
        FourthOrderIdentity(i,i) = 0.50;
    // VOIGT NOTATION AND NOT KELVIN

    Matrix IdentityCross = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            IdentityCross(i,j) = 1.0;
        }
    }

    rElastoPlasticTangentMatrix += Pressure* ( IdentityCross - 2.0 * FourthOrderIdentity);

}

void HenckyElasticPlasticUP3DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
        Vector& rStrainVector )
{

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = 0.5 * ( rRightCauchyGreen( 2, 2 ) - 1.00 );
    rStrainVector[3] = rRightCauchyGreen( 0, 1 ); // xy
    rStrainVector[4] = rRightCauchyGreen( 1, 2 ); // yz
    rStrainVector[5] = rRightCauchyGreen( 0, 2 ); // xz

}



//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************
void HenckyElasticPlasticUP3DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
        Vector& rStrainVector )
{

    // e = 0.5*(1-invFT*invF) or e = 0.5*(1-inv(b))

    //Calculating the inverse of the jacobian
    Matrix InverseLeftCauchyGreen ( 3, 3 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 2, 2 ) );
    rStrainVector[3] = - InverseLeftCauchyGreen( 0, 1 ); // xy
    rStrainVector[4] = - InverseLeftCauchyGreen( 1, 2 ); // yz
    rStrainVector[5] = - InverseLeftCauchyGreen( 0, 2 ); // xz
}
void HenckyElasticPlasticUP3DLaw::GetLawFeatures(Features& rFeatures)
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


} // namespace Kratos

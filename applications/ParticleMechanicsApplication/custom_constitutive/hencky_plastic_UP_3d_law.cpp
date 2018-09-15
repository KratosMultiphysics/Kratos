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
#include <cmath>

// External includes

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
    const double& Nu          = MaterialProperties[POISSON_RATIO];
    const double ShearModulus = Young/(2*(1 + Nu));

    // Calculate the deviatoric elastic streches eigenvalues
    Vector MainStrain      = ZeroVector(3);
    for (unsigned int i = 0; i<3; ++i)
    {
        MainStrain[i] = rNewElasticLeftCauchyGreen(i,i);
    }

    Vector DeviatoricMainStrain = ZeroVector(3);
    Vector DeviatoricPrincipalStress = ZeroVector(3);

    // First, calculate hydrostatic strain
    const double TracePrincipalStrain = MainStrain[0] + MainStrain[1] + MainStrain[2];
    const double HydrostaticStrain = TracePrincipalStrain/3.0;
    for (unsigned int i=0; i<3 ; i++)
    {
        DeviatoricMainStrain[i] = MainStrain[i] - HydrostaticStrain;
        DeviatoricPrincipalStress[i] = 2.0 * ShearModulus * DeviatoricMainStrain[i];
    }

    // We have to transform the principal deviatoric stress in cartesian stress
    Vector auxN = ZeroVector(3);
    Matrix auxM = ZeroMatrix(3,3);
    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            auxN(j) = rReturnMappingVariables.MainDirections(i,j);
        }
        auxM = MathUtils<double>::TensorProduct3(auxN, auxN);
        rStressMatrix += DeviatoricPrincipalStress(i)*auxM;
    }

    double Pressure = 0;
    GetDomainPressure( Pressure, rElasticVariables);

    for (unsigned int i = 0; i < 3; ++i)
        rStressMatrix(i,i) += Pressure * rElasticVariables.DeterminantF;

    // Now We have to apply the spectral theorem
    Matrix EigenVectors  = ZeroMatrix(3,3);
    Vector EigenValues   = ZeroVector(3);

    double tol = 1e-9;
    int iter = 100;
    SolidMechanicsMathUtilities<double>::EigenVectors(rStressMatrix, EigenVectors, EigenValues, tol, iter);

    rStressMatrix.clear();
    for(unsigned int i=0; i<3; i++)
    {
        rStressMatrix(i,i) = EigenValues(i);
    }

}

void HenckyElasticPlasticUP3DLaw::CorrectDomainPressure( Matrix& rStressMatrix, const MaterialResponseVariables & rElasticVariables)
{
    // Take out the hydrostatic term from stress matrix
    double MeanPressure = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
        MeanPressure += rStressMatrix(i,i);
    MeanPressure /=3.0;
    
    for (unsigned int i = 0; i < 3; ++i)
        rStressMatrix(i,i) -= MeanPressure;

    // Get New Pressure from interpolation and add to diagonal term of stress matrix
    double Pressure = 0;
    GetDomainPressure( Pressure, rElasticVariables);

    for (unsigned int i = 0; i < 3; ++i)
        rStressMatrix(i,i) += Pressure * rElasticVariables.DeterminantF;

}

void HenckyElasticPlasticUP3DLaw::GetDomainPressure( double& rPressure, const MaterialResponseVariables& rElasticVariables)
{
    // Interpolate Pressure from nodes to particle quadrature
    rPressure = 0.0;
    const GeometryType&  DomainGeometry =  rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues  =  rElasticVariables.GetShapeFunctionsValues();

    const unsigned int number_of_nodes  =  DomainGeometry.size();

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        rPressure += ShapeFunctionsValues[j] * DomainGeometry[j].FastGetSolutionStepValue(PRESSURE);
    }

}

void HenckyElasticPlasticUP3DLaw::CalculateElastoPlasticTangentMatrix( const MPMFlowRule::RadialReturnVariables & rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& rAlpha, Matrix& rElastoPlasticTangentMatrix, const MaterialResponseVariables& rElasticVariables )
{
    mpMPMFlowRule->ComputeElastoPlasticTangentMatrix( rReturnMappingVariables,  rNewElasticLeftCauchyGreen, rAlpha, rElastoPlasticTangentMatrix);
    
    // Obtain the domain pressure
    double Pressure;
    GetDomainPressure( Pressure, rElasticVariables);

    Pressure *= rElasticVariables.DeterminantF;

    // Material parameters
    const double Young = mpYieldCriterion->GetHardeningLaw().GetProperties()[YOUNG_MODULUS];
    const double Nu    = mpYieldCriterion->GetHardeningLaw().GetProperties()[POISSON_RATIO];

    // Bulk modulus
    double BulkModulus = Young / (3.0 * (1.0 - 2.0*Nu));

    // Check if Bulk Modulus is not NaN
    if (BulkModulus != BulkModulus)
        BulkModulus = 1.e16;

    // Subtract the Dep with Bulk Modulus to obtain Dep_deviatoric
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = 0; j < 3 ; ++j)
        {
            // TODO: Check whether this is correct or not
            rElastoPlasticTangentMatrix(i,j)  -= BulkModulus;
        }
    }

    // Adding the pressure contribution
    Matrix FourthOrderIdentity = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i)
        FourthOrderIdentity(i,i) = 1.0;

    for (unsigned int i = 3; i<6; ++i)
        FourthOrderIdentity(i,i) = 0.50;

    Matrix IdentityCross = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            IdentityCross(i,j) = 1.0;
        }
    }

    rElastoPlasticTangentMatrix += Pressure * ( IdentityCross - 2.0 * FourthOrderIdentity);
}

//***********************COMPUTE TOTAL STRAIN*****************************************
//************************************************************************************

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

void HenckyElasticPlasticUP3DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
        Vector& rStrainVector )
{
    // E = 0.5*(1-invFT*invF) or e = 0.5*(1-inv(b))
    
    // Calculating the inverse of the jacobian
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

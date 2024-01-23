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
#include "custom_constitutive/hencky_plastic_UP_3D_law.hpp"
#include "custom_utilities/mpm_math_utilities.h"
#include "mpm_application_variables.h"

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

    const Properties& material_properties  = rValues.GetMaterialProperties();
    const double& young_modulus            = material_properties[YOUNG_MODULUS];
    const double& poisson_ratio            = material_properties[POISSON_RATIO];
    const double shear_modulus             = young_modulus/(2*(1 + poisson_ratio));

    // Calculate the deviatoric elastic streches eigen_values
    Vector main_strain      = ZeroVector(3);
    for (unsigned int i = 0; i<3; ++i)
    {
        main_strain[i] = rNewElasticLeftCauchyGreen(i,i);
    }

    Vector deviatoric_main_strain      = ZeroVector(3);
    Vector deviatoric_principal_stress = ZeroVector(3);

    // First, calculate hydrostatic strain
    const double trace_main_strain = main_strain[0] + main_strain[1] + main_strain[2];
    const double hydrostatic_strain = trace_main_strain/3.0;
    for (unsigned int i=0; i<3 ; i++)
    {
        deviatoric_main_strain[i] = main_strain[i] - hydrostatic_strain;
        deviatoric_principal_stress[i] = 2.0 * shear_modulus * deviatoric_main_strain[i];
    }

    // We have to transform the principal deviatoric stress in cartesian stress
    Vector aux_N = ZeroVector(3);
    Matrix aux_M = ZeroMatrix(3,3);
    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            aux_N[j] = rReturnMappingVariables.MainDirections(i,j);
        }
        aux_M = MathUtils<double>::TensorProduct3(aux_N, aux_N);
        rStressMatrix += deviatoric_principal_stress[i]*aux_M;
    }

    double pressure = 0;
    GetDomainPressure( pressure, rElasticVariables);

    for (unsigned int i = 0; i < 3; ++i)
        rStressMatrix(i,i) += pressure * rElasticVariables.DeterminantF;

    // Now We have to apply the spectral theorem
    Matrix eigen_vectors  = ZeroMatrix(3,3);
    Vector eigen_values   = ZeroVector(3);

    double tol = 1e-9;
    int iter = 100;
    MPMMathUtilities<double>::EigenVectors(rStressMatrix, eigen_vectors, eigen_values, tol, iter);

    rStressMatrix.clear();
    for(unsigned int i=0; i<3; i++)
    {
        rStressMatrix(i,i) = eigen_values[i];
    }

}

void HenckyElasticPlasticUP3DLaw::CorrectDomainPressure( Matrix& rStressMatrix, const MaterialResponseVariables & rElasticVariables)
{
    // Take out the hydrostatic term from stress matrix
    double mean_pressure = 0.0;
    for (unsigned int i = 0; i < 3; ++i)
        mean_pressure += rStressMatrix(i,i);
    mean_pressure /=3.0;

    for (unsigned int i = 0; i < 3; ++i)
        rStressMatrix(i,i) -= mean_pressure;

    // Get New pressure from interpolation and add to diagonal term of stress matrix
    double pressure = 0;
    GetDomainPressure( pressure, rElasticVariables);

    for (unsigned int i = 0; i < 3; ++i)
        rStressMatrix(i,i) += pressure * rElasticVariables.DeterminantF;

}

void HenckyElasticPlasticUP3DLaw::GetDomainPressure( double& rPressure, const MaterialResponseVariables& rElasticVariables)
{
    // Interpolate pressure from nodes to particle quadrature
    rPressure = 0.0;
    const GeometryType&  gomain_geometry =  rElasticVariables.GetElementGeometry();
    const Vector& shape_functions  =  rElasticVariables.GetShapeFunctionsValues();

    const unsigned int number_of_nodes  =  gomain_geometry.size();

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        rPressure += shape_functions[j] * gomain_geometry[j].FastGetSolutionStepValue(PRESSURE);
    }

}

void HenckyElasticPlasticUP3DLaw::CalculateElastoPlasticTangentMatrix( const MPMFlowRule::RadialReturnVariables & rReturnMappingVariables, const Matrix& rNewElasticLeftCauchyGreen, const double& rAlpha, Matrix& rElastoPlasticTangentMatrix, const MaterialResponseVariables& rElasticVariables, const Properties& rProperties )
{
    mpMPMFlowRule->ComputeElastoPlasticTangentMatrix( rReturnMappingVariables,  rNewElasticLeftCauchyGreen, rAlpha, rElastoPlasticTangentMatrix, rProperties);

    // Obtain the domain pressure
    double pressure;
    GetDomainPressure( pressure, rElasticVariables);

    pressure *= rElasticVariables.DeterminantF;

    // Material parameters
    const double young_modulus = rProperties[YOUNG_MODULUS];
    const double poisson_ratio    = rProperties[POISSON_RATIO];

    // Bulk modulus
    double bulk_modulus = young_modulus / (3.0 * (1.0 - 2.0*poisson_ratio));

    // Check if Bulk Modulus is not NaN
    if (bulk_modulus != bulk_modulus)
        bulk_modulus = 1.e16;

    // Subtract the Dep with Bulk Modulus to obtain Dep_deviatoric
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (unsigned int j = 0; j < 3 ; ++j)
        {
            // TODO: Check whether this is correct or not
            rElastoPlasticTangentMatrix(i,j)  -= bulk_modulus;
        }
    }

    // Adding the pressure contribution
    Matrix fourth_order_identity = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i)
        fourth_order_identity(i,i) = 1.0;

    for (unsigned int i = 3; i<6; ++i)
        fourth_order_identity(i,i) = 0.50;

    Matrix identity_cross = ZeroMatrix(6,6);
    for (unsigned int i = 0; i<3; ++i)
    {
        for (unsigned int j = 0; j<3; ++j)
        {
            identity_cross(i,j) = 1.0;
        }
    }

    rElastoPlasticTangentMatrix += pressure * ( identity_cross - 2.0 * fourth_order_identity);
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
    Matrix InverseLeftCauchyGreen = ZeroMatrix( 3, 3 );
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

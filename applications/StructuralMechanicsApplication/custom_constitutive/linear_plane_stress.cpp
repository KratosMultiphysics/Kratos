// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/linear_plane_stress.h"

#include "structural_mechanics_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearPlaneStress::LinearPlaneStress()
    : ElasticIsotropic3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LinearPlaneStress::LinearPlaneStress(const LinearPlaneStress& rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearPlaneStress::Clone() const
{
    LinearPlaneStress::Pointer p_clone(new LinearPlaneStress(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LinearPlaneStress::~LinearPlaneStress()
{
}

//************************************************************************************
//************************************************************************************

bool& LinearPlaneStress::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE)
        rValue = true;
    
    return rValue;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void LinearPlaneStress::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRESS_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = 3;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 2;
}

//************************************************************************************
//************************************************************************************

void LinearPlaneStress::CalculateElasticMatrix(Matrix& C, Parameters& rValues)
{
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    const double& E = MaterialProperties[YOUNG_MODULUS];
    const double& NU = MaterialProperties[POISSON_RATIO];

    double c1 = E / (1.00 - NU * NU);
    double c2 = c1 * NU;
    double c3 = 0.5*E / (1 + NU);

    C(0, 0) = c1;
    C(0, 1) = c2;
    C(0, 2) = 0.0;
    C(1, 0) = c2;
    C(1, 1) = c1;
    C(1, 2) = 0.0;
    C(2, 0) = 0.0;
    C(2, 1) = 0.0;
    C(2, 2) = c3;
}

} // Namespace Kratos

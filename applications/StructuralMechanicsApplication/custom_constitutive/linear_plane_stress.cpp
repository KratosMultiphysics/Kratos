// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/linear_plane_stress.h"
#include "custom_utilities/constitutive_law_utilities.h"
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

/***********************************************************************************/
/***********************************************************************************/

std::string LinearPlaneStress::Info() const
{
    return "LinearPlaneStress ConstitutiveLaw instance";
}

/***********************************************************************************/
/***********************************************************************************/

void LinearPlaneStress::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;
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

void LinearPlaneStress::CalculateElasticMatrix(VoigtSizeMatrixType& rC, ConstitutiveLaw::Parameters& rValues)
{
    const auto &r_props = rValues.GetMaterialProperties();
    const double E = r_props[YOUNG_MODULUS];
    const double NU = r_props[POISSON_RATIO];
    ConstitutiveLawUtilities<3>::CalculateElasticMatrixPlaneStress(rC, E, NU);
}

//************************************************************************************
//************************************************************************************

void LinearPlaneStress::CalculatePK2Stress(
    const ConstitutiveLaw::StrainVectorType& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues
)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    ConstitutiveLawUtilities<3>::CalculatePK2StressFromStrainPlaneStress(rStressVector, rStrainVector, E, NU);
}

//************************************************************************************
//************************************************************************************

void LinearPlaneStress::CalculateCauchyGreenStrain(Parameters& rValues, Vector& rStrainVector)
{
    ConstitutiveLawUtilities<3>::CalculateCauchyGreenStrain(rValues, rStrainVector);
}

} // Namespace Kratos

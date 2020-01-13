// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Philippe Bussetta
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_advanced_constitutive/elastic_isotropic_plane_stress_uncoupled_shear.h"
#include "includes/checks.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ElasticIsotropicPlaneStressUncoupledShear::ElasticIsotropicPlaneStressUncoupledShear()
    : LinearPlaneStress()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ElasticIsotropicPlaneStressUncoupledShear::ElasticIsotropicPlaneStressUncoupledShear(const ElasticIsotropicPlaneStressUncoupledShear& rOther)
    : LinearPlaneStress(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer ElasticIsotropicPlaneStressUncoupledShear::Clone() const
{
    ElasticIsotropicPlaneStressUncoupledShear::Pointer p_clone(new ElasticIsotropicPlaneStressUncoupledShear(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

ElasticIsotropicPlaneStressUncoupledShear::~ElasticIsotropicPlaneStressUncoupledShear()
{
}

//************************************************************************************
//************************************************************************************

int ElasticIsotropicPlaneStressUncoupledShear::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0) << "YOUNG_MODULUS is null or negative." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(POISSON_RATIO);
    const double tolerance = 1.0e-12;
    const double nu_upper_bound = 0.5;
    const double nu_lower_bound = -1.0;
    const double nu = rMaterialProperties[POISSON_RATIO];
    KRATOS_ERROR_IF((nu_upper_bound - nu) < tolerance) << "POISSON_RATIO is above the upper bound 0.5." << std::endl;
    KRATOS_ERROR_IF((nu - nu_lower_bound) < tolerance) << "POISSON_RATIO is below the lower bound -1.0." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0) << "DENSITY is negative." << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(SHEAR_MODULUS);
    KRATOS_ERROR_IF(rMaterialProperties[SHEAR_MODULUS] <= 0.0) << "SHEAR_MODULUS is invalid value " << std::endl;

    KRATOS_CHECK_VARIABLE_KEY(SHEAR_MODULUS_GAMMA12);
    KRATOS_CHECK_VARIABLE_KEY(SHEAR_MODULUS_GAMMA12_2);
    KRATOS_CHECK_VARIABLE_KEY(SHEAR_MODULUS_GAMMA12_3);
    KRATOS_CHECK_VARIABLE_KEY(SHEAR_MODULUS_GAMMA12_4);

    return 0;
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropicPlaneStressUncoupledShear::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];
    const double G = r_material_properties[SHEAR_MODULUS];
    const double G1 = r_material_properties[SHEAR_MODULUS_GAMMA12];
    const double G2 = r_material_properties[SHEAR_MODULUS_GAMMA12_2];
    const double G3 = r_material_properties[SHEAR_MODULUS_GAMMA12_3];
    const double G4 = r_material_properties[SHEAR_MODULUS_GAMMA12_4];

    const Vector& r_strain_vector = rValues.GetStrainVector();
    const double abs_gamma12 = std::abs(r_strain_vector(2));

    this->CheckClearElasticMatrix(C);

    const double c1 = E / (1.00 - NU*NU);
    const double c2 = c1 * NU;
    const double c3 = G + G1 * abs_gamma12 + G2 * std::pow(abs_gamma12,2) + G3 * std::pow(abs_gamma12, 3) + G4 * std::pow(abs_gamma12, 4);

    C(0,0) = c1;
    C(0,1) = c2;
    C(1,0) = c2;
    C(1,1) = c1;
    C(2,2) = c3;
    }

//************************************************************************************
//************************************************************************************

void ElasticIsotropicPlaneStressUncoupledShear::CalculatePK2Stress(
    const Vector& rStrainVector,
    Vector& rStressVector,
    ConstitutiveLaw::Parameters& rValues
)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];
    const double G = r_material_properties[SHEAR_MODULUS];
    const double G1 = r_material_properties[SHEAR_MODULUS_GAMMA12];
    const double G2 = r_material_properties[SHEAR_MODULUS_GAMMA12_2];
    const double G3 = r_material_properties[SHEAR_MODULUS_GAMMA12_3];
    const double G4 = r_material_properties[SHEAR_MODULUS_GAMMA12_4];

    const double abs_gamma12 = std::abs(rStrainVector(2));

    const double c1 = E / (1.00 - NU * NU);
    const double c2 = c1 * NU;
    const double c3 = G + G1 * abs_gamma12 + G2 * std::pow(abs_gamma12,2) + G3 * std::pow(abs_gamma12, 3) + G4 * std::pow(abs_gamma12, 4);

    rStressVector[0] = c1 * rStrainVector[0] + c2 * rStrainVector[1];
    rStressVector[1] = c2 * rStrainVector[0] + c1 * rStrainVector[1];
    rStressVector[2] = c3 * rStrainVector[2];
}

} // Namespace Kratos

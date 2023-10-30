// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "thermal_linear_plane_stress.h"
#include "includes/checks.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void ThermalLinearPlaneStress::CalculatePK2Stress(
    const Vector& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const auto &r_geom = rValues.GetElementGeometry();
    const auto &r_N = rValues.GetShapeFunctionsValues();
    const auto &r_process_info = rValues.GetProcessInfo();
    const double E  = r_material_properties.GetValue(YOUNG_MODULUS, r_geom, r_N, r_process_info);
    const double NU = r_material_properties.GetValue(POISSON_RATIO, r_geom, r_N, r_process_info);
    ConstitutiveLawUtilities<3>::CalculatePK2StressFromStrainPlaneStress(rStressVector, rStrainVector, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void ThermalLinearPlaneStress::CalculateElasticMatrix(
    ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const auto &r_geom = rValues.GetElementGeometry();
    const auto &r_N = rValues.GetShapeFunctionsValues();
    const auto &r_process_info = rValues.GetProcessInfo();
    const double E  = r_material_properties.GetValue(YOUNG_MODULUS, r_geom, r_N, r_process_info);
    const double NU = r_material_properties.GetValue(POISSON_RATIO, r_geom, r_N, r_process_info);
    ConstitutiveLawUtilities<3>::CalculateElasticMatrixPlaneStress(rConstitutiveMatrix, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void ThermalLinearPlaneStress::SubstractThermalStrain(
    ConstitutiveLaw::StrainVectorType &rStrainVector,
    const double ReferenceTemperature,
    ConstitutiveLaw::Parameters &rParameters,
    const bool IsPlaneStrain
    )
{
    AdvancedConstitutiveLawUtilities<3>::SubstractThermalStrain(rStrainVector, ReferenceTemperature, rParameters, false);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos

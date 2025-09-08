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
//

// System includes

// External includes

// Project includes

#include "generic_anisotropic_plane_stress_2d_law.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer GenericAnisotropicPlaneStress2DLaw::Create(Kratos::Parameters NewParameters) const
{
    return Kratos::make_shared<GenericAnisotropicPlaneStress2DLaw>();
}

/***********************************************************************************/
/***********************************************************************************/

void GenericAnisotropicPlaneStress2DLaw::CalculateOrthotropicElasticMatrix(
    BoundedMatrixVoigtType& rElasticityTensor,
    const Properties& rMaterialProperties)
{
    AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateOrthotropicElasticMatrixPlaneStress(rElasticityTensor, rMaterialProperties);
}

/***********************************************************************************/
/***********************************************************************************/

} // namespace Kratos

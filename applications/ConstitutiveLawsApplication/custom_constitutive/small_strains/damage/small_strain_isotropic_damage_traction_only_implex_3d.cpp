// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi
//  Collaborators:
//

// Project includes
#include "small_strain_isotropic_damage_traction_only_implex_3d.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainIsotropicDamageTractionOnlyImplex3D::SmallStrainIsotropicDamageTractionOnlyImplex3D()
        : SmallStrainIsotropicDamageImplex3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

SmallStrainIsotropicDamageTractionOnlyImplex3D::SmallStrainIsotropicDamageTractionOnlyImplex3D(
        const SmallStrainIsotropicDamageTractionOnlyImplex3D &rOther) = default;

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainIsotropicDamageTractionOnlyImplex3D::Clone() const
{
    return Kratos::make_shared<SmallStrainIsotropicDamageTractionOnlyImplex3D>(SmallStrainIsotropicDamageTractionOnlyImplex3D(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

SmallStrainIsotropicDamageTractionOnlyImplex3D::~SmallStrainIsotropicDamageTractionOnlyImplex3D() = default;

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageTractionOnlyImplex3D::ComputePositiveStressVector(
        Vector& rStressVectorPos, Vector& rStressVector)
{
    BoundedMatrix<double, 3, 3> stress_matrix;
    BoundedMatrix<double, 3, 3> eigen_values, eigen_vectors;
    stress_matrix = MathUtils<double>::StressVectorToTensor(rStressVector);
    MathUtils<double>::GaussSeidelEigenSystem(stress_matrix, eigen_vectors, eigen_values);
    for (unsigned int i = 0; i < 3; i++)
    {
        if(eigen_values(i, i) < 0.)
        {
            eigen_values(i, i) = 0.;
        }
    }
    MathUtils<double>::BDBtProductOperation(stress_matrix, eigen_values, eigen_vectors);
    rStressVectorPos = MathUtils<double>::StressTensorToVector(stress_matrix);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageTractionOnlyImplex3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainIsotropicDamageImplex3D);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageTractionOnlyImplex3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainIsotropicDamageImplex3D);
}

} /* namespace Kratos.*/

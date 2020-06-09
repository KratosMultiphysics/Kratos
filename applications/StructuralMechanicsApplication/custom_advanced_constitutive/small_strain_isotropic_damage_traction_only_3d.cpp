// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi
//  Collaborators:
//

// Project includes
#include "small_strain_isotropic_damage_traction_only_3d.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallStrainIsotropicDamageTractionOnly3D::SmallStrainIsotropicDamageTractionOnly3D()
        : SmallStrainIsotropicDamage3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

SmallStrainIsotropicDamageTractionOnly3D::SmallStrainIsotropicDamageTractionOnly3D(
        const SmallStrainIsotropicDamageTractionOnly3D &rOther) = default;

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer SmallStrainIsotropicDamageTractionOnly3D::Clone() const
{
    return Kratos::make_shared<SmallStrainIsotropicDamageTractionOnly3D>(SmallStrainIsotropicDamageTractionOnly3D(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

SmallStrainIsotropicDamageTractionOnly3D::~SmallStrainIsotropicDamageTractionOnly3D() = default;

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageTractionOnly3D::ComputePositiveStressVector(
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

void SmallStrainIsotropicDamageTractionOnly3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainIsotropicDamage3D);
}

//************************************************************************************
//************************************************************************************

void SmallStrainIsotropicDamageTractionOnly3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainIsotropicDamage3D);
}

} /* namespace Kratos.*/

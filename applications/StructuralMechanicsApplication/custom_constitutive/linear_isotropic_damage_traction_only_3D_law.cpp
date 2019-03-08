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
#include "linear_isotropic_damage_traction_only_3D_law.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearIsotropicDamageTractionOnly3D::LinearIsotropicDamageTractionOnly3D()
        : SmallStrainIsotropicDamage3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

LinearIsotropicDamageTractionOnly3D::LinearIsotropicDamageTractionOnly3D(
        const LinearIsotropicDamageTractionOnly3D &rOther) = default;

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearIsotropicDamageTractionOnly3D::Clone() const
{
    return Kratos::make_shared<LinearIsotropicDamageTractionOnly3D>(LinearIsotropicDamageTractionOnly3D(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

LinearIsotropicDamageTractionOnly3D::~LinearIsotropicDamageTractionOnly3D() = default;

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamageTractionOnly3D::ComputePositiveStressVector(
        Vector& rStressVectorPos, Vector& rStressVector)
{
    BoundedMatrix<double, 3, 3> stress_matrix;
    BoundedMatrix<double, 3, 3> eigen_values(3, 3), eigen_vectors(3, 3), eigen_vectors_t(3, 3);
    eigen_values.clear();
    eigen_vectors.clear();
    eigen_vectors_t.clear();
    stress_matrix = MathUtils<double>::StressVectorToTensor(rStressVector);
    MathUtils<double>::EigenSystem<3> (stress_matrix, eigen_vectors_t, eigen_values);
    eigen_vectors = trans(eigen_vectors_t); // workaround, eigensystem routine returning transposed vectors
    //KRATOS_WATCH(stress_matrix)
    //KRATOS_WATCH(eigen_vectors)
    //KRATOS_WATCH(eigen_values)
    for (unsigned int i = 0; i < 3; i++)
    {
        if (eigen_values(i, i) < 0.)
        {
            eigen_values(i, i) = 0.;
        }
    }

    BoundedMatrix<double, 3, 3> aux(3, 3);
    aux = prod(eigen_values, trans(eigen_vectors));
    noalias(stress_matrix) = prod(eigen_vectors, aux);
    rStressVectorPos = MathUtils<double>::StressTensorToVector(stress_matrix);
    //KRATOS_WATCH(rStressVectorPos)
    //KRATOS_WATCH(rStressVector)
    return;
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamageTractionOnly3D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainIsotropicDamage3D);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamageTractionOnly3D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainIsotropicDamage3D);
}

} /* namespace Kratos.*/

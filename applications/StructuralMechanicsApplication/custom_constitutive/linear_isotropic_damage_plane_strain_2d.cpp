// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi

#include "linear_isotropic_damage_plane_strain_2d.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LinearIsotropicDamagePlaneStrain2D::LinearIsotropicDamagePlaneStrain2D()
    : LinearIsotropicDamage3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************
LinearIsotropicDamagePlaneStrain2D::LinearIsotropicDamagePlaneStrain2D(
    const LinearIsotropicDamagePlaneStrain2D &rOther)
    : LinearIsotropicDamage3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer LinearIsotropicDamagePlaneStrain2D::Clone() const
{
    return Kratos::make_shared<LinearIsotropicDamagePlaneStrain2D>(LinearIsotropicDamagePlaneStrain2D(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

LinearIsotropicDamagePlaneStrain2D::~LinearIsotropicDamagePlaneStrain2D() = default;

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::CalculateConstitutiveMatrix(
    Matrix &rConstitTensor,
    const Properties &rMaterialProperties
    )
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double nu = rMaterialProperties[POISSON_RATIO];
    const double Ebar = E / (1. - nu * nu);
    const double nubar = nu / (1. - nu);

    rConstitTensor.clear();

    rConstitTensor(0, 0) = 1;
    rConstitTensor(0, 1) = nubar;
    rConstitTensor(0, 2) = 0;
    rConstitTensor(1, 0) = nubar;
    rConstitTensor(1, 1) = 1;
    rConstitTensor(1, 2) = 0;
    rConstitTensor(2, 0) = 0;
    rConstitTensor(2, 1) = 0;
    rConstitTensor(2, 2) = 0.5 * (1 - nubar);

    rConstitTensor *= Ebar / (1. - nubar * nubar);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = WorkingSpaceDimension();
    rFeatures.mSpaceDimension = GetStrainSize();
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LinearIsotropicDamage3D);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LinearIsotropicDamage3D);
}

} /* namespace Kratos.*/

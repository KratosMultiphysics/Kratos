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
    const LinearIsotropicDamagePlaneStrain2D &rOther) = default;

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

void LinearIsotropicDamagePlaneStrain2D::CalculateElasticMatrix(
    const Properties &rMaterialProperties, Matrix &rElasticMatrix)
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double nu = rMaterialProperties[POISSON_RATIO];
    const double Ebar = E / (1. - nu * nu);
    const double nubar = nu / (1. - nu);

    if (rElasticMatrix.size1() != 3 || rElasticMatrix.size2() != 3)
        rElasticMatrix.resize(3, 3, false);
    rElasticMatrix.clear();

    rElasticMatrix(0, 0) = 1;
    rElasticMatrix(0, 1) = nubar;
    rElasticMatrix(0, 2) = 0;
    rElasticMatrix(1, 0) = nubar;
    rElasticMatrix(1, 1) = 1;
    rElasticMatrix(1, 2) = 0;
    rElasticMatrix(2, 0) = 0;
    rElasticMatrix(2, 1) = 0;
    rElasticMatrix(2, 2) = 0.5 * (1 - nubar);

    rElasticMatrix *= Ebar / (1. - nubar * nubar);
}

//************************************************************************************
//************************************************************************************

void LinearIsotropicDamagePlaneStrain2D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = this->WorkingSpaceDimension();
    rFeatures.mSpaceDimension = this->GetStrainSize();
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

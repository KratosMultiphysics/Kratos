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
//  Main authors:    Ruben Zorrilla
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "structural_mechanics_application_variables.h"

// Application includes
#include "custom_constitutive/hyper_elastic_simo_taylor_neo_hookean_plane_strain_2d.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

HyperElasticSimoTaylorNeoHookeanPlaneStrain2D::HyperElasticSimoTaylorNeoHookeanPlaneStrain2D()
    : BaseType()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

HyperElasticSimoTaylorNeoHookeanPlaneStrain2D::HyperElasticSimoTaylorNeoHookeanPlaneStrain2D(const HyperElasticSimoTaylorNeoHookeanPlaneStrain2D& rOther)
    : BaseType(rOther)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer HyperElasticSimoTaylorNeoHookeanPlaneStrain2D::Clone() const
{
    return Kratos::make_shared<HyperElasticSimoTaylorNeoHookeanPlaneStrain2D>(*this);
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

HyperElasticSimoTaylorNeoHookeanPlaneStrain2D::~HyperElasticSimoTaylorNeoHookeanPlaneStrain2D()
{
};

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticSimoTaylorNeoHookeanPlaneStrain2D::CalculateGreenLagrangianStrain(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector)
{
    // Get the total deformation gradient
    const auto& F = rValues.GetDeformationGradientF();

    // E = 0.5*(C - I)
    const SizeType dim = WorkingSpaceDimension();
    Matrix C_tensor(dim, dim);
    noalias(C_tensor) = prod(trans(F), F);

    rStrainVector[0] = 0.5 * (C_tensor(0,0) - 1.0);
    rStrainVector[1] = 0.5 * (C_tensor(1,1) - 1.0);
    rStrainVector[2] = C_tensor(0, 1);
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticSimoTaylorNeoHookeanPlaneStrain2D::AuxiliaryCalculateConstitutiveMatrixPK2(
    Matrix& rConstitutiveMatrix,
    const Vector& rStrain,
    const double Kappa,
    const double Mu)
{
    rConstitutiveMatrix.clear();

    //substitute_PK2_constitutive_matrix
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticSimoTaylorNeoHookeanPlaneStrain2D::AuxiliaryCalculatePK2Stress(
    Vector& rStressVector,
    const Vector& rStrain,
    const double Kappa,
    const double Mu)
{
    rStressVector.clear();

    //substitute_PK2_stress
}

} // Namespace Kratos

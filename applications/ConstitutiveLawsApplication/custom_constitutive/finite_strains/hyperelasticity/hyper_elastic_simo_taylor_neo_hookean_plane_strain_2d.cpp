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
#include "hyper_elastic_simo_taylor_neo_hookean_plane_strain_2d.h"
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
    const auto& rF = rValues.GetDeformationGradientF();

    // E = 0.5*(C - I)
    const SizeType dim = WorkingSpaceDimension();
    Matrix C_tensor(dim, dim);
    noalias(C_tensor) = prod(trans(rF), rF);

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

    const double crConstitutiveMatrix0 = 2*rStrain[1] + 1;
    const double crConstitutiveMatrix1 = std::pow(rStrain[2], 2);
    const double crConstitutiveMatrix2 = 2*rStrain[0];
    const double crConstitutiveMatrix3 = 4*rStrain[0];
    const double crConstitutiveMatrix4 = crConstitutiveMatrix0 + crConstitutiveMatrix2 + crConstitutiveMatrix3*rStrain[1];
    const double crConstitutiveMatrix5 = -crConstitutiveMatrix1 + crConstitutiveMatrix4;
    const double crConstitutiveMatrix6 = std::pow(crConstitutiveMatrix5, 6);
    const double crConstitutiveMatrix7 = 1.0/crConstitutiveMatrix6;
    const double crConstitutiveMatrix8 = Kappa*std::pow(crConstitutiveMatrix5, 4);
    const double crConstitutiveMatrix9 = crConstitutiveMatrix0*crConstitutiveMatrix8;
    const double crConstitutiveMatrix10 = rStrain[0] + rStrain[1] + 1;
    const double crConstitutiveMatrix11 = 3*crConstitutiveMatrix0;
    const double crConstitutiveMatrix12 = -crConstitutiveMatrix10*crConstitutiveMatrix11;
    const double crConstitutiveMatrix13 = 2*crConstitutiveMatrix1;
    const double crConstitutiveMatrix14 = 8*rStrain[0];
    const double crConstitutiveMatrix15 = -crConstitutiveMatrix13 + crConstitutiveMatrix14*rStrain[1] + crConstitutiveMatrix3 + 4*rStrain[1] + 2;
    const double crConstitutiveMatrix16 = Mu*std::pow(crConstitutiveMatrix5, 7.0/2.0);
    const double crConstitutiveMatrix17 = std::pow(crConstitutiveMatrix5, 7);
    const double crConstitutiveMatrix18 = 1.0/crConstitutiveMatrix17;
    const double crConstitutiveMatrix19 = Kappa*crConstitutiveMatrix17;
    const double crConstitutiveMatrix20 = Kappa*crConstitutiveMatrix6;
    const double crConstitutiveMatrix21 = crConstitutiveMatrix2 + 1;
    const double crConstitutiveMatrix22 = Kappa*std::pow(crConstitutiveMatrix5, 5);
    const double crConstitutiveMatrix23 = Mu*crConstitutiveMatrix10*std::pow(crConstitutiveMatrix5, 9.0/2.0);
    const double crConstitutiveMatrix24 = crConstitutiveMatrix18*(crConstitutiveMatrix0*crConstitutiveMatrix21*crConstitutiveMatrix22 + crConstitutiveMatrix19 - crConstitutiveMatrix20 - crConstitutiveMatrix23*(-4*crConstitutiveMatrix1 - crConstitutiveMatrix11*crConstitutiveMatrix21 + crConstitutiveMatrix14 + 16*rStrain[0]*rStrain[1] + 8*rStrain[1] + 4));
    const double crConstitutiveMatrix25 = crConstitutiveMatrix7*rStrain[2];
    const double crConstitutiveMatrix26 = -crConstitutiveMatrix25*(-crConstitutiveMatrix16*(crConstitutiveMatrix12 + crConstitutiveMatrix5) + crConstitutiveMatrix9);
    const double crConstitutiveMatrix27 = crConstitutiveMatrix21*crConstitutiveMatrix8;
    const double crConstitutiveMatrix28 = -3*crConstitutiveMatrix10*crConstitutiveMatrix21;
    const double crConstitutiveMatrix29 = -crConstitutiveMatrix25*(-crConstitutiveMatrix16*(crConstitutiveMatrix28 + crConstitutiveMatrix5) + crConstitutiveMatrix27);
    rConstitutiveMatrix(0,0)=crConstitutiveMatrix0*crConstitutiveMatrix7*(-crConstitutiveMatrix16*(crConstitutiveMatrix12 + crConstitutiveMatrix15) + crConstitutiveMatrix9);
    rConstitutiveMatrix(0,1)=crConstitutiveMatrix24;
    rConstitutiveMatrix(0,2)=crConstitutiveMatrix26;
    rConstitutiveMatrix(1,0)=crConstitutiveMatrix24;
    rConstitutiveMatrix(1,1)=crConstitutiveMatrix21*crConstitutiveMatrix7*(-crConstitutiveMatrix16*(crConstitutiveMatrix15 + crConstitutiveMatrix28) + crConstitutiveMatrix27);
    rConstitutiveMatrix(1,2)=crConstitutiveMatrix29;
    rConstitutiveMatrix(2,0)=crConstitutiveMatrix26;
    rConstitutiveMatrix(2,1)=crConstitutiveMatrix29;
    rConstitutiveMatrix(2,2)=crConstitutiveMatrix18*(crConstitutiveMatrix1*crConstitutiveMatrix22 - 1.0/2.0*crConstitutiveMatrix19 + (1.0/2.0)*crConstitutiveMatrix20 + crConstitutiveMatrix23*(crConstitutiveMatrix13 + crConstitutiveMatrix4));

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

    const double crStressVector0 = 2*rStrain[1] + 1;
    const double crStressVector1 = 2*rStrain[0];
    const double crStressVector2 = crStressVector0 + crStressVector1 + 4*rStrain[0]*rStrain[1] - std::pow(rStrain[2], 2);
    const double crStressVector3 = 1.0/crStressVector2;
    const double crStressVector4 = crStressVector1 + 1;
    const double crStressVector5 = crStressVector0*crStressVector3*crStressVector4 - 2;
    const double crStressVector6 = Mu/std::sqrt(crStressVector2);
    rStressVector[0]=-1.0/2.0*Kappa*crStressVector0*crStressVector3 + (1.0/2.0)*Kappa*crStressVector0 - 1.0/2.0*crStressVector6*(std::pow(crStressVector0, 2)*crStressVector3 + crStressVector5);
    rStressVector[1]=-1.0/2.0*Kappa*crStressVector3*crStressVector4 + (1.0/2.0)*Kappa*crStressVector4 - 1.0/2.0*crStressVector6*(crStressVector3*std::pow(crStressVector4, 2) + crStressVector5);
    rStressVector[2]=(1.0/2.0)*rStrain[2]*(Kappa*crStressVector3 - Kappa + 2*Mu*(rStrain[0] + rStrain[1] + 1)/std::pow(crStressVector2, 3.0/2.0));

}

} // Namespace Kratos

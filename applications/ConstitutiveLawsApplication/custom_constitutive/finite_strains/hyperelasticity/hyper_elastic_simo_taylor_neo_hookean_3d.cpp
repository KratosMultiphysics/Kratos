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
#include "includes/checks.h"
#include "hyper_elastic_simo_taylor_neo_hookean_3d.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

HyperElasticSimoTaylorNeoHookean3D::HyperElasticSimoTaylorNeoHookean3D()
    : BaseType()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

HyperElasticSimoTaylorNeoHookean3D::HyperElasticSimoTaylorNeoHookean3D(const HyperElasticSimoTaylorNeoHookean3D& rOther)
    : BaseType(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer HyperElasticSimoTaylorNeoHookean3D::Clone() const
{
    return Kratos::make_shared<HyperElasticSimoTaylorNeoHookean3D>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

HyperElasticSimoTaylorNeoHookean3D::~HyperElasticSimoTaylorNeoHookean3D()
{
};

/***********************************************************************************/
/***********************************************************************************/

void  HyperElasticSimoTaylorNeoHookean3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    // Get the constitutive law settings
    const auto& r_flags = rValues.GetOptions();

    // The material properties
    const auto& r_material_properties  = rValues.GetMaterialProperties();
    const double young_modulus = r_material_properties[YOUNG_MODULUS];
    const double poisson_coefficient = r_material_properties[POISSON_RATIO];

    // Calculate the required material parameters
    const double kappa = young_modulus/(3.0*(1-2.0*poisson_coefficient));
    const double lame_mu = young_modulus/(2.0*(1.0+poisson_coefficient));

    // Get the strain vector from the constitutive law data
    // This is assumed to be the Green-Lagrange strain in Voigt notation
    auto& r_strain = rValues.GetStrainVector();

    // If the strain is not provided, calculate the Green-Lagrange strain tensor calculation from the input deformation gradient
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        this->CalculateGreenLagrangianStrain(rValues, r_strain);
    }

    if (r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        auto& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        AuxiliaryCalculateConstitutiveMatrixPK2(r_constitutive_matrix, r_strain, kappa, lame_mu);
    }

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        auto& r_stress_vector = rValues.GetStressVector();
        AuxiliaryCalculatePK2Stress(r_stress_vector, r_strain, kappa, lame_mu);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

double& HyperElasticSimoTaylorNeoHookean3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    // The material properties
    const auto& r_material_properties  = rParameterValues.GetMaterialProperties();
    const double young_modulus = r_material_properties[YOUNG_MODULUS];
    const double poisson_coefficient = r_material_properties[POISSON_RATIO];

    // Calculate the required material parameters
    const double kappa = young_modulus/(3.0*(1-2.0*poisson_coefficient));
    const double lame_mu = young_modulus/(2.0*(1.0+poisson_coefficient));

    // Get the deformation gradient data
    const Matrix& r_F = rParameterValues.GetDeformationGradientF();
    const double det_F = rParameterValues.GetDeterminantF();

    if (rThisVariable == STRAIN_ENERGY) {
        // Calculate the deviatoric part of the Cauchy-Green tensor (C):
        const Matrix devC = (1.0/std::pow(det_F,2.0/3.0)) * prod(trans(r_F),r_F);

        // Get the trace of the deviatoric part of C
        double devC_trace = 0.0;
        for (IndexType i = 0; i < devC.size1();++i) {
            devC_trace += devC(i,i);
        }

        // Calculate the Helmholtz free energy
        rValue = 0.25*kappa*(std::pow(det_F,2)-1.0) - 0.5*kappa*std::log(det_F); // Volumetric contribution
        rValue += 0.5*lame_mu*(devC_trace - 3.0); // Deviatoric contribution
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticSimoTaylorNeoHookean3D::AuxiliaryCalculateConstitutiveMatrixPK2(
    Matrix& rConstitutiveMatrix,
    const Vector& rStrain,
    const double Kappa,
    const double Mu)
{
    rConstitutiveMatrix.clear();

    const double crConstitutiveMatrix0 = std::pow(rStrain[4], 2);
    const double crConstitutiveMatrix1 = (2.0/9.0)*rStrain[1];
    const double crConstitutiveMatrix2 = rStrain[1]*rStrain[2];
    const double crConstitutiveMatrix3 = (2.0/9.0)*rStrain[2];
    const double crConstitutiveMatrix4 = crConstitutiveMatrix3 + 1.0/9.0;
    const double crConstitutiveMatrix5 = 2*rStrain[0];
    const double crConstitutiveMatrix6 = std::pow(rStrain[5], 2);
    const double crConstitutiveMatrix7 = 2*rStrain[1];
    const double crConstitutiveMatrix8 = std::pow(rStrain[3], 2);
    const double crConstitutiveMatrix9 = 2*rStrain[2];
    const double crConstitutiveMatrix10 = rStrain[3]*rStrain[5];
    const double crConstitutiveMatrix11 = crConstitutiveMatrix10*rStrain[4];
    const double crConstitutiveMatrix12 = crConstitutiveMatrix2*rStrain[0];
    const double crConstitutiveMatrix13 = crConstitutiveMatrix9 + 1;
    const double crConstitutiveMatrix14 = 4*rStrain[2];
    const double crConstitutiveMatrix15 = -crConstitutiveMatrix0 + crConstitutiveMatrix14*rStrain[1];
    const double crConstitutiveMatrix16 = crConstitutiveMatrix13 + crConstitutiveMatrix15;
    const double crConstitutiveMatrix17 = 4*rStrain[1];
    const double crConstitutiveMatrix18 = crConstitutiveMatrix17*rStrain[0] - crConstitutiveMatrix8;
    const double crConstitutiveMatrix19 = crConstitutiveMatrix5 + crConstitutiveMatrix7;
    const double crConstitutiveMatrix20 = crConstitutiveMatrix14*rStrain[0] - crConstitutiveMatrix6;
    const double crConstitutiveMatrix21 = crConstitutiveMatrix19 + crConstitutiveMatrix20;
    const double crConstitutiveMatrix22 = -crConstitutiveMatrix0*crConstitutiveMatrix5 + 2*crConstitutiveMatrix11 + 8*crConstitutiveMatrix12 + crConstitutiveMatrix16 + crConstitutiveMatrix18 + crConstitutiveMatrix21 - crConstitutiveMatrix6*crConstitutiveMatrix7 - crConstitutiveMatrix8*crConstitutiveMatrix9;
    const double crConstitutiveMatrix23 = std::pow(crConstitutiveMatrix22, 7.0/3.0);
    const double crConstitutiveMatrix24 = 1.0/crConstitutiveMatrix23;
    const double crConstitutiveMatrix25 = crConstitutiveMatrix16 + crConstitutiveMatrix7;
    const double crConstitutiveMatrix26 = Kappa*std::cbrt(crConstitutiveMatrix22);
    const double crConstitutiveMatrix27 = crConstitutiveMatrix25*crConstitutiveMatrix26;
    const double crConstitutiveMatrix28 = 9*crConstitutiveMatrix27;
    const double crConstitutiveMatrix29 = crConstitutiveMatrix19 + crConstitutiveMatrix9 + 3;
    const double crConstitutiveMatrix30 = 12*rStrain[0];
    const double crConstitutiveMatrix31 = -6*crConstitutiveMatrix0*rStrain[0] - 3*crConstitutiveMatrix0 + 6*crConstitutiveMatrix11 + 24*crConstitutiveMatrix12 + 12*crConstitutiveMatrix2 + crConstitutiveMatrix30*rStrain[1] + crConstitutiveMatrix30*rStrain[2] - 6*crConstitutiveMatrix6*rStrain[1] - 3*crConstitutiveMatrix6 - 6*crConstitutiveMatrix8*rStrain[2] - 3*crConstitutiveMatrix8 + 6*rStrain[0] + 6*rStrain[1] + 6*rStrain[2] + 3;
    const double crConstitutiveMatrix32 = 4*Mu;
    const double crConstitutiveMatrix33 = Kappa*crConstitutiveMatrix13;
    const double crConstitutiveMatrix34 = crConstitutiveMatrix23*crConstitutiveMatrix33;
    const double crConstitutiveMatrix35 = std::pow(crConstitutiveMatrix22, 4.0/3.0);
    const double crConstitutiveMatrix36 = crConstitutiveMatrix33*crConstitutiveMatrix35;
    const double crConstitutiveMatrix37 = crConstitutiveMatrix13 + crConstitutiveMatrix20 + crConstitutiveMatrix5;
    const double crConstitutiveMatrix38 = 4*crConstitutiveMatrix25*crConstitutiveMatrix29;
    const double crConstitutiveMatrix39 = crConstitutiveMatrix5 + 1;
    const double crConstitutiveMatrix40 = crConstitutiveMatrix13*crConstitutiveMatrix39;
    const double crConstitutiveMatrix41 = crConstitutiveMatrix7 + 1;
    const double crConstitutiveMatrix42 = crConstitutiveMatrix13*crConstitutiveMatrix41 + crConstitutiveMatrix15 + 2;
    const double crConstitutiveMatrix43 = 3*crConstitutiveMatrix22;
    const double crConstitutiveMatrix44 = (2.0/9.0)*Mu;
    const double crConstitutiveMatrix45 = crConstitutiveMatrix24*(crConstitutiveMatrix27*crConstitutiveMatrix37 + crConstitutiveMatrix34 - crConstitutiveMatrix36 + crConstitutiveMatrix44*(crConstitutiveMatrix37*crConstitutiveMatrix38 - crConstitutiveMatrix43*(std::pow(crConstitutiveMatrix13, 2) + crConstitutiveMatrix14 + crConstitutiveMatrix21 + crConstitutiveMatrix40 + crConstitutiveMatrix42)));
    const double crConstitutiveMatrix46 = Kappa*crConstitutiveMatrix41;
    const double crConstitutiveMatrix47 = crConstitutiveMatrix23*crConstitutiveMatrix46;
    const double crConstitutiveMatrix48 = crConstitutiveMatrix35*crConstitutiveMatrix46;
    const double crConstitutiveMatrix49 = crConstitutiveMatrix18 + crConstitutiveMatrix7;
    const double crConstitutiveMatrix50 = crConstitutiveMatrix39 + crConstitutiveMatrix49;
    const double crConstitutiveMatrix51 = crConstitutiveMatrix39*crConstitutiveMatrix41 + crConstitutiveMatrix9;
    const double crConstitutiveMatrix52 = crConstitutiveMatrix24*(crConstitutiveMatrix27*crConstitutiveMatrix50 + crConstitutiveMatrix44*(crConstitutiveMatrix38*crConstitutiveMatrix50 - crConstitutiveMatrix43*(crConstitutiveMatrix17 + crConstitutiveMatrix18 + std::pow(crConstitutiveMatrix41, 2) + crConstitutiveMatrix42 + crConstitutiveMatrix5 + crConstitutiveMatrix51)) + crConstitutiveMatrix47 - crConstitutiveMatrix48);
    const double crConstitutiveMatrix53 = (1.0/9.0)*rStrain[3];
    const double crConstitutiveMatrix54 = (1.0/9.0)*rStrain[5];
    const double crConstitutiveMatrix55 = crConstitutiveMatrix3*rStrain[3] + crConstitutiveMatrix53 - crConstitutiveMatrix54*rStrain[4];
    const double crConstitutiveMatrix56 = 2*Mu;
    const double crConstitutiveMatrix57 = crConstitutiveMatrix24*(crConstitutiveMatrix28 + crConstitutiveMatrix56*(4*crConstitutiveMatrix25*crConstitutiveMatrix29 - crConstitutiveMatrix31));
    const double crConstitutiveMatrix58 = -crConstitutiveMatrix55*crConstitutiveMatrix57;
    const double crConstitutiveMatrix59 = Kappa*rStrain[4];
    const double crConstitutiveMatrix60 = crConstitutiveMatrix23*crConstitutiveMatrix59;
    const double crConstitutiveMatrix61 = crConstitutiveMatrix35*crConstitutiveMatrix59;
    const double crConstitutiveMatrix62 = -crConstitutiveMatrix10 + crConstitutiveMatrix5*rStrain[4] + rStrain[4];
    const double crConstitutiveMatrix63 = -crConstitutiveMatrix24*(crConstitutiveMatrix27*crConstitutiveMatrix62 + crConstitutiveMatrix44*(crConstitutiveMatrix38*crConstitutiveMatrix62 - crConstitutiveMatrix43*(crConstitutiveMatrix13*rStrain[4] + crConstitutiveMatrix39*rStrain[4] + crConstitutiveMatrix41*rStrain[4] + crConstitutiveMatrix62)) + crConstitutiveMatrix60 - crConstitutiveMatrix61);
    const double crConstitutiveMatrix64 = crConstitutiveMatrix1*rStrain[5] - crConstitutiveMatrix53*rStrain[4] + crConstitutiveMatrix54;
    const double crConstitutiveMatrix65 = -crConstitutiveMatrix57*crConstitutiveMatrix64;
    const double crConstitutiveMatrix66 = (2.0/9.0)*rStrain[0];
    const double crConstitutiveMatrix67 = (4.0/9.0)*rStrain[0];
    const double crConstitutiveMatrix68 = crConstitutiveMatrix26*crConstitutiveMatrix37;
    const double crConstitutiveMatrix69 = 9*crConstitutiveMatrix68;
    const double crConstitutiveMatrix70 = Kappa*crConstitutiveMatrix39;
    const double crConstitutiveMatrix71 = crConstitutiveMatrix23*crConstitutiveMatrix70;
    const double crConstitutiveMatrix72 = crConstitutiveMatrix35*crConstitutiveMatrix70;
    const double crConstitutiveMatrix73 = 4*crConstitutiveMatrix29*crConstitutiveMatrix37;
    const double crConstitutiveMatrix74 = crConstitutiveMatrix24*(crConstitutiveMatrix44*(-crConstitutiveMatrix43*(crConstitutiveMatrix20 + std::pow(crConstitutiveMatrix39, 2) + crConstitutiveMatrix40 + crConstitutiveMatrix49 + crConstitutiveMatrix51 + 4*rStrain[0] + 2) + crConstitutiveMatrix50*crConstitutiveMatrix73) + crConstitutiveMatrix50*crConstitutiveMatrix68 + crConstitutiveMatrix71 - crConstitutiveMatrix72);
    const double crConstitutiveMatrix75 = crConstitutiveMatrix24*(crConstitutiveMatrix56*(4*crConstitutiveMatrix29*crConstitutiveMatrix37 - crConstitutiveMatrix31) + crConstitutiveMatrix69);
    const double crConstitutiveMatrix76 = -crConstitutiveMatrix55*crConstitutiveMatrix75;
    const double crConstitutiveMatrix77 = -crConstitutiveMatrix53*rStrain[5] + crConstitutiveMatrix66*rStrain[4] + (1.0/9.0)*rStrain[4];
    const double crConstitutiveMatrix78 = -crConstitutiveMatrix75*crConstitutiveMatrix77;
    const double crConstitutiveMatrix79 = Kappa*rStrain[5];
    const double crConstitutiveMatrix80 = crConstitutiveMatrix23*crConstitutiveMatrix79;
    const double crConstitutiveMatrix81 = crConstitutiveMatrix35*crConstitutiveMatrix79;
    const double crConstitutiveMatrix82 = crConstitutiveMatrix7*rStrain[5] - rStrain[3]*rStrain[4] + rStrain[5];
    const double crConstitutiveMatrix83 = -crConstitutiveMatrix24*(crConstitutiveMatrix44*(-crConstitutiveMatrix43*(crConstitutiveMatrix13*rStrain[5] + crConstitutiveMatrix39*rStrain[5] + crConstitutiveMatrix41*rStrain[5] + crConstitutiveMatrix82) + crConstitutiveMatrix73*crConstitutiveMatrix82) + crConstitutiveMatrix68*crConstitutiveMatrix82 + crConstitutiveMatrix80 - crConstitutiveMatrix81);
    const double crConstitutiveMatrix84 = crConstitutiveMatrix26*crConstitutiveMatrix50;
    const double crConstitutiveMatrix85 = 9*crConstitutiveMatrix84;
    const double crConstitutiveMatrix86 = Kappa*rStrain[3];
    const double crConstitutiveMatrix87 = crConstitutiveMatrix23*crConstitutiveMatrix86;
    const double crConstitutiveMatrix88 = crConstitutiveMatrix35*crConstitutiveMatrix86;
    const double crConstitutiveMatrix89 = crConstitutiveMatrix9*rStrain[3] + rStrain[3] - rStrain[4]*rStrain[5];
    const double crConstitutiveMatrix90 = -crConstitutiveMatrix24*(crConstitutiveMatrix44*(4*crConstitutiveMatrix29*crConstitutiveMatrix50*crConstitutiveMatrix89 - crConstitutiveMatrix43*(crConstitutiveMatrix13*rStrain[3] + crConstitutiveMatrix39*rStrain[3] + crConstitutiveMatrix41*rStrain[3] + crConstitutiveMatrix89)) + crConstitutiveMatrix84*crConstitutiveMatrix89 + crConstitutiveMatrix87 - crConstitutiveMatrix88);
    const double crConstitutiveMatrix91 = crConstitutiveMatrix24*(crConstitutiveMatrix56*(4*crConstitutiveMatrix29*crConstitutiveMatrix50 - crConstitutiveMatrix31) + crConstitutiveMatrix85);
    const double crConstitutiveMatrix92 = -crConstitutiveMatrix77*crConstitutiveMatrix91;
    const double crConstitutiveMatrix93 = -crConstitutiveMatrix64*crConstitutiveMatrix91;
    const double crConstitutiveMatrix94 = std::pow(crConstitutiveMatrix89, 2);
    const double crConstitutiveMatrix95 = (1.0/9.0)*Mu*crConstitutiveMatrix29;
    const double crConstitutiveMatrix96 = crConstitutiveMatrix26*crConstitutiveMatrix89;
    const double crConstitutiveMatrix97 = 8*crConstitutiveMatrix89;
    const double crConstitutiveMatrix98 = crConstitutiveMatrix24*(crConstitutiveMatrix62*crConstitutiveMatrix96 + (1.0/2.0)*crConstitutiveMatrix80 - 1.0/2.0*crConstitutiveMatrix81 - crConstitutiveMatrix95*(crConstitutiveMatrix43*rStrain[5] - crConstitutiveMatrix62*crConstitutiveMatrix97));
    const double crConstitutiveMatrix99 = crConstitutiveMatrix24*((1.0/2.0)*crConstitutiveMatrix60 - 1.0/2.0*crConstitutiveMatrix61 + crConstitutiveMatrix82*crConstitutiveMatrix96 - crConstitutiveMatrix95*(crConstitutiveMatrix43*rStrain[4] - crConstitutiveMatrix82*crConstitutiveMatrix97));
    const double crConstitutiveMatrix100 = std::pow(crConstitutiveMatrix62, 2);
    const double crConstitutiveMatrix101 = crConstitutiveMatrix62*crConstitutiveMatrix82;
    const double crConstitutiveMatrix102 = crConstitutiveMatrix24*(crConstitutiveMatrix101*crConstitutiveMatrix26 + (1.0/2.0)*crConstitutiveMatrix87 - 1.0/2.0*crConstitutiveMatrix88 - crConstitutiveMatrix95*(-8*crConstitutiveMatrix101 + crConstitutiveMatrix43*rStrain[3]));
    const double crConstitutiveMatrix103 = std::pow(crConstitutiveMatrix82, 2);
    rConstitutiveMatrix(0,0)=crConstitutiveMatrix24*(crConstitutiveMatrix28 + crConstitutiveMatrix32*(2*crConstitutiveMatrix25*crConstitutiveMatrix29 - crConstitutiveMatrix31))*(-1.0/9.0*crConstitutiveMatrix0 + crConstitutiveMatrix1 + (4.0/9.0)*crConstitutiveMatrix2 + crConstitutiveMatrix4);
    rConstitutiveMatrix(0,1)=crConstitutiveMatrix45;
    rConstitutiveMatrix(0,2)=crConstitutiveMatrix52;
    rConstitutiveMatrix(0,3)=crConstitutiveMatrix58;
    rConstitutiveMatrix(0,4)=crConstitutiveMatrix63;
    rConstitutiveMatrix(0,5)=crConstitutiveMatrix65;
    rConstitutiveMatrix(1,0)=crConstitutiveMatrix45;
    rConstitutiveMatrix(1,1)=crConstitutiveMatrix24*(crConstitutiveMatrix32*(2*crConstitutiveMatrix29*crConstitutiveMatrix37 - crConstitutiveMatrix31) + crConstitutiveMatrix69)*(crConstitutiveMatrix4 - 1.0/9.0*crConstitutiveMatrix6 + crConstitutiveMatrix66 + crConstitutiveMatrix67*rStrain[2]);
    rConstitutiveMatrix(1,2)=crConstitutiveMatrix74;
    rConstitutiveMatrix(1,3)=crConstitutiveMatrix76;
    rConstitutiveMatrix(1,4)=crConstitutiveMatrix78;
    rConstitutiveMatrix(1,5)=crConstitutiveMatrix83;
    rConstitutiveMatrix(2,0)=crConstitutiveMatrix52;
    rConstitutiveMatrix(2,1)=crConstitutiveMatrix74;
    rConstitutiveMatrix(2,2)=crConstitutiveMatrix24*(crConstitutiveMatrix32*(2*crConstitutiveMatrix29*crConstitutiveMatrix50 - crConstitutiveMatrix31) + crConstitutiveMatrix85)*(crConstitutiveMatrix1 + crConstitutiveMatrix66 + crConstitutiveMatrix67*rStrain[1] - 1.0/9.0*crConstitutiveMatrix8 + 1.0/9.0);
    rConstitutiveMatrix(2,3)=crConstitutiveMatrix90;
    rConstitutiveMatrix(2,4)=crConstitutiveMatrix92;
    rConstitutiveMatrix(2,5)=crConstitutiveMatrix93;
    rConstitutiveMatrix(3,0)=crConstitutiveMatrix58;
    rConstitutiveMatrix(3,1)=crConstitutiveMatrix76;
    rConstitutiveMatrix(3,2)=crConstitutiveMatrix90;
    rConstitutiveMatrix(3,3)=crConstitutiveMatrix24*(crConstitutiveMatrix26*crConstitutiveMatrix94 - 1.0/2.0*crConstitutiveMatrix34 + (1.0/2.0)*crConstitutiveMatrix36 + crConstitutiveMatrix95*(crConstitutiveMatrix13*crConstitutiveMatrix43 + 8*crConstitutiveMatrix94));
    rConstitutiveMatrix(3,4)=crConstitutiveMatrix98;
    rConstitutiveMatrix(3,5)=crConstitutiveMatrix99;
    rConstitutiveMatrix(4,0)=crConstitutiveMatrix63;
    rConstitutiveMatrix(4,1)=crConstitutiveMatrix78;
    rConstitutiveMatrix(4,2)=crConstitutiveMatrix92;
    rConstitutiveMatrix(4,3)=crConstitutiveMatrix98;
    rConstitutiveMatrix(4,4)=crConstitutiveMatrix24*(crConstitutiveMatrix100*crConstitutiveMatrix26 - 1.0/2.0*crConstitutiveMatrix71 + (1.0/2.0)*crConstitutiveMatrix72 + crConstitutiveMatrix95*(8*crConstitutiveMatrix100 + crConstitutiveMatrix39*crConstitutiveMatrix43));
    rConstitutiveMatrix(4,5)=crConstitutiveMatrix102;
    rConstitutiveMatrix(5,0)=crConstitutiveMatrix65;
    rConstitutiveMatrix(5,1)=crConstitutiveMatrix83;
    rConstitutiveMatrix(5,2)=crConstitutiveMatrix93;
    rConstitutiveMatrix(5,3)=crConstitutiveMatrix99;
    rConstitutiveMatrix(5,4)=crConstitutiveMatrix102;
    rConstitutiveMatrix(5,5)=crConstitutiveMatrix24*(crConstitutiveMatrix103*crConstitutiveMatrix26 - 1.0/2.0*crConstitutiveMatrix47 + (1.0/2.0)*crConstitutiveMatrix48 + crConstitutiveMatrix95*(8*crConstitutiveMatrix103 + crConstitutiveMatrix41*crConstitutiveMatrix43));

}

/***********************************************************************************/
/***********************************************************************************/

void HyperElasticSimoTaylorNeoHookean3D::AuxiliaryCalculatePK2Stress(
    Vector& rStressVector,
    const Vector& rStrain,
    const double Kappa,
    const double Mu)
{
    rStressVector.clear();

    const double crStressVector0 = 2*rStrain[2];
    const double crStressVector1 = 2*rStrain[1];
    const double crStressVector2 = crStressVector1 + 1;
    const double crStressVector3 = std::pow(rStrain[4], 2);
    const double crStressVector4 = 4*rStrain[1];
    const double crStressVector5 = -crStressVector3 + crStressVector4*rStrain[2];
    const double crStressVector6 = crStressVector0 + crStressVector2 + crStressVector5;
    const double crStressVector7 = 2*rStrain[0];
    const double crStressVector8 = std::pow(rStrain[5], 2);
    const double crStressVector9 = std::pow(rStrain[3], 2);
    const double crStressVector10 = rStrain[4]*rStrain[5];
    const double crStressVector11 = rStrain[0]*rStrain[2];
    const double crStressVector12 = crStressVector0 + crStressVector1;
    const double crStressVector13 = crStressVector7 + 1;
    const double crStressVector14 = 4*crStressVector11 + crStressVector13 - crStressVector8;
    const double crStressVector15 = crStressVector4*rStrain[0] - crStressVector9;
    const double crStressVector16 = -crStressVector0*crStressVector9 - crStressVector1*crStressVector8 + 2*crStressVector10*rStrain[3] + 8*crStressVector11*rStrain[1] + crStressVector12 + crStressVector14 + crStressVector15 - crStressVector3*crStressVector7 + crStressVector5;
    const double crStressVector17 = 1.0/crStressVector16;
    const double crStressVector18 = 3*Kappa;
    const double crStressVector19 = crStressVector17*crStressVector6;
    const double crStressVector20 = crStressVector0 + 1;
    const double crStressVector21 = 2*Mu;
    const double crStressVector22 = crStressVector21/std::cbrt(crStressVector16);
    const double crStressVector23 = crStressVector0 + crStressVector14;
    const double crStressVector24 = crStressVector17*crStressVector23;
    const double crStressVector25 = crStressVector1 + crStressVector13 + crStressVector15;
    const double crStressVector26 = crStressVector17*crStressVector25;
    const double crStressVector27 = (1.0/6.0)*crStressVector17*crStressVector18 - 1.0/6.0*crStressVector18 + (1.0/6.0)*crStressVector21*(crStressVector12 + crStressVector7 + 3)/std::pow(crStressVector16, 4.0/3.0);
    rStressVector[0]=(1.0/2.0)*Kappa*crStressVector6 - 1.0/6.0*crStressVector17*crStressVector18*crStressVector6 - 1.0/6.0*crStressVector22*(crStressVector13*crStressVector19 + crStressVector19*crStressVector2 + crStressVector19*crStressVector20 - 3);
    rStressVector[1]=(1.0/2.0)*Kappa*crStressVector23 - 1.0/6.0*crStressVector17*crStressVector18*crStressVector23 - 1.0/6.0*crStressVector22*(crStressVector13*crStressVector24 + crStressVector2*crStressVector24 + crStressVector20*crStressVector24 - 3);
    rStressVector[2]=(1.0/2.0)*Kappa*crStressVector25 - 1.0/6.0*crStressVector17*crStressVector18*crStressVector25 - 1.0/6.0*crStressVector22*(crStressVector13*crStressVector26 + crStressVector2*crStressVector26 + crStressVector20*crStressVector26 - 3);
    rStressVector[3]=crStressVector27*(crStressVector0*rStrain[3] - crStressVector10 + rStrain[3]);
    rStressVector[4]=crStressVector27*(crStressVector7*rStrain[4] - rStrain[3]*rStrain[5] + rStrain[4]);
    rStressVector[5]=crStressVector27*(crStressVector1*rStrain[5] - rStrain[3]*rStrain[4] + rStrain[5]);

}

} // Namespace Kratos

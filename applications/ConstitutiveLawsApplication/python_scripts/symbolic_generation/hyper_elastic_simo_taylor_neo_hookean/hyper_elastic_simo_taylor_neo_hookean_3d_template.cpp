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
#include "custom_constitutive/hyper_elastic_simo_taylor_neo_hookean_3d.h"
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

    //substitute_PK2_constitutive_matrix
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

    //substitute_PK2_stress
}

} // Namespace Kratos

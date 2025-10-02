// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//
// System includes

// External includes

// Project includes
#include "reissner_mindlin_shell_elastic_constitutive_law.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ReissnerMindlinShellElasticConstitutiveLaw::ReissnerMindlinShellElasticConstitutiveLaw()
    : BeamConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ReissnerMindlinShellElasticConstitutiveLaw::ReissnerMindlinShellElasticConstitutiveLaw(const ReissnerMindlinShellElasticConstitutiveLaw& rOther)
    : BeamConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer ReissnerMindlinShellElasticConstitutiveLaw::Clone() const
{
    ReissnerMindlinShellElasticConstitutiveLaw::Pointer p_clone(new ReissnerMindlinShellElasticConstitutiveLaw(*this));
    return p_clone;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void ReissnerMindlinShellElasticConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize = 8;
    rFeatures.mSpaceDimension = 3;
}

//************************************************************************************
//************************************************************************************

int ReissnerMindlinShellElasticConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    // KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS))    << "YOUNG_MODULUS is not defined in the properties" << std::endl;
    // KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO))    << "POISSON_RATIO is not defined in the properties" << std::endl;
    // KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(CROSS_AREA))       << "CROSS_AREA is not defined in the properties"    << std::endl;
    // KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(AREA_EFFECTIVE_Y)) << "AREA_EFFECTIVE_Y is not defined in the properties" << std::endl;
    // KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(I33))              << "I33 is not defined in the properties"            << std::endl;
    // KRATOS_ERROR_IF_NOT(rMaterialProperties[YOUNG_MODULUS] > 0.0)    << "The YOUNG_MODULUS value is lower than 0.0" << std::endl;
    // KRATOS_ERROR_IF_NOT(rMaterialProperties[CROSS_AREA] > 0.0)       << "The CROSS_AREA value is lower than 0.0" << std::endl;
    // KRATOS_ERROR_IF_NOT(rMaterialProperties[I33] > 0.0)              << "The I33 value is lower than 0.0" << std::endl;
    // KRATOS_ERROR_IF_NOT(rMaterialProperties[AREA_EFFECTIVE_Y] > 0.0) << "The AREA_EFFECTIVE_Y value is lower than 0.0" << std::endl;
    // KRATOS_ERROR_IF    (rMaterialProperties[POISSON_RATIO] < 0.0)    << "The POISSON_RATIO value is lower than 0.0" << std::endl;
    // KRATOS_ERROR_IF_NOT(rMaterialProperties[POISSON_RATIO] < 0.5)    << "The POISSON_RATIO cannot be greater than or equal 0.5." << std::endl;
    return 0;
}

//************************************************************************************
//************************************************************************************

void ReissnerMindlinShellElasticConstitutiveLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ReissnerMindlinShellElasticConstitutiveLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ReissnerMindlinShellElasticConstitutiveLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ReissnerMindlinShellElasticConstitutiveLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const auto& r_cl_law_options = rValues.GetOptions();
    auto &r_material_properties = rValues.GetMaterialProperties();
    auto &r_generalized_strain_vector = rValues.GetStrainVector();
    AddInitialStrainVectorContribution(r_generalized_strain_vector);
    const auto strain_size = GetStrainSize();

    const double E  = r_material_properties[YOUNG_MODULUS];
    const double nu = r_material_properties[POISSON_RATIO];
    const double t  = r_material_properties[THICKNESS];

    if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        auto &r_generalized_stress_vector = rValues.GetStressVector();
        if (r_generalized_stress_vector.size() != strain_size)
            r_generalized_stress_vector.resize(strain_size, false);





        AddInitialStressVectorContribution(r_generalized_stress_vector);

        if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            auto &r_stress_derivatives = rValues.GetConstitutiveMatrix(); // dN_dEl, dM_dkappa, dV_dGamma_xy
            if (r_stress_derivatives.size1() != strain_size || r_stress_derivatives.size2() != strain_size)
                r_stress_derivatives.resize(strain_size, strain_size, false);
            r_stress_derivatives.clear();





        }
    }
}

} // Namespace Kratos

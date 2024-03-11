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
#include "timoshenko_beam_elastic_constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TimoshenkoBeamElasticConstitutiveLaw::TimoshenkoBeamElasticConstitutiveLaw()
    : BeamConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

TimoshenkoBeamElasticConstitutiveLaw::TimoshenkoBeamElasticConstitutiveLaw(const TimoshenkoBeamElasticConstitutiveLaw& rOther)
    : BeamConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer TimoshenkoBeamElasticConstitutiveLaw::Clone() const
{
    TimoshenkoBeamElasticConstitutiveLaw::Pointer p_clone(new TimoshenkoBeamElasticConstitutiveLaw(*this));
    return p_clone;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void TimoshenkoBeamElasticConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize = 3;
    rFeatures.mSpaceDimension = 2;
}

//************************************************************************************
//************************************************************************************

int TimoshenkoBeamElasticConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS)) << "YOUNG_MODULUS is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO)) << "POISSON_RATIO is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(CROSS_AREA))    << "CROSS_AREA is not defined in the properties"    << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(DENSITY))       << "DENSITY is not defined in the properties"       << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(IZ))            << "IZ is not defined in the properties"            << std::endl;
    return 0;
}

//************************************************************************************
//************************************************************************************

void TimoshenkoBeamElasticConstitutiveLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElasticConstitutiveLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElasticConstitutiveLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElasticConstitutiveLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const auto& r_cl_law_options = rValues.GetOptions();
    auto &r_material_properties = rValues.GetMaterialProperties();
    auto &r_strain_vector = rValues.GetStrainVector();
    auto strain_size = GetStrainSize();

    const double axial_strain = r_strain_vector[0]; // E_l
    const double curvature    = r_strain_vector[1]; // Kappa
    const double shear_strain = r_strain_vector[2]; // Gamma_xy

    const double E    = r_material_properties[YOUNG_MODULUS];
    const double nu   = r_material_properties[POISSON_RATIO];
    const double A    = r_material_properties[CROSS_AREA];
    const double I    = r_material_properties[IZ];
    const double k_s  = r_material_properties.Has(SHEAR_CORRECTION_XY) ? r_material_properties[SHEAR_CORRECTION_XY] : 5.0 / 6.0; // We assume rectangular shape

    const double G    = E / (2.0 * (1.0 + nu));
    const double A_s  = k_s * A;

    if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        auto &r_generalized_stress_vector = rValues.GetStressVector();
        if (r_generalized_stress_vector.size() != strain_size)
            r_generalized_stress_vector.resize(strain_size, false);

        const double EA  = E * A;
        const double EI  = E * I;
        const double GAs = G * A_s;

        r_generalized_stress_vector[0] = EA * axial_strain;
        r_generalized_stress_vector[1] = EI * curvature;
        r_generalized_stress_vector[2] = GAs * shear_strain;

        if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            auto &r_stress_derivatives = rValues.GetConstitutiveMatrix();
            if (r_stress_derivatives.size1() != strain_size)
                r_stress_derivatives.resize(strain_size, strain_size, false);
            r_stress_derivatives.clear();

            r_stress_derivatives(0, 0) = EA;
            r_stress_derivatives(1, 1) = EI;
            r_stress_derivatives(2, 2) = GAs;
        }
    }
}

} // Namespace Kratos

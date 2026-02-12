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
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS)) << "YOUNG_MODULUS is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO)) << "POISSON_RATIO is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(THICKNESS))     << "THICKNESS is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[YOUNG_MODULUS] > 0.0) << "The YOUNG_MODULUS value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[THICKNESS] > 0.0)     << "The THICKNESS value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF    (rMaterialProperties[POISSON_RATIO] < 0.0) << "The POISSON_RATIO value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[POISSON_RATIO] < 0.5) << "The POISSON_RATIO cannot be greater than or equal 0.5." << std::endl;
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

    if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS) || r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        auto &r_generalized_stress_vector = rValues.GetStressVector();
        if (r_generalized_stress_vector.size() != strain_size)
            r_generalized_stress_vector.resize(strain_size, false);

        const double one_minus_nu2 = 1.0 - nu * nu;
        const double half_one_minus_nu = 0.5 * (1.0 - nu);
        const double t_square = t * t;

        const double Dm_factor = E * t / one_minus_nu2;
        const double Db_factor = E * t * t_square / (12.0 * one_minus_nu2);

        const double k_shear = 5.0 / 6.0; // Shear correction factor
        const double h_max = GetMaxReferenceEdgeLength(rValues.GetElementGeometry());
        const double alpha = 0.1; // This could be an input property
        const double stenberg_stabilization = t_square / (t_square + alpha * h_max * h_max);
        const double G = ConstitutiveLawUtilities<3>::CalculateShearModulus(r_material_properties);
        const double Ds_factor = G * k_shear * t * stenberg_stabilization;

        auto &r_gen_const_matrix = rValues.GetConstitutiveMatrix();
        if (r_gen_const_matrix.size1() != strain_size || r_gen_const_matrix.size2() != strain_size)
            r_gen_const_matrix.resize(strain_size, strain_size, false);
        r_gen_const_matrix.clear();

        // Membrane constitutive matrix
        r_gen_const_matrix(0, 0) = Dm_factor;
        r_gen_const_matrix(0, 1) = Dm_factor * nu;
        r_gen_const_matrix(1, 0) = r_gen_const_matrix(0, 1);
        r_gen_const_matrix(1, 1) = Dm_factor;
        r_gen_const_matrix(2, 2) = Dm_factor * half_one_minus_nu;

        // Bending constitutive matrix
        r_gen_const_matrix(3, 3) = Db_factor;
        r_gen_const_matrix(3, 4) = Db_factor * nu;
        r_gen_const_matrix(4, 3) = r_gen_const_matrix(3, 4);
        r_gen_const_matrix(4, 4) = Db_factor;
        r_gen_const_matrix(5, 5) = Db_factor * half_one_minus_nu;

        // Shear constitutive matrix
        r_gen_const_matrix(6, 6) = Ds_factor;
        r_gen_const_matrix(7, 7) = Ds_factor;

        if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
            noalias(r_generalized_stress_vector) = prod(r_gen_const_matrix, r_generalized_strain_vector);
            AddInitialStressVectorContribution(r_generalized_stress_vector);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

double& ReissnerMindlinShellElasticConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    if (rThisVariable == VON_MISES_STRESS_TOP_SURFACE ||
        rThisVariable == VON_MISES_STRESS_BOTTOM_SURFACE ||
        rThisVariable == VON_MISES_STRESS_MIDDLE_SURFACE) {

        const auto& r_props = rValues.GetMaterialProperties();
        double z = 0.0;

        if (rThisVariable == VON_MISES_STRESS_TOP_SURFACE) {
            z = 0.5 * r_props[THICKNESS];
        } else if (rThisVariable == VON_MISES_STRESS_BOTTOM_SURFACE) {
            z = -0.5 * r_props[THICKNESS];
        }

        Vector strain_plane_stress(3), stress_plane_stress(3);
        Matrix constitutive_matrix(3, 3);

        Vector generalized_strain = rValues.GetStrainVector(); // Filled by the element

        noalias(strain_plane_stress) = boost::numeric::ublas::project(generalized_strain, boost::numeric::ublas::range(0, 3)) +
            z * boost::numeric::ublas::project(generalized_strain, boost::numeric::ublas::range(3, 6)); // gEm + z * gEb

        ConstitutiveLawUtilities<3>::CalculateElasticMatrixPlaneStress(constitutive_matrix, r_props[YOUNG_MODULUS], r_props[POISSON_RATIO]);
        noalias(stress_plane_stress) = prod(constitutive_matrix, strain_plane_stress);
        rValue = ConstitutiveLawUtilities<3>::CalculateVonMisesEquivalentStress(stress_plane_stress);
        return rValue;
    } else if (rThisVariable == VON_MISES_STRESS) {
        double top_value, mid_value, bot_value;

        CalculateValue(rValues, VON_MISES_STRESS_TOP_SURFACE, top_value);
        CalculateValue(rValues, VON_MISES_STRESS_MIDDLE_SURFACE, mid_value);
        CalculateValue(rValues, VON_MISES_STRESS_BOTTOM_SURFACE, bot_value);
        rValue = std::max({top_value, mid_value, bot_value});
        return rValue;
    }
    rValue = 0.0;
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos

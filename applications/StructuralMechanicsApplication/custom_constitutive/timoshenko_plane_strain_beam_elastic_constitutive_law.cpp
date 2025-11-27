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
#include "timoshenko_plane_strain_beam_elastic_constitutive_law.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TimoshenkoBeamPlaneStrainElasticConstitutiveLaw::TimoshenkoBeamPlaneStrainElasticConstitutiveLaw()
    : BaseType()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

TimoshenkoBeamPlaneStrainElasticConstitutiveLaw::TimoshenkoBeamPlaneStrainElasticConstitutiveLaw(const TimoshenkoBeamPlaneStrainElasticConstitutiveLaw& rOther)
    : BaseType(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer TimoshenkoBeamPlaneStrainElasticConstitutiveLaw::Clone() const
{
    TimoshenkoBeamPlaneStrainElasticConstitutiveLaw::Pointer p_clone(new TimoshenkoBeamPlaneStrainElasticConstitutiveLaw(*this));
    return p_clone;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void TimoshenkoBeamPlaneStrainElasticConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize = 5;
    rFeatures.mSpaceDimension = 2;
}

//************************************************************************************
//************************************************************************************

int TimoshenkoBeamPlaneStrainElasticConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS))         << "YOUNG_MODULUS is not defined in the properties"    << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO))         << "POISSON_RATIO is not defined in the properties"    << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(THICKNESS))             << "THICKNESS is not defined in the properties"       << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(THICKNESS_EFFECTIVE_Y)) << "THICKNESS_EFFECTIVE_Y is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[YOUNG_MODULUS] > 0.0)       << "The YOUNG_MODULUS value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[THICKNESS] > 0.0)           << "The THICKNESS value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[THICKNESS_EFFECTIVE_Y] > 0.0) << "The THICKNESS_EFFECTIVE_Y value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF    (rMaterialProperties[POISSON_RATIO] < 0.0)       << "The POISSON_RATIO value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[POISSON_RATIO] < 0.5)       << "The POISSON_RATIO cannot be greater than or equal 0.5." << std::endl;
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamPlaneStrainElasticConstitutiveLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const auto& r_cl_law_options = rValues.GetOptions();
    auto &r_material_properties = rValues.GetMaterialProperties();
    auto &r_strain_vector = rValues.GetStrainVector();
    AddInitialStrainVectorContribution(r_strain_vector);
    const auto strain_size = GetStrainSize();

    const double axial_strain = r_strain_vector[0]; // E_l
    const double curvature    = r_strain_vector[1]; // Kappa
    const double shear_strain = r_strain_vector[2]; // Gamma_xy

    const double thickness = r_material_properties[THICKNESS];
    const double E  = r_material_properties[YOUNG_MODULUS];
    const double A  = thickness; // Per unit length
    const double I  = std::pow(thickness, 3) / 12.0; // Per unit length
    const double nu = r_material_properties[POISSON_RATIO];

    const double G    = ConstitutiveLawUtilities<3>::CalculateShearModulus(r_material_properties);
    const double A_s  = r_material_properties[THICKNESS_EFFECTIVE_Y]; // Per unit length

    if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        auto &r_generalized_stress_vector = rValues.GetStressVector();
        if (r_generalized_stress_vector.size() != strain_size)
            r_generalized_stress_vector.resize(strain_size, false);

        const double one_minus_nu_squared = 1.0 - nu * nu;
        const double EA_nu  = E * A / one_minus_nu_squared;
        const double EI_nu  = E * I / one_minus_nu_squared;
        const double GAs = G * A_s;

        r_generalized_stress_vector[0] = EA_nu * axial_strain; // Nx
        r_generalized_stress_vector[1] = EI_nu * curvature;    // Mz
        r_generalized_stress_vector[2] = GAs * shear_strain;   // Vxy
        // Plane strain generalized stress components
        r_generalized_stress_vector[3] = nu * r_generalized_stress_vector[0]; // Nz
        r_generalized_stress_vector[4] = nu * r_generalized_stress_vector[1]; // Mx

        AddInitialStressVectorContribution(r_generalized_stress_vector);

        if (r_material_properties.Has(BEAM_PRESTRESS_PK2)) {
            r_generalized_stress_vector += r_material_properties[BEAM_PRESTRESS_PK2];
        }

        if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            auto &r_stress_derivatives = rValues.GetConstitutiveMatrix(); // dN_dEl, dM_dkappa, dV_dGamma_xy
            if (r_stress_derivatives.size1() != strain_size || r_stress_derivatives.size2() != strain_size)
                r_stress_derivatives.resize(strain_size, strain_size, false);
            noalias(r_stress_derivatives) = ZeroMatrix(strain_size, strain_size);

            r_stress_derivatives(0, 0) = EA_nu; // dN_dEl
            r_stress_derivatives(1, 1) = EI_nu; // dM_dkappa
            r_stress_derivatives(2, 2) = GAs;   // dV_dGamma_xy
            r_stress_derivatives(3, 3) = nu * EA_nu; // dNz_dEl
            r_stress_derivatives(4, 4) = nu * EI_nu; // dMx_dkappa
        }
    }
}

} // Namespace Kratos

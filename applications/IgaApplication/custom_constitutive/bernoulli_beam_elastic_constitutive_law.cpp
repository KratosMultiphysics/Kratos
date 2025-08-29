// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: IgaApplication/license.txt
//

// System includes

// External includes

// Project includes
#include "bernoulli_beam_elastic_constitutive_law.h"
#include "iga_application_variables.h"

namespace Kratos
{

BernoulliBeamElasticConstitutiveLaw::BernoulliBeamElasticConstitutiveLaw()
    : ConstitutiveLaw()
{
}

BernoulliBeamElasticConstitutiveLaw::BernoulliBeamElasticConstitutiveLaw(const BernoulliBeamElasticConstitutiveLaw& rOther)
    : ConstitutiveLaw(rOther)
{
}

ConstitutiveLaw::Pointer BernoulliBeamElasticConstitutiveLaw::Clone() const
{
    BernoulliBeamElasticConstitutiveLaw::Pointer p_clone(new BernoulliBeamElasticConstitutiveLaw(*this));
    return p_clone;
}

BernoulliBeamElasticConstitutiveLaw::~BernoulliBeamElasticConstitutiveLaw()
{
}

void BernoulliBeamElasticConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize = 5;         
    rFeatures.mSpaceDimension = 3;     
}

int BernoulliBeamElasticConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS)) << "YOUNG_MODULUS is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO)) << "POISSON_RATIO is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[YOUNG_MODULUS] > 0.0) << "The YOUNG_MODULUS value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF    (rMaterialProperties[POISSON_RATIO] < 0.0) << "The POISSON_RATIO value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[POISSON_RATIO] < 0.5) << "The POISSON_RATIO cannot be greater than or equal 0.5." << std::endl;
    
    // Check for IGA beam cross-section properties
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(CROSS_AREA)) << "CROSS_AREA is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(I_V)) << "I_V is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(I_N)) << "I_N is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(I_T)) << "I_T is not defined in the properties" << std::endl;
    
    // Validate property values
    KRATOS_ERROR_IF_NOT(rMaterialProperties[CROSS_AREA] > 0.0) << "CROSS_AREA must be positive" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[I_V] > 0.0) << "I_V must be positive" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[I_N] > 0.0) << "I_N must be positive" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[I_T] > 0.0) << "I_T must be positive" << std::endl;
    
    return 0;
}

void BernoulliBeamElasticConstitutiveLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void BernoulliBeamElasticConstitutiveLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void BernoulliBeamElasticConstitutiveLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

void BernoulliBeamElasticConstitutiveLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const auto& r_cl_law_options = rValues.GetOptions();
    const auto &r_material_properties = rValues.GetMaterialProperties();
    auto &r_strain_vector = rValues.GetStrainVector();

    AddInitialStrainVectorContribution(r_strain_vector);
    const auto strain_size = GetStrainSize();

    const double E = r_material_properties[YOUNG_MODULUS];
    const double nu = r_material_properties[POISSON_RATIO];
    const double G = E / (2.0 * (1.0 + nu));
    
    // Retrieve beam cross-section properties
    const double A = r_material_properties[CROSS_AREA];
    const double I_v = r_material_properties[I_V];
    const double I_n = r_material_properties[I_N];
    const double I_t = r_material_properties[I_T];

    if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        auto &r_stress_vector = rValues.GetStressVector();
        if (r_stress_vector.size() != strain_size)
            r_stress_vector.resize(strain_size, false);

        // Compute stresses according to 5x5 diagonal constitutive matrix with beam properties
        r_stress_vector[0] = E * A * r_strain_vector[0];    // N = E*A*ε (normal force)
        r_stress_vector[1] = E * I_n * r_strain_vector[1];  // M_y = E*I_y*κ_y (bending about y-axis)
        r_stress_vector[2] = E * I_n * r_strain_vector[2];  // M_z = E*I_z*κ_z (bending about z-axis)
        r_stress_vector[3] = G * I_t * r_strain_vector[3];  // T = G*I_t*γ_t1 (torsion)
        r_stress_vector[4] = G * I_t * r_strain_vector[4];  // T = G*I_t*γ_t2 (torsion)

        AddInitialStressVectorContribution(r_stress_vector);


        if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            auto &r_constitutive_matrix = rValues.GetConstitutiveMatrix();
            if (r_constitutive_matrix.size1() != strain_size || r_constitutive_matrix.size2() != strain_size)
                r_constitutive_matrix.resize(strain_size, strain_size, false);
            r_constitutive_matrix.clear();

            // 5x5 diagonal constitutive matrix with beam cross-section properties
            r_constitutive_matrix(0,0) = E * A;   // Normal stiffness
            r_constitutive_matrix(1,1) = E * I_v; // Bending stiffness about y-axis
            r_constitutive_matrix(2,2) = E * I_n; // Bending stiffness about z-axis
            r_constitutive_matrix(3,3) = 0.5 * G * I_t; // Torsional stiffness
            r_constitutive_matrix(4,4) = 0.5 * G * I_t; // Torsional stiffness
        }
    }
}

} // namespace Kratos
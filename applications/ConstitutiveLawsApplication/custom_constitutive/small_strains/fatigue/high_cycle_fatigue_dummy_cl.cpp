// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/small_strains/fatigue/high_cycle_fatigue_dummy_cl.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

HighCycleFatigueDummyCl::HighCycleFatigueDummyCl()
    : ElasticIsotropic3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

HighCycleFatigueDummyCl::HighCycleFatigueDummyCl(const HighCycleFatigueDummyCl& rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer HighCycleFatigueDummyCl::Clone() const
{
    return Kratos::make_shared<HighCycleFatigueDummyCl>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

HighCycleFatigueDummyCl::~HighCycleFatigueDummyCl()
{
}

/***********************************************************************************/
/***********************************************************************************/
void  HighCycleFatigueDummyCl::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    Flags &r_constitutive_law_options = rValues.GetOptions();
    ConstitutiveLaw::StrainVectorType &r_strain_vector = rValues.GetStrainVector();

    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        // Since we are in small strains, any strain measure works, e.g. CAUCHY_GREEN
        CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }
    AddInitialStrainVectorContribution<StrainVectorType>(r_strain_vector);

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        ConstitutiveLaw::StressVectorType &r_stress_vector = rValues.GetStressVector();
        CalculatePK2Stress(r_strain_vector, r_stress_vector, rValues);
        AddInitialStressVectorContribution<StressVectorType>(r_stress_vector);
    }

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        ConstitutiveLaw::VoigtSizeMatrixType &r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rValues);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/
void HighCycleFatigueDummyCl::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    Flags &r_constitutive_law_options = rValues.GetOptions();
    ConstitutiveLaw::StrainVectorType &r_strain_vector = rValues.GetStrainVector();
    ConstitutiveLaw::StressVectorType &r_stress_vector = rValues.GetStressVector();

    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        // Since we are in small strains, any strain measure works, e.g. CAUCHY_GREEN
        CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }
    AddInitialStrainVectorContribution<StrainVectorType>(r_strain_vector);

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        ConstitutiveLaw::StressVectorType &r_stress_vector = rValues.GetStressVector();
        CalculatePK2Stress(r_strain_vector, r_stress_vector, rValues);
        AddInitialStressVectorContribution<StressVectorType>(r_stress_vector);
    }

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        ConstitutiveLaw::VoigtSizeMatrixType &r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rValues);
    }

    double uniaxial_stress;
    uniaxial_stress = ConstitutiveLawUtilities<VoigtSize>::CalculateVonMisesEquivalentStress(r_stress_vector);

    double max_stress = mMaxStress;
    double min_stress = mMinStress;
    bool max_indicator = mMaxDetected;
    bool min_indicator = mMinDetected;

    double sign_factor = mFatigueData.CalculateTensionOrCompressionIdentifier(r_stress_vector);
    uniaxial_stress *= sign_factor;

    mFatigueData.CalculateSminAndSmax(uniaxial_stress,
                                    max_stress,
                                    min_stress,
                                    mPreviousStresses,
                                    max_indicator,
                                    min_indicator);

    mMaxStress = max_stress;
    mMinStress = min_stress;
    mMaxDetected = max_indicator;
    mMinDetected = min_indicator;

    Vector previous_stresses = ZeroVector(2);
    const Vector& r_aux_stresses = mPreviousStresses;
    previous_stresses[1] = uniaxial_stress;
    previous_stresses[0] = r_aux_stresses[1];
    mPreviousStresses = previous_stresses;

    KRATOS_WATCH(max_stress);
    KRATOS_WATCH(min_stress);
    KRATOS_WATCH(max_indicator);
    KRATOS_WATCH(min_indicator);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

// void HighCycleFatigueDummyCl::CalculatePK2Stress(
//     const ConstitutiveLaw::StrainVectorType& rStrainVector,
//     ConstitutiveLaw::StressVectorType& rStressVector,
//     ConstitutiveLaw::Parameters& rValues
//     )
// {
//     const Properties& r_material_properties = rValues.GetMaterialProperties();
//     const double E = r_material_properties[YOUNG_MODULUS];
//     const double NU = r_material_properties[POISSON_RATIO];

//     const double c1 = E / ((1.00 + NU) * (1 - 2 * NU));
//     const double c2 = c1 * (1 - NU);
//     const double c3 = c1 * NU;
//     const double c4 = c1 * 0.5 * (1 - 2 * NU);

//     rStressVector[0] = c2 * rStrainVector[0] + c3 * rStrainVector[1] + c3 * rStrainVector[2];
//     rStressVector[1] = c3 * rStrainVector[0] + c2 * rStrainVector[1] + c3 * rStrainVector[2];
//     rStressVector[2] = c3 * rStrainVector[0] + c3 * rStrainVector[1] + c2 * rStrainVector[2];
//     rStressVector[3] = c4 * rStrainVector[3];
//     rStressVector[4] = c4 * rStrainVector[4];
//     rStressVector[5] = c4 * rStrainVector[5];

//     double uniaxial_stress;
//     uniaxial_stress = ConstitutiveLawUtilities<VoigtSize>::CalculateVonMisesEquivalentStress(rStressVector);

//     double max_stress = mMaxStress;
//     double min_stress = mMinStress;
//     bool max_indicator = mMaxDetected;
//     bool min_indicator = mMinDetected;

//     double sign_factor = mFatigueData.CalculateTensionOrCompressionIdentifier(rStressVector);
//     uniaxial_stress *= sign_factor;

//     mFatigueData.CalculateSminAndSmax(uniaxial_stress,
//                                     max_stress,
//                                     min_stress,
//                                     mPreviousStresses,
//                                     max_indicator,
//                                     min_indicator);

//     mMaxStress = max_stress;
//     mMinStress = min_stress;
//     mMaxDetected = max_indicator;
//     mMinDetected = min_indicator;

//     Vector previous_stresses = ZeroVector(2);
//     const Vector& r_aux_stresses = mPreviousStresses;
//     previous_stresses[1] = uniaxial_stress;
//     previous_stresses[0] = r_aux_stresses[1];
//     mPreviousStresses = previous_stresses;

//     KRATOS_WATCH(max_stress);
//     KRATOS_WATCH(min_stress);
//     KRATOS_WATCH(max_indicator);
//     KRATOS_WATCH(min_indicator);
// }

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos

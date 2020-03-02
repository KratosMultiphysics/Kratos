// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_advanced_constitutive/wrinkling_linear_2d_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/
WrinklingLinear2DLaw::WrinklingLinear2DLaw()
    : ConstitutiveLaw()
{
}
/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/
WrinklingLinear2DLaw::WrinklingLinear2DLaw(const WrinklingLinear2DLaw& rOther)
    : ConstitutiveLaw(rOther),
      mpConstitutiveLaw(rOther.mpConstitutiveLaw)
{
}
/********************************CLONE**********************************************/
/***********************************************************************************/
ConstitutiveLaw::Pointer WrinklingLinear2DLaw::Clone() const
{
    return Kratos::make_shared<WrinklingLinear2DLaw>(*this);
}
/*******************************CONSTRUCTOR*****************************************/
/***********************************************************************************/
ConstitutiveLaw::Pointer WrinklingLinear2DLaw::Create(Kratos::Parameters NewParameters) const
{
    return Kratos::make_shared<WrinklingLinear2DLaw>();
}
std::size_t WrinklingLinear2DLaw::WorkingSpaceDimension()
{
    KRATOS_ERROR_IF(mpConstitutiveLaw->WorkingSpaceDimension()<2) << "WorkingSpaceDimension must be bigger than 1" << std::endl;
    return mpConstitutiveLaw->WorkingSpaceDimension();
}
std::size_t WrinklingLinear2DLaw::GetStrainSize()
{
    SizeType strain_size = 3;
    KRATOS_ERROR_IF_NOT(mpConstitutiveLaw->GetStrainSize()==3) << "Wrinkling law only works for 2D base laws (strain size = 3)" << std::endl;
    return strain_size;
}
void WrinklingLinear2DLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.NumberOfSubproperties()==1) << "Exactly one base claw must be given (" << rMaterialProperties.NumberOfSubproperties() << " claws are defined for baseclaw " << rMaterialProperties.Id() <<  ")" << std::endl;
    // We create the base constitutive laws
    const Properties& base_claw_prop = *(rMaterialProperties.GetSubProperties().begin());
    KRATOS_ERROR_IF_NOT(base_claw_prop.Has(CONSTITUTIVE_LAW)) << "No constitutive law set" << std::endl;
    mpConstitutiveLaw = base_claw_prop[CONSTITUTIVE_LAW]->Clone();
    mpConstitutiveLaw->InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
}
void  WrinklingLinear2DLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;
    // Do standard calculation
    Vector strain_vector = rValues.GetStrainVector();
    Vector stress_vector = ZeroVector(3);

    Matrix material_tangent_modulus = ZeroMatrix(3);
    ConstitutiveLaw::Parameters element_parameters;
    element_parameters.SetMaterialProperties(rValues.GetMaterialProperties());
    element_parameters.SetStrainVector(strain_vector);
    element_parameters.SetStressVector(stress_vector);
    element_parameters.SetConstitutiveMatrix(material_tangent_modulus);
    element_parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    element_parameters.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    element_parameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    mpConstitutiveLaw->CalculateMaterialResponsePK2(element_parameters);

    // Consider initial pre-stress state which is calculated on element level (i.e. projection)
    // This must be manually set in the element
    Vector pre_stress_vector = ZeroVector(3);
    if (rValues.IsSetStressVector()){
        pre_stress_vector = rValues.GetStressVector();
        stress_vector += pre_stress_vector;
    }

    // Check wrinkling state
    Vector wrinkling_direction = ZeroVector(2);
    WrinklingType current_wrinkling_state;
    CheckWrinklingState(current_wrinkling_state,stress_vector,strain_vector,wrinkling_direction);



    // Apply wrinkling algorithm
    // Slack
    if (current_wrinkling_state==WrinklingType::Slack){
        // We assign -pre_stress_vector because pre_stress_vector is added in the element
        // We want a zeroVector for the stress in case of slack
        stress_vector = -1.0 * pre_stress_vector;
        material_tangent_modulus = ZeroMatrix(3);
    }
    // Wrinkling
    else if (current_wrinkling_state==WrinklingType::Wrinkle){
        const double n_1 = wrinkling_direction[0];
        const double n_2 = wrinkling_direction[1];


        Vector wrinkling_operation_vector_1 = ZeroVector(3);
        wrinkling_operation_vector_1[0] = n_1*n_1;
        wrinkling_operation_vector_1[1] = n_2*n_2;
        wrinkling_operation_vector_1[2] = n_1*n_2*2.0;

        const Vector temp_vec = prod(material_tangent_modulus,wrinkling_operation_vector_1);
        const Matrix temp_mat = outer_prod(temp_vec,wrinkling_operation_vector_1);
        Matrix material_tangent_modulus_modified_1 = prod(temp_mat,material_tangent_modulus);
        const double temp_double = inner_prod(wrinkling_operation_vector_1,temp_vec);

        material_tangent_modulus_modified_1 /= temp_double;
        material_tangent_modulus_modified_1 = material_tangent_modulus - material_tangent_modulus_modified_1;

        noalias(stress_vector) = prod(material_tangent_modulus_modified_1,strain_vector);
        noalias(material_tangent_modulus) = material_tangent_modulus_modified_1;
    }
    else {
        // Else: taut, do nothing special
        // We substract the pre stress again to only obtain the material response
        // Pre-stress was only needed for wrinkling check
        // Pre-stress is added in the element
        stress_vector -= pre_stress_vector;
    }

    // Set the data on the main claw
    Flags& r_constitutive_law_options = rValues.GetOptions();
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS )) {
        Vector& r_stress_vector = rValues.GetStressVector();
        noalias(r_stress_vector) = stress_vector;
    }

    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR )) {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        noalias(r_constitutive_matrix) = material_tangent_modulus;
    }
    KRATOS_CATCH("");
}
void WrinklingLinear2DLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY;
    mpConstitutiveLaw->InitializeMaterialResponsePK2(rValues);
    KRATOS_CATCH("");
}
void WrinklingLinear2DLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY;
    mpConstitutiveLaw->FinalizeMaterialResponsePK2(rValues);
    KRATOS_CATCH("");
}
void WrinklingLinear2DLaw::ResetMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    mpConstitutiveLaw->ResetMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
}
void WrinklingLinear2DLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize = GetStrainSize();
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}
void WrinklingLinear2DLaw::PrincipalVector(Vector& rPrincipalVector, const Vector& rNonPrincipalVector)
{
    // Make sure to divide rNonPrincipalVector[2]/2 if strains are passed
    rPrincipalVector = ZeroVector(2);
    rPrincipalVector[0] = 0.50 * (rNonPrincipalVector[0]+rNonPrincipalVector[1]) + std::sqrt(0.25*(std::pow(rNonPrincipalVector[0]-rNonPrincipalVector[1],2.0)) + std::pow(rNonPrincipalVector[2],2.0));
    rPrincipalVector[1] = 0.50 * (rNonPrincipalVector[0]+rNonPrincipalVector[1]) - std::sqrt(0.25*(std::pow(rNonPrincipalVector[0]-rNonPrincipalVector[1],2.0)) + std::pow(rNonPrincipalVector[2],2.0));
}
void WrinklingLinear2DLaw::CheckWrinklingState(WrinklingType& rWrinklingState, const Vector& rStress, const Vector& rStrain, Vector& rWrinklingDirectionVector)
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();

    Vector principal_strains = ZeroVector(2);
    Vector temp_strains = ZeroVector(3);
    temp_strains = rStrain;
    temp_strains[2] /= 2.0; // Adjust voigt strain vector to calcualte principal strains
    PrincipalVector(principal_strains,temp_strains);

    Vector principal_stresses = ZeroVector(2);
    PrincipalVector(principal_stresses,rStress);

    const double min_stress = std::min(principal_stresses[0],principal_stresses[1]);
    const double max_stress = std::max(principal_stresses[0],principal_stresses[1]);
    const double max_strain = std::max(principal_strains[0],principal_strains[1]);

    rWrinklingDirectionVector = ZeroVector(2);

    // Direction check
    Vector min_stress_dir = ZeroVector(2);
    if (std::abs(rStress[2])>numerical_limit){
        min_stress_dir[0] = 1.0;
        min_stress_dir[1] = (min_stress-rStress[0]) / rStress[2];
        min_stress_dir /= MathUtils<double>::Norm(min_stress_dir);
    }
    else {
        const double stress_diff_1 = std::abs(min_stress-rStress[0]);
        const double stress_diff_2 = std::abs(min_stress-rStress[1]);
        if (stress_diff_1<=stress_diff_2) min_stress_dir[0] = 1.0;
        else min_stress_dir[1] = 1.0;
    }
    if ((min_stress > 0.0) || ((std::abs(min_stress)<numerical_limit) && (std::abs(max_stress)<numerical_limit))){
        //Second if-statement necessary for first iteration
        rWrinklingState = WrinklingType::Taut;
    } else if ((max_strain > 0.0) && (min_stress < numerical_limit)){
        rWrinklingState = WrinklingType::Wrinkle;
        rWrinklingDirectionVector[0] = min_stress_dir[0];
        rWrinklingDirectionVector[1] = min_stress_dir[1];
    } else if (max_strain<numerical_limit){
        rWrinklingState = WrinklingType::Slack;
    }
    else KRATOS_ERROR << "error in principal direction calculation 2" << std::endl;


}
int WrinklingLinear2DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR_IF_NOT(mpConstitutiveLaw) << "WrinklingLinear2DLaw is not initialized" << std::endl;
    mpConstitutiveLaw->Check(rMaterialProperties,rElementGeometry,rCurrentProcessInfo);
    return 0;
}

} // Namespace Kratos

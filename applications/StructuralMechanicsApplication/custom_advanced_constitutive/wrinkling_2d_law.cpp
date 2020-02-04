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
#include <set>

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_advanced_constitutive/wrinkling_2d_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

Wrinkling2DLaw::Wrinkling2DLaw()
    : ConstitutiveLaw()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

Wrinkling2DLaw::Wrinkling2DLaw(const Wrinkling2DLaw& rOther)
    : ConstitutiveLaw(rOther),
      mpConstitutiveLaw(rOther.mpConstitutiveLaw)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer Wrinkling2DLaw::Clone() const
{
    return Kratos::make_shared<Wrinkling2DLaw>(*this);
}

/*******************************CONSTRUCTOR*****************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer Wrinkling2DLaw::Create(Kratos::Parameters NewParameters) const
{
    return Kratos::make_shared<Wrinkling2DLaw>();
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

Wrinkling2DLaw::~Wrinkling2DLaw()
{
};



std::size_t Wrinkling2DLaw::WorkingSpaceDimension()
{
    KRATOS_ERROR_IF(mpConstitutiveLaw->WorkingSpaceDimension()<2) << "WorkingSpaceDimension must be bigger than 1" << std::endl;
    return mpConstitutiveLaw->WorkingSpaceDimension();
}


std::size_t Wrinkling2DLaw::GetStrainSize()
{
    SizeType strain_size = 3;
    KRATOS_ERROR_IF_NOT(mpConstitutiveLaw->GetStrainSize()==3) << "wrinkling law only works for 2D base laws (strain size = 3)" << std::endl;
    return strain_size;
}



template <class T>
bool Wrinkling2DLaw::THas(const Variable<T>& rTemplateVariable) const
{
    return mpConstitutiveLaw->Has(rTemplateVariable);
}


bool Wrinkling2DLaw::Has(const Variable<bool>& rThisVariable)
{
    return THas(rThisVariable);
}



bool Wrinkling2DLaw::Has(const Variable<int>& rThisVariable)
{
    return THas(rThisVariable);
}



bool Wrinkling2DLaw::Has(const Variable<double>& rThisVariable)
{
    return THas(rThisVariable);
}



bool Wrinkling2DLaw::Has(const Variable<Vector>& rThisVariable)
{
    return THas(rThisVariable);
}



bool Wrinkling2DLaw::Has(const Variable<Matrix>& rThisVariable)
{
    return THas(rThisVariable);
}



bool Wrinkling2DLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
{
    return THas(rThisVariable);
}



bool Wrinkling2DLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
{
    return THas(rThisVariable);
}



template <class T>
T& Wrinkling2DLaw::TGetValue(const Variable<T>& rTemplateVariable, T& rTemplateValue)
{
    mpConstitutiveLaw->GetValue(rTemplateVariable,rTemplateValue);
    return rTemplateValue;
}


bool& Wrinkling2DLaw::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    return TGetValue(rThisVariable,rValue);
}

int& Wrinkling2DLaw::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    return TGetValue(rThisVariable,rValue);
}

double& Wrinkling2DLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return TGetValue(rThisVariable,rValue);
}



Vector& Wrinkling2DLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return TGetValue(rThisVariable,rValue);
}



Matrix& Wrinkling2DLaw::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return TGetValue(rThisVariable,rValue);
}



array_1d<double, 3 >& Wrinkling2DLaw::GetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    return TGetValue(rThisVariable,rValue);
}



array_1d<double, 6 >& Wrinkling2DLaw::GetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    return TGetValue(rThisVariable,rValue);
}


void Wrinkling2DLaw::SetValue(
    const Variable<bool>& rThisVariable,
    const bool& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}


void Wrinkling2DLaw::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}


void Wrinkling2DLaw::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}


void Wrinkling2DLaw::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}


void Wrinkling2DLaw::SetValue(
    const Variable<Matrix>& rThisVariable,
    const Matrix& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}


void Wrinkling2DLaw::SetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    const array_1d<double, 3 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}


void Wrinkling2DLaw::SetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    const array_1d<double, 6 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}


bool& Wrinkling2DLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
    return rValue;
}


int& Wrinkling2DLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
    return rValue;
}


double& Wrinkling2DLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
    return rValue;
}


Vector& Wrinkling2DLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
    return rValue;
}


Matrix& Wrinkling2DLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
    return rValue;
}

array_1d<double, 3 >& Wrinkling2DLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
    return rValue;
}


array_1d<double, 6 >& Wrinkling2DLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
    return rValue;
}

void Wrinkling2DLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.NumberOfSubproperties()==1) << "Exactly one base claw must be given" << std::endl;
    // We create the base constitutive laws
    const Properties& base_claw_prop = *(rMaterialProperties.GetSubProperties().begin());
    KRATOS_ERROR_IF_NOT(base_claw_prop.Has(CONSTITUTIVE_LAW)) << "No constitutive law set" << std::endl;
    mpConstitutiveLaw = base_claw_prop[CONSTITUTIVE_LAW]->Clone();
    mpConstitutiveLaw->InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
}


void Wrinkling2DLaw::InitializeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const Properties& base_claw_prop = *(rMaterialProperties.GetSubProperties().begin());
    mpConstitutiveLaw->InitializeSolutionStep(base_claw_prop, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
}

void Wrinkling2DLaw::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const Properties& base_claw_prop = *(rMaterialProperties.GetSubProperties().begin());
    mpConstitutiveLaw->FinalizeSolutionStep(base_claw_prop, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
}



void Wrinkling2DLaw::InitializeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const Properties& base_claw_prop = *(rMaterialProperties.GetSubProperties().begin());
    mpConstitutiveLaw->InitializeNonLinearIteration(base_claw_prop, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
}



void Wrinkling2DLaw::FinalizeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const Properties& base_claw_prop = *(rMaterialProperties.GetSubProperties().begin());
    mpConstitutiveLaw->FinalizeNonLinearIteration(base_claw_prop, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);
}



void  Wrinkling2DLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;
    // do standard calculation
    Vector strain_vector = rValues.GetStrainVector();
    Vector stress_vector = ZeroVector(3);

    // consider initial pre-stress state which is calculated on element level (i.e. projection)
    // this must be manually set in the element
    Vector pre_stress_vector = ZeroVector(3);
    if (rValues.IsSetStressVector()){
        pre_stress_vector = rValues.GetStressVector();
        stress_vector += pre_stress_vector;
    }

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


    // check wrinkling state
    Vector wrinkling_direction = ZeroVector(2);
    WrinklingType current_wrinkling_state;
    CheckWrinklingState(current_wrinkling_state,stress_vector,strain_vector,wrinkling_direction);


    // apply wrinkling algorithm

    // slack
    if (current_wrinkling_state==WrinklingType::Slack){
        stress_vector = ZeroVector(3);
        material_tangent_modulus = ZeroMatrix(3);
    }
    // wrinkling
    else if (current_wrinkling_state==WrinklingType::Wrinkle){

        const double n_1 = wrinkling_direction[0];
        const double n_2 = wrinkling_direction[1];


        Vector wrinkling_operation_vector_1 = ZeroVector(3);
        wrinkling_operation_vector_1[0] = n_1*n_1;
        wrinkling_operation_vector_1[1] = n_2*n_2;
        wrinkling_operation_vector_1[2] = n_1*n_2*2.0;

        Vector temp_vec = prod(material_tangent_modulus,wrinkling_operation_vector_1);
        Matrix temp_mat = outer_prod(temp_vec,wrinkling_operation_vector_1);
        Matrix material_tangent_modulus_modified_1 = prod(temp_mat,material_tangent_modulus);
        double temp_double = inner_prod(wrinkling_operation_vector_1,temp_vec);

        material_tangent_modulus_modified_1 /= temp_double;
        material_tangent_modulus_modified_1 = material_tangent_modulus - material_tangent_modulus_modified_1;


        stress_vector = ZeroVector(3);
        ConstitutiveLaw::Parameters wrinkled_element_parameters;
        wrinkled_element_parameters.SetMaterialProperties(rValues.GetMaterialProperties());
        wrinkled_element_parameters.SetStrainVector(strain_vector);
        wrinkled_element_parameters.SetStressVector(stress_vector);
        wrinkled_element_parameters.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        wrinkled_element_parameters.SetConstitutiveMatrix(material_tangent_modulus_modified_1);
        mpConstitutiveLaw->CalculateMaterialResponse(wrinkled_element_parameters,ConstitutiveLaw::StressMeasure_PK2);

        stress_vector = prod(material_tangent_modulus_modified_1,strain_vector);
        material_tangent_modulus = material_tangent_modulus_modified_1;
    }
    else {
        // else: taut, do nothing special
        // we substract the pre stress again to only obtain the material reponse
        // pre stress was only needed for wrinkling check
        // pre stress is added in the element
        stress_vector -= pre_stress_vector;
    }

    // set the data on the main claw
    Flags& r_constitutive_law_options = rValues.GetOptions();
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS )) {
        rValues.SetStressVector(stress_vector);
    }

    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR )) {
        rValues.SetConstitutiveMatrix(material_tangent_modulus);
    }
    KRATOS_CATCH("");
}

void Wrinkling2DLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY;
    mpConstitutiveLaw->InitializeMaterialResponsePK2(rValues);
    KRATOS_CATCH("");
}




void Wrinkling2DLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY;
    mpConstitutiveLaw->FinalizeMaterialResponsePK2(rValues);
    KRATOS_CATCH("");
}




void Wrinkling2DLaw::ResetMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    const Properties& base_claw_prop = *(rMaterialProperties.GetSubProperties().begin());
    mpConstitutiveLaw->ResetMaterial(base_claw_prop, rElementGeometry, rShapeFunctionsValues);
}


void Wrinkling2DLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize = GetStrainSize();
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

void Wrinkling2DLaw::PrincipalVector(Vector& rPrincipalVector, const Vector& rNonPrincipalVector)
{
    // make sure to divide rNonPrincipalVector[2]/2 if strains are passed
    rPrincipalVector = ZeroVector(2);
    rPrincipalVector[0] = 0.50 * (rNonPrincipalVector[0]+rNonPrincipalVector[1]) + std::sqrt(0.25*(std::pow(rNonPrincipalVector[0]-rNonPrincipalVector[1],2.0)) + std::pow(rNonPrincipalVector[2],2.0));
    rPrincipalVector[1] = 0.50 * (rNonPrincipalVector[0]+rNonPrincipalVector[1]) - std::sqrt(0.25*(std::pow(rNonPrincipalVector[0]-rNonPrincipalVector[1],2.0)) + std::pow(rNonPrincipalVector[2],2.0));
}

void Wrinkling2DLaw::CheckWrinklingState(WrinklingType& rWrinklingState, const Vector& rStress, const Vector& rStrain, Vector& rWrinklingDirectionVector)
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();

    Vector principal_strains = ZeroVector(2);
    Vector temp_strains = ZeroVector(3);
    temp_strains = rStrain;
    temp_strains[2] /= 2.0; // adjust voigt strain vector to calcualte principal strains
    PrincipalVector(principal_strains,temp_strains);

    Vector principal_stresses = ZeroVector(2);
    PrincipalVector(principal_stresses,rStress);

    const double min_stress = std::min(principal_stresses[0],principal_stresses[1]);
    const double max_stress = std::max(principal_stresses[0],principal_stresses[1]);
    const double max_strain = std::max(principal_strains[0],principal_strains[1]);

    rWrinklingDirectionVector = ZeroVector(2);

    //direction check
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
        else if (stress_diff_2<stress_diff_1) min_stress_dir[1] = 1.0;
        else KRATOS_ERROR << "error in principal direction calculation 1"  << std::endl;
    }

    if ((min_stress > 0.0) || ((std::abs(min_stress)<numerical_limit) && (std::abs(max_stress)<numerical_limit))){
        //second if-statement necessary for first iteration
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


int Wrinkling2DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR_IF_NOT(mpConstitutiveLaw) << "claw is not initialized" << std::endl;
    const Properties& base_claw_prop = *(rMaterialProperties.GetSubProperties().begin());
    mpConstitutiveLaw->Check(base_claw_prop,rElementGeometry,rCurrentProcessInfo);
    return 0;
}

} // Namespace Kratos

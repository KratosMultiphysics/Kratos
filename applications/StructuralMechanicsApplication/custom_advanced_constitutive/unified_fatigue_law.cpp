// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:        BSD License
//                  license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//                   Sergio Jimenez
//
//
// System includes
#include <iostream>
#include <set>

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_advanced_constitutive/unified_fatigue_law.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"

// Yield surfaces
#include "custom_advanced_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/von_mises_yield_surface.h"

// Plastic potentials
#include "custom_advanced_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/von_mises_plastic_potential.h"

namespace Kratos
{

/********************************CLONE**********************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
ConstitutiveLaw::Pointer UnifiedFatigueLaw<TYieldSurfaceType>::Clone() const
{
    return Kratos::make_shared<UnifiedFatigueLaw>(*this);
}


/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::CalculateMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{
    // Auxiliar values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    const ProcessInfo& r_current_process_info = rValues.GetProcessInfo();
    const bool first_computation = (r_current_process_info[NL_ITERATION_NUMBER] == 1
                                   && r_current_process_info[STEP] == 1) ? true : false;

    if (first_computation) {
        if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ) ) {
            BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
        }
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            Vector& r_stress_vector = rValues.GetStressVector();
            if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                BaseType::CalculateElasticMatrix( r_constitutive_matrix, rValues);
                noalias(r_stress_vector) = prod( r_constitutive_matrix, r_strain_vector);
            } else {
                BaseType::CalculatePK2Stress( r_strain_vector, r_stress_vector, rValues);
            }
        }
    } else { // We check for plasticity
        // Integrate Stress plasticity
        Vector& r_integrated_stress_vector = rValues.GetStressVector();
        const double characteristic_length = 
            ConstitutiveLawUtilities<VoigtSize>::
            CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
        }

        // We compute the stress or the constitutive matrix
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {

            // We get some variables
            double threshold           = mThreshold;
            double plastic_dissipation = mPlasticDissipation;
            Vector plastic_strain      = mPlasticStrain;
            Matrix compliance_matrix   = mComplianceMatrix;

            // BoundedArrayType predictive_stress_vector;
            // // S0 = Elastic stress with strain (E-Ep)
            // Vector aux_stress = ZeroVector(VoigtSize);
            // BaseType::CalculatePK2Stress(r_strain_vector - plastic_strain, aux_stress, rValues);
            // noalias(predictive_stress_vector) = aux_stress;

            // Initialize Plastic Parameters
            double uniaxial_stress = 0.0, plastic_denominator = 0.0;
            BoundedArrayType plastic_flow = ZeroVector(VoigtSize); // DF/DS
            BoundedArrayType plastic_strain_increment = ZeroVector(VoigtSize);

            ////////////////////
            Matrix inverse_C(VoigtSize, VoigtSize);
            double aux_det_b = 0.0;
            MathUtils<double>::InvertMatrix(compliance_matrix, inverse_C, aux_det_b);
            r_integrated_stress_vector = prod(inverse_C, r_strain_vector - plastic_strain);
            noalias(r_constitutive_matrix) = inverse_C;


            TYieldSurfaceType::CalculateEquivalentStress(r_integrated_stress_vector, r_strain_vector,
                                                        uniaxial_stress, rValues);
            const double F = uniaxial_stress - mThreshold;
            if (F > machine_tolerance) {
                double J2;
                array_1d<double, VoigtSize> deviator = ZeroVector(6);
                const double I1 = r_integrated_stress_vector[0] + r_integrated_stress_vector[1] 
                                  + r_integrated_stress_vector[2];
                ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(r_integrated_stress_vector, I1, deviator, J2);
                TYieldSurfaceType::CalculateYieldSurfaceDerivative(r_integrated_stress_vector, deviator, J2, plastic_flow, rValues);

                Matrix compliance_increment(VoigtSize, VoigtSize);
                noalias(compliance_increment) = ZeroMatrix(VoigtSize, VoigtSize);
                Vector aux(VoigtSize);
                noalias(aux) = prod(inverse_C, plastic_flow);
                const double lambda_p = F / ((inner_prod(plastic_flow, aux)));

                noalias(compliance_increment) = lambda_p*(outer_prod(plastic_flow, plastic_flow)) 
                                                / (inner_prod(plastic_flow, r_integrated_stress_vector));
                
                compliance_matrix += compliance_increment;
                MathUtils<double>::InvertMatrix(compliance_matrix, inverse_C, aux_det_b);
                r_integrated_stress_vector = prod(inverse_C, r_strain_vector - plastic_strain);
                noalias(r_constitutive_matrix) = inverse_C;

                TYieldSurfaceType::CalculateEquivalentStress(r_integrated_stress_vector, r_strain_vector,
                                                            uniaxial_stress, rValues);
                KRATOS_WATCH(uniaxial_stress)
            }
            ////////////////////

            // Elastic Matrix
            // this->CalculateElasticMatrix(r_constitutive_matrix, rValues);


            // if (F <= std::abs(1.0e-4 * threshold)) { // Elastic case
            //     noalias(r_integrated_stress_vector) = predictive_stress_vector;
            // } else { // Plastic case


            //     if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            //         this->CalculateTangentTensor(rValues); // this modifies the ConstitutiveMatrix
            //     }
            // }



        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{

    // Auxiliar values
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    const ProcessInfo& r_current_process_info = rValues.GetProcessInfo();
    const bool first_computation = (r_current_process_info[NL_ITERATION_NUMBER] == 1
                                   && r_current_process_info[STEP] == 1) ? true : false;

    if (first_computation) {
        if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN ) ) {
            BaseType::CalculateCauchyGreenStrain( rValues, r_strain_vector);
        }
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            Vector& r_stress_vector = rValues.GetStressVector();
            if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                BaseType::CalculateElasticMatrix( r_constitutive_matrix, rValues);
                noalias(r_stress_vector) = prod( r_constitutive_matrix, r_strain_vector);
            } else {
                BaseType::CalculatePK2Stress( r_strain_vector, r_stress_vector, rValues);
            }
        }
    } else { // We check for plasticity
        // Integrate Stress plasticity
        Vector& r_integrated_stress_vector = rValues.GetStressVector();
        const double characteristic_length = 
            ConstitutiveLawUtilities<VoigtSize>::
            CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
        }

        // We compute the stress or the constitutive matrix
        if (r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS) ||
            r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {

            // We get some variables
            double threshold           = mThreshold;
            double plastic_dissipation = mPlasticDissipation;
            Vector plastic_strain      = mPlasticStrain;
            Matrix compliance_matrix   = mComplianceMatrix;

            // BoundedArrayType predictive_stress_vector;
            // // S0 = Elastic stress with strain (E-Ep)
            // Vector aux_stress = ZeroVector(VoigtSize);
            // BaseType::CalculatePK2Stress(r_strain_vector - plastic_strain, aux_stress, rValues);
            // noalias(predictive_stress_vector) = aux_stress;

            // Initialize Plastic Parameters
            double uniaxial_stress = 0.0, plastic_denominator = 0.0;
            BoundedArrayType plastic_flow = ZeroVector(VoigtSize); // DF/DS
            BoundedArrayType plastic_strain_increment = ZeroVector(VoigtSize);

            ////////////////////
            Matrix inverse_C(VoigtSize, VoigtSize);
            double aux_det_b = 0.0;
            MathUtils<double>::InvertMatrix(compliance_matrix, inverse_C, aux_det_b);
            r_integrated_stress_vector = prod(inverse_C, r_strain_vector - plastic_strain);
            noalias(r_constitutive_matrix) = inverse_C;


            TYieldSurfaceType::CalculateEquivalentStress(r_integrated_stress_vector, r_strain_vector,
                                                        uniaxial_stress, rValues);
            const double F = uniaxial_stress - mThreshold;
            if (F > machine_tolerance) {
                double J2;
                array_1d<double, VoigtSize> deviator = ZeroVector(6);
                const double I1 = r_integrated_stress_vector[0] + r_integrated_stress_vector[1] 
                                  + r_integrated_stress_vector[2];
                ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(r_integrated_stress_vector, I1, deviator, J2);
                TYieldSurfaceType::CalculateYieldSurfaceDerivative(r_integrated_stress_vector, deviator, J2, plastic_flow, rValues);

                Matrix compliance_increment(VoigtSize, VoigtSize);
                noalias(compliance_increment) = ZeroMatrix(VoigtSize, VoigtSize);
                Vector aux(VoigtSize);
                noalias(aux) = prod(inverse_C, plastic_flow);
                const double lambda_p = F / ((inner_prod(plastic_flow, aux)));

                noalias(compliance_increment) = lambda_p*(outer_prod(plastic_flow, plastic_flow)) 
                                                / (inner_prod(plastic_flow, r_integrated_stress_vector));
                
                compliance_matrix += compliance_increment;
                noalias(mComplianceMatrix) = compliance_matrix;
                MathUtils<double>::InvertMatrix(compliance_matrix, inverse_C, aux_det_b);
                r_integrated_stress_vector = prod(inverse_C, r_strain_vector - plastic_strain);
                noalias(r_constitutive_matrix) = inverse_C;
            }
            ////////////////////

            // Elastic Matrix
            // this->CalculateElasticMatrix(r_constitutive_matrix, rValues);


            // if (F <= std::abs(1.0e-4 * threshold)) { // Elastic case
            //     noalias(r_integrated_stress_vector) = predictive_stress_vector;
            // } else { // Plastic case


            //     if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            //         this->CalculateTangentTensor(rValues); // this modifies the ConstitutiveMatrix
            //     }
            // }



        }
    }





}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::CalculateElasticComplianceMatrix(
    Matrix& rC,
    const Properties& rMaterialProperties
    )
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double NU = rMaterialProperties[POISSON_RATIO];

    this->CheckClearElasticMatrix(rC);

    const double G = E / (2.0 * (1.0 + NU));
    const double c1 = 1.0 / E;
    const double c2 = -NU / E;
    const double c3 = 1.0 / G;

    rC(0,0) = c1; rC(0,1) = c2; rC(0,2) = c2;
    rC(1,0) = c2; rC(1,1) = c1; rC(1,2) = c2;
    rC(2,0) = c2; rC(2,1) = c2; rC(2,2) = c1;

    rC(3,3) = c3;
    rC(4,4) = c3;
    rC(5,5) = c3;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(
    const Variable<bool>& rThisVariable
    )
{
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(
    const Variable<int>& rThisVariable
    )
{
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(
    const Variable<double>& rThisVariable
    )
{
    // At least one layer should have the value
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(
    const Variable<Vector>& rThisVariable
    )
{
    // At least one layer should have the value
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(
    const Variable<Matrix>& rThisVariable
    )
{
    // At least one layer should have the value
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(
    const Variable<array_1d<double, 3 > >& rThisVariable
    )
{
    // At least one layer should have the value
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(
    const Variable<array_1d<double, 6 > >& rThisVariable
    )
{
    // At least one layer should have the value
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool& UnifiedFatigueLaw<TYieldSurfaceType>::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    // At least one layer should have the value
    rValue = false;
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
int& UnifiedFatigueLaw<TYieldSurfaceType>::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    // At least one layer should have the value
    rValue = 0;
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
double& UnifiedFatigueLaw<TYieldSurfaceType>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    // We combine the values of the layers
    rValue = 0.0;
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
Vector& UnifiedFatigueLaw<TYieldSurfaceType>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    // We combine the values of the layers
    rValue.clear();
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
Matrix& UnifiedFatigueLaw<TYieldSurfaceType>::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    // We combine the values of the layers
    rValue.clear();
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
array_1d<double, 3 >& UnifiedFatigueLaw<TYieldSurfaceType>::GetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    // We combine the values of the layers
    rValue = ZeroVector(3);
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
array_1d<double, 6 >& UnifiedFatigueLaw<TYieldSurfaceType>::GetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    // We combine the values of the layers
    rValue = ZeroVector(6);
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::SetValue(
    const Variable<bool>& rThisVariable,
    const bool& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::SetValue(
    const Variable<Matrix>& rThisVariable,
    const Matrix& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::SetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    const array_1d<double, 3 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::SetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    const array_1d<double, 6 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}


/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
double& UnifiedFatigueLaw<TYieldSurfaceType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
Vector& UnifiedFatigueLaw<TYieldSurfaceType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
Matrix& UnifiedFatigueLaw<TYieldSurfaceType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return rValue;
}
/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
array_1d<double, 3 >& UnifiedFatigueLaw<TYieldSurfaceType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
array_1d<double, 6 >& UnifiedFatigueLaw<TYieldSurfaceType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We construct the CL parameters
    ProcessInfo dummy_process_info;
    ConstitutiveLaw::Parameters aux_param(rElementGeometry, rMaterialProperties, 
                                          dummy_process_info);

    // We call the integrator
    double initial_threshold;
    TYieldSurfaceType::GetInitialUniaxialThreshold(aux_param, initial_threshold);
    mThreshold = initial_threshold;

    Matrix initial_compliance(VoigtSize, VoigtSize);
    CalculateElasticComplianceMatrix(initial_compliance, rMaterialProperties);
    noalias(mComplianceMatrix) = initial_compliance;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::CalculateMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void  UnifiedFatigueLaw<TYieldSurfaceType>::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::CalculateMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::InitializeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::InitializeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::InitializeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::InitializeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::FinalizeMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues
    )
{
    FinalizeMaterialResponseCauchy(rValues);
}

template<class TYieldSurfaceType>
/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueLaw<TYieldSurfaceType>::FinalizeMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::FinalizeMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues
    )
{
    FinalizeMaterialResponseCauchy(rValues);
}


/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
int UnifiedFatigueLaw<TYieldSurfaceType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // The auxiliar output
    int aux_out = 0;
    return aux_out;
}


/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::CalculateTangentTensor(
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        KRATOS_ERROR << "Analytic solution not available" << std::endl;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class UnifiedFatigueLaw<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>;
// template class UnifiedFatigueLaw<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>;

} // Namespace Kratos

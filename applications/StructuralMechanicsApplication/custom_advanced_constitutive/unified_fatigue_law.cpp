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
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
UnifiedFatigueLaw<TYieldSurfaceType>::UnifiedFatigueLaw()
    : ConstitutiveLaw()
{
}

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
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(const Variable<bool>& rThisVariable)
{
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(const Variable<int>& rThisVariable)
{
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(const Variable<double>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(const Variable<Vector>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(const Variable<Matrix>& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
{
    // At least one layer should have the value
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
bool UnifiedFatigueLaw<TYieldSurfaceType>::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
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
bool UnifiedFatigueLaw<TYieldSurfaceType>::ValidateInput(const Properties& rMaterialProperties)
{
    // We check it layer by layer
    bool valid_input = true;
    return valid_input;
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
    ConstitutiveLaw::Parameters aux_param(rElementGeometry, rMaterialProperties, dummy_process_info);

    // We call the integrator
    double initial_threshold;
    TConstLawIntegratorType::GetInitialUniaxialThreshold(aux_param, initial_threshold);
    this->SetThreshold(initial_threshold);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void  UnifiedFatigueLaw<TYieldSurfaceType>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
{




}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::InitializeMaterialResponsePK1(Parameters& rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::InitializeMaterialResponsePK2(Parameters& rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::InitializeMaterialResponseKirchhoff(Parameters& rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::InitializeMaterialResponseCauchy(Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

template<class TYieldSurfaceType>
/***********************************************************************************/
/***********************************************************************************/

void UnifiedFatigueLaw<TYieldSurfaceType>::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TYieldSurfaceType>
void UnifiedFatigueLaw<TYieldSurfaceType>::FinalizeMaterialResponseCauchy(Parameters& rValues)
{







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
void UnifiedFatigueLaw<TYieldSurfaceType>::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues, const ConstitutiveLaw::StressMeasure& rStressMeasure)
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
template class UnifiedFatigueLaw<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>;

} // Namespace Kratos

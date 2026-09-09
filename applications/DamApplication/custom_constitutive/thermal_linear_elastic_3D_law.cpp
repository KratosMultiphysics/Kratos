//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

/* Project includes */
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"

namespace Kratos
{

//Default Constructor
ThermalLinearElastic3DLaw::ThermalLinearElastic3DLaw() : BaseType() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalLinearElastic3DLaw::ThermalLinearElastic3DLaw(const ThermalLinearElastic3DLaw& rOther) : BaseType(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
ThermalLinearElastic3DLaw::~ThermalLinearElastic3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer ThermalLinearElastic3DLaw::Clone() const
{
    ThermalLinearElastic3DLaw::Pointer p_clone(new ThermalLinearElastic3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic3DLaw::CalculateDamThermalStrain(Vector& rThermalStrain, Parameters& rValues) const
{
    const double alpha =
        rValues.GetMaterialProperties()[THERMAL_EXPANSION];
    const double temperature =
        AdvancedConstitutiveLawUtilities<6>::CalculateInGaussPoint(TEMPERATURE, rValues);
    const double reference_temperature =
        AdvancedConstitutiveLawUtilities<6>::CalculateInGaussPoint(NODAL_REFERENCE_TEMPERATURE, rValues);

    rThermalStrain = ZeroVector(this->GetStrainSize());
    const double delta_temperature = temperature - reference_temperature;
    // 3D isotropic thermal strain: epsilon_th = alpha * (T - T_ref) * [1,1,1,0,0,0].
    for (SizeType i = 0; i < 3; ++i) {
        rThermalStrain[i] = alpha * delta_temperature;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic3DLaw::SubstractThermalStrain(
    ConstitutiveLaw::StrainVectorType& rStrainVector,
    const double ReferenceTemperature,
    ConstitutiveLaw::Parameters& rValues,
    const bool IsPlaneStrain)
{
    // Dam thermal strain: THERMAL_EXPANSION coefficient and shape-function
    // interpolated NODAL_REFERENCE_TEMPERATURE. ReferenceTemperature (the
    // generic scalar member) is intentionally not used.
    Vector thermal_strain(this->GetStrainSize());
    this->CalculateDamThermalStrain(thermal_strain, rValues);
    for (SizeType i = 0; i < thermal_strain.size(); ++i) {
        rStrainVector[i] -= thermal_strain[i];
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

int ThermalLinearElastic3DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
    // Standard elastic material checks (shared with the CLA elastic base).
    const double tolerance = 1.0e-12;
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] < 0.0)
        << "YOUNG_MODULUS is negative." << std::endl;
    KRATOS_ERROR_IF((0.5 - rMaterialProperties[POISSON_RATIO]) < tolerance)
        << "POISSON_RATIO is above the upper bound 0.5." << std::endl;
    KRATOS_ERROR_IF((rMaterialProperties[POISSON_RATIO] + 1.0) < tolerance)
        << "POISSON_RATIO is below the lower bound -1.0." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0)
        << "DENSITY is negative." << std::endl;
    return 0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Vector& ThermalLinearElastic3DLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    KRATOS_TRY

    if (rThisVariable == THERMAL_STRAIN_VECTOR ||
        rThisVariable == THERMAL_STRESS_VECTOR ||
        rThisVariable == MECHANICAL_STRESS_VECTOR) {

        // These outputs are computed from the current state carried by the
        // Parameters (element-provided total strain, shape functions, geometry),
        // so they reach the law through CalculateOnConstitutiveLaw. Has() is
        // intentionally left false so that the parameter-dependent CalculateValue
        // path is used instead of GetValue.
        const SizeType strain_size = this->GetStrainSize();

        // Constitutive matrix (reused from the inherited CLA law).
        ConstitutiveLaw::VoigtSizeMatrixType constitutive_matrix(strain_size, strain_size);
        noalias(constitutive_matrix) = ZeroMatrix(strain_size, strain_size);
        this->CalculateElasticMatrix(constitutive_matrix, rParameterValues);

        // Dam thermal strain (THERMAL_EXPANSION + NODAL_REFERENCE_TEMPERATURE).
        Vector thermal_strain_vector(strain_size);
        this->CalculateDamThermalStrain(thermal_strain_vector, rParameterValues);

        if (rThisVariable == MECHANICAL_STRESS_VECTOR) {
            // MECHANICAL_STRESS_VECTOR = C * epsilon
            const Vector& r_strain = rParameterValues.GetStrainVector();
            if (rValue.size() != strain_size)
                rValue.resize(strain_size, false);
            noalias(rValue) = prod(constitutive_matrix, r_strain);
            return rValue;
        }

        if (rThisVariable == THERMAL_STRAIN_VECTOR) {
            // THERMAL_STRAIN_VECTOR = epsilon_th
            if (rValue.size() != strain_size)
                rValue.resize(strain_size, false);
            noalias(rValue) = thermal_strain_vector;
            return rValue;
        }

        // THERMAL_STRESS_VECTOR = C * epsilon_th
        if (rValue.size() != strain_size)
            rValue.resize(strain_size, false);
        noalias(rValue) = prod(constitutive_matrix, thermal_strain_vector);
        return rValue;
    }

    // Not one of the specialized outputs: keep the base behaviour.
    return rValue;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Matrix& ThermalLinearElastic3DLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue)
{
    KRATOS_TRY

    const SizeType strain_size = this->GetStrainSize();
    const std::size_t dimension = (strain_size == 6) ? 3 : 2;

    if (rThisVariable == THERMAL_STRAIN_TENSOR) {
        Vector strain_vector = ZeroVector(strain_size);
        this->CalculateValue(rParameterValues, THERMAL_STRAIN_VECTOR, strain_vector);
        if (rValue.size1() != dimension || rValue.size2() != dimension)
            rValue.resize(dimension, dimension, false);
        noalias(rValue) = MathUtils<double>::StrainVectorToTensor(strain_vector);
        return rValue;
    }

    if (rThisVariable == THERMAL_STRESS_TENSOR) {
        Vector stress_vector = ZeroVector(strain_size);
        this->CalculateValue(rParameterValues, THERMAL_STRESS_VECTOR, stress_vector);
        if (rValue.size1() != dimension || rValue.size2() != dimension)
            rValue.resize(dimension, dimension, false);
        noalias(rValue) = MathUtils<double>::StressVectorToTensor(stress_vector);
        return rValue;
    }

    if (rThisVariable == MECHANICAL_STRESS_TENSOR) {
        Vector stress_vector = ZeroVector(strain_size);
        this->CalculateValue(rParameterValues, MECHANICAL_STRESS_VECTOR, stress_vector);
        if (rValue.size1() != dimension || rValue.size2() != dimension)
            rValue.resize(dimension, dimension, false);
        noalias(rValue) = MathUtils<double>::StressVectorToTensor(stress_vector);
        return rValue;
    }

    // Not one of the specialized outputs: keep the base behaviour.
    return rValue;

    KRATOS_CATCH( "" )
}

} // Namespace Kratos
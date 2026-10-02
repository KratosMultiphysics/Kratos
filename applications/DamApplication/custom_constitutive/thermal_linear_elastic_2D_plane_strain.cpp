
//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

/* Project includes */
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"

namespace Kratos
{

ThermalLinearElastic2DPlaneStrain::ThermalLinearElastic2DPlaneStrain() : BaseType() {}

ThermalLinearElastic2DPlaneStrain::ThermalLinearElastic2DPlaneStrain(const ThermalLinearElastic2DPlaneStrain& rOther) : BaseType(rOther) {}

ThermalLinearElastic2DPlaneStrain::~ThermalLinearElastic2DPlaneStrain() {}

ConstitutiveLaw::Pointer ThermalLinearElastic2DPlaneStrain::Clone() const
{
    ThermalLinearElastic2DPlaneStrain::Pointer p_clone(new ThermalLinearElastic2DPlaneStrain(*this));
    return p_clone;
}

void ThermalLinearElastic2DPlaneStrain::CalculateDamThermalStrain(Vector& rThermalStrain, Parameters& rValues) const
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double alpha = r_material_properties[THERMAL_EXPANSION];
    const double poisson = r_material_properties[POISSON_RATIO];
    const double temperature =
        AdvancedConstitutiveLawUtilities<3>::CalculateInGaussPoint(TEMPERATURE, rValues);
    const double reference_temperature =
        AdvancedConstitutiveLawUtilities<3>::CalculateInGaussPoint(NODAL_REFERENCE_TEMPERATURE, rValues);

    rThermalStrain = ZeroVector(this->GetStrainSize());
    const double delta_temperature = temperature - reference_temperature;
    // Plane-strain thermal strain: epsilon_th = alpha*(1+nu)*(T - T_ref)*[1,1,0].
    const double factor = alpha * (1.0 + poisson) * delta_temperature;
    rThermalStrain[0] = factor;
    rThermalStrain[1] = factor;
}

void ThermalLinearElastic2DPlaneStrain::SubstractThermalStrain(
    ConstitutiveLaw::StrainVectorType& rStrainVector,
    const double ReferenceTemperature,
    ConstitutiveLaw::Parameters& rValues,
    const bool IsPlaneStrain)
{
    Vector thermal_strain(this->GetStrainSize());
    this->CalculateDamThermalStrain(thermal_strain, rValues);
    for (SizeType i = 0; i < thermal_strain.size(); ++i) {
        rStrainVector[i] -= thermal_strain[i];
    }
}

int ThermalLinearElastic2DPlaneStrain::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
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

Vector& ThermalLinearElastic2DPlaneStrain::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    KRATOS_TRY

    if (rThisVariable == THERMAL_STRAIN_VECTOR ||
        rThisVariable == THERMAL_STRESS_VECTOR ||
        rThisVariable == MECHANICAL_STRESS_VECTOR) {
        const SizeType strain_size = this->GetStrainSize();
        ConstitutiveLaw::VoigtSizeMatrixType constitutive_matrix(strain_size, strain_size);
        noalias(constitutive_matrix) = ZeroMatrix(strain_size, strain_size);
        this->CalculateElasticMatrix(constitutive_matrix, rParameterValues);
        Vector thermal_strain_vector(strain_size);
        this->CalculateDamThermalStrain(thermal_strain_vector, rParameterValues);
        if (rThisVariable == MECHANICAL_STRESS_VECTOR) {
            const Vector& r_strain = rParameterValues.GetStrainVector();
            if (rValue.size() != strain_size) rValue.resize(strain_size, false);
            noalias(rValue) = prod(constitutive_matrix, r_strain);
            return rValue;
        }
        if (rThisVariable == THERMAL_STRAIN_VECTOR) {
            if (rValue.size() != strain_size) rValue.resize(strain_size, false);
            noalias(rValue) = thermal_strain_vector;
            return rValue;
        }
        if (rValue.size() != strain_size) rValue.resize(strain_size, false);
        noalias(rValue) = prod(constitutive_matrix, thermal_strain_vector);
        return rValue;
    }
    return rValue;
    KRATOS_CATCH( "" )
}

Matrix& ThermalLinearElastic2DPlaneStrain::CalculateValue(
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
        if (rValue.size1() != dimension || rValue.size2() != dimension) rValue.resize(dimension, dimension, false);
        noalias(rValue) = MathUtils<double>::StrainVectorToTensor(strain_vector);
        return rValue;
    }
    if (rThisVariable == THERMAL_STRESS_TENSOR) {
        Vector stress_vector = ZeroVector(strain_size);
        this->CalculateValue(rParameterValues, THERMAL_STRESS_VECTOR, stress_vector);
        if (rValue.size1() != dimension || rValue.size2() != dimension) rValue.resize(dimension, dimension, false);
        noalias(rValue) = MathUtils<double>::StressVectorToTensor(stress_vector);
        return rValue;
    }
    if (rThisVariable == MECHANICAL_STRESS_TENSOR) {
        Vector stress_vector = ZeroVector(strain_size);
        this->CalculateValue(rParameterValues, MECHANICAL_STRESS_VECTOR, stress_vector);
        if (rValue.size1() != dimension || rValue.size2() != dimension) rValue.resize(dimension, dimension, false);
        noalias(rValue) = MathUtils<double>::StressVectorToTensor(stress_vector);
        return rValue;
    }
    return rValue;
    KRATOS_CATCH( "" )
}

void ThermalLinearElastic2DPlaneStrain::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);
    rFeatures.mStrainSize = GetStrainSize();
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

} // Namespace Kratos

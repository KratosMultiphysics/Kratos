
//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

/* Project includes */
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress.hpp"
#include "includes/checks.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/thermal_output_utilities.hpp"

namespace Kratos
{

ThermalLinearElastic2DPlaneStress::ThermalLinearElastic2DPlaneStress() : BaseType() {}

ThermalLinearElastic2DPlaneStress::ThermalLinearElastic2DPlaneStress(const ThermalLinearElastic2DPlaneStress& rOther) : BaseType(rOther) {}

ThermalLinearElastic2DPlaneStress::~ThermalLinearElastic2DPlaneStress() {}

ConstitutiveLaw::Pointer ThermalLinearElastic2DPlaneStress::Clone() const
{
    ThermalLinearElastic2DPlaneStress::Pointer p_clone(new ThermalLinearElastic2DPlaneStress(*this));
    return p_clone;
}

void ThermalLinearElastic2DPlaneStress::CalculateDamThermalStrain(Vector& rThermalStrain, Parameters& rValues) const
{
    const double alpha = rValues.GetMaterialProperties()[THERMAL_EXPANSION];
    const double temperature =
        AdvancedConstitutiveLawUtilities<3>::CalculateInGaussPoint(TEMPERATURE, rValues);
    const double reference_temperature =
        AdvancedConstitutiveLawUtilities<3>::CalculateInGaussPoint(NODAL_REFERENCE_TEMPERATURE, rValues);

    rThermalStrain = ZeroVector(this->GetStrainSize());
    const double delta_temperature = temperature - reference_temperature;
    // Plane-stress thermal strain: epsilon_th = alpha*(T - T_ref)*[1,1,0].
    const double factor = alpha * delta_temperature;
    rThermalStrain[0] = factor;
    rThermalStrain[1] = factor;
}

void ThermalLinearElastic2DPlaneStress::SubstractThermalStrain(
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

int ThermalLinearElastic2DPlaneStress::Check(
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

Vector& ThermalLinearElastic2DPlaneStress::CalculateValue(
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
        // Undamaged thermoelastic assembly (damage factor 1.0), shared with the
        // thermal damage families.
        Vector thermal_strain_out, thermal_stress_out, mechanical_stress_out;
        ThermalOutputUtilities::AssembleOutputs(
            thermal_strain_out, thermal_stress_out, mechanical_stress_out,
            rParameterValues.GetStrainVector(), thermal_strain_vector,
            constitutive_matrix, 1.0);
        if (rThisVariable == THERMAL_STRAIN_VECTOR) {
            if (rValue.size() != strain_size) rValue.resize(strain_size, false);
            noalias(rValue) = thermal_strain_out;
            return rValue;
        }
        if (rThisVariable == THERMAL_STRESS_VECTOR) {
            if (rValue.size() != strain_size) rValue.resize(strain_size, false);
            noalias(rValue) = thermal_stress_out;
            return rValue;
        }
        // MECHANICAL_STRESS_VECTOR = C * epsilon
        if (rValue.size() != strain_size) rValue.resize(strain_size, false);
        noalias(rValue) = mechanical_stress_out;
        return rValue;
    }
    // Not one of the specialized outputs: delegate to the base-class behaviour.
    // Qualified through ElasticIsotropic3D, which owns these overloads: the CLA
    // thermal base only declares the Variable<double> overload.
    return ElasticIsotropic3D::CalculateValue(rParameterValues, rThisVariable, rValue);
    KRATOS_CATCH( "" )
}

Matrix& ThermalLinearElastic2DPlaneStress::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue)
{
    KRATOS_TRY
    const SizeType strain_size = this->GetStrainSize();
    if (rThisVariable == THERMAL_STRAIN_TENSOR) {
        Vector strain_vector = ZeroVector(strain_size);
        this->CalculateValue(rParameterValues, THERMAL_STRAIN_VECTOR, strain_vector);
        ThermalOutputUtilities::AssignStrainTensor(rValue, strain_vector);
        return rValue;
    }
    if (rThisVariable == THERMAL_STRESS_TENSOR) {
        Vector stress_vector = ZeroVector(strain_size);
        this->CalculateValue(rParameterValues, THERMAL_STRESS_VECTOR, stress_vector);
        ThermalOutputUtilities::AssignStressTensor(rValue, stress_vector);
        return rValue;
    }
    if (rThisVariable == MECHANICAL_STRESS_TENSOR) {
        Vector stress_vector = ZeroVector(strain_size);
        this->CalculateValue(rParameterValues, MECHANICAL_STRESS_VECTOR, stress_vector);
        ThermalOutputUtilities::AssignStressTensor(rValue, stress_vector);
        return rValue;
    }
    // Not one of the specialized outputs: delegate to the base-class behaviour.
    // Qualified through ElasticIsotropic3D, which owns these overloads: the CLA
    // thermal base only declares the Variable<double> overload.
    return ElasticIsotropic3D::CalculateValue(rParameterValues, rThisVariable, rValue);
    KRATOS_CATCH( "" )
}

void ThermalLinearElastic2DPlaneStress::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set( PLANE_STRESS_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);
    rFeatures.mStrainSize = GetStrainSize();
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

} // Namespace Kratos

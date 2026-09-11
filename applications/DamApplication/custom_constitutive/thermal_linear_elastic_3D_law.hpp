//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

#if !defined (KRATOS_THERMAL_LINEAR_ELASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_THERMAL_LINEAR_ELASTIC_3D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// ConstitutiveLawsApplication standard thermal elastic 3D law (base).
#include "custom_constitutive/thermal/small_strains/elastic/thermal_elastic_isotropic_3d.h"

#include "dam_application_variables.h"

namespace Kratos
{

/**
 * @brief Thin Dam compatibility adapter over ConstitutiveLawsApplication
 *        ThermalElasticIsotropic3D.
 *
 * The generic small-strain thermoelastic kernel (constitutive matrix, stress,
 * thermal strain, material response for every stress measure, stateless
 * lifecycle) is inherited from the ConstitutiveLawsApplication law. This class
 * only retains the Dam-specific behavior:
 *   - the spatial reference temperature read from the interpolated
 *     NODAL_REFERENCE_TEMPERATURE (shape-function evaluation) instead of the
 *     generic scalar REFERENCE_TEMPERATURE;
 *   - the historical Dam material coefficient THERMAL_EXPANSION;
 *   - the specialized Dam thermo-mechanical outputs resolved through the
 *     parameter-aware CalculateValue path (Has() == false).
 */
class KRATOS_API(DAM_APPLICATION) ThermalLinearElastic3DLaw : public ThermalElasticIsotropic3D
{

public:

    /// The ConstitutiveLawsApplication base law.
    using BaseType = ThermalElasticIsotropic3D;

    KRATOS_CLASS_POINTER_DEFINITION(ThermalLinearElastic3DLaw);

    // Bring base-class overloads of CalculateValue into scope to avoid hiding warnings.
    using BaseType::CalculateValue;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ThermalLinearElastic3DLaw();

    // Copy Constructor
    ThermalLinearElastic3DLaw (const ThermalLinearElastic3DLaw& rOther);

    // Destructor
    ~ThermalLinearElastic3DLaw() override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ConstitutiveLaw::Pointer Clone() const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief Dam thermal-strain subtraction for the inherited response kernel.
     * @details The generic ConstitutiveLawsApplication response applies the
     * thermal strain using the scalar REFERENCE_TEMPERATURE and
     * THERMAL_EXPANSION_COEFFICIENT. Dam uses instead the shape-function
     * interpolated NODAL_REFERENCE_TEMPERATURE and the historical
     * THERMAL_EXPANSION coefficient. Only this Dam-specific thermal strain is
     * overridden; the inherited C/stress/response machinery is reused.
     */
    void SubstractThermalStrain(
        ConstitutiveLaw::StrainVectorType& rStrainVector,
        const double ReferenceTemperature,
        ConstitutiveLaw::Parameters& rValues,
        const bool IsPlaneStrain = false) override;

    /**
     * @brief Performs the checks of the law with the Dam material contract
     * (TEMPERATURE available, THERMAL_EXPANSION set), in addition to the
     * standard elastic material checks.
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * Computes the specialized thermo-mechanical vector outputs from the current
     * state carried by the Parameters:
     *   THERMAL_STRAIN_VECTOR    = epsilon_th
     *   THERMAL_STRESS_VECTOR    = C * epsilon_th
     *   MECHANICAL_STRESS_VECTOR = C * epsilon
     * so that the total constitutive stress satisfies
     *   stress = MECHANICAL_STRESS_VECTOR - THERMAL_STRESS_VECTOR.
     * The output is read-only with respect to the constitutive state. The
     * constitutive matrix is reused from the inherited CLA law.
     */
    Vector& CalculateValue(Parameters& rParameterValues, const Variable<Vector>& rThisVariable, Vector& rValue) override;

    /**
     * Computes the specialized thermo-mechanical tensor outputs
     * (THERMAL_STRAIN_TENSOR, THERMAL_STRESS_TENSOR, MECHANICAL_STRESS_TENSOR)
     * as the tensor representations of the corresponding vector outputs.
     */
    Matrix& CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /**
     * @brief Dam thermal strain (3D): epsilon_th = alpha*(T - T_ref)*[1,1,1,0,0,0]
     * with alpha = THERMAL_EXPANSION and T_ref from the shape-function
     * interpolated NODAL_REFERENCE_TEMPERATURE.
     */
    void CalculateDamThermalStrain(Vector& rThermalStrain, Parameters& rValues) const;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
    }

};

} // namespace Kratos

#endif // KRATOS_THERMAL_LINEAR_ELASTIC_3D_LAW_H_INCLUDED defined
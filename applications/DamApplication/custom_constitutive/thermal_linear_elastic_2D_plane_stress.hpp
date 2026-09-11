//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

#if !defined (KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRESS_H_INCLUDED)
#define  KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRESS_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// ConstitutiveLawsApplication standard plane-stress thermal elastic law (base).
#include "custom_constitutive/thermal/small_strains/elastic/thermal_linear_plane_stress.h"

#include "dam_application_variables.h"

namespace Kratos
{

/**
 * @brief Thin Dam compatibility adapter over ConstitutiveLawsApplication
 *        ThermalLinearPlaneStress.
 * @details The generic plane-stress thermoelastic kernel (constitutive matrix,
 * stress, thermal strain, response, stateless lifecycle) is inherited. Only the
 * Dam-specific behavior is retained: NODAL_REFERENCE_TEMPERATURE interpolation,
 * the historical THERMAL_EXPANSION coefficient, the specialized Dam outputs and
 * the historical GetLawFeatures (PLANE_STRESS_LAW).
 */
class KRATOS_API(DAM_APPLICATION) ThermalLinearElastic2DPlaneStress : public ThermalLinearPlaneStress
{

public:

    /// The ConstitutiveLawsApplication base law.
    using BaseType = ThermalLinearPlaneStress;

    KRATOS_CLASS_POINTER_DEFINITION(ThermalLinearElastic2DPlaneStress);

    // Bring base-class overloads of CalculateValue into scope to avoid hiding warnings.
    using BaseType::CalculateValue;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ThermalLinearElastic2DPlaneStress();

    // Copy Constructor
    ThermalLinearElastic2DPlaneStress (const ThermalLinearElastic2DPlaneStress& rOther);

    // Destructor
    ~ThermalLinearElastic2DPlaneStress() override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ConstitutiveLaw::Pointer Clone() const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief Dam plane-stress thermal-strain subtraction for the inherited
     * response kernel (NODAL_REFERENCE_TEMPERATURE + THERMAL_EXPANSION).
     */
    void SubstractThermalStrain(
        ConstitutiveLaw::StrainVectorType& rStrainVector,
        const double ReferenceTemperature,
        ConstitutiveLaw::Parameters& rValues,
        const bool IsPlaneStrain = false) override;

    /**
     * @brief Performs the checks of the law with the Dam material contract.
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * Computes the specialized thermo-mechanical vector outputs
     * (THERMAL_STRAIN_VECTOR, THERMAL_STRESS_VECTOR, MECHANICAL_STRESS_VECTOR)
     * from the current state carried by the Parameters.
     */
    Vector& CalculateValue(Parameters& rParameterValues, const Variable<Vector>& rThisVariable, Vector& rValue) override;

    /**
     * Computes the specialized thermo-mechanical tensor outputs.
     */
    Matrix& CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue) override;

    /**
     * @brief Historical Dam law features (plane stress).
     */
    void GetLawFeatures(Features& rFeatures) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /**
     * @brief Dam plane-stress thermal strain:
     * epsilon_th = alpha*(T - T_ref)*[1,1,0] with alpha = THERMAL_EXPANSION and
     * T_ref from the shape-function interpolated NODAL_REFERENCE_TEMPERATURE.
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

#endif // KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRESS_H_INCLUDED defined
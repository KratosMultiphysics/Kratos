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
#include "custom_constitutive/continuum_laws/linear_elastic_3D_law.hpp"

#include "dam_application_variables.h"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) ThermalLinearElastic3DLaw : public LinearElastic3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ThermalLinearElastic3DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ThermalLinearElastic3DLaw();

    // Copy Constructor
    ThermalLinearElastic3DLaw (const ThermalLinearElastic3DLaw& rOther);

    // Destructor
    virtual ~ThermalLinearElastic3DLaw();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ConstitutiveLaw::Pointer Clone() const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (Parameters & rValues) override;

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues) override;

    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters & rValues) override;

    /**
     * This law has no evolving internal state, so the material-response
     * initialization and finalization callbacks are not required.
     */
    bool RequiresInitializeMaterialResponse() override;

    /**
     * This law has no evolving internal state, so the material-response
     * initialization and finalization callbacks are not required.
     */
    bool RequiresFinalizeMaterialResponse() override;

    /**
     * Computes the specialized thermo-mechanical vector outputs from the current
     * state carried by the Parameters:
     *   THERMAL_STRAIN_VECTOR    = epsilon_th
     *   THERMAL_STRESS_VECTOR    = C * epsilon_th
     *   MECHANICAL_STRESS_VECTOR = C * epsilon
     * so that the total constitutive stress satisfies
     *   stress = MECHANICAL_STRESS_VECTOR - THERMAL_STRESS_VECTOR.
     * The output is read-only with respect to the constitutive state.
     */
    Vector& CalculateValue(Parameters& rParameterValues, const Variable<Vector>& rThisVariable, Vector& rValue) override;

    /**
     * Computes the specialized thermo-mechanical tensor outputs
     * (THERMAL_STRAIN_TENSOR, THERMAL_STRESS_TENSOR, MECHANICAL_STRESS_TENSOR)
     * as the tensor representations of the corresponding vector outputs,
     * obtained by reusing CalculateValue(const Variable<Vector>&, ...) and the
     * standard MathUtils Voigt-to-tensor conversions. The output is read-only
     * with respect to the constitutive state.
     */
    Matrix& CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * Common infinitesimal thermo-elastic response shared by the PK2, Kirchhoff
     * and Cauchy stress measures:
     *     epsilon_th = alpha * (T - T_ref) * [1,1,1,0,0,0]
     *     stress     = C * (epsilon - epsilon_th)
     * @param rValues
     */
    void CalculateThermoElasticResponse(Parameters& rValues);

    double& CalculateDomainTemperature ( const MaterialResponseVariables & rElasticVariables, double & rTemperature) override;
    
    double& CalculateNodalReferenceTemperature ( const MaterialResponseVariables & rElasticVariables, double & rNodalReferenceTemperature);

    virtual void CalculateThermalStrain( Vector& rThermalStrainVector, const MaterialResponseVariables & rElasticVariables, double & rTemperature, double & rNodalReferenceTemperature);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LinearElastic3DLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LinearElastic3DLaw )
    }

}; // Class ThermalLinearElastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_THERMAL_LINEAR_ELASTIC_3D_LAW_H_INCLUDED  defined

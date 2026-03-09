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
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

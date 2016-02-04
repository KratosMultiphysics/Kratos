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
#include "custom_constitutive/linear_elastic_3D_law.hpp"

namespace Kratos
{

class ThermalLinearElastic3DLaw : public LinearElastic3DLaw
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
   
    ConstitutiveLaw::Pointer Clone() const;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues);
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double& CalculateDomainTemperature ( const MaterialResponseVariables & rElasticVariables, double & rTemperature);

    virtual void CalculateThermalStrain( Vector& rThermalStrainVector, const MaterialResponseVariables & rElasticVariables, double & rTemperature);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LinearElastic3DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LinearElastic3DLaw )
    }

}; // Class ThermalLinearElastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_THERMAL_LINEAR_ELASTIC_3D_LAW_H_INCLUDED  defined 

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
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp" 

namespace Kratos
{

class ThermalLinearElastic2DPlaneStress : public ThermalLinearElastic3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ThermalLinearElastic2DPlaneStress);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ThermalLinearElastic2DPlaneStress();

    // Copy Constructor
    ThermalLinearElastic2DPlaneStress (const ThermalLinearElastic2DPlaneStress& rOther);

    // Destructor
    virtual ~ThermalLinearElastic2DPlaneStress();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   
    ConstitutiveLaw::Pointer Clone() const;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures);

    void CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix, const double &rYoungModulus, const double &rPoissonCoefficient );

    void CalculateThermalStrain( Vector& rThermalStrainVector, const MaterialResponseVariables & rElasticVariables, double & rTemperature);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ThermalLinearElastic3DLaw)
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ThermalLinearElastic3DLaw)
    }

}; // Class ThermalLinearElastic2DPlaneStress
}  // namespace Kratos.
#endif // KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRESS_H_INCLUDED  defined 

//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_LOCAL_DAMAGE_PLANE_STRESS_2D_LAW_H_INCLUDED)
#define  KRATOS_LOCAL_DAMAGE_PLANE_STRESS_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/local_damage_plane_strain_2D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) LocalDamagePlaneStress2DLaw : public LocalDamagePlaneStrain2DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(LocalDamagePlaneStress2DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    LocalDamagePlaneStress2DLaw();
    
    /// Second Constructor
    LocalDamagePlaneStress2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 
    
    /// Copy Constructor
    LocalDamagePlaneStress2DLaw (const LocalDamagePlaneStress2DLaw& rOther);

    /// Destructor
    virtual ~LocalDamagePlaneStress2DLaw();

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void GetLawFeatures(Features& rFeatures);
    
    ConstitutiveLaw::Pointer Clone() const;
        
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
protected:

    /// Member Variables
        
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    void CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,const double& YoungModulus,const double& PoissonCoefficient );

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    /// Serialization
    
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class LocalDamagePlaneStress2DLaw
}  // namespace Kratos.
#endif // KRATOS_LOCAL_DAMAGE_PLANE_STRESS_2D_LAW_H_INCLUDED  defined 

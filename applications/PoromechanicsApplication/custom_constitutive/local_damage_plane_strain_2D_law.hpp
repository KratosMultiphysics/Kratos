//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_LOCAL_DAMAGE_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define  KRATOS_LOCAL_DAMAGE_PLANE_STRAIN_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/local_damage_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) LocalDamagePlaneStrain2DLaw : public LocalDamage3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(LocalDamagePlaneStrain2DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    LocalDamagePlaneStrain2DLaw();
    
    /// Second Constructor
    LocalDamagePlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 
    
    /// Copy Constructor
    LocalDamagePlaneStrain2DLaw (const LocalDamagePlaneStrain2DLaw& rOther);

    /// Destructor
    ~LocalDamagePlaneStrain2DLaw() override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void GetLawFeatures(Features& rFeatures) override;
    
    ConstitutiveLaw::Pointer Clone() const override;
        
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    SizeType WorkingSpaceDimension() override
    {
        return 2;
    }
    
    SizeType GetStrainSize() override
    {
        return 3;
    }
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
        
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,const double& YoungModulus,const double& PoissonCoefficient ) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    /// Serialization
    
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class LocalDamagePlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_LOCAL_DAMAGE_PLANE_STRAIN_2D_LAW_H_INCLUDED  defined 
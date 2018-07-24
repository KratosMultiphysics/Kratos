//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_NONLOCAL_DAMAGE_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define  KRATOS_NONLOCAL_DAMAGE_PLANE_STRAIN_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/nonlocal_damage_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) NonlocalDamagePlaneStrain2DLaw : public NonlocalDamage3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(NonlocalDamagePlaneStrain2DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    NonlocalDamagePlaneStrain2DLaw();
    
    /// Second Constructor
    NonlocalDamagePlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 
    
    /// Copy Constructor
    NonlocalDamagePlaneStrain2DLaw (const NonlocalDamagePlaneStrain2DLaw& rOther);

    /// Destructor
    ~NonlocalDamagePlaneStrain2DLaw() override;

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

}; // Class NonlocalDamagePlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_NONLOCAL_DAMAGE_PLANE_STRAIN_2D_LAW_H_INCLUDED  defined 
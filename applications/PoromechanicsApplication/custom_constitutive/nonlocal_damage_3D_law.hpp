//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED)
#define  KRATOS_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/local_damage_3D_law.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) NonlocalDamage3DLaw : public LocalDamage3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(NonlocalDamage3DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    NonlocalDamage3DLaw();
    
    /// Second Constructor
    NonlocalDamage3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 
    
    /// Copy Constructor
    NonlocalDamage3DLaw (const NonlocalDamage3DLaw& rOther);

    /// Destructor
    ~NonlocalDamage3DLaw() override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;
    
    ConstitutiveLaw::Pointer Clone() const override;

    void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues ) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateMaterialResponseCauchy (Parameters & rValues) override;
    
    void FinalizeMaterialResponseCauchy (Parameters & rValues) override;
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) override;

    void SetValue( const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo ) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    
    double mNonlocalEquivalentStrain;
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalReturnMapping( FlowRule::RadialReturnVariables& rReturnMappingVariables, Matrix& rStressMatrix, 
                                        Vector& rStressVector, const Matrix& LinearElasticMatrix, const Vector& StrainVector );

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

}; // Class NonlocalDamage3DLaw
}  // namespace Kratos.
#endif // KRATOS_NONLOCAL_DAMAGE_3D_LAW_H_INCLUDED  defined 

//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_RESTORE_SIMO_JU_NONLOCAL_3D_LAW_H_INCLUDED)
#define  KRATOS_RESTORE_SIMO_JU_NONLOCAL_3D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/restore_simo_ju_3D_law.hpp"
#include "custom_constitutive/custom_flow_rules/restore_nonlocal_damage_flow_rule.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) RestoreSimoJuNonlocal3DLaw : public RestoreSimoJu3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(RestoreSimoJuNonlocal3DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    RestoreSimoJuNonlocal3DLaw();
    
    /// Second Constructor
    RestoreSimoJuNonlocal3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 
    
    /// Copy Constructor
    RestoreSimoJuNonlocal3DLaw (const RestoreSimoJuNonlocal3DLaw& rOther);

    /// Destructor
    virtual ~RestoreSimoJuNonlocal3DLaw();

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);
    
    ConstitutiveLaw::Pointer Clone() const;

    void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues );

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateMaterialResponseCauchy (Parameters & rValues);
    
    void FinalizeMaterialResponseCauchy (Parameters & rValues);
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double& GetValue( const Variable<double>& rThisVariable, double& rValue );

    void SetValue( const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo );

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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class RestoreSimoJuNonlocal3DLaw
}  // namespace Kratos.
#endif // KRATOS_RESTORE_SIMO_JU_NONLOCAL_3D_LAW_H_INCLUDED  defined 

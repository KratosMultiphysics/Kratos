//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

#if !defined (KRATOS_EXPONENTIAL_COHESIVE_3D_LAW_H_INCLUDED)
#define  KRATOS_EXPONENTIAL_COHESIVE_3D_LAW_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/serializer.h"
// #include "includes/constitutive_law.h"
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) ExponentialCohesive3DLaw : public BilinearCohesive3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ExponentialCohesive3DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ExponentialCohesive3DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<ExponentialCohesive3DLaw>(ExponentialCohesive3DLaw(*this));
    }

    // Copy Constructor
    ExponentialCohesive3DLaw (const ExponentialCohesive3DLaw& rOther) : BilinearCohesive3DLaw(rOther)
    {
    }

    // Destructor
    ~ExponentialCohesive3DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;
    
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;
        
    void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMaterialResponseCauchy (Parameters & rValues) override;
    
    void FinalizeMaterialResponseCauchy (Parameters & rValues) override;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) override;
    
    void SetValue( const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo ) override;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables
    
    double mStateVariable;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void ComputeEquivalentStrain(double& rEquivalentStrain,const Vector& StrainVector,const double& CriticalDisplacement);
    
    virtual void ComputeEquivalentStrainContact(double& rEquivalentStrain,const Vector& StrainVector,const double& CriticalDisplacement);
    
    
    virtual void ComputeConstitutiveMatrixLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& JointStrength,
                                                        const double& DamageThreshold,const double& CriticalDisplacement);

    virtual void ComputeConstitutiveMatrixContactLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                            const double& JointStrength,const double& DamageThreshold,const double& CriticalDisplacement);
                                      
                                                            
    virtual void ComputeConstitutiveMatrixUnloading(Matrix& rConstitutiveMatrix,const double& JointStrength,
                                                        const double& DamageThreshold,const double& CriticalDisplacement);

    virtual void ComputeConstitutiveMatrixContactUnloading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                            const double& JointStrength,const double& DamageThreshold,const double& CriticalDisplacement);
                                                            
                                                            
    virtual void ComputeStressVector(Vector& rStressVector,const Vector& StrainVector,const double& JointStrength,
                                                const double& DamageThreshold,const double& CriticalDisplacement);
    
    virtual void ComputeStressVectorContact(Vector& rStressVector,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                        const double& JointStrength,const double& DamageThreshold,const double& CriticalDisplacement);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class ExponentialCohesive3DLaw
}  // namespace Kratos.
#endif // KRATOS_EXPONENTIAL_COHESIVE_3D_LAW_H_INCLUDED  defined 

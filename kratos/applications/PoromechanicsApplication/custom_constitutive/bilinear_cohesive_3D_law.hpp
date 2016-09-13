//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_BILINEAR_COHESIVE_3D_LAW_H_INCLUDED)
#define  KRATOS_BILINEAR_COHESIVE_3D_LAW_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/serializer.h"
#include "includes/constitutive_law.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) BilinearCohesive3DLaw : public ConstitutiveLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(BilinearCohesive3DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    BilinearCohesive3DLaw();

    // Copy Constructor
    BilinearCohesive3DLaw (const BilinearCohesive3DLaw& rOther);

    // Destructor
    virtual ~BilinearCohesive3DLaw();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures);
    
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);
    
    ConstitutiveLaw::Pointer Clone() const;
    
    void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMaterialResponseCauchy (Parameters & rValues);
    
    void FinalizeMaterialResponseCauchy (Parameters & rValues);
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double& GetValue( const Variable<double>& rThisVariable, double& rValue );
    
    void SetValue( const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables
    
    double mStateVariable;
    double mStateVariableEquilibrium;
    
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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class BilinearCohesive3DLaw
}  // namespace Kratos.
#endif // KRATOS_BILINEAR_COHESIVE_3D_LAW_H_INCLUDED  defined 

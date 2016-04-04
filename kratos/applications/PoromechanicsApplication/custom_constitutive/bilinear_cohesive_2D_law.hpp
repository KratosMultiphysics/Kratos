//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_BILINEAR_COHESIVE_2D_LAW_H_INCLUDED)
#define  KRATOS_BILINEAR_COHESIVE_2D_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) BilinearCohesive2DLaw : public BilinearCohesive3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(BilinearCohesive2DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    BilinearCohesive2DLaw();

    // Copy Constructor
    BilinearCohesive2DLaw (const BilinearCohesive2DLaw& rOther);

    // Destructor
    virtual ~BilinearCohesive2DLaw();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures);
        
    ConstitutiveLaw::Pointer Clone() const;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables
        
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ComputeEffectiveDisplacement(double& rEffectiveDisplacement,const Vector& StrainVector,const double& CriticalDisplacement);
    
    void ComputeEffectiveDisplacementContact(double& rEffectiveDisplacement,const Vector& StrainVector,const double& CriticalDisplacement);
    
    
    void ComputeConstitutiveMatrixLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& JointStrength,
                                                        const double& DamageThreshold,const double& CriticalDisplacement);

    void ComputeConstitutiveMatrixContactLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                            const double& JointStrength,const double& DamageThreshold,const double& CriticalDisplacement);
                                      
                                                            
    void ComputeConstitutiveMatrixUnloading(Matrix& rConstitutiveMatrix,const double& JointStrength,
                                                        const double& DamageThreshold,const double& CriticalDisplacement);

    void ComputeConstitutiveMatrixContactUnloading(Matrix& rConstitutiveMatrix,const double& YoungModulus,const double& FrictionCoefficient,
                                                            const double& JointStrength,const double& DamageThreshold,const double& CriticalDisplacement);
                                                            
                                                            
    void ComputeStressVector(Vector& rStressVector,const Vector& StrainVector,const double& JointStrength,
                                                const double& DamageThreshold,const double& CriticalDisplacement);
    
    void ComputeStressVectorContact(Vector& rStressVector,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                        const double& JointStrength,const double& DamageThreshold,const double& CriticalDisplacement);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BilinearCohesive3DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BilinearCohesive3DLaw )
    }

}; // Class BilinearCohesive2DLaw
}  // namespace Kratos.
#endif // KRATOS_BILINEAR_COHESIVE_2D_LAW_H_INCLUDED  defined 

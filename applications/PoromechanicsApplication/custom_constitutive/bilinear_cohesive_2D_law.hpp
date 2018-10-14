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
    ~BilinearCohesive2DLaw() override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;
        
    ConstitutiveLaw::Pointer Clone() const override;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables
        
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ComputeEquivalentStrain(double& rEquivalentStrain,const Vector& StrainVector,const double& CriticalDisplacement) override;
    
    void ComputeEquivalentStrainContact(double& rEquivalentStrain,const Vector& StrainVector,const double& CriticalDisplacement) override;
    
    
    void ComputeConstitutiveMatrixLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& JointStrength,
                                                        const double& DamageThreshold,const double& CriticalDisplacement) override;

    void ComputeConstitutiveMatrixContactLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                            const double& JointStrength,const double& DamageThreshold,const double& CriticalDisplacement) override;
                                      
                                                            
    void ComputeConstitutiveMatrixUnloading(Matrix& rConstitutiveMatrix,const double& JointStrength,
                                                        const double& DamageThreshold,const double& CriticalDisplacement) override;

    void ComputeConstitutiveMatrixContactUnloading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                            const double& JointStrength,const double& DamageThreshold,const double& CriticalDisplacement) override;
                                                            
                                                            
    void ComputeStressVector(Vector& rStressVector,const Vector& StrainVector,const double& JointStrength,
                                                const double& DamageThreshold,const double& CriticalDisplacement) override;
    
    void ComputeStressVectorContact(Vector& rStressVector,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                        const double& JointStrength,const double& DamageThreshold,const double& CriticalDisplacement) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BilinearCohesive3DLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BilinearCohesive3DLaw )
    }

}; // Class BilinearCohesive2DLaw
}  // namespace Kratos.
#endif // KRATOS_BILINEAR_COHESIVE_2D_LAW_H_INCLUDED  defined 

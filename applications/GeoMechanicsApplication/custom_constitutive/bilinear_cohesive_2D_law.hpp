// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) BilinearCohesive2DLaw : public BilinearCohesive3DLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BilinearCohesive2DLaw);

    BilinearCohesive2DLaw();

    void GetLawFeatures(Features& rFeatures) override;

    ConstitutiveLaw::Pointer Clone() const override;

    SizeType GetStrainSize() const override { return 2; }

    SizeType WorkingSpaceDimension() override { return 2; }

protected:
    void ComputeEquivalentStrain(double&       rEquivalentStrain,
                                 const Vector& StrainVector,
                                 const double& CriticalDisplacement) override;

    void ComputeEquivalentStrainContact(double&       rEquivalentStrain,
                                        const Vector& StrainVector,
                                        const double& CriticalDisplacement) override;

    void ComputeConstitutiveMatrixLoading(Matrix&       rConstitutiveMatrix,
                                          const Vector& StrainVector,
                                          const double& JointStrength,
                                          const double& DamageThreshold,
                                          const double& CriticalDisplacement) override;

    void ComputeConstitutiveMatrixContactLoading(Matrix&       rConstitutiveMatrix,
                                                 const Vector& StrainVector,
                                                 const double& YoungModulus,
                                                 const double& FrictionCoefficient,
                                                 const double& JointStrength,
                                                 const double& DamageThreshold,
                                                 const double& CriticalDisplacement) override;

    void ComputeConstitutiveMatrixUnloading(Matrix&       rConstitutiveMatrix,
                                            const double& JointStrength,
                                            const double& DamageThreshold,
                                            const double& CriticalDisplacement) override;

    void ComputeConstitutiveMatrixContactUnloading(Matrix&       rConstitutiveMatrix,
                                                   const Vector& StrainVector,
                                                   const double& YoungModulus,
                                                   const double& FrictionCoefficient,
                                                   const double& JointStrength,
                                                   const double& DamageThreshold,
                                                   const double& CriticalDisplacement) override;

    void ComputeStressVector(Vector&       rStressVector,
                             const Vector& StrainVector,
                             const double& JointStrength,
                             const double& DamageThreshold,
                             const double& CriticalDisplacement) override;

    void ComputeStressVectorContact(Vector&       rStressVector,
                                    const Vector& StrainVector,
                                    const double& YoungModulus,
                                    const double& FrictionCoefficient,
                                    const double& JointStrength,
                                    const double& DamageThreshold,
                                    const double& CriticalDisplacement) override;

private:
    // Serialization
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BilinearCohesive3DLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BilinearCohesive3DLaw)
    }

}; // Class BilinearCohesive2DLaw

} // namespace Kratos
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

// Application includes
#include "custom_constitutive/bilinear_cohesive_2D_law.hpp"
#include "geo_mechanics_application_constants.h"
#include "utilities/math_utils.h"

namespace Kratos
{

BilinearCohesive2DLaw::BilinearCohesive2DLaw() = default;

void BilinearCohesive2DLaw::GetLawFeatures(Features& rFeatures)
{
    // Set the type of law
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);

    // Set strain measure required by the constitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

    // Set the space dimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();

    // Set the strain size
    rFeatures.mStrainSize = GetStrainSize();
}

ConstitutiveLaw::Pointer BilinearCohesive2DLaw::Clone() const
{
    return Kratos::make_shared<BilinearCohesive2DLaw>(*this);
}

void BilinearCohesive2DLaw::ComputeEquivalentStrain(double&       rEquivalentStrain,
                                                    const Vector& StrainVector,
                                                    const double& CriticalDisplacement)
{
    KRATOS_TRY

    rEquivalentStrain = MathUtils<>::Norm(StrainVector) / CriticalDisplacement;

    KRATOS_CATCH("")
}

void BilinearCohesive2DLaw::ComputeEquivalentStrainContact(double&       rEquivalentStrain,
                                                           const Vector& StrainVector,
                                                           const double& CriticalDisplacement)
{
    KRATOS_TRY

    rEquivalentStrain = std::abs(StrainVector[0]) / CriticalDisplacement;

    KRATOS_CATCH("")
}

void BilinearCohesive2DLaw::ComputeConstitutiveMatrixLoading(Matrix&       rConstitutiveMatrix,
                                                             const Vector& StrainVector,
                                                             const double& YieldStress,
                                                             const double& DamageThreshold,
                                                             const double& CriticalDisplacement)
{
    KRATOS_TRY

    rConstitutiveMatrix(0, 0) =
        YieldStress / ((1.0 - DamageThreshold) * CriticalDisplacement) *
        ((1.0 - mStateVariable) / mStateVariable -
         StrainVector[0] * StrainVector[0] /
             (CriticalDisplacement * CriticalDisplacement * mStateVariable * mStateVariable * mStateVariable));
    rConstitutiveMatrix(1, 1) =
        YieldStress / ((1.0 - DamageThreshold) * CriticalDisplacement) *
        ((1.0 - mStateVariable) / mStateVariable -
         StrainVector[1] * StrainVector[1] /
             (CriticalDisplacement * CriticalDisplacement * mStateVariable * mStateVariable * mStateVariable));

    rConstitutiveMatrix(0, 1) = -YieldStress * StrainVector[0] * StrainVector[1] /
                                ((1.0 - DamageThreshold) * CriticalDisplacement * CriticalDisplacement *
                                 CriticalDisplacement * mStateVariable * mStateVariable * mStateVariable);
    rConstitutiveMatrix(1, 0) = rConstitutiveMatrix(0, 1);

    KRATOS_CATCH("")
}

void BilinearCohesive2DLaw::ComputeConstitutiveMatrixContactLoading(Matrix& rConstitutiveMatrix,
                                                                    const Vector& StrainVector,
                                                                    const double& YoungModulus,
                                                                    const double& FrictionCoefficient,
                                                                    const double& YieldStress,
                                                                    const double& DamageThreshold,
                                                                    const double& CriticalDisplacement)
{
    KRATOS_TRY

    rConstitutiveMatrix(0, 0) =
        YieldStress / ((1.0 - DamageThreshold) * CriticalDisplacement) *
        ((1.0 - mStateVariable) / mStateVariable -
         StrainVector[0] * StrainVector[0] /
             (CriticalDisplacement * CriticalDisplacement * mStateVariable * mStateVariable * mStateVariable));
    rConstitutiveMatrix(1, 1) = YoungModulus / (DamageThreshold * CriticalDisplacement);

    if (std::abs(StrainVector[0]) <= 1.e-20) {
        rConstitutiveMatrix(0, 1) = 0.0;
    } else {
        rConstitutiveMatrix(0, 1) =
            -YieldStress * StrainVector[0] * StrainVector[1] /
                ((1.0 - DamageThreshold) * CriticalDisplacement * CriticalDisplacement *
                 CriticalDisplacement * mStateVariable * mStateVariable * mStateVariable) -
            std::copysign(1.0, StrainVector[0]) * YoungModulus * FrictionCoefficient /
                (DamageThreshold * CriticalDisplacement);
    }

    rConstitutiveMatrix(1, 0) = 0.0;

    KRATOS_CATCH("")
}

void BilinearCohesive2DLaw::ComputeConstitutiveMatrixUnloading(Matrix&       rConstitutiveMatrix,
                                                               const double& YieldStress,
                                                               const double& DamageThreshold,
                                                               const double& CriticalDisplacement)
{
    KRATOS_TRY

    rConstitutiveMatrix(0, 0) = YieldStress / (CriticalDisplacement * mStateVariable) *
                                (1.0 - mStateVariable) / (1.0 - DamageThreshold);
    rConstitutiveMatrix(1, 1) = rConstitutiveMatrix(0, 0);

    rConstitutiveMatrix(0, 1) = 0.0;
    rConstitutiveMatrix(1, 0) = 0.0;

    KRATOS_CATCH("")
}

void BilinearCohesive2DLaw::ComputeConstitutiveMatrixContactUnloading(Matrix& rConstitutiveMatrix,
                                                                      const Vector& StrainVector,
                                                                      const double& YoungModulus,
                                                                      const double& FrictionCoefficient,
                                                                      const double& YieldStress,
                                                                      const double& DamageThreshold,
                                                                      const double& CriticalDisplacement)
{
    KRATOS_TRY

    rConstitutiveMatrix(0, 0) = YieldStress / (CriticalDisplacement * mStateVariable) *
                                (1.0 - mStateVariable) / (1.0 - DamageThreshold);
    rConstitutiveMatrix(1, 1) = YoungModulus / (DamageThreshold * CriticalDisplacement);

    if (std::abs(StrainVector[0]) <= 1.e-20) {
        rConstitutiveMatrix(0, 1) = 0.0;
    } else {
        rConstitutiveMatrix(0, 1) = -std::copysign(1.0, StrainVector[0]) * YoungModulus *
                                    FrictionCoefficient / (DamageThreshold * CriticalDisplacement);
    }

    rConstitutiveMatrix(1, 0) = 0.0;

    KRATOS_CATCH("")
}

void BilinearCohesive2DLaw::ComputeStressVector(Vector&       rStressVector,
                                                const Vector& StrainVector,
                                                const double& YieldStress,
                                                const double& DamageThreshold,
                                                const double& CriticalDisplacement)
{
    KRATOS_TRY

    rStressVector[0] = YieldStress / (CriticalDisplacement * mStateVariable) *
                       (1.0 - mStateVariable) / (1.0 - DamageThreshold) * StrainVector[0];
    rStressVector[1] = YieldStress / (CriticalDisplacement * mStateVariable) *
                       (1.0 - mStateVariable) / (1.0 - DamageThreshold) * StrainVector[1];

    KRATOS_CATCH("")
}

void BilinearCohesive2DLaw::ComputeStressVectorContact(Vector&       rStressVector,
                                                       const Vector& StrainVector,
                                                       const double& YoungModulus,
                                                       const double& FrictionCoefficient,
                                                       const double& YieldStress,
                                                       const double& DamageThreshold,
                                                       const double& CriticalDisplacement)
{
    KRATOS_TRY

    // Note: StrainVector[1] < 0.0
    rStressVector[indexStress2DInterface::INDEX_2D_INTERFACE_ZZ] =
        YoungModulus / (DamageThreshold * CriticalDisplacement) *
        StrainVector[indexStress2DInterface::INDEX_2D_INTERFACE_ZZ];

    if (std::abs(StrainVector[indexStress2DInterface::INDEX_2D_INTERFACE_XZ]) <= 1.e-20) {
        rStressVector[indexStress2DInterface::INDEX_2D_INTERFACE_XZ] = 0.0;
    } else {
        rStressVector[indexStress2DInterface::INDEX_2D_INTERFACE_XZ] =
            YieldStress / (CriticalDisplacement * mStateVariable) * (1.0 - mStateVariable) /
                (1.0 - DamageThreshold) * StrainVector[indexStress2DInterface::INDEX_2D_INTERFACE_XZ] -
            std::copysign(1.0, StrainVector[indexStress2DInterface::INDEX_2D_INTERFACE_XZ]) *
                FrictionCoefficient * rStressVector[indexStress2DInterface::INDEX_2D_INTERFACE_ZZ];
    }

    KRATOS_CATCH("")
}

} // Namespace Kratos
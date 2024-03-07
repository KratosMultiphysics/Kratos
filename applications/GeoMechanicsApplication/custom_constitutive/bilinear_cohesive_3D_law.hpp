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

// System includes
#include <cmath>

// Project includes
#include "includes/constitutive_law.h"
#include "includes/serializer.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) BilinearCohesive3DLaw : public ConstitutiveLaw
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(BilinearCohesive3DLaw);

    BilinearCohesive3DLaw();

    void GetLawFeatures(Features& rFeatures) override;

    int Check(const Properties&   rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo&  rCurrentProcessInfo) const override;

    ConstitutiveLaw::Pointer Clone() const override;

    void InitializeMaterial(const Properties&   rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector&       rShapeFunctionsValues) override;

    SizeType GetStrainSize() const override { return 3; }

    SizeType WorkingSpaceDimension() override { return 3; }

    void CalculateMaterialResponseCauchy(Parameters& rValues) override;

    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;

    using ConstitutiveLaw::GetValue;
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

    using ConstitutiveLaw::SetValue;
    void SetValue(const Variable<double>& rVariable, const double& rValue, const ProcessInfo& rCurrentProcessInfo) override;

protected:
    // Member Variables
    double mStateVariable;

    virtual void ComputeEquivalentStrain(double&       rEquivalentStrain,
                                         const Vector& StrainVector,
                                         const double& CriticalDisplacement);

    virtual void ComputeEquivalentStrainContact(double&       rEquivalentStrain,
                                                const Vector& StrainVector,
                                                const double& CriticalDisplacement);

    virtual void ComputeConstitutiveMatrixLoading(Matrix&       rConstitutiveMatrix,
                                                  const Vector& StrainVector,
                                                  const double& JointStrength,
                                                  const double& DamageThreshold,
                                                  const double& CriticalDisplacement);

    virtual void ComputeConstitutiveMatrixContactLoading(Matrix&       rConstitutiveMatrix,
                                                         const Vector& StrainVector,
                                                         const double& YoungModulus,
                                                         const double& FrictionCoefficient,
                                                         const double& JointStrength,
                                                         const double& DamageThreshold,
                                                         const double& CriticalDisplacement);

    virtual void ComputeConstitutiveMatrixUnloading(Matrix&       rConstitutiveMatrix,
                                                    const double& JointStrength,
                                                    const double& DamageThreshold,
                                                    const double& CriticalDisplacement);

    virtual void ComputeConstitutiveMatrixContactUnloading(Matrix&       rConstitutiveMatrix,
                                                           const Vector& StrainVector,
                                                           const double& YoungModulus,
                                                           const double& FrictionCoefficient,
                                                           const double& JointStrength,
                                                           const double& DamageThreshold,
                                                           const double& CriticalDisplacement);

    virtual void ComputeStressVector(Vector&       rStressVector,
                                     const Vector& StrainVector,
                                     const double& JointStrength,
                                     const double& DamageThreshold,
                                     const double& CriticalDisplacement);

    virtual void ComputeStressVectorContact(Vector&       rStressVector,
                                            const Vector& StrainVector,
                                            const double& YoungModulus,
                                            const double& FrictionCoefficient,
                                            const double& JointStrength,
                                            const double& DamageThreshold,
                                            const double& CriticalDisplacement);

private:
    // Serialization
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

}; // Class BilinearCohesive3DLaw

} // namespace Kratos
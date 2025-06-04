// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

#include "custom_constitutive/small_strain_udsm_3D_law.hpp"
#include "includes/define.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUDSM3DInterfaceLaw : public SmallStrainUDSMLaw
{
public:
    // The base class ConstitutiveLaw type definition
    using BaseType = ConstitutiveLaw;

    // The size type definition
    using SizeType = std::size_t;

    // Pointer definition of SmallStrainUDSM3DInterfaceLaw
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainUDSM3DInterfaceLaw);

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue) override;
    using SmallStrainUDSMLaw::GetValue;

    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    using SmallStrainUDSMLaw::SetValue;

    SizeType WorkingSpaceDimension() override;

    [[nodiscard]] SizeType GetStrainSize() const override;

    [[nodiscard]] std::string Info() const override;
    void                      PrintInfo(std::ostream& rOStream) const override;
    void                      PrintData(std::ostream& rOStream) const override;

protected:
    void UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters& rValues) override;
    void SetExternalStressVector(Vector& rStressVector) override;
    void SetInternalStressVector(const Vector& rStressVector) override;
    void SetInternalStrainVector(const Vector& rStrainVector) override;
    void CopyConstitutiveMatrix(ConstitutiveLaw::Parameters& rValues, Matrix& rConstitutiveMatrix) override;

private:
    [[nodiscard]] indexStress3D getIndex3D(indexStress3DInterface index3D) const;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
}; // Class SmallStrainUDSM3DLaw

} // namespace Kratos
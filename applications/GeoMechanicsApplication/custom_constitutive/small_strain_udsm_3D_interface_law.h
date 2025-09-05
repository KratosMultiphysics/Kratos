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

#include "custom_constitutive/small_strain_udsm_law.h"
#include "includes/define.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUDSM3DInterfaceLaw : public SmallStrainUDSMLaw
{
public:
    using SizeType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainUDSM3DInterfaceLaw);

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue) override;
    using SmallStrainUDSMLaw::GetValue;

    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;
    using SmallStrainUDSMLaw::SetValue;

    SizeType WorkingSpaceDimension() override;

    [[nodiscard]] SizeType GetStrainSize() const override;

    [[nodiscard]] std::string Info() const override;
    void                      PrintData(std::ostream& rOStream) const override;

protected:
    void UpdateInternalDeltaStrainVector(Parameters& rValues) override;
    void SetExternalStressVector(Vector& rStressVector) override;
    void SetInternalStressVector(const Vector& rStressVector) override;
    void SetInternalStrainVector(const Vector& rStrainVector) override;
    void CopyConstitutiveMatrix(Parameters& rValues, Matrix& rConstitutiveMatrix) override;

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos
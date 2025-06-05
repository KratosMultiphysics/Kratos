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

#include "includes/define.h"
#include "small_strain_udsm_law.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUDSM2DInterfaceLaw : public SmallStrainUDSMLaw
{
public:
    // The base class ConstitutiveLaw type definition
    using BaseType = ConstitutiveLaw;

    // The size type definition
    using SizeType = std::size_t;

    // Pointer definition of SmallStrainUDSM2DInterfaceLaw
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainUDSM2DInterfaceLaw);

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    using SmallStrainUDSMLaw::GetValue;
    Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue) override;

    using SmallStrainUDSMLaw::SetValue;
    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override;

    SizeType WorkingSpaceDimension() override;

    [[nodiscard]] SizeType GetStrainSize() const override;

    [[nodiscard]] std::string Info() const override;
    void                      PrintInfo(std::ostream& rOStream) const override;
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
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

class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUDSM2DPlaneStrainLaw : public SmallStrainUDSM3DLaw
{
public:
    using BaseType = ConstitutiveLaw;
    using SizeType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainUDSM2DPlaneStrainLaw);

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;

    SizeType WorkingSpaceDimension() override;

    [[nodiscard]] SizeType GetStrainSize() const override;

    [[nodiscard]] std::string Info() const override;

    void PrintInfo(std::ostream& rOStream) const override;
    void PrintData(std::ostream& rOStream) const override;

private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
}; // Class SmallStrainUDSM3DLaw

} // namespace Kratos
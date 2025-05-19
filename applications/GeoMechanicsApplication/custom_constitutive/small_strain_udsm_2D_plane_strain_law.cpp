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

// External includes
#include "custom_constitutive/small_strain_udsm_2D_plane_strain_law.hpp"

namespace Kratos
{

ConstitutiveLaw::Pointer SmallStrainUDSM2DPlaneStrainLaw::Clone() const
{
    KRATOS_TRY

    return Kratos::make_shared<SmallStrainUDSM2DPlaneStrainLaw>(*this);

    KRATOS_CATCH("")
}

void SmallStrainUDSM2DPlaneStrainLaw::SetExternalStressVector(Vector& rStressVector)
{
    KRATOS_TRY
    for (unsigned int i = 0; i < GetStrainSize(); ++i) {
        rStressVector(i) = mStressVector[i];
    }
    KRATOS_CATCH("")
}

void SmallStrainUDSM2DPlaneStrainLaw::SetValue(const Variable<Vector>& rVariable,
                                               const Vector&           rValue,
                                               const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == STATE_VARIABLES) {
        SmallStrainUDSM3DLaw::SetValue(rVariable, rValue, rCurrentProcessInfo);
    } else if ((rVariable == CAUCHY_STRESS_VECTOR) && (rValue.size() == GetStrainSize())) {
        this->SetInternalStressVector(rValue);
    }
}

SizeType SmallStrainUDSM2DPlaneStrainLaw::GetStrainSize() const
{
    return VOIGT_SIZE_2D_PLANE_STRAIN;
}

} // namespace Kratos
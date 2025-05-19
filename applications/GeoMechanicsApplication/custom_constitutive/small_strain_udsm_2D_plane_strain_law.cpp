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
    return std::make_shared<SmallStrainUDSM2DPlaneStrainLaw>(*this);
}

SizeType SmallStrainUDSM2DPlaneStrainLaw::GetStrainSize() const
{
    return VOIGT_SIZE_2D_PLANE_STRAIN;
}

} // namespace Kratos
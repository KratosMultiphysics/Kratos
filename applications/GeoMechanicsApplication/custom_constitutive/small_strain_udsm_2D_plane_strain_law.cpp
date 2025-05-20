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

#include "custom_constitutive/small_strain_udsm_2D_plane_strain_law.hpp"

namespace Kratos
{

ConstitutiveLaw::Pointer SmallStrainUDSM2DPlaneStrainLaw::Clone() const
{
    return std::make_shared<SmallStrainUDSM2DPlaneStrainLaw>(*this);
}

SizeType SmallStrainUDSM2DPlaneStrainLaw::WorkingSpaceDimension() { return N_DIM_2D; }

SizeType SmallStrainUDSM2DPlaneStrainLaw::GetStrainSize() const
{
    return VOIGT_SIZE_2D_PLANE_STRAIN;
}

std::string SmallStrainUDSM2DPlaneStrainLaw::Info() const
{
    return "SmallStrainUDSM2DPlaneStrainLaw";
}

void SmallStrainUDSM2DPlaneStrainLaw::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

void SmallStrainUDSM2DPlaneStrainLaw::PrintData(std::ostream& rOStream) const
{
    rOStream << "SmallStrainUDSM2DPlaneStrainLaw Data";
}

void SmallStrainUDSM2DPlaneStrainLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
}

void SmallStrainUDSM2DPlaneStrainLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
}

} // namespace Kratos
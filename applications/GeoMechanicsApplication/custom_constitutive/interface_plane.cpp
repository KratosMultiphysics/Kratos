// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "interface_plane.h"

#include "geo_mechanics_application_constants.h"
#include "includes/constitutive_law.h"
#include "includes/exception.h"

namespace Kratos
{

Matrix InterfacePlane::CalculateElasticMatrix(double YoungsModulus, double PoissonsRatio) const
{
    KRATOS_ERROR << "InterfacePlane::CalculateElasticMatrix is not yet implemented";
    return ZeroMatrix(4, 4);
}

std::unique_ptr<ConstitutiveLawDimension> InterfacePlane::Clone() const
{
    return std::make_unique<InterfacePlane>();
}

std::size_t InterfacePlane::GetStrainSize() const { return VOIGT_SIZE_2D_INTERFACE; }

std::size_t InterfacePlane::GetDimension() const { return N_DIM_2D; }

std::size_t InterfacePlane::GetNumberOfNormalComponents() const
{
    KRATOS_ERROR << "InterfacePlane::GetNumberOfNormalComponents is not yet implemented";
    return 3;
}

Flags InterfacePlane::GetSpatialType() const
{
    KRATOS_ERROR << "InterfacePlane::GetSpatialType is not yet implemented";
    return ConstitutiveLaw::PLANE_STRAIN_LAW;
}

void InterfacePlane::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void InterfacePlane::load(Serializer&)
{
    // No data members to be saved (yet)
}

} // namespace Kratos

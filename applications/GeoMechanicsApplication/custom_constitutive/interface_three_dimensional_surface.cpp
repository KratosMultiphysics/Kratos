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

#include "interface_three_dimensional_surface.h"

#include "geo_mechanics_application_constants.h"
#include "includes/constitutive_law.h"
#include "includes/exception.h"

namespace Kratos
{

Matrix InterfaceThreeDimensionalSurface::CalculateElasticMatrix(double YoungsModulus, double PoissonsRatio) const
{
    KRATOS_ERROR << "InterfaceThreeDimensionalSurface::CalculateElasticMatrix is not yet implemented";
    return ZeroMatrix(4, 4);// NOSONAR: required to satisfy return type
}

std::unique_ptr<ConstitutiveLawDimension> InterfaceThreeDimensionalSurface::Clone() const
{
    return std::make_unique<InterfaceThreeDimensionalSurface>();
}

std::size_t InterfaceThreeDimensionalSurface::GetStrainSize() const { return VOIGT_SIZE_3D_INTERFACE; }

std::size_t InterfaceThreeDimensionalSurface::GetDimension() const { return N_DIM_3D; }

std::size_t InterfaceThreeDimensionalSurface::GetNumberOfNormalComponents() const
{
    KRATOS_ERROR << "InterfaceThreeDimensionalSurface::GetNumberOfNormalComponents is not yet implemented";
    return 3;// NOSONAR: required to satisfy return type
}

Flags InterfaceThreeDimensionalSurface::GetSpatialType() const
{
    KRATOS_ERROR << "InterfaceThreeDimensionalSurface::GetSpatialType is not yet implemented";
    return ConstitutiveLaw::PLANE_STRAIN_LAW;// NOSONAR: required to satisfy return type
}

void InterfaceThreeDimensionalSurface::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void InterfaceThreeDimensionalSurface::load(Serializer&)
{
    // No data members to be saved (yet)
}

} // namespace Kratos

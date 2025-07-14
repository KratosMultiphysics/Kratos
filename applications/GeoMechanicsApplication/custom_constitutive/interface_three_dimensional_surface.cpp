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
//                   Gennady Markelov
//

#include "interface_three_dimensional_surface.h"

#include "custom_utilities/constitutive_law_utilities.h"
#include "geo_mechanics_application_constants.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

Matrix InterfaceThreeDimensionalSurface::CalculateElasticMatrix(const Properties& rProperties) const
{
    auto result = ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(
        rProperties[INTERFACE_NORMAL_STIFFNESS], rProperties[INTERFACE_SHEAR_STIFFNESS],
        GetStrainSize(), GetNumberOfNormalComponents());

    return result;
}

std::unique_ptr<ConstitutiveLawDimension> InterfaceThreeDimensionalSurface::Clone() const
{
    return std::make_unique<InterfaceThreeDimensionalSurface>();
}

std::size_t InterfaceThreeDimensionalSurface::GetStrainSize() const
{
    return VOIGT_SIZE_3D_INTERFACE;
}

std::size_t InterfaceThreeDimensionalSurface::GetDimension() const { return N_DIM_3D; }

std::size_t InterfaceThreeDimensionalSurface::GetNumberOfNormalComponents() const { return 1; }

Flags InterfaceThreeDimensionalSurface::GetSpatialType() const
{
    return ConstitutiveLaw::THREE_DIMENSIONAL_LAW;
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

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

namespace Kratos
{

Matrix InterfaceThreeDimensionalSurface::MakeInterfaceConstitutiveMatrix(double NormalStiffness,
                                                                         double ShearStiffness,
                                                                         std::size_t TractionSize) const
{
    auto result = ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(
        NormalStiffness, ShearStiffness, TractionSize);
    result(2, 2) = result(1, 1);

    return result;
}

std::unique_ptr<InterfaceConstitutiveLawDimension> InterfaceThreeDimensionalSurface::Clone() const
{
    return std::make_unique<InterfaceThreeDimensionalSurface>();
}

std::size_t InterfaceThreeDimensionalSurface::GetStrainSize() const
{
    return VOIGT_SIZE_3D_INTERFACE;
}

std::size_t InterfaceThreeDimensionalSurface::GetDimension() const { return N_DIM_3D; }

void InterfaceThreeDimensionalSurface::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void InterfaceThreeDimensionalSurface::load(Serializer&)
{
    // No data members to be saved (yet)
}

} // namespace Kratos

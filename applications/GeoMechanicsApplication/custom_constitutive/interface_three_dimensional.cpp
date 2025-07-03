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

#include "interface_three_dimensional.h"

#include "geo_mechanics_application_constants.h"
#include "includes/constitutive_law.h"
#include "includes/exception.h"

namespace Kratos
{

Matrix InterfaceThreeDimensional::CalculateElasticMatrix(double YoungsModulus, double PoissonsRatio) const
{
    KRATOS_ERROR << "InterfaceThreeDimensional::CalculateElasticMatrix is not yet implemented";
    return ZeroMatrix(4, 4);
}

std::unique_ptr<ConstitutiveLawDimension> InterfaceThreeDimensional::Clone() const
{
    return std::make_unique<InterfaceThreeDimensional>();
}

std::size_t InterfaceThreeDimensional::GetStrainSize() const { return VOIGT_SIZE_3D_INTERFACE; }

std::size_t InterfaceThreeDimensional::GetDimension() const { return N_DIM_3D; }

std::size_t InterfaceThreeDimensional::GetNumberOfNormalComponents() const
{
    KRATOS_ERROR << "InterfaceThreeDimensional::GetNumberOfNormalComponents is not yet implemented";
    return 3;
}

Flags InterfaceThreeDimensional::GetSpatialType() const
{
    KRATOS_ERROR << "InterfaceThreeDimensional::GetSpatialType is not yet implemented";
    return ConstitutiveLaw::PLANE_STRAIN_LAW;
}

void InterfaceThreeDimensional::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void InterfaceThreeDimensional::load(Serializer&)
{
    // No data members to be saved (yet)
}

} // namespace Kratos

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

#include "interface_linear_strain.h"

#include "geo_mechanics_application_constants.h"
#include "includes/constitutive_law.h"
#include "includes/exception.h"

namespace Kratos
{

Matrix InterfaceLinearStrain::CalculateElasticMatrix(double YoungsModulus, double PoissonsRatio) const
{
    KRATOS_ERROR << "not yet implemented";
    return ZeroMatrix(4, 4);
}

std::unique_ptr<ConstitutiveLawDimension> InterfaceLinearStrain::Clone() const
{
    return std::make_unique<InterfaceLinearStrain>();
}

std::size_t InterfaceLinearStrain::GetStrainSize() const { return VOIGT_SIZE_2D_INTERFACE; }

std::size_t InterfaceLinearStrain::GetDimension() const { return N_DIM_2D; }

std::size_t InterfaceLinearStrain::GetNumberOfNormalComponents() const
{
    KRATOS_ERROR << "not yet implemented";
    return 3;
}

Flags InterfaceLinearStrain::GetSpatialType() const
{
    KRATOS_ERROR << "not yet implemented";
    return ConstitutiveLaw::PLANE_STRAIN_LAW;
}

void InterfaceLinearStrain::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void InterfaceLinearStrain::load(Serializer&)
{
    // No data members to be saved (yet)
}

} // namespace Kratos

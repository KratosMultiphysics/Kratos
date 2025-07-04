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

#include "interface_plane_strain.h"

#include "geo_mechanics_application_constants.h"
#include "includes/constitutive_law.h"
#include "includes/exception.h"

namespace Kratos
{

Matrix InterfacePlaneStrain::CalculateElasticMatrix(double YoungsModulus, double PoissonsRatio) const
{
    KRATOS_ERROR << "InterfacePlaneStrain::CalculateElasticMatrix is not yet implemented";
    return ZeroMatrix(4, 4);// NOSONAR: required to satisfy return type
}

std::unique_ptr<ConstitutiveLawDimension> InterfacePlaneStrain::Clone() const
{
    return std::make_unique<InterfacePlaneStrain>();
}

std::size_t InterfacePlaneStrain::GetStrainSize() const { return VOIGT_SIZE_2D_INTERFACE; }

std::size_t InterfacePlaneStrain::GetDimension() const { return N_DIM_2D; }

std::size_t InterfacePlaneStrain::GetNumberOfNormalComponents() const
{
    KRATOS_ERROR << "InterfacePlaneStrain::GetNumberOfNormalComponents is not yet implemented";
    return 3;// NOSONAR: required to satisfy return type
}

Flags InterfacePlaneStrain::GetSpatialType() const
{
    KRATOS_ERROR << "InterfacePlaneStrain::GetSpatialType is not yet implemented";
    return ConstitutiveLaw::PLANE_STRAIN_LAW;// NOSONAR: required to satisfy return type
}

void InterfacePlaneStrain::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void InterfacePlaneStrain::load(Serializer&)
{
    // No data members to be saved (yet)
}

} // namespace Kratos

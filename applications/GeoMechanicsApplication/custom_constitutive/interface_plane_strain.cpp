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

#include "interface_plane_strain.h"

#include "custom_utilities/constitutive_law_utilities.h"
#include "geo_mechanics_application_constants.h"

namespace Kratos
{

Matrix InterfacePlaneStrain::CalculateElasticMatrix(double NormalStiffness, double ShearStiffness) const
{
    return ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(
        NormalStiffness, ShearStiffness, GetStrainSize());
}

std::unique_ptr<ConstitutiveLawDimension> InterfacePlaneStrain::Clone() const
{
    return std::make_unique<InterfacePlaneStrain>();
}

std::size_t InterfacePlaneStrain::GetStrainSize() const { return VOIGT_SIZE_2D_INTERFACE; }

std::size_t InterfacePlaneStrain::GetDimension() const { return N_DIM_2D; }

std::size_t InterfacePlaneStrain::GetNumberOfNormalComponents() const { return 1; }

Flags InterfacePlaneStrain::GetSpatialType() const { return ConstitutiveLaw::PLANE_STRAIN_LAW; }

void InterfacePlaneStrain::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void InterfacePlaneStrain::load(Serializer&)
{
    // No data members to be saved (yet)
}

} // namespace Kratos

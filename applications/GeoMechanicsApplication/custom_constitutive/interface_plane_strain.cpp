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

#include "custom_utilities/constitutive_law_utilities.h"
#include "geo_mechanics_application_constants.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

Matrix InterfacePlaneStrain::MakeInterfaceConstitutiveMatrix(double      NormalStiffness,
                                                             double      ShearStiffness,
                                                             std::size_t TractionSize) const
{
    return ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(NormalStiffness, ShearStiffness, TractionSize);
}

std::unique_ptr<InterfaceConstitutiveLawDimension> InterfacePlaneStrain::Clone() const
{
    return std::make_unique<InterfacePlaneStrain>();
}

std::size_t InterfacePlaneStrain::GetStrainSize() const { return VOIGT_SIZE_2D_INTERFACE; }

std::size_t InterfacePlaneStrain::GetDimension() const { return N_DIM_2D; }

void InterfacePlaneStrain::save(Serializer&) const
{
    // No data members to be saved (yet)
}

void InterfacePlaneStrain::load(Serializer&)
{
    // No data members to be saved (yet)
}

} // namespace Kratos

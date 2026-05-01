// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "custom_elements/axisymmetric_stress_state.h"
#include "custom_elements/interface_stress_state.h"
#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/three_dimensional_stress_state.h"

#include "custom_constitutive/interface_plane_strain.h"
#include "custom_constitutive/interface_three_dimensional_surface.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_constitutive/thermal_filter_law.h"
#include "custom_constitutive/three_dimensional.h"
#include "custom_retention/saturated_below_phreatic_level_law.h"
#include "custom_retention/saturated_law.h"
#include "custom_retention/van_genuchten_law.h"

#include "custom_utilities/registration_utilities.hpp"
#include "includes/serializer.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{
template <typename TConcrete, typename TBase>
void CheckSerializable(const std::string& rName)
{
    const auto scoped_registration = ScopedSerializerRegistration{std::make_pair(rName, TConcrete{})};

    const auto p_obj      = std::unique_ptr<TBase>{std::make_unique<TConcrete>()};
    auto       serializer = StreamSerializer{};

    serializer.save("test_tag", p_obj);
    auto p_loaded = std::unique_ptr<TBase>{};
    serializer.load("test_tag", p_loaded);

    KRATOS_EXPECT_NE(p_loaded, nullptr);
    KRATOS_EXPECT_NE(dynamic_cast<TConcrete*>(p_loaded.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(GeoMechanics_CheckSerializationRegistration, KratosGeoMechanicsFastSuite)
{
    // Stress state policies
    CheckSerializable<PlaneStrainStressState, StressStatePolicy>("PlaneStrainStressState");
    CheckSerializable<ThreeDimensionalStressState, StressStatePolicy>(
        "ThreeDimensionalStressState");
    CheckSerializable<AxisymmetricStressState, StressStatePolicy>("AxisymmetricStressState");
    CheckSerializable<Line2DInterfaceStressState, StressStatePolicy>("Line2DInterfaceStressState");
    CheckSerializable<SurfaceInterfaceStressState, StressStatePolicy>(
        "SurfaceInterfaceStressState");

    // Constitutive law dimensions already covered elsewhere but include smoke checks
    CheckSerializable<PlaneStrain, ConstitutiveLawDimension>("PlaneStrain");
    CheckSerializable<ThreeDimensional, ConstitutiveLawDimension>("ThreeDimensional");

    // Retention laws
    CheckSerializable<SaturatedBelowPhreaticLevelLaw, RetentionLaw>(
        "SaturatedBelowPhreaticLevelLaw");
    CheckSerializable<SaturatedLaw, RetentionLaw>("SaturatedLaw");
    CheckSerializable<VanGenuchtenLaw, RetentionLaw>("VanGenuchtenLaw");

    // Constitutive law dimensions
    CheckSerializable<InterfacePlaneStrain, ConstitutiveLawDimension>("InterfacePlaneStrain");
    CheckSerializable<InterfaceThreeDimensionalSurface, ConstitutiveLawDimension>(
        "InterfaceThreeDimensionalSurface");

    // Thermal law
    CheckSerializable<GeoThermalFilterLaw, GeoThermalLaw>("GeoThermalFilterLaw");
}

} // namespace Kratos::Testing

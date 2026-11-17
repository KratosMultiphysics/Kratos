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

#include "custom_constitutive/coulomb_yield_surface.h"
#include "custom_constitutive/interface_plane_strain.h"
#include "custom_constitutive/interface_three_dimensional_surface.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_constitutive/tension_cutoff.h"
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
void CheckSerializable()
{
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
    CheckSerializable<PlaneStrainStressState, StressStatePolicy>();
    CheckSerializable<ThreeDimensionalStressState, StressStatePolicy>();
    CheckSerializable<AxisymmetricStressState, StressStatePolicy>();
    CheckSerializable<Line2DInterfaceStressState, StressStatePolicy>();
    CheckSerializable<SurfaceInterfaceStressState, StressStatePolicy>();

    // Constitutive law dimensions already covered elsewhere but include smoke checks
    CheckSerializable<PlaneStrain, ConstitutiveLawDimension>();
    CheckSerializable<ThreeDimensional, ConstitutiveLawDimension>();

    // Retention laws
    CheckSerializable<SaturatedBelowPhreaticLevelLaw, RetentionLaw>();
    CheckSerializable<SaturatedLaw, RetentionLaw>();
    CheckSerializable<VanGenuchtenLaw, RetentionLaw>();

    // Constitutive law dimensions
    CheckSerializable<InterfacePlaneStrain, ConstitutiveLawDimension>();
    CheckSerializable<InterfaceThreeDimensionalSurface, ConstitutiveLawDimension>();

    // Thermal law
    CheckSerializable<GeoThermalFilterLaw, GeoThermalLaw>();

    // Yield surfaces
    CheckSerializable<CoulombYieldSurface, CoulombYieldSurface>();
    CheckSerializable<TensionCutoff, TensionCutoff>();
}

} // namespace Kratos::Testing

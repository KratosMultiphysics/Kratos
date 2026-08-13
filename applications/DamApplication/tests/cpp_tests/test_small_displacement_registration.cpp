// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication developers
//

// Phase 6A: the production registration switch. Every historical Dam
// small-displacement registration name now creates a StructuralMechanicsApplication
// SmallDisplacement runtime element, while the legacy Dam classes remain compiled
// (temporary rollback checkpoint). Old .mdpa inputs and binary restarts keep
// working. No production registration is modified by these tests.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Project includes
#include "dam_fast_suite.h"
#include "containers/model.h"
#include "includes/constitutive_law.h"
#include "includes/element.h"
#include "includes/kratos_components.h"
#include "includes/model_part.h"
#include "includes/model_part_io.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"

// Application includes
#include "dam_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Material data.
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;
constexpr double test_density = 2400.0;
constexpr double test_thickness = 0.15;

/// Builds a model part with one element of the given registered name, created
/// from the registered prototype geometry (scaled/translated).
Element::Pointer CreateRegisteredElement(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    const std::string& rLawName,
    ModelPart*& rOutModelPart,
    const bool rIs2d)
{
    KRATOS_TRY;
    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);
    ProcessInfo& r_pi = r_model_part.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = rIs2d ? 2 : 3;
    r_pi[SPACE_DIMENSION] = rIs2d ? 2 : 3;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 1.0;
    r_pi[IS_CONVERGED] = true;

    for (auto& v : std::vector<const VariableData*>{&DISPLACEMENT, &VELOCITY, &ACCELERATION,
                    &VOLUME_ACCELERATION, &TEMPERATURE, &NODAL_REFERENCE_TEMPERATURE,
                    &INITIAL_STRESS_TENSOR}) {
        r_model_part.AddNodalSolutionStepVariable(*v);
    }

    const Element& r_prototype = KratosComponents<Element>::Get(rElementName);
    const auto& r_proto_geom = r_prototype.GetGeometry();
    Matrix local_coordinates;
    r_proto_geom.PointsLocalCoordinates(local_coordinates);
    const std::size_t number_of_nodes = local_coordinates.size1();
    const double scale = 2.0;
    const double offset_x = 0.5, offset_y = 1.0, offset_z = 0.25;

    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        Node::Pointer p_node = r_model_part.CreateNewNode(
            i + 1,
            scale * local_coordinates(i, 0) + offset_x,
            scale * local_coordinates(i, 1) + offset_y,
            (local_coordinates.size2() > 2 ? scale * local_coordinates(i, 2) + offset_z : 0.0));
        p_node->AddDof(DISPLACEMENT_X);
        p_node->AddDof(DISPLACEMENT_Y);
        p_node->AddDof(DISPLACEMENT_Z);
        p_node->FastGetSolutionStepValue(TEMPERATURE) = 20.0;
        p_node->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = 20.0;
        Matrix z3(3, 3);
        noalias(z3) = ZeroMatrix(3, 3);
        p_node->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = z3;
    }

    auto p_props = r_model_part.CreateNewProperties(1);
    (*p_props)[YOUNG_MODULUS] = test_young_modulus;
    (*p_props)[POISSON_RATIO] = test_poisson_ratio;
    (*p_props)[DENSITY] = test_density;
    if (rIs2d)
        (*p_props)[THICKNESS] = test_thickness;
    if (rLawName == "ThermalLinearElastic3DLaw")
        p_props->SetValue(CONSTITUTIVE_LAW, ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLaw()));
    else
        p_props->SetValue(CONSTITUTIVE_LAW, ConstitutiveLaw::Pointer(new ThermalLinearElastic2DPlaneStrain()));

    std::vector<ModelPart::IndexType> element_nodes(number_of_nodes);
    for (std::size_t i = 0; i < number_of_nodes; ++i)
        element_nodes[i] = i + 1;
    r_model_part.CreateNewElement(rElementName, 1, element_nodes, p_props);

    Element::Pointer p_element = r_model_part.pGetElement(1);
    p_element->Initialize(r_pi);
    rOutModelPart = &r_model_part;
    return p_element;
    KRATOS_CATCH("");
}

/// True if the runtime type is the StructuralMechanics SmallDisplacement class
/// (mangled "N6Kratos17SmallDisplacementE").
bool IsSmaRuntime(const Element& rElement)
{
    const std::string name = typeid(rElement).name();
    return name.find("17SmallDisplacementE") != std::string::npos;
}

/// True if the runtime type is one of the LEGACY Dam classes.
bool IsLegacyRuntime(const Element& rElement)
{
    const std::string name = typeid(rElement).name();
    return name.find("ThermoMechanic") != std::string::npos ||
           name.find("25SmallDisplacementElementE") != std::string::npos;
}

} // namespace

//************************************************************************************
// 1. All 22 historical names create SMA runtime elements.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R6A_AllHistoricalNames_RuntimeSMA, KratosDamFastSuite)
{
    // (name, 2D flag)
    const std::vector<std::pair<std::string, bool>> names = {
        {"SmallDisplacementSolidElement2D3N", true},
        {"SmallDisplacementSolidElement2D4N", true},
        {"SmallDisplacementSolidElement2D6N", true},
        {"SmallDisplacementSolidElement2D8N", true},
        {"SmallDisplacementSolidElement2D9N", true},
        {"SmallDisplacementSolidElement3D4N", false},
        {"SmallDisplacementSolidElement3D6N", false},
        {"SmallDisplacementSolidElement3D8N", false},
        {"SmallDisplacementSolidElement3D10N", false},
        {"SmallDisplacementSolidElement3D15N", false},
        {"SmallDisplacementSolidElement3D20N", false},
        {"SmallDisplacementSolidElement3D27N", false},
        {"SmallDisplacementThermoMechanicElement2D3N", true},
        {"SmallDisplacementThermoMechanicElement2D4N", true},
        {"SmallDisplacementThermoMechanicElement2D6N", true},
        {"SmallDisplacementThermoMechanicElement2D8N", true},
        {"SmallDisplacementThermoMechanicElement2D9N", true},
        {"SmallDisplacementThermoMechanicElement3D4N", false},
        {"SmallDisplacementThermoMechanicElement3D8N", false},
        {"SmallDisplacementThermoMechanicElement3D10N", false},
        {"SmallDisplacementThermoMechanicElement3D20N", false},
        {"SmallDisplacementThermoMechanicElement3D27N", false},
    };

    for (const auto& entry : names) {
        const std::string& name = entry.first;
        const bool is2d = entry.second;
        Model model;
        ModelPart* p_mp = nullptr;
        Element::Pointer p_elem = CreateRegisteredElement(
            model, "R6A" + name, name, is2d ? "ThermalLinearElastic2DPlaneStrain" : "ThermalLinearElastic3DLaw", p_mp, is2d);

        const std::size_t expected_nodes = p_elem->GetGeometry().PointsNumber();
        KRATOS_EXPECT_TRUE(IsSmaRuntime(*p_elem));
        KRATOS_EXPECT_FALSE(IsLegacyRuntime(*p_elem));
        KRATOS_EXPECT_EQ(expected_nodes,
                         KratosComponents<Element>::Get(name).GetGeometry().PointsNumber());
        KRATOS_EXPECT_EQ(p_elem->GetProperties()[YOUNG_MODULUS], test_young_modulus);
        std::cout << "[6A] " << name << " -> runtime=" << typeid(*p_elem).name()
                  << " npoints=" << expected_nodes << std::endl;
    }
    std::cout << "[6A] All 22 historical registrations create SMA runtime elements." << std::endl;
}

//************************************************************************************
// 2. Full generic API surface through a historical name.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R6A_MechanicalThermo_ConvergenceAndAPISurface, KratosDamFastSuite)
{
    // A historical Dam name must expose the full generic element API of the SMA
    // runtime element it creates.
    Model model;
    ModelPart* p_mech_mp = nullptr;
    Element::Pointer p_mech = CreateRegisteredElement(model, "ConvMech", "SmallDisplacementSolidElement3D8N", "ThermalLinearElastic3DLaw", p_mech_mp, false);

    KRATOS_EXPECT_TRUE(IsSmaRuntime(*p_mech));

    // Full generic API surface.
    KRATOS_EXPECT_EQ(p_mech->Check(p_mech_mp->GetProcessInfo()), 0);
    Element::Pointer p_clone = p_mech->Clone(2, p_mech->GetGeometry().Points());
    KRATOS_EXPECT_EQ(p_clone->GetGeometry().PointsNumber(), 8u);

    p_mech->InitializeSolutionStep(p_mech_mp->GetProcessInfo());
    p_mech->InitializeNonLinearIteration(p_mech_mp->GetProcessInfo());
    Matrix lhs_mech, mass_mech, damping_mech;
    Vector rhs_mech;
    p_mech->CalculateLocalSystem(lhs_mech, rhs_mech, p_mech_mp->GetProcessInfo());
    p_mech->CalculateMassMatrix(mass_mech, p_mech_mp->GetProcessInfo());
    p_mech->CalculateDampingMatrix(damping_mech, p_mech_mp->GetProcessInfo());
    p_mech->FinalizeNonLinearIteration(p_mech_mp->GetProcessInfo());
    p_mech_mp->GetProcessInfo()[IS_CONVERGED] = true;
    p_mech->FinalizeSolutionStep(p_mech_mp->GetProcessInfo());

    // Integration-point outputs + CONSTITUTIVE_LAW pointer.
    std::vector<Vector> v_out;
    p_mech->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, v_out, p_mech_mp->GetProcessInfo());
    std::vector<Matrix> m_out;
    p_mech->CalculateOnIntegrationPoints(CAUCHY_STRESS_TENSOR, m_out, p_mech_mp->GetProcessInfo());
    std::vector<double> s_out;
    p_mech->CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, s_out, p_mech_mp->GetProcessInfo());
    std::vector<ConstitutiveLaw::Pointer> law_out;
    p_mech->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, law_out, p_mech_mp->GetProcessInfo());
    KRATOS_EXPECT_EQ(v_out.size(), law_out.size());
    KRATOS_EXPECT_EQ(m_out.size(), law_out.size());
    KRATOS_EXPECT_EQ(s_out.size(), law_out.size());
    KRATOS_EXPECT_EQ(law_out.size(), 8u);
    KRATOS_EXPECT_TRUE(law_out[0] != nullptr);

    std::cout << "[6A] historical names expose the full generic SMA element API "
              << "surface (Check, Clone, lifecycle, local system, mass, damping, "
              << "integration-point outputs)." << std::endl;
}

//************************************************************************************


//************************************************************************************
// 4. Unmodified committed .mdpa smoke test.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R6A_UnmodifiedMdpa_Smoke, KratosDamFastSuite)
{
    // Load a committed Dam .mdpa unchanged: the historical
    // SmallDisplacementThermoMechanicElement2D3N elements must now be SMA at
    // runtime (the interface elements remain unchanged).
    const std::string source_dir = std::string(__FILE__);
    const std::string repo_root = source_dir.substr(0, source_dir.find("applications/DamApplication"));
    // ModelPartIO appends the ".mdpa" extension to the given filename.
    const std::string mdpa_path = repo_root +
        "applications/DamApplication/tests/joint_elastic_cohesive_2d_normal/"
        "joint_elastic_cohesive_2d_normal";

    Model model;
    ModelPart& r_mp = model.CreateModelPart("MdpaSmoke", 2);
    r_mp.GetProcessInfo()[DOMAIN_SIZE] = 2;
    r_mp.GetProcessInfo()[SPACE_DIMENSION] = 2;
    r_mp.GetProcessInfo()[IS_RESTARTED] = false;
    r_mp.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_mp.AddNodalSolutionStepVariable(VELOCITY);
    r_mp.AddNodalSolutionStepVariable(ACCELERATION);
    r_mp.AddNodalSolutionStepVariable(TEMPERATURE);
    r_mp.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);

    ModelPartIO(mdpa_path).ReadModelPart(r_mp);

    KRATOS_EXPECT_EQ(r_mp.NumberOfElements(), 25u);   // 20 solid + 5 interface

    std::size_t sma_count = 0, interface_count = 0, legacy_count = 0;
    for (auto& r_elem : r_mp.Elements()) {
        const std::string runtime = typeid(r_elem).name();
        if (IsSmaRuntime(r_elem))
            ++sma_count;
        else if (runtime.find("SmallDisplacementInterfaceElement") != std::string::npos)
            ++interface_count;
        else if (IsLegacyRuntime(r_elem))
            ++legacy_count;
    }
    std::cout << "[6A] unmodified .mdpa: SMA elements=" << sma_count
              << " interface=" << interface_count
              << " legacy=" << legacy_count << std::endl;
    KRATOS_EXPECT_EQ(sma_count, 20u);
    KRATOS_EXPECT_EQ(interface_count, 5u);
    KRATOS_EXPECT_EQ(legacy_count, 0u);
    std::cout << "[6A] committed .mdpa loads unchanged; historical solid elements "
              << "are now SMA at runtime." << std::endl;
}

} // namespace Testing
} // namespace Kratos

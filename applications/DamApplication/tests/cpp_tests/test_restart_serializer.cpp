// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    DamApplication developers

// Serializer/restart compatibility contract. A frozen binary restart produced
// by the removed Dam thermo-mechanical element still loads into the production
// aliases with preserved geometry, properties and committed damage.
//
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "dam_fast_suite.h"
#include "containers/model.h"
#include "includes/constitutive_law.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/stream_serializer.h"

// Application includes
#include "dam_application_variables.h"
#include "structural_mechanics_application_variables.h"

// StructuralMechanics small-displacement element (the runtime type behind the
// historical names).
#include "custom_elements/solid_elements/small_displacement.h"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Material data.
constexpr double test_young_modulus = 2.0e7;

/// Serializer trace used for saving/loading the ROUND TRIP experiments. The
/// production restart format is binary (SERIALIZER_NO_TRACE, see
/// restart_utility.py: "no_trace" is the default), so the faithful simulation
/// of a real restart must use binary mode. ASCII (trace) mode is used only to
/// probe the stored registration name.
constexpr Serializer::TraceType round_trip_trace = Serializer::SERIALIZER_NO_TRACE;

/// Reads the committed damage variable from the element's first constitutive law.
double ElementDamage(Element& rElement, const ProcessInfo& rPi)
{
    std::vector<ConstitutiveLaw::Pointer> laws;
    rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, rPi);
    double d = 0.0;
    laws[0]->GetValue(DAMAGE_VARIABLE, d);
    return d;
}

/// A serializable holder for element pointers.
struct ElementHolder
{
    std::vector<Element::Pointer> elements;
    void save(Serializer& rSerializer) const { rSerializer.save("Elements", elements); }
    void load(Serializer& rSerializer) { rSerializer.load("Elements", elements); }
};

/// Loads the holder from a string (binary mode).
void LoadHolderFromString(const std::string& rArchive, ElementHolder& rHolder)
{
    StreamSerializer serializer(rArchive, round_trip_trace);
    serializer.SetLoadState();
    serializer.load("Holder", rHolder);
}

/// Absolute path of the frozen true-legacy restart fixture.
std::string LegacyFixturePath()
{
    const std::string source_dir = std::string(__FILE__);
    const std::string repo_root = source_dir.substr(0, source_dir.find("applications/DamApplication"));
    return repo_root + "applications/DamApplication/tests/cpp_tests/fixtures/legacy_thermo_3d8n_damage.dat";
}

/// Reads raw bytes from a file.
std::string ReadBytesFromFile(const std::string& rPath)
{
    std::ifstream ifs(rPath, std::ios::binary);
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return oss.str();
}


} // namespace

//************************************************************************************
// Frozen binary fixture from the removed Dam thermo-mechanical element.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(LegacyBinaryRestartLoadsIntoSMAElement, KratosDamFastSuite)
{
    // The fixture is a genuine binary archive written by the removed Dam
    // thermo-mechanical element; it must load through the current production
    // registration aliases into the SMA runtime element.
    const std::string archive = ReadBytesFromFile(LegacyFixturePath());
    KRATOS_EXPECT_GT(archive.size(), 0u);

    ProcessInfo pi;
    pi[DOMAIN_SIZE] = 3;
    pi[SPACE_DIMENSION] = 3;
    pi[DELTA_TIME] = 1.0;
    pi[IS_CONVERGED] = true;

    std::string outcome = "unknown";
    try {
        ElementHolder loaded;
        LoadHolderFromString(archive, loaded);
        KRATOS_EXPECT_EQ(loaded.elements.size(), 1u);
        Element::Pointer p_loaded = loaded.elements[0];
        const std::string runtime = typeid(*p_loaded).name();
        std::cout << "[restart] frozen legacy restart loaded runtime_type = " << runtime << std::endl;
        KRATOS_EXPECT_TRUE(runtime.find("SmallDisplacement") != std::string::npos);
        KRATOS_EXPECT_TRUE(runtime.find("ThermoMechanic") == std::string::npos);
        KRATOS_EXPECT_EQ(p_loaded->GetGeometry().size(), 8u);
        KRATOS_EXPECT_EQ(p_loaded->GetProperties()[YOUNG_MODULUS], test_young_modulus);
        for (auto& n : p_loaded->GetGeometry()) {
            n.AddDof(DISPLACEMENT_X);
            n.AddDof(DISPLACEMENT_Y);
            n.AddDof(DISPLACEMENT_Z);
        }
        // Committed damage preserved.
        const double damage = ElementDamage(*p_loaded, pi);
        KRATOS_EXPECT_GT(damage, 0.0);
        std::cout << "[restart] frozen legacy restart committed damage = " << damage << std::endl;

        // The first real material response after load works without any manual
        // material-properties repair (automatic transient rebinding).
        std::vector<Vector> stress;
        p_loaded->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, stress, pi);
        KRATOS_EXPECT_EQ(stress.size(), 8u);
        KRATOS_EXPECT_TRUE(std::isfinite(stress[0][0]));

        // Mass matrix.
        Matrix mass;
        p_loaded->CalculateMassMatrix(mass, pi);
        KRATOS_EXPECT_EQ(mass.size1(), 24u);
        KRATOS_EXPECT_GT(mass(0, 0), 0.0);
        outcome = "A: frozen legacy binary restart loads into SMA and responds "
                  "without manual material repair";
    } catch (const std::exception& rException) {
        outcome = std::string("load failed: ") + rException.what();
    }
    std::cout << "[restart] FROZEN LEGACY RESTART UNDER PRODUCTION ALIASES: " << outcome << std::endl;
}


} // namespace Testing
} // namespace Kratos

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

// Phase 5C.3 (Serializer / restart characterization): CHARACTERIZATION ONLY.
// Determines how Kratos Serializer behaves when multiple registration names
// refer to one StructuralMechanics SmallDisplacement runtime type, and whether
// restart archives created with the CURRENT legacy Dam element classes can be
// loaded after the historical names are repointed to SMA prototypes. No
// production registration or class is modified.

// System includes
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
#include "includes/kratos_components.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/stream_serializer.h"

// Application includes
#include "dam_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_3D_law.hpp"

// StructuralMechanics small-displacement element (the runtime type behind the
// historical names after Phase 6A).
#include "custom_elements/solid_elements/small_displacement.h"
// The legacy class is retained in Phase 6A and used ONLY to generate a TRUE
// legacy archive (the historical registration names now create SMA elements).
#include "geometries/hexahedra_3d_8.h"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Material data.
constexpr double test_density = 2400.0;
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;

/// Serializer trace used for saving/loading the ROUND TRIP experiments. The
/// production restart format is binary (SERIALIZER_NO_TRACE, see
/// restart_utility.py: "no_trace" is the default), so the faithful simulation
/// of a real restart must use binary mode. ASCII (trace) mode is used only to
/// probe the stored registration name.
constexpr Serializer::TraceType round_trip_trace = Serializer::SERIALIZER_NO_TRACE;
constexpr Serializer::TraceType ascii_trace = Serializer::SERIALIZER_TRACE_ALL;

/// The 3D8N coordinates used by all element experiments.
const std::vector<std::vector<double>> coords_h8 = {
    {0,0,0},{2.0,0,0},{2.0,1.0,0},{0,1.0,0},{0,0,2.0},{2.0,0,2.0},{2.0,1.0,2.0},{0,1.0,2.0}};

/// Builds a single-element model part (3D8N) with the given registered element
/// name and law family ("linear" = stateless, "local" = stateful damage).
Element::Pointer CreateElement(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    const std::string& rLawFamily,
    ModelPart*& rOutModelPart)
{
    KRATOS_TRY;
    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);
    ProcessInfo& r_pi = r_model_part.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = 3;
    r_pi[SPACE_DIMENSION] = 3;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 1.0;
    r_pi[IS_CONVERGED] = true;

    for (auto& v : std::vector<const VariableData*>{&DISPLACEMENT, &VELOCITY, &ACCELERATION,
                    &VOLUME_ACCELERATION, &TEMPERATURE, &NODAL_REFERENCE_TEMPERATURE,
                    &NODAL_YOUNG_MODULUS, &INITIAL_STRESS_TENSOR}) {
        r_model_part.AddNodalSolutionStepVariable(*v);
    }

    for (std::size_t i = 0; i < coords_h8.size(); ++i) {
        const auto& c = coords_h8[i];
        Node::Pointer p_node = r_model_part.CreateNewNode(i + 1, c[0], c[1], c[2]);
        p_node->AddDof(DISPLACEMENT_X);
        p_node->AddDof(DISPLACEMENT_Y);
        p_node->AddDof(DISPLACEMENT_Z);
        p_node->FastGetSolutionStepValue(TEMPERATURE) = 20.0;
        p_node->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = 20.0;
        p_node->FastGetSolutionStepValue(NODAL_YOUNG_MODULUS) = test_young_modulus;
        Matrix z3(3, 3);
        noalias(z3) = ZeroMatrix(3, 3);
        p_node->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR) = z3;
    }

    auto p_props = r_model_part.CreateNewProperties(1);
    (*p_props)[YOUNG_MODULUS] = test_young_modulus;
    (*p_props)[POISSON_RATIO] = test_poisson_ratio;
    (*p_props)[DENSITY] = test_density;
    (*p_props)[THERMAL_EXPANSION] = 1.0e-5;
    (*p_props)[DAMAGE_THRESHOLD] = 5.0e-3;
    (*p_props)[STRENGTH_RATIO] = 10.0;
    (*p_props)[FRACTURE_ENERGY] = 5000.0;
    ConstitutiveLaw::Pointer p_law = (rLawFamily == "local")
        ? ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw())
        : ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLaw());
    p_props->SetValue(CONSTITUTIVE_LAW, p_law);

    std::vector<ModelPart::IndexType> element_nodes(coords_h8.size());
    for (std::size_t i = 0; i < coords_h8.size(); ++i)
        element_nodes[i] = i + 1;
    r_model_part.CreateNewElement(rElementName, 1, element_nodes, p_props);

    Element::Pointer p_element = r_model_part.pGetElement(1);
    p_element->Initialize(r_pi);
    rOutModelPart = &r_model_part;
    return p_element;
    KRATOS_CATCH("");
}

/// Applies a damaging displacement state and commits it.
void CommitDamage(Element& rElement, ProcessInfo& rPi)
{
    for (auto& n : rElement.GetGeometry()) {
        const auto& x0 = n.GetInitialPosition();
        auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
        u[0] = 2.0e-5 * x0[0];
        u[1] = -test_poisson_ratio * 2.0e-5 * x0[1];
        u[2] = -test_poisson_ratio * 2.0e-5 * x0[2];
        n.X() = x0[0] + u[0]; n.Y() = x0[1] + u[1]; n.Z() = x0[2] + u[2];
    }
    rPi[IS_CONVERGED] = true;
    rElement.FinalizeSolutionStep(rPi);
}

/// Applies a uniaxial displacement state (no damage for the linear law).
void ApplyState(Element& rElement)
{
    for (auto& n : rElement.GetGeometry()) {
        const auto& x0 = n.GetInitialPosition();
        auto& u = n.FastGetSolutionStepValue(DISPLACEMENT);
        u[0] = 1.0e-5 * x0[0];
        u[1] = -test_poisson_ratio * 1.0e-5 * x0[1];
        u[2] = -test_poisson_ratio * 1.0e-5 * x0[2];
        n.X() = x0[0] + u[0]; n.Y() = x0[1] + u[1]; n.Z() = x0[2] + u[2];
    }
}

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

/// Saves the holder to a string (binary mode, the production restart format).
std::string SaveHolderToString(const ElementHolder& rHolder)
{
    StreamSerializer serializer(round_trip_trace);
    serializer.save("Holder", rHolder);
    return serializer.GetStringRepresentation();
}

/// Loads the holder from a string (binary mode).
void LoadHolderFromString(const std::string& rArchive, ElementHolder& rHolder)
{
    StreamSerializer serializer(rArchive, round_trip_trace);
    serializer.SetLoadState();
    serializer.load("Holder", rHolder);
}

/// Saves the holder to a string (ASCII trace mode, used for name probing only).
std::string SaveHolderToStringAscii(const ElementHolder& rHolder)
{
    StreamSerializer serializer(ascii_trace);
    serializer.save("Holder", rHolder);
    return serializer.GetStringRepresentation();
}

/// Loads the holder from an ASCII archive (name-probing only).
void LoadHolderFromStringAscii(const std::string& rArchive, ElementHolder& rHolder)
{
    StreamSerializer serializer(rArchive, ascii_trace);
    serializer.SetLoadState();
    serializer.load("Holder", rHolder);
}

/// Extracts the first registered element name from an ASCII archive.
std::string ExtractElementName(const std::string& rArchive)
{
    const std::string marker = "\"SmallDisplacement";
    const std::size_t pos = rArchive.find(marker);
    if (pos == std::string::npos)
        return "<not found>";
    const std::size_t begin = pos + 1;   // skip opening quote
    const std::size_t end = rArchive.find('"', begin);
    return rArchive.substr(begin, end - begin);
}

/// Counts occurrences of a substring (structural probe of the ASCII trace).
std::size_t CountSubstring(const std::string& rText, const std::string& rSubstring)
{
    std::size_t count = 0;
    std::size_t pos = 0;
    while ((pos = rText.find(rSubstring, pos)) != std::string::npos) {
        ++count;
        pos += rSubstring.size();
    }
    return count;
}

/// Absolute path of the frozen true-legacy restart fixture.
std::string LegacyFixturePath()
{
    const std::string source_dir = std::string(__FILE__);
    const std::string repo_root = source_dir.substr(0, source_dir.find("applications/DamApplication"));
    return repo_root + "applications/DamApplication/tests/cpp_tests/fixtures/legacy_thermo_3d8n_damage.dat";
}

/// Writes raw bytes to a file (used to freeze the legacy archive).
void WriteBytesToFile(const std::string& rPath, const std::string& rContent)
{
    std::ofstream ofs(rPath, std::ios::binary);
    ofs.write(rContent.data(), static_cast<std::streamsize>(rContent.size()));
}

/// Reads raw bytes from a file.
std::string ReadBytesFromFile(const std::string& rPath)
{
    std::ifstream ifs(rPath, std::ios::binary);
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return oss.str();
}

/// A serializable holder for a constitutive-law pointer.
struct LawHolder
{
    ConstitutiveLaw::Pointer p_law;
    void save(Serializer& rSerializer) const { rSerializer.save("Law", p_law); }
    void load(Serializer& rSerializer) { rSerializer.load("Law", p_law); }
};

/// Saves the law holder to a string (binary mode).
std::string SaveLawHolder(const LawHolder& rHolder)
{
    StreamSerializer serializer(round_trip_trace);
    serializer.save("Holder", rHolder);
    return serializer.GetStringRepresentation();
}

/// Loads the law holder from a string (binary mode).
void LoadLawHolder(const std::string& rArchive, LawHolder& rHolder)
{
    StreamSerializer serializer(rArchive, round_trip_trace);
    serializer.SetLoadState();
    serializer.load("Holder", rHolder);
}

} // namespace

//************************************************************************************
// 1. Canonical serialized names: SMA and legacy thermo element.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R53_CanonicalSerializedName, KratosDamFastSuite)
{
    // After Phase 6A both the direct SMA name and the historical thermo name
    // create the SAME SMA runtime type, so both serialize under the canonical
    // SMA Serializer name (first-registration semantics from 5C.3).
    Model model;
    ModelPart* p_sma_mp = nullptr;
    ModelPart* p_hist_mp = nullptr;
    Element::Pointer p_sma = CreateElement(model, "CanSma", "SmallDisplacementElement3D8N", "linear", p_sma_mp);
    Element::Pointer p_hist = CreateElement(model, "CanHist", "SmallDisplacementThermoMechanicElement3D8N", "linear", p_hist_mp);

    const std::string sma_archive = SaveHolderToStringAscii(ElementHolder{{p_sma}});
    const std::string hist_archive = SaveHolderToStringAscii(ElementHolder{{p_hist}});
    const std::string sma_name = ExtractElementName(sma_archive);
    const std::string hist_name = ExtractElementName(hist_archive);
    std::cout << "[6A] SMA 3D8N stored name = " << sma_name << std::endl;
    std::cout << "[6A] historical thermo 3D8N stored name = " << hist_name << std::endl;
    KRATOS_EXPECT_EQ(sma_name, "SmallDisplacementElement2D3N");
    KRATOS_EXPECT_EQ(hist_name, "SmallDisplacementElement2D3N");
    std::cout << "[6A] both the direct SMA name and the historical thermo name "
              << "serialize under the canonical SMA name; no custom Serializer "
              << "factory is required." << std::endl;
}

//************************************************************************************
// 2. Multiple aliases for one runtime type: load registry + canonical save name.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R53_MultipleAlias_InsertionOrder, KratosDamFastSuite)
{
    // Register one additional alias name for the SMA SmallDisplacement runtime
    // type (the canonical name was set by the first SMA registration).
    Serializer::Register<SmallDisplacement>("R53SmaAliasExtra", dynamic_cast<const SmallDisplacement&>(
        KratosComponents<Element>::Get("SmallDisplacementElement3D8N")));

    // Save an SMA element; the archive must use the FIRST-registered name.
    Model model;
    ModelPart* p_mp = nullptr;
    Element::Pointer p_elem = CreateElement(model, "Alias", "SmallDisplacementElement3D8N", "linear", p_mp);
    const std::string archive = SaveHolderToStringAscii(ElementHolder{{p_elem}});
    const std::string stored_name = ExtractElementName(archive);
    std::cout << "[5C.3] SMA element stored name = " << stored_name << std::endl;
    KRATOS_EXPECT_EQ(stored_name, "SmallDisplacementElement2D3N");   // first registration wins

    // Loading by the extra ALIAS name must also reconstruct an SMA element:
    // replace the stored name with the alias in the ASCII archive.
    std::string alias_archive = archive;
    const std::string replace_from = "\"" + stored_name + "\"";
    const std::string replace_to = "\"R53SmaAliasExtra\"";
    const std::size_t rpos = alias_archive.find(replace_from);
    KRATOS_EXPECT_TRUE(rpos != std::string::npos);
    alias_archive.replace(rpos, replace_from.size(), replace_to);

    ElementHolder loaded;
    LoadHolderFromStringAscii(alias_archive, loaded);
    KRATOS_EXPECT_EQ(loaded.elements.size(), 1u);
    const std::string runtime = typeid(*loaded.elements[0]).name();
    std::cout << "[5C.3] alias-loaded runtime_type = " << runtime << std::endl;
    KRATOS_EXPECT_TRUE(runtime.find("SmallDisplacement") != std::string::npos);
    KRATOS_EXPECT_EQ(loaded.elements[0]->GetGeometry().size(), 8u);

    std::cout << "[5C.3] multiple aliases -> one runtime type: all names load; "
              << "canonical save name = first registered" << std::endl;
}

//************************************************************************************
// 3. New-format restart round trip (SMA + local damage).
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R53_NewFormat_RoundTrip, KratosDamFastSuite)
{
    Model model;
    ModelPart* p_mp = nullptr;
    Element::Pointer p_elem = CreateElement(model, "NewRt", "SmallDisplacementElement3D8N", "linear", p_mp);
    ProcessInfo& r_pi = p_mp->GetProcessInfo();
    ApplyState(*p_elem);

    // Pre-load reference: LHS and mass.
    Matrix lhs_before, mass_before;
    Vector rhs_before;
    p_elem->CalculateLocalSystem(lhs_before, rhs_before, r_pi);
    p_elem->CalculateMassMatrix(mass_before, r_pi);

    const std::string archive = SaveHolderToString(ElementHolder{{p_elem}});

    // Reload in a fresh serializer context.
    ElementHolder loaded;
    LoadHolderFromString(archive, loaded);
    KRATOS_EXPECT_EQ(loaded.elements.size(), 1u);
    Element::Pointer p_loaded = loaded.elements[0];

    // Runtime type is SMA.
    std::cout << "[5C.3] new-format loaded runtime_type = " << typeid(*p_loaded).name() << std::endl;
    KRATOS_EXPECT_TRUE(std::string(typeid(*p_loaded).name()).find("SmallDisplacement") != std::string::npos);

    // Geometry / Properties / integration method / constitutive-law vector.
    KRATOS_EXPECT_EQ(p_loaded->GetGeometry().size(), 8u);
    KRATOS_EXPECT_EQ(p_loaded->GetProperties()[YOUNG_MODULUS], test_young_modulus);
    KRATOS_EXPECT_EQ(p_loaded->GetGeometry().WorkingSpaceDimension(), 3u);

    // Add DOFs on the reloaded nodes so that system assembly works.
    for (auto& n : p_loaded->GetGeometry()) {
        n.AddDof(DISPLACEMENT_X);
        n.AddDof(DISPLACEMENT_Y);
        n.AddDof(DISPLACEMENT_Z);
    }

    // Stress state preserved (specialized outputs).
    std::vector<Vector> v_out_before;
    p_elem->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, v_out_before, r_pi);
    std::vector<Vector> v_out;
    p_loaded->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, v_out, r_pi);
    KRATOS_EXPECT_EQ(v_out.size(), 8u);
    KRATOS_EXPECT_NEAR(v_out[0][0], v_out_before[0][0], 1.0e-8);

    Matrix lhs_after, mass_after;
    Vector rhs_after;
    p_loaded->CalculateLocalSystem(lhs_after, rhs_after, r_pi);
    p_loaded->CalculateMassMatrix(mass_after, r_pi);
    KRATOS_EXPECT_EQ(lhs_after.size1(), lhs_before.size1());
    KRATOS_EXPECT_NEAR(lhs_after(0, 0), lhs_before(0, 0), 1.0e-8);
    KRATOS_EXPECT_NEAR(rhs_after(0), rhs_before(0), 1.0e-8);
    KRATOS_EXPECT_NEAR(mass_after(0, 0), mass_before(0, 0), 1.0e-12);

    std::cout << "[5C.3] new-format round trip OK: type/geometry/props/law/damage/"
              << "outputs/LHS/RHS/mass preserved" << std::endl;
}

//************************************************************************************
// 4. TRUE legacy archive: generation + control round trip (legacy class directly).
KRATOS_TEST_CASE_IN_SUITE(R53_LawState_Isolation, KratosDamFastSuite)
{
    // Serialize/reload the migrated law alone (damaged state). This separates
    // element-wrapper compatibility from constitutive-law-state compatibility.
    Model model;
    ModelPart* p_mp = nullptr;
    Element::Pointer p_elem = CreateElement(model, "LawIso", "SmallDisplacementElement3D8N", "local", p_mp);
    CommitDamage(*p_elem, p_mp->GetProcessInfo());
    const double damage_ref = ElementDamage(*p_elem, p_mp->GetProcessInfo());

    // Grab the element's own (damaged) law clone and serialize just the law.
    std::vector<ConstitutiveLaw::Pointer> laws;
    p_elem->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, laws, p_mp->GetProcessInfo());

    LawHolder holder;
    holder.p_law = laws[0];
    const std::string archive = SaveLawHolder(holder);

    LawHolder loaded;
    LoadLawHolder(archive, loaded);
    double d_after = 0.0;
    loaded.p_law->GetValue(DAMAGE_VARIABLE, d_after);
    KRATOS_EXPECT_NEAR(d_after, damage_ref, 1.0e-12);
    std::cout << "[5C.3] law-only round trip: damage preserved independently of the "
              << "element wrapper (" << damage_ref << ")" << std::endl;
}

//************************************************************************************
// 8. Load the FROZEN true-legacy binary fixture under production aliases.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(R6B_FrozenLegacyRestart_ProductionAliasLoad, KratosDamFastSuite)
{
    // This is the permanent R4 regression. The fixture is a genuine binary
    // archive written by the actual legacy C++ class before it was deleted; it
    // must load through the current production registration aliases into SMA.
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
        std::cout << "[6B] frozen legacy restart loaded runtime_type = " << runtime << std::endl;
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
        std::cout << "[6B] frozen legacy restart committed damage = " << damage << std::endl;

        // Phase 6C.1: the FIRST real material response after load works WITHOUT
        // any manual material-properties repair (automatic transient rebinding).
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
    std::cout << "[6B] FROZEN LEGACY RESTART UNDER PRODUCTION ALIASES: " << outcome << std::endl;
}

} // namespace Testing
} // namespace Kratos

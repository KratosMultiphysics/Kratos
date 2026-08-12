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

// Phase 5C.1 (MASS-01): CHARACTERIZATION ONLY. Determine whether the mass-matrix
// behavior of the legacy Dam SmallDisplacementElement /
// SmallDisplacementThermoMechanicElement can be reproduced by the
// StructuralMechanicsApplication SmallDisplacement element without modifying
// StructuralMechanicsApplication. No production code is modified.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
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
#include "includes/variables.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

// Application includes
#include "dam_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "custom_constitutive/thermal_simo_ju_local_damage_3D_law.hpp"
#include "custom_constitutive/thermal_simo_ju_nonlocal_damage_3D_law.hpp"

namespace Kratos
{
namespace Testing
{

namespace
{

/// Tight comparison tolerance for identical-quadrature / analytical references.
constexpr double exact_tolerance = 1.0e-12;

/// Material data.
constexpr double test_density = 2400.0;
constexpr double test_thickness = 0.15;
constexpr double test_young_modulus = 2.0e7;
constexpr double test_poisson_ratio = 0.2;

/// Constructs the constitutive law directly (registered-law lookup is avoided
/// in the characterization tests).
ConstitutiveLaw::Pointer CreateLaw(const std::string& rName)
{
    if (rName == "ThermalLinearElastic2DPlaneStrain")
        return ConstitutiveLaw::Pointer(new ThermalLinearElastic2DPlaneStrain());
    if (rName == "ThermalLinearElastic3DLaw")
        return ConstitutiveLaw::Pointer(new ThermalLinearElastic3DLaw());
    if (rName == "ThermalSimoJuLocalDamage3DLaw")
        return ConstitutiveLaw::Pointer(new ThermalSimoJuLocalDamage3DLaw());
    if (rName == "ThermalSimoJuNonlocalDamage3DLaw")
        return ConstitutiveLaw::Pointer(new ThermalSimoJuNonlocalDamage3DLaw());
    KRATOS_ERROR << "Unknown law name " << rName << std::endl;
}

/// Builds a 2D or 3D model part with one element of the given registered name,
/// created from the given nodal coordinates. The model part pointer is stored in
/// rOutModelPart so the caller can access the ProcessInfo.
Element::Pointer CreateElement(
    Model& rModel,
    const std::string& rModelPartName,
    const std::string& rElementName,
    const std::vector<std::vector<double>>& rCoords,
    const bool rIs2d,
    const double rDensity,
    const bool rSetThickness,
    const double rThickness,
    const bool rSetLumpedFlag,
    ModelPart*& rOutModelPart)
{
    KRATOS_TRY;
    ModelPart& r_model_part = rModel.CreateModelPart(rModelPartName, 2);
    ProcessInfo& r_pi = r_model_part.GetProcessInfo();
    r_pi[DOMAIN_SIZE] = rIs2d ? 2 : 3;
    r_pi[SPACE_DIMENSION] = rIs2d ? 2 : 3;
    r_pi[IS_RESTARTED] = false;
    r_pi[DELTA_TIME] = 1.0;
    r_pi[IS_CONVERGED] = true;
    if (rSetLumpedFlag)
        r_pi[COMPUTE_LUMPED_MASS_MATRIX] = true;

    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    r_model_part.AddNodalSolutionStepVariable(NODAL_REFERENCE_TEMPERATURE);
    for (std::size_t i = 0; i < rCoords.size(); ++i) {
        const auto& c = rCoords[i];
        Node::Pointer p_node = r_model_part.CreateNewNode(i + 1, c[0], c[1], c[2]);
        p_node->AddDof(DISPLACEMENT_X);
        p_node->AddDof(DISPLACEMENT_Y);
        p_node->AddDof(DISPLACEMENT_Z);
        p_node->FastGetSolutionStepValue(TEMPERATURE) = 0.0;
        p_node->FastGetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE) = 0.0;
    }

    auto p_props = r_model_part.CreateNewProperties(1);
    (*p_props)[YOUNG_MODULUS] = test_young_modulus;
    (*p_props)[POISSON_RATIO] = test_poisson_ratio;
    (*p_props)[DENSITY] = rDensity;
    if (rSetThickness)
        (*p_props)[THICKNESS] = rThickness;

    const std::string law_name = rIs2d ? "ThermalLinearElastic2DPlaneStrain" : "ThermalLinearElastic3DLaw";
    p_props->SetValue(CONSTITUTIVE_LAW, CreateLaw(law_name));

    Geometry<Node>::PointsArrayType pts;
    for (std::size_t i = 0; i < rCoords.size(); ++i)
        pts.push_back(r_model_part.pGetNode(i + 1));

    Element::Pointer p_element = KratosComponents<Element>::Get(rElementName).Create(1, pts, p_props);
    p_element->Initialize(r_pi);
    rOutModelPart = &r_model_part;
    return p_element;
    KRATOS_CATCH("");
}

/// Returns the element mass matrix.
Matrix ElementMass(Element& rElement, const ProcessInfo& rPi)
{
    Matrix M;
    rElement.CalculateMassMatrix(M, rPi);
    return M;
}

/// Sum of the first-direction (x) translational block of the mass matrix. Each
/// displacement component carries the same nodal block, so this equals the
/// physical scalar element mass.
double TotalDirectionalMass(const Matrix& rMassMatrix, const std::size_t rDimension)
{
    const std::size_t number_of_nodes = rMassMatrix.size1() / rDimension;
    double total = 0.0;
    for (std::size_t i = 0; i < number_of_nodes; ++i)
        for (std::size_t j = 0; j < number_of_nodes; ++j)
            total += rMassMatrix(i * rDimension, j * rDimension);
    return total;
}

/// Max absolute entry difference between two matrices.
double MaxAbsDiff(const Matrix& rA, const Matrix& rB)
{
    double max_diff = 0.0;
    for (std::size_t i = 0; i < rA.size1(); ++i)
        for (std::size_t j = 0; j < rA.size2(); ++j)
            max_diff = std::max(max_diff, std::abs(rA(i, j) - rB(i, j)));
    return max_diff;
}

/// Scale-aware tolerance for the Dam+flag == SMA assertions on large matrices
/// (roundoff accumulation grows with the matrix size).
double MassTolerance(const Matrix& rReference)
{
    double max_abs_entry = 0.0;
    for (std::size_t i = 0; i < rReference.size1(); ++i)
        for (std::size_t j = 0; j < rReference.size2(); ++j)
            max_abs_entry = std::max(max_abs_entry, std::abs(rReference(i, j)));
    return 1.0e-9 * (1.0 + max_abs_entry);
}

/// Relative Frobenius norm difference ||A-B||_F / ||B||_F.
double FrobeniusRelDiff(const Matrix& rA, const Matrix& rB)
{
    double diff_norm = 0.0, ref_norm = 0.0;
    for (std::size_t i = 0; i < rA.size1(); ++i)
        for (std::size_t j = 0; j < rA.size2(); ++j) {
            diff_norm += std::pow(rA(i, j) - rB(i, j), 2);
            ref_norm += std::pow(rB(i, j), 2);
        }
    return std::sqrt(diff_norm / ref_norm);
}

/// Independent consistent-mass matrix by explicit numerical integration with
/// the given quadrature (used as a non-production oracle).
Matrix IndependentConsistentMass(
    const Geometry<Node>& rGeometry,
    const double rDensity,
    const double rThickness,
    const GeometryData::IntegrationMethod rMethod,
    const bool rIs2d)
{
    const auto& r_integration_points = rGeometry.IntegrationPoints(rMethod);
    const Matrix& r_N = rGeometry.ShapeFunctionsValues(rMethod);
    const std::size_t number_of_nodes = rGeometry.size();
    const std::size_t dimension = rIs2d ? 2 : 3;

    Matrix M(number_of_nodes * dimension, number_of_nodes * dimension);
    noalias(M) = ZeroMatrix(number_of_nodes * dimension, number_of_nodes * dimension);

    for (std::size_t p = 0; p < r_integration_points.size(); ++p) {
        Matrix J0(dimension, dimension);
        GeometryUtils::JacobianOnInitialConfiguration(rGeometry, r_integration_points[p], J0);
        const double detJ0 = MathUtils<double>::Det(J0);
        double weight = r_integration_points[p].Weight() * detJ0;
        if (rIs2d)
            weight *= rThickness;

        for (std::size_t i = 0; i < number_of_nodes; ++i) {
            const std::size_t index_i = i * dimension;
            for (std::size_t j = 0; j < number_of_nodes; ++j) {
                const std::size_t index_j = j * dimension;
                const double value = r_N(p, i) * r_N(p, j) * weight * rDensity;
                for (std::size_t k = 0; k < dimension; ++k)
                    M(index_i + k, index_j + k) += value;
            }
        }
    }
    return M;
}

/// Registered element name pairs used in the comparison.
struct ElementNamePair
{
    std::string dam;
    std::string sma;
};

/// Prints a comparison metrics line for the report.
void PrintMassMetrics(
    const std::string& rWhat,
    const Matrix& rDam,
    const Matrix& rSma,
    const std::size_t rDimension)
{
    std::cout << "[5C.1] " << rWhat << ": max_abs_diff=" << MaxAbsDiff(rDam, rSma)
              << " frob_rel_diff=" << FrobeniusRelDiff(rDam, rSma)
              << " total_mass_dam=" << TotalDirectionalMass(rDam, rDimension)
              << " total_mass_sma=" << TotalDirectionalMass(rSma, rDimension) << std::endl;
}

} // namespace

//************************************************************************************
// 1. Integration-rule probe: which quadrature does each implementation use?
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(Mass01_IntegrationRuleProbe, KratosDamFastSuite)
{
    // For the reference geometries, report the number of Gauss points of the
    // geometry-default rule vs the elevated (geometry_default + 1) rule used by
    // the SMA consistent mass.
    struct GeoInfo { std::string name; std::string dam_name; };
    const std::vector<GeoInfo> geos = {
        {"2D3N", "SmallDisplacementSolidElement2D3N"},
        {"2D4N", "SmallDisplacementSolidElement2D4N"},
        {"3D4N", "SmallDisplacementSolidElement3D4N"},
        {"3D8N", "SmallDisplacementSolidElement3D8N"},
        {"2D6N", "SmallDisplacementSolidElement2D6N"},
        {"2D9N", "SmallDisplacementSolidElement2D9N"},
        {"3D10N", "SmallDisplacementSolidElement3D10N"},
        {"3D27N", "SmallDisplacementSolidElement3D27N"}};

    for (const auto& g : geos) {
        const auto& proto_geom = KratosComponents<Element>::Get(g.dam_name).GetGeometry();
        const auto default_method = proto_geom.GetDefaultIntegrationMethod();
        const auto elevated_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(proto_geom);
        std::cout << "[5C.1] " << g.name
                  << ": geometry_default=" << int(default_method)
                  << " (npoints=" << proto_geom.IntegrationPoints(default_method).size()
                  << ") elevated=" << int(elevated_method)
                  << " (npoints=" << proto_geom.IntegrationPoints(elevated_method).size()
                  << ")" << std::endl;
    }
    std::cout << "[6C.2] SMA consistent mass uses the elevated (geometry_default+1) rule "
              << "for every geometry (permanent M1 behavior)." << std::endl;
}

//************************************************************************************
// 2. T3 consistent mass: default, explicit-consistent, lumped.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(Mass01_T3_Consistent, KratosDamFastSuite)
{
    const ElementNamePair names{"SmallDisplacementSolidElement2D3N", "SmallDisplacementElement2D3N"};
    const std::vector<std::vector<double>> coords = {{0,0,0},{2.0,0,0},{0,1.0,0}};
    const double area = 1.0;         // right triangle
    const double physical_mass = test_density * test_thickness * area;

    Model model;
    ModelPart* p_dam_mp = nullptr;
    ModelPart* p_sma_mp = nullptr;
    Element::Pointer p_dam = CreateElement(model, "T3Dam", names.dam, coords, true, test_density, true, test_thickness, false, p_dam_mp);
    Element::Pointer p_sma = CreateElement(model, "T3Sma", names.sma, coords, true, test_density, true, test_thickness, false, p_sma_mp);
    const Matrix M_dam_default = ElementMass(*p_dam, p_dam_mp->GetProcessInfo());
    const Matrix M_sma_default = ElementMass(*p_sma, p_sma_mp->GetProcessInfo());

    // Analytical exact simplex consistent mass: diag = A*t*rho/6, off = A*t*rho/12.
    Matrix M_exact(6, 6);
    noalias(M_exact) = ZeroMatrix(6, 6);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            const double value = (i == j) ? physical_mass / 6.0 : physical_mass / 12.0;
            M_exact(i * 2, j * 2) = value;
            M_exact(i * 2 + 1, j * 2 + 1) = value;
        }

    // SMA default (3-point) == analytical exact.
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_sma_default, M_exact), 0.0, exact_tolerance);
    // Phase 6A (M1): the HISTORICAL Dam name now creates an SMA runtime element,
    // so its default consistent mass is the SMA exact simplex mass. The legacy
    // no-flag uniform 1/9 subintegrated mass is INTENTIONALLY replaced (M1).
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam_default, M_sma_default), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam_default, M_exact), 0.0, exact_tolerance);

    // Total mass conserved in both.
    KRATOS_EXPECT_NEAR(TotalDirectionalMass(M_dam_default, 2), physical_mass, exact_tolerance);
    KRATOS_EXPECT_NEAR(TotalDirectionalMass(M_sma_default, 2), physical_mass, exact_tolerance);

    // Lumped: Dam == SMA (both use geometry LumpingFactors * total mass).
    ModelPart* p_dl_mp = nullptr;
    ModelPart* p_sl_mp = nullptr;
    Element::Pointer p_dam_lump = CreateElement(model, "T3DamLump", names.dam, coords, true, test_density, true, test_thickness, true, p_dl_mp);
    Element::Pointer p_sma_lump = CreateElement(model, "T3SmaLump", names.sma, coords, true, test_density, true, test_thickness, true, p_sl_mp);
    const Matrix M_dam_lumped = ElementMass(*p_dam_lump, p_dl_mp->GetProcessInfo());
    const Matrix M_sma_lumped = ElementMass(*p_sma_lump, p_sl_mp->GetProcessInfo());
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam_lumped, M_sma_lumped), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(TotalDirectionalMass(M_dam_lumped, 2), physical_mass, exact_tolerance);
    // Lumped is diagonal.
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            if (i != j)
                KRATOS_EXPECT_NEAR(M_dam_lumped(i, j), 0.0, exact_tolerance);

    PrintMassMetrics("T3 default", M_dam_default, M_sma_default, 2);
    std::cout << "[6A][M1] T3: historical Dam name -> SMA exact consistent mass "
              << "(legacy no-flag uniform 1/9 intentionally replaced); "
              << "lumped remains compatible" << std::endl;
}

//************************************************************************************
// 3. T4 consistent mass: default, explicit-consistent, lumped.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(Mass01_T4_Consistent, KratosDamFastSuite)
{
    const ElementNamePair names{"SmallDisplacementSolidElement3D4N", "SmallDisplacementElement3D4N"};
    const std::vector<std::vector<double>> coords = {{0,0,0},{2.0,0,0},{0,1.0,0},{0,0,1.0}};
    const double volume = 1.0 / 3.0;
    const double physical_mass = test_density * volume;

    Model model;
    ModelPart* p_dam_mp = nullptr;
    ModelPart* p_sma_mp = nullptr;
    Element::Pointer p_dam = CreateElement(model, "T4Dam", names.dam, coords, false, test_density, false, 0.0, false, p_dam_mp);
    Element::Pointer p_sma = CreateElement(model, "T4Sma", names.sma, coords, false, test_density, false, 0.0, false, p_sma_mp);
    const Matrix M_dam_default = ElementMass(*p_dam, p_dam_mp->GetProcessInfo());
    const Matrix M_sma_default = ElementMass(*p_sma, p_sma_mp->GetProcessInfo());

    // Analytical exact tetra consistent mass: diag = V*rho/10, off = V*rho/20.
    Matrix M_exact(12, 12);
    noalias(M_exact) = ZeroMatrix(12, 12);
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) {
            const double value = (i == j) ? physical_mass / 10.0 : physical_mass / 20.0;
            for (std::size_t k = 0; k < 3; ++k)
                M_exact(i * 3 + k, j * 3 + k) = value;
        }
    // SMA default (elevated 4-point rule) == independent integration with the
    // SAME rule (exact match), and == the analytical simplex consistent mass up
    // to the quadrature-point truncation of the rule (the tetra 4-point rule
    // stores truncated barycentric coordinates).
    const auto& r_geom_sma = p_sma->GetGeometry();
    const auto elevated_t4 = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom_sma);
    const Matrix M_ind_t4 = IndependentConsistentMass(r_geom_sma, test_density, 0.0, elevated_t4, false);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_sma_default, M_ind_t4), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_sma_default, M_exact), 0.0, 1.0e-5);

    // Phase 6A (M1): the HISTORICAL Dam name now creates an SMA runtime element,
    // so its default consistent mass is the SMA exact simplex mass. The legacy
    // no-flag uniform 1/16 subintegrated mass is INTENTIONALLY replaced (M1).
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam_default, M_sma_default), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam_default, M_exact), 0.0, 1.0e-5);

    KRATOS_EXPECT_NEAR(TotalDirectionalMass(M_dam_default, 3), physical_mass, exact_tolerance);
    KRATOS_EXPECT_NEAR(TotalDirectionalMass(M_sma_default, 3), physical_mass, exact_tolerance);

    // Lumped Dam == SMA.
    ModelPart* p_dl_mp = nullptr;
    ModelPart* p_sl_mp = nullptr;
    Element::Pointer p_dam_lump = CreateElement(model, "T4DamLump", names.dam, coords, false, test_density, false, 0.0, true, p_dl_mp);
    Element::Pointer p_sma_lump = CreateElement(model, "T4SmaLump", names.sma, coords, false, test_density, false, 0.0, true, p_sl_mp);
    const Matrix M_dam_lumped = ElementMass(*p_dam_lump, p_dl_mp->GetProcessInfo());
    const Matrix M_sma_lumped = ElementMass(*p_sma_lump, p_sl_mp->GetProcessInfo());
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam_lumped, M_sma_lumped), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(TotalDirectionalMass(M_dam_lumped, 3), physical_mass, exact_tolerance);

    PrintMassMetrics("T4 default", M_dam_default, M_sma_default, 3);
    std::cout << "[6A][M1] T4: historical Dam name -> SMA exact consistent mass "
              << "(legacy no-flag uniform 1/16 intentionally replaced)" << std::endl;
}

//************************************************************************************
// 4. Affine Q4 and H8 consistent mass: default equivalence.
//************************************************************************************

namespace
{
void VerifyAffineConsistent(
    const std::string& rLabel,
    const ElementNamePair& rNames,
    const std::vector<std::vector<double>>& rCoords,
    const bool rIs2d,
    const double rPhysicalMass)
{
    Model model;
    ModelPart* p_dam_mp = nullptr;
    ModelPart* p_sma_mp = nullptr;
    Element::Pointer p_dam = CreateElement(model, rLabel + "Dam", rNames.dam, rCoords, rIs2d, test_density, rIs2d, test_thickness, false, p_dam_mp);
    Element::Pointer p_sma = CreateElement(model, rLabel + "Sma", rNames.sma, rCoords, rIs2d, test_density, rIs2d, test_thickness, false, p_sma_mp);
    const Matrix M_dam = ElementMass(*p_dam, p_dam_mp->GetProcessInfo());
    const Matrix M_sma = ElementMass(*p_sma, p_sma_mp->GetProcessInfo());

    const auto& r_geom = p_sma->GetGeometry();
    const auto default_method = r_geom.GetDefaultIntegrationMethod();
    const auto elevated_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom);
    const Matrix M_ind_default = IndependentConsistentMass(r_geom, test_density, test_thickness, default_method, rIs2d);
    const Matrix M_ind_elevated = IndependentConsistentMass(r_geom, test_density, test_thickness, elevated_method, rIs2d);

    // For affine geometry the two quadrature rules agree (shape-function
    // products integrated exactly by both).
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_ind_default, M_ind_elevated), 0.0, exact_tolerance);
    // Dam default and SMA default both equal the (rule-independent) reference.
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam, M_ind_default), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_sma, M_ind_elevated), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam, M_sma), 0.0, exact_tolerance);

    KRATOS_EXPECT_NEAR(TotalDirectionalMass(M_dam, rIs2d ? 2 : 3), rPhysicalMass, 1.0e-9);
    KRATOS_EXPECT_NEAR(TotalDirectionalMass(M_sma, rIs2d ? 2 : 3), rPhysicalMass, 1.0e-9);

    // Matrix structure: symmetric, zero coupling between orthogonal components.
    const std::size_t dim = rIs2d ? 2 : 3;
    for (std::size_t i = 0; i < M_dam.size1(); ++i)
        for (std::size_t j = 0; j < M_dam.size2(); ++j) {
            KRATOS_EXPECT_NEAR(M_dam(i, j), M_dam(j, i), exact_tolerance);
            if (i % dim != j % dim)
                KRATOS_EXPECT_NEAR(M_dam(i, j), 0.0, exact_tolerance);
        }

    PrintMassMetrics(rLabel + " affine", M_dam, M_sma, rIs2d ? 2 : 3);
    std::cout << "[5C.1] " << rLabel << " (affine): Dam default == SMA default == "
              << "rule-independent analytical reference" << std::endl;
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(Mass01_Q4_AffineConsistent, KratosDamFastSuite)
{
    const std::vector<std::vector<double>> coords = {{0,0,0},{2.0,0,0},{2.0,1.0,0},{0,1.0,0}};
    VerifyAffineConsistent("Q4",
        {"SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N"},
        coords, true, test_density * test_thickness * 2.0);
}

KRATOS_TEST_CASE_IN_SUITE(Mass01_H8_AffineConsistent, KratosDamFastSuite)
{
    const std::vector<std::vector<double>> coords = {
        {0,0,0},{2.0,0,0},{2.0,1.0,0},{0,1.0,0},{0,0,2.0},{2.0,0,2.0},{2.0,1.0,2.0},{0,1.0,2.0}};
    VerifyAffineConsistent("H8",
        {"SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N"},
        coords, false, test_density * 4.0);
}

//************************************************************************************
// 5. Distorted Q4 and H8: quadrature sensitivity.
//************************************************************************************

namespace
{
/// Checks that the geometry Jacobian determinant is strictly positive on the
/// reference configuration at the given quadrature points.
bool JacobianPositive(const Geometry<Node>& rGeometry, const GeometryData::IntegrationMethod rMethod)
{
    const auto& r_ips = rGeometry.IntegrationPoints(rMethod);
    const std::size_t dimension = rGeometry.WorkingSpaceDimension();
    for (const auto& r_gp : r_ips) {
        Matrix J0(dimension, dimension);
        GeometryUtils::JacobianOnInitialConfiguration(rGeometry, r_gp, J0);
        if (MathUtils<double>::Det(J0) <= 0.0)
            return false;
    }
    return true;
}

void VerifyDistortedSensitivity(
    const std::string& rLabel,
    const ElementNamePair& rNames,
    const std::vector<std::vector<double>>& rCoords,
    const bool rIs2d)
{
    Model model;
    ModelPart* p_dam_mp = nullptr;
    ModelPart* p_sma_mp = nullptr;
    Element::Pointer p_dam = CreateElement(model, rLabel + "Dam", rNames.dam, rCoords, rIs2d, test_density, rIs2d, test_thickness, false, p_dam_mp);
    Element::Pointer p_sma = CreateElement(model, rLabel + "Sma", rNames.sma, rCoords, rIs2d, test_density, rIs2d, test_thickness, false, p_sma_mp);

    // Non-degenerate check on the actual (distorted) geometry.
    const auto& r_geom = p_sma->GetGeometry();
    const auto elevated = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom);
    KRATOS_EXPECT_TRUE(JacobianPositive(r_geom, elevated));

    const Matrix M_dam = ElementMass(*p_dam, p_dam_mp->GetProcessInfo());
    const Matrix M_sma = ElementMass(*p_sma, p_sma_mp->GetProcessInfo());

    const auto default_method = r_geom.GetDefaultIntegrationMethod();
    const Matrix M_ind_default = IndependentConsistentMass(r_geom, test_density, test_thickness, default_method, rIs2d);
    const Matrix M_ind_elevated = IndependentConsistentMass(r_geom, test_density, test_thickness, elevated, rIs2d);

    // Both the historical Dam alias and the direct SMA element use the SMA
    // runtime: their masses coincide with the independent references.
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam, M_ind_default), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_sma, M_ind_elevated), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam, M_sma), 0.0, exact_tolerance);

    // For the quad/hex families the bilinear/trilinear Jacobian makes the mass
    // integrand of total degree <= 3 per direction, which the DEFAULT rule
    // already integrates exactly. Hence Dam default == SMA default even on this
    // distorted geometry: the elevated quadrature is a numerical no-op here.
    const double max_abs = MaxAbsDiff(M_dam, M_sma);
    const double frob_rel = FrobeniusRelDiff(M_dam, M_sma);
    const double total_diff = std::abs(TotalDirectionalMass(M_dam, rIs2d ? 2 : 3) -
                                       TotalDirectionalMass(M_sma, rIs2d ? 2 : 3));
    KRATOS_EXPECT_NEAR(max_abs, 0.0, exact_tolerance);

    std::cout << "[6C.2] " << rLabel << " distorted: historical-alias vs SMA "
              << "max_abs_diff=" << max_abs
              << " frob_rel_diff=" << frob_rel
              << " total_mass_diff=" << total_diff
              << " -> default quadrature already exact for this family" << std::endl;
    PrintMassMetrics(rLabel + " distorted default", M_dam, M_sma, rIs2d ? 2 : 3);
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(Mass01_Q4_Distorted, KratosDamFastSuite)
{
    // Non-affine quadrilateral: opposite edges are not parallel.
    const std::vector<std::vector<double>> coords = {{0,0,0},{2.0,0,0},{2.3,1.4,0},{0.2,1.1,0}};
    VerifyDistortedSensitivity("Q4",
        {"SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N"},
        coords, true);
}

KRATOS_TEST_CASE_IN_SUITE(Mass01_H8_Distorted, KratosDamFastSuite)
{
    // Non-affine hexahedron: one corner displaced out of the affine position.
    const std::vector<std::vector<double>> coords = {
        {0,0,0},{2.0,0,0},{2.0,1.0,0},{0,1.0,0},
        {0,0,2.0},{2.0,0,2.0},{2.4,1.2,2.2},{0,1.0,2.0}};
    VerifyDistortedSensitivity("H8",
        {"SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N"},
        coords, false);
}

//************************************************************************************
// 6. Thickness behavior (2D): linear scaling, and both implementations.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(Mass01_Thickness, KratosDamFastSuite)
{
    const ElementNamePair names{"SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N"};
    const std::vector<std::vector<double>> coords = {{0,0,0},{2.0,0,0},{2.0,1.0,0},{0,1.0,0}};
    const double t1 = 0.1, t2 = 0.3;

    Model model;
    ModelPart* p_mp1 = nullptr; ModelPart* p_mp2 = nullptr; ModelPart* p_mp3 = nullptr;
    ModelPart* p_mp4 = nullptr; ModelPart* p_mp5 = nullptr; ModelPart* p_mp6 = nullptr;
    Element::Pointer p_dam_t1 = CreateElement(model, "ThDamT1", names.dam, coords, true, test_density, true, t1, false, p_mp1);
    Element::Pointer p_dam_t2 = CreateElement(model, "ThDamT2", names.dam, coords, true, test_density, true, t2, false, p_mp2);
    Element::Pointer p_sma_t1 = CreateElement(model, "ThSmaT1", names.sma, coords, true, test_density, true, t1, false, p_mp3);
    Element::Pointer p_sma_t2 = CreateElement(model, "ThSmaT2", names.sma, coords, true, test_density, true, t2, false, p_mp4);
    const Matrix M_dam_1 = ElementMass(*p_dam_t1, p_mp1->GetProcessInfo());
    const Matrix M_dam_2 = ElementMass(*p_dam_t2, p_mp2->GetProcessInfo());
    const Matrix M_sma_1 = ElementMass(*p_sma_t1, p_mp3->GetProcessInfo());
    const Matrix M_sma_2 = ElementMass(*p_sma_t2, p_mp4->GetProcessInfo());

    // Mass scales linearly with thickness, identically for Dam and SMA.
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam_2, (t2 / t1) * M_dam_1), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_sma_2, (t2 / t1) * M_sma_1), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam_1, M_sma_1), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam_2, M_sma_2), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(TotalDirectionalMass(M_dam_2, 2), test_density * t2 * 2.0, 1.0e-9);

    // Absent THICKNESS: both behave as unit thickness.
    Element::Pointer p_dam_no = CreateElement(model, "ThDamNo", names.dam, coords, true, test_density, false, 0.0, false, p_mp5);
    Element::Pointer p_sma_no = CreateElement(model, "ThSmaNo", names.sma, coords, true, test_density, false, 0.0, false, p_mp6);
    const Matrix M_dam_no = ElementMass(*p_dam_no, p_mp5->GetProcessInfo());
    const Matrix M_sma_no = ElementMass(*p_sma_no, p_mp6->GetProcessInfo());
    KRATOS_EXPECT_NEAR(TotalDirectionalMass(M_dam_no, 2), test_density * 2.0, 1.0e-9);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam_no, M_sma_no), 0.0, exact_tolerance);

    std::cout << "[5C.1] thickness: linear scaling t2/t1=" << (t2 / t1)
              << " identical for Dam and SMA; absent THICKNESS -> unit thickness on both" << std::endl;
}

//************************************************************************************
// 7. Density behavior.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(Mass01_Density, KratosDamFastSuite)
{
    const ElementNamePair names{"SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N"};
    const std::vector<std::vector<double>> coords = {
        {0,0,0},{2.0,0,0},{2.0,1.0,0},{0,1.0,0},{0,0,2.0},{2.0,0,2.0},{2.0,1.0,2.0},{0,1.0,2.0}};

    const double rho_nonunit = 3141.5926;
    Model model;
    ModelPart* p_mp1 = nullptr; ModelPart* p_mp2 = nullptr; ModelPart* p_mp3 = nullptr; ModelPart* p_mp4 = nullptr;
    Element::Pointer p_dam = CreateElement(model, "DensDam", names.dam, coords, false, rho_nonunit, false, 0.0, false, p_mp1);
    Element::Pointer p_sma = CreateElement(model, "DensSma", names.sma, coords, false, rho_nonunit, false, 0.0, false, p_mp2);
    const Matrix M_dam = ElementMass(*p_dam, p_mp1->GetProcessInfo());
    const Matrix M_sma = ElementMass(*p_sma, p_mp2->GetProcessInfo());
    KRATOS_EXPECT_NEAR(TotalDirectionalMass(M_dam, 3), rho_nonunit * 4.0, 1.0e-9);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam, M_sma), 0.0, exact_tolerance);

    // Zero density: zero mass matrix on both.
    Element::Pointer p_dam_z = CreateElement(model, "DensDamZ", names.dam, coords, false, 0.0, false, 0.0, false, p_mp3);
    Element::Pointer p_sma_z = CreateElement(model, "DensSmaZ", names.sma, coords, false, 0.0, false, 0.0, false, p_mp4);
    const Matrix M_dam_z = ElementMass(*p_dam_z, p_mp3->GetProcessInfo());
    const Matrix M_sma_z = ElementMass(*p_sma_z, p_mp4->GetProcessInfo());
    const Matrix zero = ZeroMatrix(M_dam_z.size1(), M_dam_z.size2());
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam_z, zero), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_sma_z, zero), 0.0, exact_tolerance);

    std::cout << "[5C.1] density: non-unit and zero density behave identically on Dam and SMA "
              << "(both read Properties[DENSITY]; missing DENSITY is an error on both)" << std::endl;
}

//************************************************************************************
// 9. Higher-order geometry contract: Dam+flag == SMA default, for every direct
//    counterpart geometry.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(Mass01_HigherOrder_FlagContract, KratosDamFastSuite)
{
    const std::vector<std::pair<ElementNamePair, bool>> cases = {
        {{"SmallDisplacementSolidElement2D3N", "SmallDisplacementElement2D3N"}, true},
        {{"SmallDisplacementSolidElement2D4N", "SmallDisplacementElement2D4N"}, true},
        {{"SmallDisplacementSolidElement2D6N", "SmallDisplacementElement2D6N"}, true},
        {{"SmallDisplacementSolidElement2D8N", "SmallDisplacementElement2D8N"}, true},
        {{"SmallDisplacementSolidElement2D9N", "SmallDisplacementElement2D9N"}, true},
        {{"SmallDisplacementSolidElement3D4N", "SmallDisplacementElement3D4N"}, false},
        {{"SmallDisplacementSolidElement3D6N", "SmallDisplacementElement3D6N"}, false},
        {{"SmallDisplacementSolidElement3D8N", "SmallDisplacementElement3D8N"}, false},
        {{"SmallDisplacementSolidElement3D10N", "SmallDisplacementElement3D10N"}, false},
        {{"SmallDisplacementSolidElement3D15N", "SmallDisplacementElement3D15N"}, false},
        {{"SmallDisplacementSolidElement3D20N", "SmallDisplacementElement3D20N"}, false},
        {{"SmallDisplacementSolidElement3D27N", "SmallDisplacementElement3D27N"}, false},
    };

    for (const auto& entry : cases) {
        const ElementNamePair& names = entry.first;
        const bool is2d = entry.second;

        // Build coordinates from the registered Dam prototype geometry (scaled).
        const Element& r_proto = KratosComponents<Element>::Get(names.dam);
        Matrix local_coordinates;
        r_proto.GetGeometry().PointsLocalCoordinates(local_coordinates);
        const double scale = 2.0;
        const double offset_x = 0.5, offset_y = 1.0, offset_z = 0.25;
        std::vector<std::vector<double>> coords(local_coordinates.size1(), std::vector<double>(3, 0.0));
        for (std::size_t i = 0; i < coords.size(); ++i) {
            coords[i][0] = scale * local_coordinates(i, 0) + offset_x;
            coords[i][1] = scale * local_coordinates(i, 1) + offset_y;
            if (local_coordinates.size2() > 2)
                coords[i][2] = scale * local_coordinates(i, 2) + offset_z;
        }

        Model model;
        ModelPart* p_dam_mp = nullptr;
        ModelPart* p_sma_mp = nullptr;
        Element::Pointer p_dam = CreateElement(model, "HO" + names.dam, names.dam, coords, is2d, test_density, is2d, test_thickness, false, p_dam_mp);
        Element::Pointer p_sma = CreateElement(model, "HO" + names.sma, names.sma, coords, is2d, test_density, is2d, test_thickness, false, p_sma_mp);
        const Matrix M_dam = ElementMass(*p_dam, p_dam_mp->GetProcessInfo());
        const Matrix M_sma = ElementMass(*p_sma, p_sma_mp->GetProcessInfo());

        // The historical Dam alias and the direct SMA name produce the SAME
        // (SMA) mass matrix for every higher-order geometry.
        KRATOS_EXPECT_NEAR(MaxAbsDiff(M_dam, M_sma), 0.0, MassTolerance(M_sma));

        const double max_abs = MaxAbsDiff(M_dam, M_sma);
        const double frob_rel = FrobeniusRelDiff(M_dam, M_sma);

        std::cout << "[6C.2] " << names.dam.substr(28)
                  << " (" << (is2d ? "2D" : "3D") << ")"
                  << " npoints=" << M_sma.size1()
                  << ": historical-alias vs SMA max_abs_diff=" << max_abs
                  << " frob_rel_diff=" << frob_rel << std::endl;
    }
    std::cout << "[6C.2] Higher-order contract: every historical Dam alias yields the "
              << "same SMA mass matrix as the direct SMA name." << std::endl;
}

//************************************************************************************
// 10. Rayleigh damping (alpha*M + beta*K) and inertial force (M*a, C*v).
//************************************************************************************

namespace
{
void VerifyDampingAndInertia(const std::string& rLabel, const ElementNamePair& rNames,
                             const std::vector<std::vector<double>>& rCoords, const bool rIs2d)
{
    Model model;
    ModelPart* p_dam_mp = nullptr;
    ModelPart* p_sma_mp = nullptr;
    Element::Pointer p_dam = CreateElement(model, rLabel + "Dam", rNames.dam, rCoords, rIs2d, test_density, rIs2d, test_thickness, false, p_dam_mp);
    Element::Pointer p_sma = CreateElement(model, rLabel + "Sma", rNames.sma, rCoords, rIs2d, test_density, rIs2d, test_thickness, false, p_sma_mp);
    ProcessInfo& r_dam_pi = p_dam_mp->GetProcessInfo();
    ProcessInfo& r_sma_pi = p_sma_mp->GetProcessInfo();

    const double alpha = 0.02, beta = 0.03;
    for (auto& p_prop : std::vector<Properties::Pointer>{p_dam->pGetProperties(), p_sma->pGetProperties()}) {
        (*p_prop)[RAYLEIGH_ALPHA] = alpha;
        (*p_prop)[RAYLEIGH_BETA] = beta;
    }

    Matrix C_dam, C_sma, K_dam, K_sma, M_dam, M_sma;
    p_dam->CalculateDampingMatrix(C_dam, r_dam_pi);
    p_sma->CalculateDampingMatrix(C_sma, r_sma_pi);
    p_dam->CalculateLeftHandSide(K_dam, r_dam_pi);
    p_sma->CalculateLeftHandSide(K_sma, r_sma_pi);
    p_dam->CalculateMassMatrix(M_dam, r_dam_pi);
    p_sma->CalculateMassMatrix(M_sma, r_sma_pi);

    // Stiffness identical (validated in earlier phases).
    KRATOS_EXPECT_NEAR(MaxAbsDiff(K_dam, K_sma), 0.0, exact_tolerance);

    // Damping is alpha*M + beta*K; under M1 the historical alias and the direct
    // SMA path share the same mass, so the damping matrices coincide.
    const double damping_diff = MaxAbsDiff(C_dam, C_sma);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(C_dam, alpha * M_dam + beta * K_dam), 0.0, 1.0e-9);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(C_sma, alpha * M_sma + beta * K_sma), 0.0, 1.0e-9);
    KRATOS_EXPECT_NEAR(damping_diff, 0.0, 1.0e-9);

    // Inertial force M*a with a prescribed, node-varying acceleration vector.
    Vector a(M_dam.size1());
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] = 0.5 + 0.01 * i;
    const Vector f_dam = prod(M_dam, a);
    const Vector f_sma = prod(M_sma, a);
    // M1: the historical Dam name creates the SMA element, so the inertial force
    // is identical to the direct SMA path.
    KRATOS_EXPECT_NEAR(norm_2(f_dam - f_sma), 0.0, 1.0e-9);

    // Damping force C*v.
    const Vector v = a;
    const Vector cv_dam = prod(C_dam, v);
    const Vector cv_sma = prod(C_sma, v);
    KRATOS_EXPECT_NEAR(norm_2(cv_dam - cv_sma), 0.0, 1.0e-9);

    std::cout << "[6C.2] " << rLabel << " damping: historical-alias vs SMA "
              << "damping_diff=" << damping_diff
              << " |M*a|=" << norm_2(f_dam) << std::endl;
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(Mass01_DampingAndInertia_T3, KratosDamFastSuite)
{
    // T3 (simplex): the default-mode mass genuinely differs, so the dynamic
    // quantities differ unless Dam uses the consistent-mass flag.
    const std::vector<std::vector<double>> coords = {{0,0,0},{2.0,0,0},{0,1.0,0}};
    VerifyDampingAndInertia("T3", {"SmallDisplacementSolidElement2D3N", "SmallDisplacementElement2D3N"}, coords, true);
}

KRATOS_TEST_CASE_IN_SUITE(Mass01_DampingAndInertia_T4, KratosDamFastSuite)
{
    const std::vector<std::vector<double>> coords = {{0,0,0},{2.0,0,0},{0,1.0,0},{0,0,1.0}};
    VerifyDampingAndInertia("T4", {"SmallDisplacementSolidElement3D4N", "SmallDisplacementElement3D4N"}, coords, false);
}

//************************************************************************************
// 11. Minimal dynamic system: characteristic frequency via Rayleigh quotient.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(Mass01_RayleighQuotient, KratosDamFastSuite)
{
    // Single tetra (simplex: default mass differs). The stiffness K is shared;
    // only the mass matrix varies. omega^2 = phi^T K phi / phi^T M phi.
    const std::vector<std::vector<double>> coords = {{0,0,0},{2.0,0,0},{0,1.0,0},{0,0,1.0}};
    Model model;
    ModelPart* p_dam_mp = nullptr;
    ModelPart* p_sma_mp = nullptr;
    Element::Pointer p_dam = CreateElement(model, "RQDam", "SmallDisplacementSolidElement3D4N", coords, false, test_density, false, 0.0, false, p_dam_mp);
    Element::Pointer p_sma = CreateElement(model, "RQSma", "SmallDisplacementElement3D4N", coords, false, test_density, false, 0.0, false, p_sma_mp);

    Matrix K, M_dam, M_sma;
    p_sma->CalculateLeftHandSide(K, p_sma_mp->GetProcessInfo());
    p_dam->CalculateMassMatrix(M_dam, p_dam_mp->GetProcessInfo());
    p_sma->CalculateMassMatrix(M_sma, p_sma_mp->GetProcessInfo());

    // Mode: affine expansion phi = (x, y, z) per node (a non-rigid mode).
    Vector phi(M_dam.size1());
    for (std::size_t n = 0; n < 4; ++n) {
        phi[n * 3 + 0] = coords[n][0];
        phi[n * 3 + 1] = coords[n][1];
        phi[n * 3 + 2] = coords[n][2];
    }

    auto rayleigh = [&](const Matrix& rM) {
        const double num = inner_prod(phi, prod(K, phi));
        const double den = inner_prod(phi, prod(rM, phi));
        return std::sqrt(num / den);  // omega
    };

    const double omega_dam = rayleigh(M_dam);
    const double omega_sma = rayleigh(M_sma);
    // M1: the historical Dam name creates the SMA element, so the characteristic
    // frequency equals the SMA one (the legacy no-flag simplex mass is gone).
    KRATOS_EXPECT_NEAR(omega_dam, omega_sma, 1.0e-9);

    std::cout << "[6C.2][M1] Rayleigh quotient (T4): historical-name omega=" << omega_dam
              << " == SMA omega=" << omega_sma << std::endl;
}

//************************************************************************************
// 12. Thermo-mechanical element mass neutrality.
//************************************************************************************

namespace
{
std::string ThermoCounterpart(const std::string& rSolidName)
{
    // "SmallDisplacementSolidElement2D4N" -> "SmallDisplacementThermoMechanicElement2D4N"
    const std::string marker = "Solid";
    const std::size_t pos = rSolidName.find(marker);
    return rSolidName.substr(0, pos) + "ThermoMechanic" + rSolidName.substr(pos + marker.size());
}
} // namespace

KRATOS_TEST_CASE_IN_SUITE(Mass01_ThermoMechanicalNeutrality, KratosDamFastSuite)
{
    // Dam SmallDisplacementSolid vs Dam SmallDisplacementThermoMechanic on
    // identical geometry must produce identical mass matrices.
    const std::vector<std::pair<std::string, std::vector<std::vector<double>>>> cases = {
        {"SmallDisplacementSolidElement2D4N", {{0,0,0},{2.0,0,0},{2.0,1.0,0},{0,1.0,0}}},
        {"SmallDisplacementSolidElement3D8N", {{0,0,0},{2.0,0,0},{2.0,1.0,0},{0,1.0,0},{0,0,2.0},{2.0,0,2.0},{2.0,1.0,2.0},{0,1.0,2.0}}}};
    for (const auto& c : cases) {
        const std::string& solid_name = c.first;
        const bool is2d = (c.second.size() <= 4);
        const std::string thermo_name = ThermoCounterpart(solid_name);
        Model model;
        ModelPart* p_mp1 = nullptr;
        ModelPart* p_mp2 = nullptr;
        Element::Pointer p_mech = CreateElement(model, "NeutMech" + solid_name, solid_name, c.second, is2d, test_density, is2d, test_thickness, false, p_mp1);
        Element::Pointer p_thermo = CreateElement(model, "NeutThermo" + thermo_name, thermo_name, c.second, is2d, test_density, is2d, test_thickness, false, p_mp2);
        const Matrix M_mech = ElementMass(*p_mech, p_mp1->GetProcessInfo());
        const Matrix M_thermo = ElementMass(*p_thermo, p_mp2->GetProcessInfo());
        KRATOS_EXPECT_EQ(M_mech.size1(), M_thermo.size1());
        KRATOS_EXPECT_NEAR(MaxAbsDiff(M_mech, M_thermo), 0.0, exact_tolerance);
        std::cout << "[5C.1] " << solid_name << " vs " << thermo_name
                  << ": mass identical (thermo element is mass-neutral)" << std::endl;
    }
}

//************************************************************************************
// 13. Constitutive-law independence.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(Mass01_ConstitutiveIndependence, KratosDamFastSuite)
{
    // Mass must not depend on the constitutive law: compare the H8 mass matrix
    // computed with linear, local-damage and nonlocal-damage laws. The mass path
    // never invokes the constitutive law, so the matrices must be identical.
    const std::vector<std::vector<double>> coords = {
        {0,0,0},{2.0,0,0},{2.0,1.0,0},{0,1.0,0},{0,0,2.0},{2.0,0,2.0},{2.0,1.0,2.0},{0,1.0,2.0}};

    auto make_with_law = [&](const std::string& rLawName, Model& rModel, ModelPart*& rOutMp) {
        ModelPart& r_mp = rModel.CreateModelPart("CLI" + rLawName, 2);
        r_mp.GetProcessInfo()[DOMAIN_SIZE] = 3;
        r_mp.GetProcessInfo()[SPACE_DIMENSION] = 3;
        r_mp.GetProcessInfo()[IS_RESTARTED] = false;
        r_mp.GetProcessInfo()[IS_CONVERGED] = true;
        r_mp.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_mp.AddNodalSolutionStepVariable(VELOCITY);
        r_mp.AddNodalSolutionStepVariable(ACCELERATION);
        r_mp.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
        for (std::size_t i = 0; i < coords.size(); ++i) {
            Node::Pointer p_node = r_mp.CreateNewNode(i + 1, coords[i][0], coords[i][1], coords[i][2]);
            p_node->AddDof(DISPLACEMENT_X); p_node->AddDof(DISPLACEMENT_Y); p_node->AddDof(DISPLACEMENT_Z);
        }
        auto p_props = r_mp.CreateNewProperties(1);
        (*p_props)[YOUNG_MODULUS] = test_young_modulus;
        (*p_props)[POISSON_RATIO] = test_poisson_ratio;
        (*p_props)[DENSITY] = test_density;
        (*p_props)[DAMAGE_THRESHOLD] = 5.0e-3;
        (*p_props)[STRENGTH_RATIO] = 10.0;
        (*p_props)[FRACTURE_ENERGY] = 5000.0;
        p_props->SetValue(CONSTITUTIVE_LAW, CreateLaw(rLawName));
        Geometry<Node>::PointsArrayType pts;
        for (std::size_t i = 0; i < coords.size(); ++i) pts.push_back(r_mp.pGetNode(i + 1));
        Element::Pointer p_elem = KratosComponents<Element>::Get("SmallDisplacementElement3D8N").Create(1, pts, p_props);
        p_elem->Initialize(r_mp.GetProcessInfo());
        rOutMp = &r_mp;
        return p_elem;
    };

    Model model;
    ModelPart* p_mp1 = nullptr; ModelPart* p_mp2 = nullptr; ModelPart* p_mp3 = nullptr;
    Element::Pointer p_linear = make_with_law("ThermalLinearElastic3DLaw", model, p_mp1);
    Element::Pointer p_local = make_with_law("ThermalSimoJuLocalDamage3DLaw", model, p_mp2);
    Element::Pointer p_nonlocal = make_with_law("ThermalSimoJuNonlocalDamage3DLaw", model, p_mp3);

    const Matrix M_linear = ElementMass(*p_linear, p_mp1->GetProcessInfo());
    const Matrix M_local = ElementMass(*p_local, p_mp2->GetProcessInfo());
    const Matrix M_nonlocal = ElementMass(*p_nonlocal, p_mp3->GetProcessInfo());

    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_linear, M_local), 0.0, exact_tolerance);
    KRATOS_EXPECT_NEAR(MaxAbsDiff(M_linear, M_nonlocal), 0.0, exact_tolerance);
    std::cout << "[5C.1] constitutive-law independence: mass identical for linear / "
              << "local-damage / nonlocal-damage laws (mass never invokes the law)" << std::endl;
}

//************************************************************************************
// 14. Repeat/thread determinism.
//************************************************************************************

KRATOS_TEST_CASE_IN_SUITE(Mass01_Determinism, KratosDamFastSuite)
{
    const std::vector<std::vector<double>> coords = {
        {0,0,0},{2.0,0,0},{2.0,1.0,0},{0,1.0,0},{0,0,2.0},{2.0,0,2.0},{2.4,1.2,2.2},{0,1.0,2.0}};
    Model model;
    ModelPart* p_mp1 = nullptr; ModelPart* p_mp2 = nullptr;
    Element::Pointer p_dam = CreateElement(model, "DetDam", "SmallDisplacementSolidElement3D8N", coords, false, test_density, false, 0.0, false, p_mp1);
    Element::Pointer p_sma = CreateElement(model, "DetSma", "SmallDisplacementElement3D8N", coords, false, test_density, false, 0.0, false, p_mp2);

    const Matrix M_dam_first = ElementMass(*p_dam, p_mp1->GetProcessInfo());
    const Matrix M_sma_first = ElementMass(*p_sma, p_mp2->GetProcessInfo());
    for (std::size_t repeat = 0; repeat < 5; ++repeat) {
        KRATOS_EXPECT_NEAR(MaxAbsDiff(ElementMass(*p_dam, p_mp1->GetProcessInfo()), M_dam_first), 0.0, exact_tolerance);
        KRATOS_EXPECT_NEAR(MaxAbsDiff(ElementMass(*p_sma, p_mp2->GetProcessInfo()), M_sma_first), 0.0, exact_tolerance);
    }
    std::cout << "[5C.1] determinism: 5 repeated mass computations identical on Dam and SMA" << std::endl;
}

} // namespace Testing
} // namespace Kratos

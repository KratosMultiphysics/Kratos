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

// Permanent M1 mass-policy regression: the historical Dam element names resolve to
// the StructuralMechanics SmallDisplacement element, whose exact/elevated consistent
// mass is the adopted Dam mass policy (M1). These tests lock the analytical simplex
// (T3/T4) exact consistent mass and the unchanged non-simplex (affine Q4) behavior
// through the historical registrations. No production code is modified.

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



//************************************************************************************
// 5. Distorted Q4 and H8: quadrature sensitivity.
//************************************************************************************

namespace
{
/// Checks that the geometry Jacobian determinant is strictly positive on the



} // namespace



//************************************************************************************
// 6. Thickness behavior (2D): linear scaling, and both implementations.
//************************************************************************************


//************************************************************************************
// 7. Density behavior.
//************************************************************************************


//************************************************************************************
// 9. Higher-order geometry contract: Dam+flag == SMA default, for every direct
//    counterpart geometry.
//************************************************************************************


//************************************************************************************
// 10. Rayleigh damping (alpha*M + beta*K) and inertial force (M*a, C*v).
//************************************************************************************

namespace
{

} // namespace



//************************************************************************************
// 11. Minimal dynamic system: characteristic frequency via Rayleigh quotient.
//************************************************************************************


//************************************************************************************
// 12. Thermo-mechanical element mass neutrality.
//************************************************************************************

namespace
{

} // namespace


//************************************************************************************
// 13. Constitutive-law independence.
//************************************************************************************


//************************************************************************************
// 14. Repeat/thread determinism.
//************************************************************************************


} // namespace Testing
} // namespace Kratos

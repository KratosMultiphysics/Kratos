//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Gorgi 

#include "containers/model.h"
#include "testing/testing.h"
#include "test_creation_utility.h"
#include "custom_elements/solid_iga_element.h"
#include "iga_structural_mechanics_fast_suite.h"
#include "iga_fast_suite.h"

namespace Kratos::Testing
{

namespace
{
    typedef std::size_t SizeType;

    typename Element::Pointer GetSolidIGAElement(
        ModelPart& rModelPart, SizeType PolynomialDegree, IntegrationPoint<3> IntegrationPoint)
    {
        auto p_elem_prop = rModelPart.CreateNewProperties(0);
        p_elem_prop->SetValue(THICKNESS, 1.0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 100);
        p_elem_prop->SetValue(POISSON_RATIO, 0.3);

        // Usa il KratosComponents per istanziare la legge costitutiva
        const auto& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStrain2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        auto p_quadrature_point = TestCreationUtility::GetQuadraturePointGeometry(
            rModelPart, PolynomialDegree, IntegrationPoint);

        // body force
        Vector bf(3, 0.0);
        bf[0] = 30.0;
        bf[1] = 0.0;
        p_quadrature_point->SetValue(BODY_FORCE, bf);

        return Kratos::make_intrusive<SolidIGAElement>(1, p_quadrature_point, p_elem_prop);
    }

}

// Full test for SolidIGAElement with p=3
KRATOS_TEST_CASE_IN_SUITE(SolidIGAElementP3, KratosIgaSMFastSuite)
{

    Model model;
    auto &r_model_part = model.CreateModelPart("ModelPart");
    r_model_part.SetBufferSize(2);
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    // Define integration point
    IntegrationPoint<3> integration_point(0.0694318442029737, 0.211324865405187, 0.0, 0.086963711284364);

    // Create the solid element
    auto p_solid_element = GetSolidIGAElement(r_model_part, 3, integration_point);

    // Add DOFs
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
    }

    // Initialize constitutive law
    p_solid_element->Initialize(r_process_info);
    p_solid_element->InitializeSolutionStep(r_process_info);

    // Calculate LHS and RHS
    Matrix lhs;
    Vector rhs;
    p_solid_element->CalculateLocalSystem(lhs, rhs, r_process_info);

    const double tolerance = 1.0e-4;

    // Expected LHS and RHS (inizialmente mettiamo 0, poi da aggiornare)
    const std::vector<std::vector<double>> expected_LHS = { 
        {26.634, 13.8059, 0.680681, -2.84413, -0.34324, -0.654986, -0.0183364, -0.0273023, -20.4029, -6.80383, -5.98196, -3.11307, -0.551908, -0.350916, -0.0163522, -0.0116783},
        {13.8059, 77.423, -5.81133, 15.8213, -1.09777, 1.06788, -0.0438207, 0.0237591, -3.30279, -75.6427, -3.12446, -17.3359, -0.411087, -1.32364, -0.0146501, -0.0336702},
        {0.680681, -5.81133, 4.64531, -2.62913, 0.681826, -0.295276, 0.0252949, -0.00980869, -5.98196, 7.37865, -0.135104, 1.29569, 0.0797439, 0.0701177, 0.00421729, 0.0010834},
        {-2.84413, 15.8213, -2.62913, 4.82509, -0.344831, 0.455791, -0.0122736, 0.013718, 5.19511, -17.3359, 0.628968, -3.53646, 0.00709391, -0.238199, -0.000814292, -0.00528591},
        {-0.34324, -1.09777, 0.681826, -0.344831, 0.107478, -0.0331235, 0.00408086, -0.00100772, -0.551908, 1.21471, 0.0797439, 0.24534, 0.0211172, 0.0163239, 0.000902422, 0.000356709},
        {-0.654986, 1.06788, -0.295276, 0.455791, -0.0331235, 0.0501805, -0.00109968, 0.00165026, 0.830399, -1.32364, 0.146039, -0.238199, 0.00792415, -0.0134391, 0.00012316, -0.000226463},
        {-0.0183364, -0.0438207, 0.0252949, -0.0122736, 0.00408086, -0.00109968, 0.00015605, -3.1924e-05, -0.0163522, 0.0467291, 0.00421729, 0.00979922, 0.000902422, 0.000681862, 3.7062e-05, 1.57328e-05},
        {-0.0273023, 0.0237591, -0.00980869, 0.013718, -0.00100772, 0.00165026, -3.1924e-05, 5.66304e-05, 0.0316649, -0.0336702, 0.00609706, -0.00528591, 0.000380995, -0.000226463, 7.6372e-06, -1.45571e-06},
        {-20.4029, -3.30279, -5.98196, 5.19511, -0.551908, 0.830399, -0.0163522, 0.0316649, 22.0725, -3.69928, 4.56149, 0.762082, 0.312053, 0.175503, 0.00705742, 0.00731562},
        {-6.80383, -75.6427, 7.37865, -17.3359, 1.21471, -1.32364, 0.0467291, -0.0336702, -3.69928, 76.1198, 1.55714, 16.9301, 0.294145, 1.25511, 0.0117417, 0.0310145},
        {-5.98196, -3.12446, -0.135104, 0.628968, 0.0797439, 0.146039, 0.00421729, 0.00609706, 4.56149, 1.55714, 1.34361, 0.704473, 0.124318, 0.0791189, 0.00369049, 0.00262823},
        {-3.11307, -17.3359, 1.29569, -3.53646, 0.24534, -0.238199, 0.00979922, -0.00528591, 0.762082, 16.9301, 0.704473, 3.88174, 0.0923972, 0.296503, 0.00328871, 0.00754528},
        {-0.551908, -0.411087, 0.0797439, 0.00709391, 0.0211172, 0.00792415, 0.000902422, 0.000380995, 0.312053, 0.294145, 0.124318, 0.0923972, 0.0133397, 0.00887541, 0.000432846, 0.000270018},
        {-0.350916, -1.32364, 0.0701177, -0.238199, 0.0163239, -0.0134391, 0.000681862, -0.000226463, 0.175503, 1.25511, 0.0791189, 0.296503, 0.00887541, 0.0232839, 0.000294658, 0.000607968},
        {-0.0163522, -0.0146501, 0.00421729, -0.000814292, 0.000902422, 0.00012316, 3.7062e-05, 7.6372e-06, 0.00705742, 0.0117417, 0.00369049, 0.00328871, 0.000432846, 0.000294658, 1.46821e-05, 8.55402e-06},
        {-0.0116783, -0.0336702, 0.0010834, -0.00528591, 0.000356709, -0.000226463, 1.57328e-05, -1.45571e-06, 0.00731562, 0.0310145, 0.00262823, 0.00754528, 0.000270018, 0.000607968, 8.55402e-06, 1.62397e-05}
    };

    const std::vector<double> expected_RHS = {
        0.165807,
        0,
        0.0371137,
        0,
        0.00276914,
        0,
        6.88706e-05,
        0,
        0.0444278,
        0,
        0.00994458,
        0,
        0.000741988,
        0,
        1.84538e-05,
        0
    };

    for (SizeType i = 0; i < lhs.size1(); ++i) {
        for (SizeType j = 0; j < lhs.size2(); ++j) {
            KRATOS_EXPECT_NEAR(lhs(i,j), expected_LHS[i][j], tolerance);
        }
    }

    for (SizeType i = 0; i < rhs.size(); ++i) {
        KRATOS_EXPECT_NEAR(rhs(i), expected_RHS[i], tolerance);
    }

}

} // namespace Kratos::Testing

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nunez and Inigo Lopez
//
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "compressible_potential_flow_application_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "custom_elements/embedded_transonic_perturbation_potential_flow_element.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "tests/cpp_tests/test_utilities.h"
#include "processes/find_global_nodal_entity_neighbours_process.h"

namespace Kratos {
namespace Testing {

typedef ModelPart::IndexType IndexType;
typedef ModelPart::NodeIterator NodeIteratorType;

void GenerateEmbeddedTransonicPerturbationElement(ModelPart& rModelPart) {
    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
    rModelPart.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);
    rModelPart.AddNodalSolutionStepVariable(GEOMETRY_DISTANCE);

    // Set the element properties
    Properties::Pointer pElemProp = rModelPart.CreateNewProperties(0);
    rModelPart.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.225;
    rModelPart.GetProcessInfo()[FREE_STREAM_MACH] = 0.6;
    rModelPart.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    rModelPart.GetProcessInfo()[SOUND_VELOCITY] = 340.3;
    rModelPart.GetProcessInfo()[MACH_LIMIT] = 1.73205080756887729;
    rModelPart.GetProcessInfo()[CRITICAL_MACH] = 0.99;
    rModelPart.GetProcessInfo()[UPWIND_FACTOR_CONSTANT] = 1.0;
    rModelPart.GetProcessInfo()[PENALTY_COEFFICIENT] = 100.0;

    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = rModelPart.GetProcessInfo().GetValue(FREE_STREAM_MACH) * rModelPart.GetProcessInfo().GetValue(SOUND_VELOCITY);
    rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3};
    rModelPart.CreateNewElement("EmbeddedTransonicPerturbationPotentialFlowElement2D3N", 1, elemNodes, pElemProp);
}

void GenerateEmbeddedTransonicPerturbationUpwindElement(ModelPart& rModelPart) {
    // Variables addition
    // Set the element properties
    Properties::Pointer pElemProp = rModelPart.CreateNewProperties(1);

    // Geometry creation
    rModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 3, 4};
    rModelPart.CreateNewElement("EmbeddedTransonicPerturbationPotentialFlowElement2D3N", 2, elemNodes, pElemProp);
}

void AssignPotentialsToNormalEmbeddedTransonicPerturbationElement(Element::Pointer pElement)
{
    std::array<double, 3> potential{1.0, 100.0, 150.0};

    for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
}

void AssignDistancesEmbeddedTransonicPerturbationElement(Element::Pointer pElement)
{
    std::array<double, 3> distances{-1.0, 1.0, 1.0};
    for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = distances[i];
}

void AssignPerturbationPotentialsToEmbeddedTransonicElement(Element& rElement, const std::array<double, 3> rPotential) {
    for (unsigned int i = 0; i < 3; i++){
        rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i];
    }
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the RHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(EmbeddedTransonicPerturbationPotentialFlowElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateEmbeddedTransonicPerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);

    AssignPotentialsToNormalEmbeddedTransonicPerturbationElement(p_element);
    AssignDistancesEmbeddedTransonicPerturbationElement(p_element);

    // Compute RHS
    Vector RHS = ZeroVector(3);

    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{109.69824459475,-91.606971325612,-18.091273269139};

    KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-12);
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the inlet RHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(EmbeddedTransonicPerturbationPotentialFlowInletElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateEmbeddedTransonicPerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->SetFlags(INLET);
    p_element->Initialize(r_current_process_info);

    AssignDistancesEmbeddedTransonicPerturbationElement(p_element);
    AssignPotentialsToNormalEmbeddedTransonicPerturbationElement(p_element);

    // Compute RHS
    Vector RHS = ZeroVector(3);

    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{109.69824459475,-91.606971325612,-18.091273269139};

    KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-12);
}

/** Checks the TransonicPerturbationPotentialFlowElement.
 * Checks the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(EmbeddedTransonicPerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateEmbeddedTransonicPerturbationElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);

    AssignDistancesEmbeddedTransonicPerturbationElement(p_element);
    AssignPotentialsToNormalEmbeddedTransonicPerturbationElement(p_element);

    // Compute LHS
    Matrix LHS = ZeroMatrix(4, 4);

    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::array<double, 9> reference{ 0.045857088483312,-0.097966128805804,0.052109040322493,
                                     -0.097966128805804,0.50330688816856,-0.40534075936275,
                                      0.052109040322493,-0.40534075936275,0.35323171904026,
                                      };

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_EXPECT_RELATIVE_NEAR(LHS(i, j), reference[i * 3 + j], 1e-12);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(EmbeddedTransonicPerturbationPotentialFlowSupersonicElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateEmbeddedTransonicPerturbationElement(model_part);
    GenerateEmbeddedTransonicPerturbationUpwindElement(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    Element::Pointer p_upwind_element = model_part.pGetElement(2);

    FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType> find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    p_element -> Set(INLET, false);
    p_upwind_element -> Set(ACTIVE, true);
    AssignDistancesEmbeddedTransonicPerturbationElement(p_element);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);

    std::array<double, 3> high_potential{1.0, 200.0, 100.0};  // node id order 23 74 55
    std::array<double, 3> low_potential{1.0, 100.0, 150.0};   // node id order 23 55 67
    // mach number 1.92516
    AssignPerturbationPotentialsToEmbeddedTransonicElement(*p_element, high_potential);
    // mach number 0.3999
    AssignPerturbationPotentialsToEmbeddedTransonicElement(*p_upwind_element, low_potential);


    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    p_element->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    p_upwind_element->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55};
    std::vector<int> upwind_ids{23, 55, 67};
    for (int i = 0; i < 3; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[i]->SetEquationId(upwind_ids[i]);
    }

    p_element->Initialize(r_current_process_info);
    p_upwind_element->SetFlags(INLET);

    // Compute RHS
    Vector RHS = ZeroVector(4);

    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{138.942250054897102,-173.403842905459442,34.4615928505623046,0.0};

    KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-15);
}

/** Checks the EmbeddedTransonicPerturbationPotentialFlowElement
 * Tests the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(PingEmbeddedTransonicPerturbationPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateEmbeddedTransonicPerturbationElement(model_part);
    GenerateEmbeddedTransonicPerturbationUpwindElement(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    Element::Pointer p_upwind_element = model_part.pGetElement(2);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();
    p_element -> Set(INLET, false);
    p_upwind_element -> Set(ACTIVE, true);
    AssignDistancesEmbeddedTransonicPerturbationElement(p_element);

    std::array<double, 3> potential{1.0, 100.0, 150.0};
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);

    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    PotentialFlowTestUtilities::ComputeElementalSensitivities<3>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

/** Checks the EmbeddedTransonicPerturbationPotentialFlowElement when it's INLET
 * Tests the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(PingEmbeddedTransonicPerturbationInletPotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateEmbeddedTransonicPerturbationElement(model_part);
    GenerateEmbeddedTransonicPerturbationUpwindElement(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    Element::Pointer p_upwind_element = model_part.pGetElement(2);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();
    p_element -> Set(INLET, true);
    p_upwind_element -> Set(ACTIVE, true);
    AssignDistancesEmbeddedTransonicPerturbationElement(p_element);

    std::array<double, 3> potential{1.0, 100.0, 150.0};
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->Initialize(r_current_process_info);

    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    PotentialFlowTestUtilities::ComputeElementalSensitivities<3>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}


/** Checks the EmbeddedTransonicPerturbationPotentialFlowElement when it's a WAKE element.
 * Tests the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(PingEmbeddedTransonicPerturbationWakePotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {

    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateEmbeddedTransonicPerturbationElement(model_part);
    GenerateEmbeddedTransonicPerturbationUpwindElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    Element::Pointer p_upwind_element = model_part.pGetElement(2);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();
    p_upwind_element -> Set(ACTIVE, true);
    p_element -> Set(INLET, true);

    std::array<double, 6> potential{1.39572, 110.69275, 121.1549827, 104.284736, 2.39572, 46.69275};

    Matrix LHS_finite_diference = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<3>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}



/** Checks the EmbeddedTransonicPerturbationPotentialFlowElement.
 * Tests the LHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(PingEmbeddedTransonicPerturbationPotentialFlowSupersonicAccElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateEmbeddedTransonicPerturbationElement(model_part);
    GenerateEmbeddedTransonicPerturbationUpwindElement(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    AssignDistancesEmbeddedTransonicPerturbationElement(p_element);
    Element::Pointer p_upwind_element = model_part.pGetElement(2);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType> find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element -> Set(INLET, false);
    p_upwind_element -> Set(INLET, true);
    p_upwind_element -> Set(ACTIVE, true);
    p_element->Initialize(r_current_process_info);

    std::array<double, 3> high_potential{1.0, 200.0, 100.0};  // node id order 23 74 55
    std::array<double, 3> low_potential{1.0, 100.0, 150.0};   // node id order 23 55 67
    // mach number 1.92516
    AssignPerturbationPotentialsToEmbeddedTransonicElement(*p_element, high_potential);
    // mach number 0.39943
    AssignPerturbationPotentialsToEmbeddedTransonicElement(*p_upwind_element, low_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    p_element->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    p_upwind_element->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55}; // 1 2 3
    std::vector<int> upwind_ids{23, 55, 67};  // 1 3 4
    for (int i = 0; i < 3; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[i]->SetEquationId(upwind_ids[i]);
    }

    Vector RHS_original = ZeroVector(number_of_nodes + 1);
    Matrix LHS_original = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);

    // Compute original RHS and LHS
    p_element->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);

    double delta = 1e-3;
    for(unsigned int i = 0; i < 4; i++){
        // Pinging
        if (i < 3) {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
        }
        else {
            p_upwind_element->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
        }

        Vector RHS_pinged = ZeroVector(number_of_nodes + 1);
        Matrix LHS_pinged = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
        // Compute pinged LHS and RHS
        p_element->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

        for(unsigned int k = 0; k < number_of_nodes + 1; k++){
            // Compute the finite difference estimate of the sensitivity
            LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
            // Compute the average of the original and pinged analytic sensitivities
            LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
        }

        // Unpinging
        if (i < 3) {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
        else {
            p_upwind_element->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
    }
    KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingEmbeddedTransonicPerturbationPotentialFlowSupersonicDecElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateEmbeddedTransonicPerturbationElement(model_part);
    GenerateEmbeddedTransonicPerturbationUpwindElement(model_part);

    Element::Pointer p_element = model_part.pGetElement(1);
    AssignDistancesEmbeddedTransonicPerturbationElement(p_element);
    Element::Pointer p_upwind_element = model_part.pGetElement(2);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType> find_nodal_neighbours_process(model_part);
    find_nodal_neighbours_process.Execute();

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element -> Set(INLET, false);
    p_upwind_element -> Set(INLET, true);
    p_upwind_element -> Set(ACTIVE, true);
    p_element->Initialize(r_current_process_info);

    std::array<double, 3> high_potential{1.0, 200.0, 100.0};  // node id order 23 74 55
    std::array<double, 3> low_potential{1.0, 100.0, 150.0};   // node id order 23 55 67
    // mach number 0.39943
    AssignPerturbationPotentialsToEmbeddedTransonicElement(*p_element, low_potential);
    // mach number 1.92516
    AssignPerturbationPotentialsToEmbeddedTransonicElement(*p_upwind_element, high_potential);

    for (auto& r_node : model_part.Nodes()){
        r_node.AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType CurrentElementalDofList;
    p_element->GetDofList(CurrentElementalDofList, r_current_process_info);

    Element::DofsVectorType UpwindElementalDofList;
    p_upwind_element->GetDofList(UpwindElementalDofList, r_current_process_info);

    std::vector<int> current_ids{23, 74, 55}; // 1 2 3
    std::vector<int> upwind_ids{23, 55, 67};  // 1 3 4
    for (int i = 0; i < 3; i++) {
        CurrentElementalDofList[i]->SetEquationId(current_ids[i]);
        UpwindElementalDofList[i]->SetEquationId(upwind_ids[i]);
    }

    Vector RHS_original = ZeroVector(number_of_nodes + 1);
    Matrix LHS_original = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);

    // Compute original RHS and LHS
    p_element->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);

    double delta = 1e-3;
    for(unsigned int i = 0; i < 4; i++){
        // Pinging
        if (i < 3) {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
        }
        else {
            p_upwind_element->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
        }

        Vector RHS_pinged = ZeroVector(number_of_nodes + 1);
        Matrix LHS_pinged = ZeroMatrix(number_of_nodes + 1, number_of_nodes + 1);
        // Compute pinged LHS and RHS
        p_element->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

        for(unsigned int k = 0; k < number_of_nodes + 1; k++){
            // Compute the finite difference estimate of the sensitivity
            LHS_finite_diference( k, i) = -(RHS_pinged(k)-RHS_original(k)) / delta;
            // Compute the average of the original and pinged analytic sensitivities
            LHS_analytical( k, i) = 0.5 * (LHS_original(k,i) + LHS_pinged(k,i));
        }

        // Unpinging
        if (i < 3) {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
        else {
            p_upwind_element->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
    }
    KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

} // namespace Testing
} // namespace Kratos.

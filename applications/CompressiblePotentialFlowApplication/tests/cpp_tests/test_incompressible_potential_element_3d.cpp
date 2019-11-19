//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//
//

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_elements/incompressible_potential_flow_element.h"
#include "custom_elements/incompressible_potential_flow_pressure_element.h"
#include "custom_elements/embedded_incompressible_potential_flow_element.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos {
  namespace Testing {

    typedef ModelPart::IndexType IndexType;
    typedef ModelPart::NodeIterator NodeIteratorType;

    void Generate3DPressureElement(ModelPart& rModelPart)
    {
        // Variables addition
        rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
        rModelPart.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);

        // Set the element properties
        rModelPart.CreateNewProperties(0);
        Properties::Pointer pElemProp = rModelPart.pGetProperties(0);
        BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
        free_stream_velocity(0) = 10.0;

        rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;
        rModelPart.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.0;

        // Geometry creation
        rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
        rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
        rModelPart.CreateNewNode(4, 0.0, 0.0, 1.0);
        std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3, 4 };
        rModelPart.CreateNewElement("IncompressiblePotentialFlowPressureElement3D4N", 1, elemNodes, pElemProp);
    }

    void AssignPotentialsToWakeElement3D(Element::Pointer pElement, const array_1d<double, 4>& rDistances)
    {
        // Define the nodal values
        Vector potential(4);
        potential(0) = 0.0;
        potential(1) = 1.0;
        potential(2) = 1.0;
        potential(3) = 1.0;

        for (unsigned int i = 0; i < 4; i++){
          if (rDistances(i) > 0.0)
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential(i);
          else
            pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential(i);
        }
        for (unsigned int i = 0; i < 4; i++){
          if (rDistances(i) < 0.0)
            pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential(i)+5+i;
          else
            pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential(i)+5+i;
        }
    }

    BoundedVector<double,4> AssignDistances3D()
    {
        BoundedVector<double,4> distances;
        distances(0) = 1e-12;
        distances(1) = -1e-12;
        distances(2) = 1e-12;
        distances(3) = 1.0;
        return distances;
    }

    KRATOS_TEST_CASE_IN_SUITE(IncompressiblePotentialFlowPressureElementCalculateLocalSystemWake, CompressiblePotentialApplicationFastSuite)
    {
        Model this_model;
        ModelPart& model_part = this_model.CreateModelPart("Main", 3);

        Generate3DPressureElement(model_part);
        Element::Pointer pElement = model_part.pGetElement(1);

        BoundedVector<double,4> distances = AssignDistances3D();

        pElement->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
        pElement->GetValue(WAKE) = true;
        pElement->GetValue(WING_TIP) = true;
        pElement->GetValue(WING_TIP_ELEMENT) = true;
        pElement->Set(STRUCTURE);

        pElement->GetGeometry()[0].SetValue(TRAILING_EDGE, true);
        pElement->GetGeometry()[2].SetValue(TRAILING_EDGE, true);

        AssignPotentialsToWakeElement3D(pElement, distances);

        // Compute RHS and LHS
        Vector RHS = ZeroVector(8);
        Matrix LHS = ZeroMatrix(8, 8);

        pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

        // std::cout.precision(5);
        // std::cout << std::scientific;
        // std::cout << std::showpos;

        const double delta = 1e-4;
        Vector sensitivity_fd = ZeroVector(8);
        Vector sensitivity_an = ZeroVector(8);
        Vector diff = ZeroVector(8);
        Vector LHS_or = ZeroVector(8);
        Vector LHS_perturbed_vec = ZeroVector(8);
        for(unsigned int i = 0; i < 4; i++){
            // Perturbing
            if(distances(i)>0.0){
                pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
            }
            else{
                pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;
            }

            // Compute RHS and LHS
            Vector RHS_perturbed = ZeroVector(8);
            Matrix LHS_perturbed = ZeroMatrix(8, 8);

            pElement->CalculateLocalSystem(LHS_perturbed, RHS_perturbed, model_part.GetProcessInfo());
            sensitivity_fd(i) = -(RHS_perturbed(8) - RHS(8)) / delta;
            sensitivity_an(i) = 0.5 * (LHS_perturbed(8,i) + LHS(8,i));
            diff(i) = std::abs(sensitivity_fd(i) - sensitivity_an(i));
            LHS_or(i) = LHS(8,i);
            LHS_perturbed_vec(i) = LHS_perturbed(8,i);

            // Returning to original value
            if(distances(i)>0.0){
                pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
            }
            else{
                pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
            }
        }

        for(unsigned int i = 0; i < 4; i++){
            // Perturbing
            if(distances(i)>0.0){
                pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;
            }
            else{
                pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
            }

            // Compute RHS and LHS
            Vector RHS_perturbed = ZeroVector(8);
            Matrix LHS_perturbed = ZeroMatrix(8, 8);

            pElement->CalculateLocalSystem(LHS_perturbed, RHS_perturbed, model_part.GetProcessInfo());
            sensitivity_fd(i+4) = -(RHS_perturbed(8) - RHS(8)) / delta;
            sensitivity_an(i+4) = 0.5 * (LHS_perturbed(8,i+4) + LHS(8,i+4));
            diff(i+4) = std::abs(sensitivity_fd(i+4) - sensitivity_an(i+4));
            LHS_or(i+4) = LHS(8,i+4);
            LHS_perturbed_vec(i+4) = LHS_perturbed(8,i+4);

            // Returning to original value
            if(distances(i)>0.0){
                pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
            }
            else{
                pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
            }
        }

        std::cout << std::endl;
        KRATOS_WATCH(diff)
        KRATOS_WATCH(sensitivity_fd)
        KRATOS_WATCH(LHS_or)

        // Check the RHS values (the RHS is computed as the LHS x previous_solution,
        // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
        // std::array<double, 6> reference({0.5, 0.0, 0.0, 0.0, 0.0, -0.5});

        // for (unsigned int i = 0; i < RHS.size(); i++) {
        //   KRATOS_CHECK_NEAR(RHS(i), reference[i], 1e-6);
        // }
    }
  } // namespace Testing
}  // namespace Kratos.

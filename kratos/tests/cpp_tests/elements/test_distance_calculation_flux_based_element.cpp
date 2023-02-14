//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pablo Becker
//

// System includes
#include <limits>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "elements/distance_calculation_flux_based_element.h"

// Utility includes
#include "utilities/geometry_utilities.h"

namespace Kratos
{
namespace Testing
{
    KRATOS_TEST_CASE_IN_SUITE(TestDistanceCalculationFluxBasedElement, KratosCoreFastSuite)
    {
        typedef GlobalPointersVector<Element> ElementPointerVector;

        Model current_model;
        ModelPart& model_part = current_model.CreateModelPart("Main");
        model_part.SetBufferSize(1);

        // Variables addition
        model_part.AddNodalSolutionStepVariable(DISTANCE);

        // Process info creation
        auto& r_process_info = model_part.GetProcessInfo();

        r_process_info.SetValue(CHARACTERISTIC_LENGTH, 1.0);

        // Geometry creation
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
        std::vector<ModelPart::IndexType> elem_nodes {1, 2, 3, 4};
        auto p_elem_prop = model_part.CreateNewProperties(0);
        auto p_element = model_part.CreateNewElement("DistanceCalculationFluxBasedElement3D4N", 1, elem_nodes, p_elem_prop);

        // Define the geom and nodal values, adding 4 fake neighbours so that convective term is added
        auto& r_geom = p_element->GetGeometry();
        ElementPointerVector empty_vector;
        p_element->SetValue(NEIGHBOUR_ELEMENTS, empty_vector);
        ElementPointerVector& rElemNeighbours = p_element->GetValue(NEIGHBOUR_ELEMENTS);
        for (unsigned int i = 0; i < 4; i++) {
            r_geom[i].AddDof(DISTANCE);
            r_geom[i].FastGetSolutionStepValue(DISTANCE) = i;
            rElemNeighbours.push_back(GlobalPointer<Element>(p_element));
        }

        // Compute RHS and LHS
        Vector rhs = ZeroVector(4);
        Matrix lhs = ZeroMatrix(4,4);

        //Checking the potential system
        p_element->Check(r_process_info);
        r_process_info.SetValue(REDISTANCE_STEP,1);
        p_element->CalculateLocalSystem(lhs, rhs, r_process_info);
        const double tolerance = 1.0e-5;
        Vector rhs_reference = ZeroVector(4);
        Matrix lhs_reference = ZeroMatrix(4,4);
        rhs_reference[0] = 1.04167;
        rhs_reference[1] = -0.166667;
        rhs_reference[2] = -0.375;
        rhs_reference[3] = -0.583333;
        lhs_reference(0,0) = 0.541667;
        lhs_reference(0,1) = -0.166667;
        lhs_reference(0,2) = -0.166667;
        lhs_reference(0,3) = -0.166667;
        lhs_reference(1,0) = -0.166667;
        lhs_reference(1,1) = 0.208333;
        lhs_reference(1,2) = 0.0;
        lhs_reference(1,3) = 0.0;
        lhs_reference(2,0) = -0.166667;
        lhs_reference(2,1) = 0.0;
        lhs_reference(2,2) = 0.208333;
        lhs_reference(2,3) = 0.0;
        lhs_reference(3,0) = -0.166667;
        lhs_reference(3,1) = 0.0;
        lhs_reference(3,2) = 0.0;
        lhs_reference(3,3) = 0.208333;
        KRATOS_CHECK_VECTOR_NEAR(rhs, rhs_reference, tolerance);
        KRATOS_CHECK_MATRIX_NEAR(lhs, lhs_reference, tolerance);

        //checking the explicit assembly
        p_element->AddExplicitContribution(r_process_info);
        for (unsigned int i = 0; i < 4; i++) {
            KRATOS_CHECK_NEAR(r_geom[i].GetValue(NODAL_VOLUME),1.0/24.0,tolerance);
            KRATOS_WATCH(r_geom[i].GetValue(POTENTIAL_GRADIENT));
        }

        //Checking the distance system
        r_process_info.SetValue(REDISTANCE_STEP,2);
        p_element->CalculateLocalSystem(lhs, rhs, r_process_info);
        rhs_reference[0] = 0.0535326;
        rhs_reference[1] = -0.0237635;
        rhs_reference[2] = -0.0297173;
        rhs_reference[3] = -0.0356711;
        lhs_reference(0,0) = 0.0309465;
        lhs_reference(0,1) = -0.00551524;
        lhs_reference(0,2) = -0.0103155;
        lhs_reference(0,3) = -0.0151158;
        lhs_reference(1,0) = -0.0141958;
        lhs_reference(1,1) = 0.00296177;
        lhs_reference(1,2) = 0.00449361;
        lhs_reference(1,3) = 0.00674041;
        lhs_reference(2,0) = -0.01726;
        lhs_reference(2,1) = 0.0027575;
        lhs_reference(2,2) = 0.00622996;
        lhs_reference(2,3) = 0.00827249;
        lhs_reference(3,0) = -0.0203241;
        lhs_reference(3,1) = 0.00326819;
        lhs_reference(3,2) = 0.00653638;
        lhs_reference(3,3) = 0.0105195;
        KRATOS_CHECK_VECTOR_NEAR(rhs, rhs_reference, tolerance);
        KRATOS_CHECK_MATRIX_NEAR(lhs, lhs_reference, tolerance);

    }

} // namespace Testing.
} // namespace Kratos.

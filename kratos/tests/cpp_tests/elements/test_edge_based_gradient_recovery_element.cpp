//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//

// System includes
#include <limits>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "elements/levelset_convection_element_simplex.h"

// Utility includes
#include "utilities/geometry_utilities.h"

namespace Kratos
{
namespace Testing
{
    KRATOS_TEST_CASE_IN_SUITE(EdgeBasedGradientRecoveryElement2D, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& model_part = current_model.CreateModelPart("Main");
        model_part.SetBufferSize(1);

        // Variables addition
        model_part.AddNodalSolutionStepVariable(NODAL_MAUX);
        model_part.AddNodalSolutionStepVariable(NODAL_VAUX);

        // Process info creation
        model_part.GetProcessInfo().SetValue(GRADIENT_PENALTY_COEFFICIENT, 1.0e-6);

        // Geometry creation
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        model_part.CreateNewNode(2, 2.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> elem_nodes {1, 2};
        auto p_elem_prop = model_part.CreateNewProperties(0);
        auto p_element = model_part.CreateNewElement("EdgeBasedGradientRecoveryElement2D2N", 1, elem_nodes, p_elem_prop);

        // Define the nodal values
        auto& r_geom = p_element->GetGeometry();
        r_geom[0].SetValue(NODAL_MAUX, 1.0);
        r_geom[1].SetValue(NODAL_MAUX, 2.0);

        // Compute RHS and LHS
        Vector rhs = ZeroVector(4);
        Matrix lhs = ZeroMatrix(4,4);
        const auto& r_process_info = model_part.GetProcessInfo();
        p_element->CalculateLocalSystem(lhs, rhs, r_process_info);

        // Check values
        const double tolerance = 1.0e-8;
        Vector rhs_reference = ZeroVector(4);
        Matrix lhs_reference = ZeroMatrix(4,4);
        rhs_reference[0] = 0.8;
        rhs_reference[1] = 0.4;
        rhs_reference[2] = 0.8;
        rhs_reference[3] = 0.4;
        lhs_reference(0,0) = 0.800002236068;
        lhs_reference(0,1) = 0.4;
        lhs_reference(0,2) = 0.799997763932;
        lhs_reference(0,3) = 0.4;
        lhs_reference(1,0) = 0.4;
        lhs_reference(1,1) = 0.200002236068;
        lhs_reference(1,2) = 0.4;
        lhs_reference(1,3) = 0.199997763932;
        lhs_reference(2,0) = 0.799997763932;
        lhs_reference(2,1) = 0.4;
        lhs_reference(2,2) = 0.800002236068;
        lhs_reference(2,3) = 0.4;
        lhs_reference(3,0) = 0.4;
        lhs_reference(3,1) = 0.199997763932;
        lhs_reference(3,2) = 0.4;
        lhs_reference(3,3) = 0.200002236068;
        KRATOS_EXPECT_VECTOR_NEAR(rhs, rhs_reference, tolerance);
        KRATOS_EXPECT_MATRIX_NEAR(lhs, lhs_reference, tolerance);
    }

} // namespace Testing.
} // namespace Kratos.

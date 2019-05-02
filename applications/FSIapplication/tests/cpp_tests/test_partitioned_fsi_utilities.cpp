//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/global_variables.h"
#include "containers/model.h"
#include "spaces/ublas_space.h"
#include "testing/testing.h"

// Application includes
#include "custom_utilities/partitioned_fsi_utilities.hpp"

namespace Kratos {
namespace Testing {

    typedef UblasSpace<double, Matrix, Vector> SpaceType;

    void SetTestInterface(ModelPart &rTestModelPart)
    {
        rTestModelPart.AddNodalSolutionStepVariable(TEMPERATURE);
        rTestModelPart.AddNodalSolutionStepVariable(SCALAR_PROJECTED);
        rTestModelPart.AddNodalSolutionStepVariable(SCALAR_INTERFACE_RESIDUAL);
        rTestModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
        rTestModelPart.AddNodalSolutionStepVariable(MESH_DISPLACEMENT);
        rTestModelPart.AddNodalSolutionStepVariable(FSI_INTERFACE_RESIDUAL);

        rTestModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        rTestModelPart.CreateNewNode(2, 0.0, 1.0, 0.0);
        rTestModelPart.CreateNewNode(3, 0.0, 2.0, 0.0);
        rTestModelPart.CreateNewNode(4, 0.0, 3.0, 0.0);
    }

    void SetTestInterfaceConditions2D(ModelPart &rTestModelPart)
    {
        Properties::Pointer p_cond_prop = Kratos::make_shared<Properties>();
        rTestModelPart.AddProperties(p_cond_prop, 0);
        rTestModelPart.CreateNewCondition("Condition2D2N", 1, {{1,2}}, p_cond_prop);
        rTestModelPart.CreateNewCondition("Condition2D2N", 2, {{2,3}}, p_cond_prop);
        rTestModelPart.CreateNewCondition("Condition2D2N", 3, {{3,4}}, p_cond_prop);
    }

    void SetTestDoubleValues(ModelPart &rTestModelPart)
    {
        for (auto &it_node : rTestModelPart.Nodes()) {
            it_node.FastGetSolutionStepValue(TEMPERATURE) = it_node.Y();
            it_node.FastGetSolutionStepValue(SCALAR_PROJECTED) = it_node.Y() + 1.0;
        }
    }

    void SetTestArrayValues(ModelPart &rTestModelPart)
    {
        for (auto &it_node : rTestModelPart.Nodes()) {
            array_1d<double,3> aux_disp = ZeroVector(3);
            aux_disp[0] = it_node.Y();
            aux_disp[1] = 2.0 * it_node.Y();
            aux_disp[2] = 3.0 * it_node.Y();
            array_1d<double,3> aux_mesh_disp = ZeroVector(3);
            aux_mesh_disp[0] = 2.0 * it_node.Y();
            aux_mesh_disp[1] = 3.0 * it_node.Y();
            aux_mesh_disp[2] = 4.0 * it_node.Y();
            it_node.FastGetSolutionStepValue(DISPLACEMENT) = aux_disp;
            it_node.FastGetSolutionStepValue(MESH_DISPLACEMENT) = aux_mesh_disp;
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesDoubleGetInterfaceResidualsize, FSIApplicationFastSuite)
    {
        // Set test model part
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetTestInterface(r_test_model_part);

        // Test GetInterfaceResidualsize()
        PartitionedFSIUtilities<SpaceType, double, 2> part_fsi_utils_double_2D;
        const int residual_size = part_fsi_utils_double_2D.GetInterfaceResidualSize(r_test_model_part);
        KRATOS_CHECK_EQUAL(residual_size, 4);
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesArray2DGetInterfaceResidualsize, FSIApplicationFastSuite)
    {
        // Set test model part
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetTestInterface(r_test_model_part);

        // Test GetInterfaceResidualsize()
        PartitionedFSIUtilities<SpaceType, array_1d<double,3>, 2> part_fsi_utils_array_2D;
        const int residual_size = part_fsi_utils_array_2D.GetInterfaceResidualSize(r_test_model_part);
        KRATOS_CHECK_EQUAL(residual_size, 8);
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesArray3DGetInterfaceResidualsize, FSIApplicationFastSuite)
    {
        // Set test model part
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetTestInterface(r_test_model_part);

        // Test GetInterfaceResidualsize()
        PartitionedFSIUtilities<SpaceType, array_1d<double,3>, 3> part_fsi_utils_array_3D;
        const int residual_size = part_fsi_utils_array_3D.GetInterfaceResidualSize(r_test_model_part);
        KRATOS_CHECK_EQUAL(residual_size, 12);
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesDoubleSetUpInterfaceVector, FSIApplicationFastSuite)
    {
        // Set test model part
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetTestInterface(r_test_model_part);

        // Test SetUpInterfaceVector()
        PartitionedFSIUtilities<SpaceType, double, 2> part_fsi_utils_array_2D;
        auto p_interface_vector = part_fsi_utils_array_2D.SetUpInterfaceVector(r_test_model_part);
        KRATOS_CHECK_EQUAL(p_interface_vector->size(), 4);
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesArray2DSetUpInterfaceVector, FSIApplicationFastSuite)
    {
        // Set test model part
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetTestInterface(r_test_model_part);

        // Test SetUpInterfaceVector()
        PartitionedFSIUtilities<SpaceType, array_1d<double,3>, 2> part_fsi_utils_array_2D;
        auto p_interface_vector = part_fsi_utils_array_2D.SetUpInterfaceVector(r_test_model_part);
        KRATOS_CHECK_EQUAL(p_interface_vector->size(), 8);
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesArray3DSetUpInterfaceVector, FSIApplicationFastSuite)
    {
        // Set test model part
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetTestInterface(r_test_model_part);

        // Test SetUpInterfaceVector()
        PartitionedFSIUtilities<SpaceType, array_1d<double,3>, 3> part_fsi_utils_array_3D;
        auto p_interface_vector = part_fsi_utils_array_3D.SetUpInterfaceVector(r_test_model_part);
        KRATOS_CHECK_EQUAL(p_interface_vector->size(), 12);
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesDouble2DComputeInterfaceResidualVector, FSIApplicationFastSuite)
    {
        // Set test model part
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetTestInterface(r_test_model_part);
        SetTestDoubleValues(r_test_model_part);

        // Test ComputeInterfaceResidualVector()
        PartitionedFSIUtilities<SpaceType, double, 2> part_fsi_utils_array_2D;
        auto p_interface_vector = part_fsi_utils_array_2D.SetUpInterfaceVector(r_test_model_part);
        part_fsi_utils_array_2D.ComputeInterfaceResidualVector(
            r_test_model_part,
            TEMPERATURE,
            SCALAR_PROJECTED,
            SCALAR_INTERFACE_RESIDUAL,
            *p_interface_vector);

        const double tolerance = 1.0e-8;
        std::vector<double> expected_values = {-1.0, -1.0, -1.0, -1.0};
        for (unsigned int i = 0; i < expected_values.size(); ++i) {
            KRATOS_CHECK_NEAR(expected_values[i], (*p_interface_vector)[i], tolerance);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesArray2DComputeInterfaceResidualVectorNodal, FSIApplicationFastSuite)
    {
        // Set test model part
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetTestInterface(r_test_model_part);
        SetTestArrayValues(r_test_model_part);

        // Test ComputeInterfaceResidualVector()
        PartitionedFSIUtilities<SpaceType, array_1d<double,3>, 2> part_fsi_utils_array_2D;
        auto p_interface_vector = part_fsi_utils_array_2D.SetUpInterfaceVector(r_test_model_part);

        part_fsi_utils_array_2D.ComputeInterfaceResidualVector(
            r_test_model_part,
            DISPLACEMENT,
            MESH_DISPLACEMENT,
            FSI_INTERFACE_RESIDUAL,
            *p_interface_vector);

        const double tolerance = 1.0e-8;
        std::vector<double> expected_values = {0.0,0.0,-1.0,-1.0,-2.0,-2.0,-3.0,-3.0};
        for (unsigned int i = 0; i < expected_values.size(); ++i) {
            KRATOS_CHECK_NEAR(expected_values[i], (*p_interface_vector)[i], tolerance);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesArray2DComputeInterfaceResidualVectorConsistent, FSIApplicationFastSuite)
    {
        // Set test model part
        Model model;
        auto &r_test_model_part = model.CreateModelPart("TestModelPart");
        SetTestInterface(r_test_model_part);
        SetTestArrayValues(r_test_model_part);
        SetTestInterfaceConditions2D(r_test_model_part);

        // Test ComputeInterfaceResidualVector()
        PartitionedFSIUtilities<SpaceType, array_1d<double,3>, 2> part_fsi_utils_array_2D;
        auto p_interface_vector = part_fsi_utils_array_2D.SetUpInterfaceVector(r_test_model_part);

        part_fsi_utils_array_2D.ComputeInterfaceResidualVector(
            r_test_model_part,
            DISPLACEMENT,
            MESH_DISPLACEMENT,
            FSI_INTERFACE_RESIDUAL,
            *p_interface_vector,
            "consistent");

        const double tolerance = 1.0e-8;
        std::vector<double> expected_values = {-0.166666667,-0.166666667,-1.0,-1.0,-2.0,-2.0,-1.333333333,-1.333333333};
        for (unsigned int i = 0; i < expected_values.size(); ++i) {
            KRATOS_CHECK_NEAR(expected_values[i], (*p_interface_vector)[i], tolerance);
        }
    }

} // namespace Testing
}  // namespace Kratos.

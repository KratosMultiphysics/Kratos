//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
#include "includes/define.h"
#include "includes/model_part_io.h"
#include "includes/global_variables.h"
#include "processes/structured_mesh_generator_process.h"
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
        rTestModelPart.CreateNewCondition("LineCondition2D2N", 1, {{1,2}}, p_cond_prop);
        rTestModelPart.CreateNewCondition("LineCondition2D2N", 2, {{2,3}}, p_cond_prop);
        rTestModelPart.CreateNewCondition("LineCondition2D2N", 3, {{3,4}}, p_cond_prop);
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

    void GenerateTestSkinModelPart(ModelPart &rTestSkinModelPart)
    {
        rTestSkinModelPart.SetBufferSize(1);
        Properties::Pointer p_properties_1(new Properties(1));
        rTestSkinModelPart.AddProperties(p_properties_1);
        rTestSkinModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        rTestSkinModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
        rTestSkinModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
        rTestSkinModelPart.CreateNewNode(4, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> cond_nodes_1 = {1,2};
        std::vector<ModelPart::IndexType> cond_nodes_2 = {2,3};
        std::vector<ModelPart::IndexType> cond_nodes_3 = {3,4};
        rTestSkinModelPart.CreateNewCondition("LineCondition2D2N", 1, cond_nodes_1, p_properties_1);
        rTestSkinModelPart.CreateNewCondition("LineCondition2D2N", 2, cond_nodes_2, p_properties_1);
        rTestSkinModelPart.CreateNewCondition("LineCondition2D2N", 3, cond_nodes_3, p_properties_1);
    }

    /**
     * @brief Helper accessing class
     * This helper class has the unique purpose to give access to the protectedmethods
     * of PartitionedFSIUtilities. It aims to make it possible its standalone testing
     * @tparam TSpace Linear algebra space type
     * @tparam TValueType Residual and unknown value type (array_1d<double, 3> or double)
     * @tparam TDim Problem domain size
     */
    template<class TSpace, class TValueType, unsigned int TDim>
    class PartitionedFSIUtilitiesAccessor : public PartitionedFSIUtilities<TSpace, TValueType, TDim>
    {
    public:

        PartitionedFSIUtilitiesAccessor()
        {
        }

        void AccessComputeConsistentResidual(
            ModelPart &rModelPart,
            Variable<TValueType> &rOriginalVariable,
            Variable<TValueType> &rModifiedVariable,
            Variable<TValueType> &rResidualVariable)
        {
            this->ComputeConsistentResidual(rModelPart, rOriginalVariable, rModifiedVariable, rResidualVariable);
        }

        void AccessComputeNodeByNodeResidual(
            ModelPart &rModelPart,
            Variable<TValueType> &rOriginalVariable,
            Variable<TValueType> &rModifiedVariable,
            Variable<TValueType> &rResidualVariable)
        {
            this->ComputeNodeByNodeResidual(rModelPart, rOriginalVariable, rModifiedVariable, rResidualVariable);
        }
    };

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesComputeConsistentResidual, FSIApplicationFastSuite)
    {
        // Set the partitioned FSI utilities
        PartitionedFSIUtilities<SpaceType,array_1d<double,3>,2> partitioned_fsi_utilities;

        // Set the model part in where the error is to be computed
        Model model;
        ModelPart &main_model_part = model.CreateModelPart("OriginModelPart");
        main_model_part.AddNodalSolutionStepVariable(VELOCITY);
        main_model_part.AddNodalSolutionStepVariable(VECTOR_PROJECTED);
        main_model_part.AddNodalSolutionStepVariable(FSI_INTERFACE_RESIDUAL);
        GenerateTestSkinModelPart(main_model_part);

        // Set a fake field to compute the error
        array_1d<double,3> unit_val;
        unit_val[0] = 1.0;
        unit_val[1] = 1.0;
        unit_val[2] = 1.0;
        for (auto &it_node : main_model_part.Nodes()) {
            it_node.FastGetSolutionStepValue(VELOCITY) = unit_val;
            it_node.FastGetSolutionStepValue(VECTOR_PROJECTED) = 2.0 * unit_val;
        }

        // Compute the consistent residual
        PartitionedFSIUtilitiesAccessor<SpaceType,array_1d<double,3>,2> aux_partitioned_fsi_utilities_accessor;
        aux_partitioned_fsi_utilities_accessor.AccessComputeConsistentResidual(
            main_model_part,
            VECTOR_PROJECTED,
            VELOCITY,
            FSI_INTERFACE_RESIDUAL);

        // Check results
        double tol = 1.0e-10;
        unsigned int aux_count = 0;
        std::array<double, 12> expected_values = {0.5,0.5,0.5,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.5,0.5};
        for (auto &r_node : main_model_part.Nodes()) {
            const auto &r_fsi_int_res = r_node.FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL);
            for (unsigned int i = 0; i < 3; ++i) {
                KRATOS_CHECK_NEAR(r_fsi_int_res[i], expected_values[3 * aux_count + i], tol);
            }
            aux_count++;
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesComputeNodeByNodeResidual, FSIApplicationFastSuite)
    {
        // Set the partitioned FSI utilities
        PartitionedFSIUtilities<SpaceType,array_1d<double,3>,2> partitioned_fsi_utilities;

        // Set the model part in where the error is to be computed
        Model model;
        ModelPart &main_model_part = model.CreateModelPart("OriginModelPart");
        main_model_part.AddNodalSolutionStepVariable(VELOCITY);
        main_model_part.AddNodalSolutionStepVariable(VECTOR_PROJECTED);
        main_model_part.AddNodalSolutionStepVariable(FSI_INTERFACE_RESIDUAL);
        GenerateTestSkinModelPart(main_model_part);

        // Set a fake field to compute the error
        array_1d<double,3> unit_val;
        unit_val[0] = 1.0;
        unit_val[1] = 1.0;
        unit_val[2] = 1.0;
        for (auto &it_node : main_model_part.Nodes()) {
            it_node.FastGetSolutionStepValue(VELOCITY) = unit_val;
            it_node.FastGetSolutionStepValue(VECTOR_PROJECTED) = 2.0 * unit_val;
        }

        // Compute the node by node residual
        PartitionedFSIUtilitiesAccessor<SpaceType,array_1d<double,3>,2> aux_partitioned_fsi_utilities_accessor;
        aux_partitioned_fsi_utilities_accessor.AccessComputeNodeByNodeResidual(
            main_model_part,
            VECTOR_PROJECTED,
            VELOCITY,
            FSI_INTERFACE_RESIDUAL);

        // Check results
        double tol = 1.0e-10;
        unsigned int aux_count = 0;
        std::array<double, 12> expected_values = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
        for (auto &r_node : main_model_part.Nodes()) {
            const auto &r_fsi_int_res = r_node.FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL);
            for (unsigned int i = 0; i < 3; ++i) {
                KRATOS_CHECK_NEAR(r_fsi_int_res[i], expected_values[3 * aux_count + i], tol);
            }
            aux_count++;
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesCopySkinToElements, FSIApplicationFastSuite)
    {
        // Set the partitioned FSI utilities
        PartitionedFSIUtilities<SpaceType,array_1d<double,3>,2> partitioned_fsi_utilities;

        // Set the model part that contains the condition based skin
        Model model;
        ModelPart &main_model_part = model.CreateModelPart("OriginModelPart");
        GenerateTestSkinModelPart(main_model_part);

        // Set the new skin model part
        ModelPart &element_based_skin = model.CreateModelPart("ElementBasedSkin");
        partitioned_fsi_utilities.CreateCouplingElementBasedSkin(main_model_part, element_based_skin);

        // Check results
        KRATOS_CHECK_EQUAL(element_based_skin.NumberOfNodes(), 4);
        KRATOS_CHECK_EQUAL(element_based_skin.NumberOfElements(), 3);
        KRATOS_CHECK_EQUAL(element_based_skin.NumberOfConditions(), 0);
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

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesArray2DCreateCouplingElementBasedSkin, FSIApplicationFastSuite)
    {
        // Set the partitioned FSI utilities
        PartitionedFSIUtilities<SpaceType,array_1d<double,3>,2> partitioned_fsi_utilities;

        // Set the model part containing the origin skin
        Model model;
        ModelPart &r_skin_model_part = model.CreateModelPart("OriginModelPart");
        GenerateTestSkinModelPart(r_skin_model_part);

        // Create the destination submodelpart
        ModelPart &r_destination_model_part = model.CreateModelPart("DestinationModelPart");

        // Generate the auxiliary element based skin
        partitioned_fsi_utilities.CreateCouplingElementBasedSkin(
            r_skin_model_part,
            r_destination_model_part);

        // Check results
        KRATOS_CHECK_EQUAL(r_destination_model_part.NumberOfNodes(), 4);
        KRATOS_CHECK_EQUAL(r_destination_model_part.NumberOfElements(), 3);
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesArray2DInitializeInterfaceVector, FSIApplicationFastSuite)
    {
        // Set the partitioned FSI utilities
        PartitionedFSIUtilities<SpaceType,array_1d<double,3>,2> partitioned_fsi_utilities;

        // Set the model part containing the origin skin
        Model model;
        ModelPart &r_skin_model_part = model.CreateModelPart("OriginModelPart");
        r_skin_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        GenerateTestSkinModelPart(r_skin_model_part);

        // Set a fake field to fill the interface vector
        for (auto &r_node : r_skin_model_part.Nodes()) {
            auto &r_disp = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            r_disp[0] = r_node.Id();
            r_disp[1] = 2.0 * r_node.Id();
            r_disp[2] = 3.0 * r_node.Id();
        }

        // Generate the auxiliary element based skin
        Vector interface_vector;
        interface_vector.resize(partitioned_fsi_utilities.GetInterfaceResidualSize(r_skin_model_part), false);
        partitioned_fsi_utilities.InitializeInterfaceVector(
            r_skin_model_part,
            DISPLACEMENT,
            interface_vector);

        // Check results
        const double tolerance = 1.0e-8;
        std::array<double, 8> expected_results = {1.0,2.0,2.0,4.0,3.0,6.0,4.0,8.0};
        for (unsigned int i = 0; i < 8; ++i) {
            KRATOS_CHECK_NEAR(interface_vector(i), expected_results[i], tolerance);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesArray2DEmbeddedToPositiveFacePressureInterpolator, FSIApplicationFastSuite)
    {
        // Set the partitioned FSI utilities
        PartitionedFSIUtilities<SpaceType,array_1d<double,3>,2> partitioned_fsi_utilities;

        // Set the model part containing the origin skin
        Model model;
        ModelPart &r_skin_model_part = model.CreateModelPart("OriginModelPart");
        r_skin_model_part.AddNodalSolutionStepVariable(PRESSURE);
        r_skin_model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
        GenerateTestSkinModelPart(r_skin_model_part);

        // Create a background mesh
        ModelPart &r_background_model_part = model.CreateModelPart("BackgroundModelPart");
        r_background_model_part.AddNodalSolutionStepVariable(PRESSURE);

        auto p_point_1 = Kratos::make_intrusive<Node<3>>(1, -2.0, -2.0, 0.0);
        auto p_point_2 = Kratos::make_intrusive<Node<3>>(2, -2.0,  3.0, 0.0);
        auto p_point_3 = Kratos::make_intrusive<Node<3>>(3,  3.0,  3.0, 0.0);
        auto p_point_4 = Kratos::make_intrusive<Node<3>>(4,  3.0, -3.0, 0.0);

        Quadrilateral2D4<Node<3>> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

        Parameters mesher_parameters(R"(
        {
            "number_of_divisions": 7,
            "element_name": "Element2D3N"
        })");
        StructuredMeshGeneratorProcess(geometry, r_background_model_part, mesher_parameters).Execute();

        // Set a fake pressure field in the background mesh
        for (auto &r_node : r_background_model_part.Nodes()) {
            r_node.FastGetSolutionStepValue(PRESSURE) = r_node.Id();
        }

        // Call the embedded pressure interpolation method
        partitioned_fsi_utilities.EmbeddedPressureToPositiveFacePressureInterpolator(
            r_background_model_part,
            r_skin_model_part);

        // Copy the obtained values from PRESSURE to POSITIVE_FACE_PRESSURE
        for (auto &r_node : r_skin_model_part.Nodes()) {
            r_node.FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) = r_node.FastGetSolutionStepValue(PRESSURE);
        }

        // Check results
        unsigned int i = 0;
        const double tolerance = 1.0e-4;
        std::array<double, 8> expected_results = {26.5105,37.8462,39.0974,27.8053};
        for (const auto &r_node : r_skin_model_part.Nodes()) {
            KRATOS_CHECK_NEAR(r_node.FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE), expected_results[i], tolerance);
            i++;
        }
    }


} // namespace Testing
}  // namespace Kratos.

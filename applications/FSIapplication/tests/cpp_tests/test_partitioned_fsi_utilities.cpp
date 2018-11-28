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
#include "includes/define.h"
#include "includes/model_part_io.h"
#include "spaces/ublas_space.h"
#include "testing/testing.h"

// Application includes
#include "custom_utilities/partitioned_fsi_utilities.hpp"

namespace Kratos {
namespace Testing {

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
        rTestSkinModelPart.CreateNewCondition("Condition2D2N", 1, cond_nodes_1, p_properties_1);
        rTestSkinModelPart.CreateNewCondition("Condition2D2N", 2, cond_nodes_2, p_properties_1);
        rTestSkinModelPart.CreateNewCondition("Condition2D2N", 3, cond_nodes_3, p_properties_1);
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesCopySkinToElements, FSIApplicationFastSuite)
    {
        // Set the partitioned FSI utilities
        typedef UblasSpace<double, Matrix, Vector > TSpace;
        PartitionedFSIUtilities<TSpace, 2> partitioned_fsi_utilities;

        // Set the model part that contains the condition based skin
        Model model;
        ModelPart &main_model_part = model.CreateModelPart("OriginModelPart");
        GenerateTestSkinModelPart(main_model_part);

        // Set the new skin model part
        ModelPart &element_based_skin = model.CreateModelPart("ElementBasedSkin");
        partitioned_fsi_utilities.CopySkinToElements(main_model_part, element_based_skin);

        // Check results
        KRATOS_CHECK_EQUAL(element_based_skin.NumberOfNodes(), 4);
        KRATOS_CHECK_EQUAL(element_based_skin.NumberOfElements(), 3);
        KRATOS_CHECK_EQUAL(element_based_skin.NumberOfConditions(), 0);
    }

    KRATOS_TEST_CASE_IN_SUITE(PartitionedFSIUtilitiesComputeConsistentResidual, FSIApplicationFastSuite)
    {
        // Set the partitioned FSI utilities
        typedef UblasSpace<double, Matrix, Vector > TSpace;
        PartitionedFSIUtilities<TSpace, 2> partitioned_fsi_utilities;

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
        
        // Compute the consistent residual vector
        partitioned_fsi_utilities.ComputeConsistentResidual(
            main_model_part,
            VECTOR_PROJECTED,
            VELOCITY,
            FSI_INTERFACE_RESIDUAL);

        KRATOS_WATCH(main_model_part.GetNode(1).FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL));
        KRATOS_WATCH(main_model_part.GetNode(2).FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL));
        KRATOS_WATCH(main_model_part.GetNode(3).FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL));
        KRATOS_WATCH(main_model_part.GetNode(4).FastGetSolutionStepValue(FSI_INTERFACE_RESIDUAL));

    }

} // namespace Testing
} // namespace Kratos.

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_parameters.h"
#include "containers/model.h"
#include "utilities/model_part_utils.h"

namespace Kratos::Testing {
    using PropertiesType = Properties;
    using ElementType = Element;
    using ConditionType = Condition;
    using MeshType = Mesh<Node, PropertiesType, ElementType, ConditionType>;
    using NodesContainerType = MeshType::NodesContainerType;

    KRATOS_TEST_CASE_IN_SUITE(AddNodesFromOrderedContainer, KratosCoreFastSuite)
    {
        Model model;

        auto& r_model_part = model.CreateModelPart("Main");

        auto& r_smp = r_model_part.CreateSubModelPart("Inlet1");
        auto& r_ssmp = r_model_part.CreateSubModelPart("Inlet1.sub_inlet");
        auto& r_smp2 = r_model_part.CreateSubModelPart("Inlet2");
        auto& r_ssmp2 = r_model_part.CreateSubModelPart("Inlet2.sub_inlet");
        NodesContainerType aux, aux2, aux3; //auxiliary containers

        // Now we create some nodes and store some of them in a container
        for(std::size_t id=1; id<25; id++) {
             auto p_node = r_model_part.CreateNewNode(id, 0.00,0.00,0.00);
             if (id%2!=0) {
                aux.push_back(p_node);
             }
            if (id<15) {
                aux2.push_back(p_node);
             }
             aux3.push_back(p_node);
        }
        // We check the first one
        r_ssmp.AddNodes(aux.begin(), aux.end());
        KRATOS_EXPECT_EQ(r_ssmp.NumberOfNodes(), 12);
        KRATOS_EXPECT_EQ(r_smp.NumberOfNodes(), 12);

        ModelPartUtils::AddNodesFromOrderedContainer(r_ssmp2, aux.begin(), aux.end());
        KRATOS_EXPECT_EQ(r_ssmp2.NumberOfNodes(), 12);
        KRATOS_EXPECT_EQ(r_smp2.NumberOfNodes(), 12);

        // We check the second one
        r_ssmp.AddNodes(aux2.begin(), aux2.end());
        KRATOS_EXPECT_EQ(r_ssmp.NumberOfNodes(), 19);
        KRATOS_EXPECT_EQ(r_smp.NumberOfNodes(), 19);

        ModelPartUtils::AddNodesFromOrderedContainer(r_ssmp2, aux2.begin(), aux2.end());
        KRATOS_EXPECT_EQ(r_ssmp2.NumberOfNodes(), 19);
        KRATOS_EXPECT_EQ(r_smp2.NumberOfNodes(), 19);

        // Now with aux3, we start with ordered because it needs to be unique
        ModelPartUtils::AddNodesFromOrderedContainer(r_ssmp2, aux3.begin(), aux3.end());
        KRATOS_EXPECT_EQ(r_ssmp2.NumberOfNodes(), 24);
        KRATOS_EXPECT_EQ(r_smp2.NumberOfNodes(), 24);

        // Here we can go a bit further. No need to be Unique
        for(auto it=aux.begin();it!=aux.end(); it++){
            aux3.push_back(*(it.base()));
        }
        r_ssmp.AddNodes(aux3.begin(), aux3.end());
        KRATOS_EXPECT_EQ(r_ssmp.NumberOfNodes(), 24);
        KRATOS_EXPECT_EQ(r_smp.NumberOfNodes(), 24);
    }

	} // namespace Testing


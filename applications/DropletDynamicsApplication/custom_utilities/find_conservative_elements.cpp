//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors: 
//
//

// System includes
#include <limits>

// External includes

// Project includes
#include "custom_utilities/find_conservative_elements.h"
#include "utilities/variable_utils.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/rbf_shape_functions_utility.h"
#include "includes/model_part_io.h"


namespace Kratos
{

template<bool THistorical>
void FindConservativeElementsProcess<THistorical>::FINDCUTELEMENTNODES()
{

    //detect all the nodes that belong to cons elements
    ModelPart::NodesContainerType cons_elem_nodes;
    ModelPart::ElementsContainerType cons_elem;
    cons_elem_nodes.clear();
    cons_elem.clear();

    // Create a new ModelPart to store cons elements and nodes
    if (!mrModelPart.HasSubModelPart("ConsModelPart")) {
    ModelPart& consModelPart = mrModelPart.CreateSubModelPart("ConsModelPart");
    // Additional initialization if needed
    }

    // Get the reference to the sub-model part
    ModelPart& consModelPart = mrModelPart.GetSubModelPart("ConsModelPart");

    // Clear the existing elements and nodes in the sub-model part
    consModelPart.Elements().clear();
    consModelPart.Nodes().clear();

    for(const auto& r_elem : mrModelPart.Elements()){
        const auto& r_geom = r_elem.GetGeometry();

        double dist;
        unsigned int num = 0;
        
        for (const auto& r_node : r_geom) {
            dist = r_node.FastGetSolutionStepValue(DISTANCE);
            if (dist > -0.4999 && dist < 0.4999){++num;}
        }

        if(num > 0){ //element is inside the conservative bound
            cons_elem.push_back(mrModelPart.pGetElement(r_elem.Id()));
            for(unsigned int i = 0; i < r_geom.size(); ++i) {
                cons_elem_nodes.push_back(r_geom(i));
            }
        }
    }

    cons_elem_nodes.Unique();
    cons_elem.Unique();

    // Add cons element nodes to the new ModelPart
    for (const auto& pNode : cons_elem_nodes) {
        consModelPart.Nodes().push_back(Kratos::intrusive_ptr<Node>(const_cast<Node*>(&pNode))); 
    }

    // Add cons elements to the new ModelPart
    for (const auto& r_elem : cons_elem) {
        consModelPart.Elements().push_back(Kratos::intrusive_ptr<Element>(const_cast<Element*>(&r_elem)));
    }
    //ModelPartIO consModelPartIO("ConsModelPart", IO::WRITE);
    //consModelPartIO.WriteModelPart(consModelPart); 
    
    // Open the mdpa file for writing
    std::ofstream mdpa_file("ConsModelPart.mdpa");

    mdpa_file << "Begin ModelPartData\n";
    mdpa_file << "End ModelPartData\n";
    mdpa_file << "  Begin Properties 0\n"; // You may need to adjust the property ID
    mdpa_file << "  End Properties\n";
    mdpa_file << "  Begin Nodes\n";

    for (const auto& pNode : cons_elem_nodes) {
        mdpa_file << "    " << pNode.Id() << "    " << pNode.X() << "    " << pNode.Y() << "    " << pNode.Z() << "\n";
    }

    mdpa_file << "  End Nodes\n";
    mdpa_file << "  Begin Elements Element2D3N\n";

    for (const auto& r_elem : cons_elem) {
        mdpa_file << "    " << r_elem.Id() << "    0    " << r_elem.GetGeometry()[0].Id() << "    " << r_elem.GetGeometry()[1].Id() << "    " << r_elem.GetGeometry()[2].Id() << "\n";
    }

    mdpa_file << "  End Elements\n";

    // Write the nodes of the "ConsModelPart" sub-model part
    mdpa_file << "Begin SubModelPart ConsSubModelPart\n";
    mdpa_file << "  Begin SubModelPartNodes" << std::endl;
    //mdpa_file << cons_elem_nodes.size() << std::endl;
    for (const auto& pNode : cons_elem_nodes) {
        mdpa_file << pNode.Id() << std::endl;
    }
    mdpa_file << "  End SubModelPartNodes\n" << std::endl;
    mdpa_file << "  Begin SubModelPartElements" << std::endl;
    //mdpa_file << cons_elem.size() << std::endl;
    for (const auto& r_elem : cons_elem) {
        mdpa_file << r_elem.Id() << std::endl;
    }
    mdpa_file << "  End SubModelPartElements\n" << std::endl;
    mdpa_file << "End SubModelPart ConsSubModelPart";
    // Close the mdpa file
    mdpa_file.close();

    // Print cons element nodes to the output
    KRATOS_INFO("CutElementNodes") << "Cut Element Nodes: ";
    for (const auto& pNode : cons_elem_nodes) {
        KRATOS_INFO("CutElementNodes") << "Node ID: " << pNode.Id() << ", Coordinates: " << pNode.Coordinates() << std::endl;
    }
    // Print cons element IDs to the output
    KRATOS_INFO("CutElementIDs") << "Cut Element IDs: ";
    for (const auto& r_elem : cons_elem) {
        KRATOS_INFO("CutElementIDs") << "Element ID: " << r_elem.Id() << std::endl;
    }
}

template<bool THistorical>
void FindConservativeElementsProcess<THistorical>::Execute()
{
    KRATOS_TRY
    FINDCUTELEMENTNODES();

    // Check if variables are available
    if (THistorical) {
        KRATOS_ERROR_IF_NOT(mrModelPart.NodesBegin()->SolutionStepsDataHas( DISTANCE )) << "Variable DISTANCE not in the model part!" << std::endl;
    }


    KRATOS_CATCH("")
}

// /***********************************************************************************/
// /***********************************************************************************/

template<>
double& FindConservativeElementsProcess<true>::GetDISValue(NodeType& rNode)
{
    return rNode.FastGetSolutionStepValue(DISTANCE);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double& FindConservativeElementsProcess<false>::GetDISValue(NodeType& rNode)
{
    return rNode.GetValue(DISTANCE);
}

/***********************************************************************************/
/***********************************************************************************/

template class FindConservativeElementsProcess<true>;
template class FindConservativeElementsProcess<false>;

} // namespace Kratos

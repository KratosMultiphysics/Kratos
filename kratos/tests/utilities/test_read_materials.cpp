// System includes

// External includes

// Project includes
//#include "includes/define.h"
#include "testing/testing.h"
#include "utilities/read_materials_utility.hpp"
#include "geometries/triangle_2d_3.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ReadMaterialsUtility, KratosCoreFastSuite)
{
    Parameters parameters(R"(
        {
        "properties": [{
            "model_part_name": "MATERIAL_HOMOGENEOUS",
            "properties_id": 1,
            "Material": {
                "name": "homog",
                "constitutive_law": {
                    "name": "LinearElastic3DLaw"
                },
                "Variables": {
                    "DENSITY": 7850.0,
                    "YOUNG_MODULUS": 210000000000.0,
                    "POISSON_RATIO": 0.3
                    },
                 "Tables": {
                    "Table1": {
        		    "input_variable": "TEMPERATURE",
        		    "output_variable": "YOUNG_MODULUS",
        		    "data": [ [0, 100], [1, 200], [2, 300] ]
                    },
        		        "Table2": {
                        "input_variable": "TEMPERATURE",
        	                "output_variable": "YOUNG_MODULUS",
        	                "data": [ [3, 400], [4, 500], [5, 600] ]
                        }
        	            }
                }
            }, {
            "model_part_name": "MATERIAL_MULTISCALE",
            "properties_id": 2,
            "Material": {
                "name": "multiscale",
                "constitutive_law": {
                    "name": "LinearElastic3DLaw"
                },
                "Variables": {
                    "DENSITY": 7.0,
                    "YOUNG_MODULUS": 2100.0,
                    "POISSON_RATIO": 0.3
                    },
                "Tables":  {}
                }
            }]
        })");

    ModelPart model_part("dummy");
    //Properties::Pointer p_elem_prop = model_part.pGetProperties(0);
    //NodeType::Pointer p_node_1 = model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    //NodeType::Pointer p_node_2 = model_part.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    //NodeType::Pointer p_node_3 = model_part.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    //std::vector<NodeType::Pointer> element_nodes_0 (3);
    //element_nodes_0[0] = p_node_1;
    //element_nodes_0[1] = p_node_2;
    //element_nodes_0[2] = p_node_3;
    //Triangle2D3 <NodeType> triangle_0( element_nodes_0 );
    //Element::Pointer p_elem_0 = model_part.CreateNewElement("Element2D3N", 1, triangle_0, p_elem_prop);
    //p_sub_modelpart_1->AddElement(p_elem_0);

    ReadMaterialsUtility(model_part, parameters);

    Properties::Pointer p_prop = model_part.pGetProperties(1, 0);

    //p_prop->
    //auto variable = KratosComponents<Variable<double>>().Get(iter.name());
    KRATOS_CHECK_EQUAL(p_prop["Variables"]["DENSITY"].GetDouble(), 1000);
}
}
} // namespace Kratos.

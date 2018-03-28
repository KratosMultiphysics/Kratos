// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
// 

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h" 
#include "custom_utilities/interface_preprocess.h" 

/* Geometries */
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

namespace Kratos
{
template<>
void InterfacePreprocessCondition::GenerateInterfacePart<2>(
    ModelPart& rOriginPart,
    ModelPart& rInterfacePart,
    Parameters ThisParameters
    )
{
    KRATOS_TRY;
    
    Parameters default_parameters = Parameters(R"(
    {
        "simplify_geometry"                    : false
    })" );
    
    ThisParameters.ValidateAndAssignDefaults(default_parameters);
    
    const bool simplest_geometry = ThisParameters["simplify_geometry"].GetBool();
    
    IndexType cond_counter = 0;
    
    // We reorder the conditions
    IndexType cond_id = ReorderConditions();
    
    // Store new properties in a map
    std::unordered_map<IndexType, Properties::Pointer> new_properties = CreateNewProperties(rOriginPart);

    // Generate Conditions from original the edges that can be considered interface
    for (auto it_elem = rOriginPart.ElementsBegin(); it_elem != rOriginPart.ElementsEnd(); ++it_elem) {
        GeometryType& this_geometry = it_elem->GetGeometry();
        Properties::Pointer p_prop = new_properties[it_elem->pGetProperties()->Id()];
        KRATOS_DEBUG_ERROR_IF(p_prop == nullptr) << "ERROR:: Property not well initialized" << std::endl;
        
        for (IndexType i_edge = 0; i_edge < this_geometry.EdgesNumber(); ++i_edge)
            GenerateEdgeCondition(rInterfacePart, p_prop, this_geometry.Edges()[i_edge], simplest_geometry, cond_counter, cond_id);
    }
    
    // NOTE: Reorder ID if parallellization
    
    const IndexType num_nodes = mrMainModelPart.Nodes().size();
    PrintNodesAndConditions(num_nodes, cond_counter);
    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void InterfacePreprocessCondition::GenerateInterfacePart<3>(
    ModelPart& rOriginPart,
    ModelPart& rInterfacePart,
    Parameters ThisParameters
    )
{
    KRATOS_TRY;
    
    Parameters default_parameters = Parameters(R"(
    {
        "simplify_geometry"                    : false
    })" );
    
    ThisParameters.ValidateAndAssignDefaults(default_parameters);
    
    const bool simplest_geometry = ThisParameters["simplify_geometry"].GetBool();
    
    IndexType cond_counter = 0;
    
    // We reorder the conditions
    IndexType cond_id = ReorderConditions();
    
    // Store new properties in a map
    std::unordered_map<IndexType, Properties::Pointer> new_properties = CreateNewProperties(rOriginPart);

    // Generate Conditions from original the faces that can be considered interface
    for (auto it_elem = rOriginPart.ElementsBegin(); it_elem != rOriginPart.ElementsEnd(); ++it_elem) {          
        GeometryType& this_geometry = it_elem->GetGeometry();
        Properties::Pointer p_prop = new_properties[it_elem->pGetProperties()->Id()];
        KRATOS_DEBUG_ERROR_IF(p_prop == nullptr) << "ERROR:: Property not well initialized" << std::endl;
        
        if (this_geometry.LocalSpaceDimension() == 3) {
            for (IndexType i_face = 0; i_face < this_geometry.FacesNumber(); ++i_face)
                GenerateFaceCondition(rInterfacePart, p_prop, this_geometry.Faces()[i_face], simplest_geometry, cond_counter, cond_id);
        } else
            GenerateFaceCondition(rInterfacePart, p_prop, this_geometry, simplest_geometry, cond_counter, cond_id);
    }
    
    // NOTE: Reorder ID if parallellization
    
    const IndexType num_nodes = mrMainModelPart.Nodes().size();
    PrintNodesAndConditions(num_nodes, cond_counter);
    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

std::unordered_map<IndexType, Properties::Pointer> InterfacePreprocessCondition::CreateNewProperties(ModelPart& rOriginPart)
{
    const SizeType number_properties = mrMainModelPart.NumberOfProperties();

    std::unordered_map<IndexType, Properties::Pointer> new_properties;
    // Reorder id Properties
    IndexType count = 0;
    for (auto it_prop = mrMainModelPart.PropertiesBegin(); it_prop < mrMainModelPart.PropertiesEnd(); ++it_prop) {
        ++count;
        it_prop->SetId(count);
    }
    // Create new properties
    count = 0;
    for (auto it_prop = rOriginPart.PropertiesBegin(); it_prop < rOriginPart.PropertiesEnd(); ++it_prop) {
        Properties::Pointer p_original_prop = *(it_prop.base());
        ++count;
        Properties new_prop(number_properties + count + 1);
        Properties::Pointer p_new_prop = Kratos::make_shared<Properties>(new_prop);
        new_properties.insert({p_original_prop->Id(), p_new_prop});

        // Now we copy (an remove) the properties we have interest
        CopyProperties(p_original_prop, p_new_prop, FRICTION_COEFFICIENT);
        CopyProperties(p_original_prop, p_new_prop, THICKNESS, false);
        CopyProperties(p_original_prop, p_new_prop, YOUNG_MODULUS);
    }

    return new_properties;
}

/***********************************************************************************/
/***********************************************************************************/

void InterfacePreprocessCondition::CreateNewCondition(
    Properties::Pointer pThisProperties,
    GeometryType& rGeometry,
    const IndexType CondId,
    Condition const& rCondition
    )
{
    KRATOS_TRY;

    Condition::Pointer p_cond = Condition::Pointer(rCondition.Create(CondId, rGeometry, pThisProperties));
    mrMainModelPart.AddCondition(p_cond);

    // We set the condition as master or slave (master by default)
    GeometryType& this_geometry = p_cond->GetGeometry();
    bool is_slave = true;
    for (IndexType it_node = 0; it_node < this_geometry.size(); ++it_node)
        if (this_geometry[it_node].Is(SLAVE) == false) is_slave = false;
    if (is_slave == true)  p_cond->Set(SLAVE, true);
    else  p_cond->Set(MASTER, true);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void InterfacePreprocessCondition::PrintNodesAndConditions(
    const IndexType NodesCounter,
    const IndexType CondCounter
    )
{
    KRATOS_INFO("Nodes found") << "\t" << NodesCounter << " nodes ";
    KRATOS_INFO("Conditions found") << "and " << CondCounter <<  " conditions found." << std::endl;

    // Check that we actually found something
    KRATOS_ERROR_IF(NodesCounter == 0) << "No interface nodes found. Please check that nodes on both sides of the interface have been assigned Is(INTERFACE) = true." << std::endl;
    KRATOS_ERROR_IF(CondCounter == 0) << "No interface conditions found. Please check that nodes on both sides of the interface have been assigned Is(INTERFACE) = true and that the contact surfaces have been assigned conditions." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

IndexType InterfacePreprocessCondition::ReorderConditions()
{
    // We reorder the conditions
    ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
    const IndexType num_conditions = static_cast<int>(conditions_array.size());
    
    for(IndexType i = 0; i < num_conditions; ++i) {
        auto it_condition = conditions_array.begin() + i;
        it_condition->SetId(i + 1);
    }
    
    return num_conditions;
}

/***********************************************************************************/
/***********************************************************************************/

inline void InterfacePreprocessCondition::GenerateEdgeCondition(
    ModelPart& rInterfacePart,
    Properties::Pointer pThisProperties,
    GeometryType& EdgeGeometry,
    const bool SimplestGeometry,
    IndexType& CondCounter,
    IndexType& CondId
    )
{    
    IndexType count = 0;
    const IndexType number_of_points = EdgeGeometry.PointsNumber();
    for (IndexType it_node = 0; it_node < number_of_points; ++it_node) {
        if (EdgeGeometry[it_node].IsDefined(INTERFACE) == true)
            if (EdgeGeometry[it_node].Is(INTERFACE) == true) ++count;
    }
 
    const std::string condition_name = (number_of_points == 2 || SimplestGeometry) ? "Condition2D2N" : "Condition2D3N";
 
    Condition const& r_condition =  KratosComponents<Condition>::Get(condition_name);
    
    if (count == number_of_points) {
        ++CondId; // NOTE: To paralellize be careful with this ID
        if (number_of_points == 2) {
            // We initialize a vector for the IDs
            std::vector<std::size_t> condition_ids(1);
            
            CreateNewCondition(pThisProperties, EdgeGeometry, CondId, r_condition);
            condition_ids[0] = CondId;
            ++CondCounter;
        
            rInterfacePart.AddConditions(condition_ids);
        } else {                            
            if (SimplestGeometry == false) {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);
                
                CreateNewCondition(pThisProperties, EdgeGeometry, CondId, r_condition);
                condition_ids[0] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            } else {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(2);

                Line2D2< NodeType > lin_1(EdgeGeometry(0), EdgeGeometry(1));
                CreateNewCondition(pThisProperties, lin_1, CondId, r_condition);
                condition_ids[0] = CondId;
                ++CondCounter;
                ++CondId;
                Line2D2< NodeType > lin_2(EdgeGeometry(1), EdgeGeometry(2));
                CreateNewCondition(pThisProperties, lin_2, CondId, r_condition);
                condition_ids[1] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

inline void InterfacePreprocessCondition::GenerateFaceCondition(
    ModelPart& rInterfacePart,
    Properties::Pointer pThisProperties,
    GeometryType& FaceGeometry,
    const bool SimplestGeometry,
    IndexType& CondCounter,
    IndexType& CondId
    )
{
    IndexType count = 0;
    const IndexType number_of_points = FaceGeometry.PointsNumber();
    for (IndexType it_node = 0; it_node < number_of_points; ++it_node) {
        if (FaceGeometry[it_node].IsDefined(INTERFACE) == true)  
            if (FaceGeometry[it_node].Is(INTERFACE) == true) ++count;
    }
    
    const std::string condition_name = (number_of_points == 3 || SimplestGeometry) ? "Condition3D" : (number_of_points == 4) ? "Condition3D4N" : (number_of_points == 6) ? "Condition3D6N" : (number_of_points == 8) ? "Condition3D8N" : "Condition3D9N";
 
    Condition const& r_condition =  KratosComponents<Condition>::Get(condition_name);
    
    if (count == number_of_points) {
        ++CondId;
        if (number_of_points == 3) {
            // We initialize a vector for the IDs
            std::vector<std::size_t> condition_ids(1);
            
            CreateNewCondition(pThisProperties, FaceGeometry, CondId, r_condition);
            condition_ids[0] = CondId;
            ++CondCounter;
            
            rInterfacePart.AddConditions(condition_ids);
        } else if (number_of_points == 4) {
            if (SimplestGeometry == false) {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);
            
                CreateNewCondition(pThisProperties, FaceGeometry, CondId, r_condition);
                condition_ids[0] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            } else {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(2);
                
                Triangle3D3< NodeType > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(2));
                CreateNewCondition(pThisProperties, tri_1, CondId, r_condition);
                condition_ids[0] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_2(FaceGeometry(2), FaceGeometry(3), FaceGeometry(0));
                CreateNewCondition(pThisProperties, tri_2, CondId, r_condition);
                condition_ids[1] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
        } else if (number_of_points == 6) {
            if (SimplestGeometry == false) {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);
                
                CreateNewCondition(pThisProperties, FaceGeometry, CondId, r_condition);
                condition_ids[0] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            } else {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(4);
                
                Triangle3D3< NodeType > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(5));
                CreateNewCondition(pThisProperties, tri_1, CondId, r_condition);
                condition_ids[0] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_2(FaceGeometry(1), FaceGeometry(2), FaceGeometry(3));
                CreateNewCondition(pThisProperties, tri_2, CondId, r_condition);
                condition_ids[1] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_3(FaceGeometry(1), FaceGeometry(3), FaceGeometry(5));
                CreateNewCondition(pThisProperties, tri_3, CondId, r_condition);
                condition_ids[2] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_4(FaceGeometry(3), FaceGeometry(4), FaceGeometry(5));
                CreateNewCondition(pThisProperties, tri_4, CondId, r_condition);
                condition_ids[3] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
        } else if (number_of_points == 8) {
            if (SimplestGeometry == false) {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);
                
                CreateNewCondition(pThisProperties, FaceGeometry, CondId, r_condition);
                condition_ids[0] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            } else {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(6);

                Triangle3D3< NodeType > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(7));
                CreateNewCondition(pThisProperties, tri_1, CondId, r_condition);
                condition_ids[0] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_2(FaceGeometry(1), FaceGeometry(5), FaceGeometry(7));
                CreateNewCondition(pThisProperties, tri_2, CondId, r_condition);
                condition_ids[1] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_3(FaceGeometry(1), FaceGeometry(3), FaceGeometry(5));
                CreateNewCondition(pThisProperties, tri_3, CondId, r_condition);
                condition_ids[2] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_4(FaceGeometry(1), FaceGeometry(2), FaceGeometry(3));
                CreateNewCondition(pThisProperties, tri_4, CondId, r_condition);
                condition_ids[3] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_5(FaceGeometry(3), FaceGeometry(4), FaceGeometry(5));
                CreateNewCondition(pThisProperties, tri_5, CondId, r_condition);
                condition_ids[4] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_6(FaceGeometry(5), FaceGeometry(6), FaceGeometry(7));
                CreateNewCondition(pThisProperties, tri_6, CondId, r_condition);
                condition_ids[5] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
        } else { // Assuming it will not be a very weird geometry
            if (SimplestGeometry == false) {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);
                
                CreateNewCondition(pThisProperties, FaceGeometry, CondId, r_condition);
                condition_ids[0] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            } else {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(8);

                Triangle3D3< NodeType > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(8));
                CreateNewCondition(pThisProperties, tri_1, CondId, r_condition);
                condition_ids[0] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_2(FaceGeometry(1), FaceGeometry(2), FaceGeometry(3));
                CreateNewCondition(pThisProperties, tri_2, CondId, r_condition);
                condition_ids[1] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_3(FaceGeometry(1), FaceGeometry(3), FaceGeometry(8));
                CreateNewCondition(pThisProperties, tri_3, CondId, r_condition);
                condition_ids[2] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_4(FaceGeometry(8), FaceGeometry(3), FaceGeometry(4));
                CreateNewCondition(pThisProperties, tri_4, CondId, r_condition);
                condition_ids[3] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_5(FaceGeometry(8), FaceGeometry(4), FaceGeometry(5));
                CreateNewCondition(pThisProperties, tri_5, CondId, r_condition);
                condition_ids[4] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_6(FaceGeometry(5), FaceGeometry(6), FaceGeometry(7));
                CreateNewCondition(pThisProperties, tri_6, CondId, r_condition);
                condition_ids[5] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_7(FaceGeometry(5), FaceGeometry(7), FaceGeometry(8));
                CreateNewCondition(pThisProperties, tri_7, CondId, r_condition);
                condition_ids[6] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< NodeType > tri_8(FaceGeometry(0), FaceGeometry(8), FaceGeometry(7));
                CreateNewCondition(pThisProperties, tri_8, CondId, r_condition);
                condition_ids[7] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
        }
    }
}
}  // namespace Kratos.

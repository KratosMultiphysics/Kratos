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
    ModelPart& rInterfacePart,
    Parameters ThisParameters
    )
{
    KRATOS_TRY;
    
    Parameters default_parameters = Parameters(R"(
    {
        "simplify_geometry"                    : false,
        "contact_property_id"                  : 0
    })" );
    
    ThisParameters.ValidateAndAssignDefaults(default_parameters);
    
    const bool simplest_geometry = ThisParameters["simplify_geometry"].GetBool();
    const int contact_property_id = ThisParameters["contact_property_id"].GetInt();
    
    IndexType cond_counter = 0;

    // Generate Conditions from original the edges that can be considered interface
    if (rInterfacePart.Conditions().size() > 0) { // We use the already existant conditions geometry (recommended)
        cond_counter = rInterfacePart.Conditions().size();
        // Check and creates the properties
        CheckAndCreateProperties(rInterfacePart);
    } else if (rInterfacePart.Nodes().size() > 0) { // Only in case we have assigned the flag directly to nodes (no conditions)
        // We reorder the conditions
        IndexType cond_id = ReorderConditions();

        // Store new properties in a map
        std::unordered_map<IndexType, Properties::Pointer> new_properties;
        if (contact_property_id == 0) new_properties = CreateNewProperties();

        // We iterate over the elements and check the nodes on the interface
        for (auto it_elem = mrMainModelPart.ElementsBegin(); it_elem != mrMainModelPart.ElementsEnd(); ++it_elem) {
            GeometryType& this_geometry = it_elem->GetGeometry();
            Properties::Pointer p_prop = (contact_property_id == 0) ? new_properties[it_elem->pGetProperties()->Id()] : mrMainModelPart.pGetProperties(contact_property_id);
            KRATOS_DEBUG_ERROR_IF(p_prop == nullptr) << "ERROR:: Property not well initialized" << std::endl;

            if (this_geometry.LocalSpaceDimension() == 2) {
                for (IndexType i_edge = 0; i_edge < this_geometry.EdgesNumber(); ++i_edge)
                    GenerateEdgeCondition(rInterfacePart, p_prop, this_geometry.Edges()[i_edge], simplest_geometry, cond_counter, cond_id);
            } else {
                GenerateEdgeCondition(rInterfacePart, p_prop, this_geometry, simplest_geometry, cond_counter, cond_id);
            }
        }
    } else {
        KRATOS_ERROR << "ERROR:: Nor conditions or nodes on the interface. Check your flags" << std::endl;
    }
    
    // NOTE: Reorder ID if parallellization
    
    const IndexType num_nodes = rInterfacePart.Nodes().size();
    PrintNodesAndConditions(num_nodes, cond_counter);
    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void InterfacePreprocessCondition::GenerateInterfacePart<3>(
    ModelPart& rInterfacePart,
    Parameters ThisParameters
    )
{
    KRATOS_TRY;
    
    Parameters default_parameters = Parameters(R"(
    {
        "simplify_geometry"                    : false,
        "contact_property_id"                  : 0
    })" );
    
    ThisParameters.ValidateAndAssignDefaults(default_parameters);
    
    const bool simplest_geometry = ThisParameters["simplify_geometry"].GetBool();
    const int contact_property_id = ThisParameters["contact_property_id"].GetInt();
    
    IndexType cond_counter = 0;
    
    // Generate Conditions from original the edges that can be considered interface
    if (rInterfacePart.Conditions().size() > 0) { // We use the already existant conditions geometry (recommended)
        cond_counter = rInterfacePart.Conditions().size();
        // Check and creates the properties
        CheckAndCreateProperties(rInterfacePart);
    } else if (rInterfacePart.Nodes().size() > 0) { // Only in case we have assigned the flag directly to nodes (no conditions)
        // We reorder the conditions
        IndexType cond_id = ReorderConditions();

        // Store new properties in a map
        std::unordered_map<IndexType, Properties::Pointer> new_properties;
        if (contact_property_id == 0) new_properties = CreateNewProperties();

        // Generate Conditions from original the faces that can be considered interface
        for (auto it_elem = mrMainModelPart.ElementsBegin(); it_elem != mrMainModelPart.ElementsEnd(); ++it_elem) {
            GeometryType& this_geometry = it_elem->GetGeometry();
            Properties::Pointer p_prop = (contact_property_id == 0) ? new_properties[it_elem->pGetProperties()->Id()] : mrMainModelPart.pGetProperties(contact_property_id);
            KRATOS_DEBUG_ERROR_IF(p_prop == nullptr) << "ERROR:: Property not well initialized" << std::endl;

            if (this_geometry.LocalSpaceDimension() == 3) {
                for (IndexType i_face = 0; i_face < this_geometry.FacesNumber(); ++i_face)
                    GenerateFaceCondition(rInterfacePart, p_prop, this_geometry.Faces()[i_face], simplest_geometry, cond_counter, cond_id);
            } else
                GenerateFaceCondition(rInterfacePart, p_prop, this_geometry, simplest_geometry, cond_counter, cond_id);
        }
    } else {
        KRATOS_ERROR << "ERROR:: Nor conditions or nodes on the interface. Check your flags" << std::endl;
    }
    
    // NOTE: Reorder ID if parallellization
    
    const IndexType num_nodes = rInterfacePart.Nodes().size();
    PrintNodesAndConditions(num_nodes, cond_counter);
    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void InterfacePreprocessCondition::CheckAndCreateProperties(ModelPart& rInterfacePart)
{
    // We check that the properties define what must be defined
    Properties::Pointer p_prop_old = rInterfacePart.Conditions().begin()->pGetProperties();
    if (!(p_prop_old->Has(YOUNG_MODULUS))) {
        // Store new properties in a map
        const std::size_t number_properties = mrMainModelPart.NumberOfProperties();
        Properties::Pointer p_prop_new = mrMainModelPart.pGetProperties(number_properties + 1);

        GeometryType& this_geometry_cond = rInterfacePart.Conditions().begin()->GetGeometry();
        const std::size_t number_of_nodes = this_geometry_cond.size();
        std::vector<IndexType> index_vector(number_of_nodes);
        for (IndexType i_node = 0; i_node < number_of_nodes; i_node++) {
            index_vector[i_node] = this_geometry_cond[i_node].Id();
        }
        std::sort(index_vector.begin(), index_vector.end());

        IndexType counter = 0;
        for (auto it_elem = mrMainModelPart.ElementsBegin(); it_elem != mrMainModelPart.ElementsEnd(); ++it_elem) {
            GeometryType& this_geometry = it_elem->GetGeometry();

            const bool is_on_the_face = CheckOnTheFace(index_vector, this_geometry);

            if (is_on_the_face) {
                Properties::Pointer p_prop = it_elem->pGetProperties();

                // Now we copy (an remove) the properties we have interest
                CopyProperties(p_prop, p_prop_new, FRICTION_COEFFICIENT);
                CopyProperties(p_prop, p_prop_new, THICKNESS, false);
                CopyProperties(p_prop, p_prop_new, YOUNG_MODULUS);

                counter++;
                break;
            }
        }

        // Now we iterate over the conditions
        if (counter > 0) {
            ConditionsArrayType& conditions_array = rInterfacePart.Conditions();

            #pragma omp parallel for
            for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i) {
                auto it_cond = conditions_array.begin() + i;
                it_cond->SetProperties(p_prop_new);
            }
        } else {
            KRATOS_ERROR << "It was not possible to add a property" << std::endl;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool InterfacePreprocessCondition::CheckOnTheFace(
    const std::vector<std::size_t>& rIndexVector,
    GeometryType& rElementGeometry
    )
{
    if (rElementGeometry.WorkingSpaceDimension() == 2) {
        for (IndexType i_edge = 0; i_edge < rElementGeometry.EdgesNumber(); ++i_edge) {
            const IndexType number_of_nodes = rElementGeometry.Edges()[i_edge].size();
            std::vector<IndexType> index_vector_face(number_of_nodes);
            if (number_of_nodes == rIndexVector.size()) {
                for (IndexType i_node = 0; i_node < number_of_nodes; i_node++) {
                    index_vector_face[i_node] = rElementGeometry.Edges()[i_edge][i_node].Id();
                }
                std::sort(index_vector_face.begin(), index_vector_face.end());

                bool is_here = true;
                for (IndexType i_node = 0; i_node < number_of_nodes; i_node++) {
                    if (index_vector_face[i_node] != rIndexVector[i_node]) {
                        is_here = false;
                        break;
                    }
                }
                if (is_here) return true;
            }
        }
    } else {
        if (rElementGeometry.LocalSpaceDimension() == 3) {
            for (IndexType i_face = 0; i_face < rElementGeometry.FacesNumber(); ++i_face) {
                const IndexType number_of_nodes = rElementGeometry.Faces()[i_face].size();
                if (number_of_nodes == rIndexVector.size()) {
                    std::vector<IndexType> index_vector_face(number_of_nodes);
                    for (IndexType i_node = 0; i_node < number_of_nodes; i_node++) {
                        index_vector_face[i_node] = rElementGeometry.Faces()[i_face][i_node].Id();
                    }

                    std::sort(index_vector_face.begin(), index_vector_face.end());

                    bool is_here = true;
                    for (IndexType i_node = 0; i_node < number_of_nodes; i_node++) {
                        if (index_vector_face[i_node] != rIndexVector[i_node]) {
                            is_here = false;
                            break;
                        }
                    }
                    if (is_here) return true;
                }
            }
        } else {
            const IndexType number_of_nodes = rElementGeometry.size();
            std::vector<IndexType> index_vector_face(number_of_nodes);
            for (IndexType i_node = 0; i_node < number_of_nodes; i_node++) {
                index_vector_face[i_node] = rElementGeometry[i_node].Id();
            }
            std::sort(index_vector_face.begin(), index_vector_face.end());

            for (IndexType i_node = 0; i_node < number_of_nodes; i_node++) {
                if (index_vector_face[i_node] != rIndexVector[i_node]) return false;
            }
            return true;
        }
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

std::unordered_map<IndexType, Properties::Pointer> InterfacePreprocessCondition::CreateNewProperties()
{
    const SizeType number_properties = mrMainModelPart.NumberOfProperties();

    std::unordered_map<IndexType, Properties::Pointer> new_properties;
    // Reorder id Properties
    IndexType count = 0;
    for (auto it_prop = mrMainModelPart.PropertiesBegin(); it_prop < mrMainModelPart.PropertiesEnd(); ++it_prop) {
        ++count;
        it_prop->SetId(count);
    }

    // Create a vector with the ids of the old properties
    const SizeType number_properties_origin = mrMainModelPart.NumberOfProperties();
    std::vector<IndexType> index_vector(number_properties_origin);
    count = 0;
    for (auto it_prop = mrMainModelPart.PropertiesBegin(); it_prop < mrMainModelPart.PropertiesEnd(); ++it_prop) {
        index_vector[count] = it_prop->Id();
        ++count;
    }

    // Copy to the new properties
    count = 0;
    for (auto& i_prop : index_vector) {
        Properties::Pointer p_original_prop = mrMainModelPart.pGetProperties(i_prop);
        ++count;
        Properties::Pointer p_new_prop = mrMainModelPart.pGetProperties(number_properties + count + 1);
        new_properties.insert({i_prop, p_new_prop});

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
    AssignMasterSlaveCondition(p_cond);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void InterfacePreprocessCondition::AssignMasterSlaveCondition(Condition::Pointer pCond)
{
    // We set the condition as master or slave (master by default)
    GeometryType& this_geometry = pCond->GetGeometry();
    bool is_slave = true;
    for (IndexType i_node = 0; i_node < this_geometry.size(); ++i_node) {
        if (this_geometry[i_node].Is(SLAVE) == false) {
            is_slave = false;
            break;
        }
    }
    if (is_slave == true) pCond->Set(SLAVE, true);
    else pCond->Set(MASTER, true);
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
    
    const std::string condition_name = (number_of_points == 3 || SimplestGeometry) ? "SurfaceCondition3D3N" : (number_of_points == 4) ? "SurfaceCondition3D4N" : (number_of_points == 6) ? "SurfaceCondition3D6N" : (number_of_points == 8) ? "SurfaceCondition3D8N" : "SurfaceCondition3D9N";
 
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

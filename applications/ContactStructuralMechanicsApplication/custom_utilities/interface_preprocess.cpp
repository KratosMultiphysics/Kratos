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
void InterfacePreprocessCondition::GenerateInterfacePart(
    ModelPart& rInterfacePart,
    Parameters ThisParameters
    )
{
    KRATOS_TRY;

    Parameters default_parameters = Parameters(R"(
    {
        "simplify_geometry"                    : false,
        "contact_property_id"                  : 0
    })");

    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    // Getting some parameters
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
        for (auto& r_elem : mrMainModelPart.Elements()) {
            GeometryType& r_geometry = r_elem.GetGeometry();
            const std::size_t working_space_dimension = r_geometry.WorkingSpaceDimension();
            const std::size_t local_space_dimension = r_geometry.LocalSpaceDimension();

            Properties::Pointer p_prop = (contact_property_id == 0) ? new_properties[r_elem.pGetProperties()->Id()] : mrMainModelPart.CreateNewProperties(contact_property_id);
            KRATOS_DEBUG_ERROR_IF(p_prop == nullptr) << "ERROR:: Property not well initialized" << std::endl;

            if (working_space_dimension == local_space_dimension) {
                const auto faces = r_geometry.GenerateBoundariesEntities();
                for (IndexType i_face = 0; i_face < faces.size(); ++i_face) {
                    if (working_space_dimension == 2) {
                        GenerateEdgeCondition(rInterfacePart, p_prop, faces[i_face], simplest_geometry, cond_counter, cond_id);
                    } else {
                        GenerateFaceCondition(rInterfacePart, p_prop, faces[i_face], simplest_geometry, cond_counter, cond_id);
                    }
                }
            } else {
                if (working_space_dimension == 2) {
                    GenerateEdgeCondition(rInterfacePart, p_prop, r_geometry, simplest_geometry, cond_counter, cond_id);
                } else {
                    GenerateFaceCondition(rInterfacePart, p_prop, r_geometry, simplest_geometry, cond_counter, cond_id);
                }
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

void InterfacePreprocessCondition::CheckAndCreateProperties(ModelPart& rInterfacePart)
{
    KRATOS_TRY;

    // We check that the properties define what must be defined
    Properties::Pointer p_prop_old = rInterfacePart.Conditions().begin()->pGetProperties();
    if (!(p_prop_old->Has(YOUNG_MODULUS))) {
        // Store new properties in a map
        const std::size_t number_properties = mrMainModelPart.NumberOfProperties();
        Properties::Pointer p_prop_new = mrMainModelPart.CreateNewProperties(number_properties + 1);

        GeometryType& this_geometry_cond = rInterfacePart.Conditions().begin()->GetGeometry();
        const std::size_t number_of_nodes = this_geometry_cond.size();
        std::vector<IndexType> index_vector(number_of_nodes);
        for (IndexType i_node = 0; i_node < number_of_nodes; i_node++) {
            index_vector[i_node] = this_geometry_cond[i_node].Id();
        }
        std::sort(index_vector.begin(), index_vector.end());

        IndexType counter = 0;
        for (auto& r_elem : mrMainModelPart.Elements()) {
            GeometryType& r_geometry = r_elem.GetGeometry();

            const bool is_on_the_face = CheckOnTheFace(index_vector, r_geometry);

            if (is_on_the_face) {
                Properties::Pointer p_prop = r_elem.pGetProperties();

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

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

bool InterfacePreprocessCondition::CheckOnTheFace(
    const std::vector<std::size_t>& rIndexVector,
    GeometryType& rElementGeometry
    )
{
    KRATOS_TRY;

    const std::size_t working_space_dimension = rElementGeometry.WorkingSpaceDimension();
    const std::size_t local_space_dimension = rElementGeometry.LocalSpaceDimension();
    if (working_space_dimension == local_space_dimension) {
        const auto faces  = rElementGeometry.GenerateBoundariesEntities();
        for (IndexType i_face = 0; i_face < faces.size(); ++i_face) {
            const IndexType number_of_nodes = faces[i_face].size();
            std::vector<IndexType> index_vector_face(number_of_nodes);
            if (number_of_nodes == rIndexVector.size()) {
                for (IndexType i_node = 0; i_node < number_of_nodes; i_node++) {
                    index_vector_face[i_node] = faces[i_face][i_node].Id();
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

    return false;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

std::unordered_map<IndexType, Properties::Pointer> InterfacePreprocessCondition::CreateNewProperties()
{
    KRATOS_TRY;

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
        Properties::Pointer p_new_prop = mrMainModelPart.CreateNewProperties(number_properties + count + 1);
        new_properties.insert({i_prop, p_new_prop});

        // Now we copy (an remove) the properties we have interest
        CopyProperties(p_original_prop, p_new_prop, FRICTION_COEFFICIENT);
        CopyProperties(p_original_prop, p_new_prop, THICKNESS, false);
        CopyProperties(p_original_prop, p_new_prop, YOUNG_MODULUS);
    }

    return new_properties;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void InterfacePreprocessCondition::CreateNewCondition(
    Properties::Pointer pThisProperties,
    const GeometryType& rGeometry,
    const IndexType ConditionId,
    Condition const& rCondition
    )
{
    KRATOS_TRY;

    Condition::Pointer p_cond = Condition::Pointer(rCondition.Create(ConditionId, rGeometry.Points(), pThisProperties));
    mrMainModelPart.AddCondition(p_cond);
    AssignMasterSlaveCondition(p_cond);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void InterfacePreprocessCondition::AssignMasterSlaveCondition(Condition::Pointer pCond)
{
    KRATOS_TRY;

    // We set the condition as master or slave (master by default)
    GeometryType& this_geometry = pCond->GetGeometry();
    bool is_slave = true;
    for (IndexType i_node = 0; i_node < this_geometry.size(); ++i_node) {
        if (this_geometry[i_node].Is(SLAVE) == false) {
            is_slave = false;
            break;
        }
    }
    if (is_slave) pCond->Set(SLAVE, true);
    else pCond->Set(MASTER, true);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void InterfacePreprocessCondition::PrintNodesAndConditions(
    const IndexType NodesCounter,
    const IndexType CondCounter
    )
{
    KRATOS_TRY;

    KRATOS_INFO("Nodes found") << "\t" << NodesCounter << " nodes ";
    KRATOS_INFO("Conditions found") << "and " << CondCounter <<  " conditions found." << std::endl;

    // Check that we actually found something
    KRATOS_ERROR_IF(NodesCounter == 0) << "No interface nodes found. Please check that nodes on both sides of the interface have been assigned Is(INTERFACE) = true." << std::endl;
    KRATOS_ERROR_IF(CondCounter == 0) << "No interface conditions found. Please check that nodes on both sides of the interface have been assigned Is(INTERFACE) = true and that the contact surfaces have been assigned conditions." << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

IndexType InterfacePreprocessCondition::ReorderConditions()
{
    KRATOS_TRY;

    // We reorder the conditions
    ConditionsArrayType& r_conditions_array = mrMainModelPart.GetRootModelPart().Conditions();
    const IndexType num_conditions = static_cast<int>(r_conditions_array.size());
    const auto it_cond_begin = r_conditions_array.begin();

    for(IndexType i = 0; i < num_conditions; ++i) {
        auto it_condition = it_cond_begin + i;
        it_condition->SetId(i + 1);
    }

    return num_conditions;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

inline void InterfacePreprocessCondition::GenerateEdgeCondition(
    ModelPart& rInterfacePart,
    Properties::Pointer pThisProperties,
    const GeometryType& rEdgeGeometry,
    const bool SimplestGeometry,
    IndexType& rCondCounter,
    IndexType& rConditionId
    )
{
    KRATOS_TRY;

    // First we check if the nodes belong to the interface modelpart
    const auto& r_nodes_array = rInterfacePart.Nodes();
    for (auto& r_node : rEdgeGeometry) {
        auto it_found = r_nodes_array.find(r_node.Id());
        if(it_found == r_nodes_array.end()) { // Node does not exist in the top model part
            return void();
        }
    }

    IndexType count = 0;
    const IndexType number_of_points = rEdgeGeometry.PointsNumber();
    for (IndexType it_node = 0; it_node < number_of_points; ++it_node) {
        if (rEdgeGeometry[it_node].IsDefined(INTERFACE))
            if (rEdgeGeometry[it_node].Is(INTERFACE)) ++count;
    }

    const std::string condition_name = (number_of_points == 2 || SimplestGeometry) ? "LineCondition2D2N" : "LineCondition2D3N";

    Condition const& r_condition =  KratosComponents<Condition>::Get(condition_name);

    if (count == number_of_points) {
        ++rConditionId; // NOTE: To paralellize be careful with this ID
        if (number_of_points == 2) {
            // We initialize a vector for the IDs
            std::vector<std::size_t> condition_ids(1);

            CreateNewCondition(pThisProperties, rEdgeGeometry, rConditionId, r_condition);
            condition_ids[0] = rConditionId;
            ++rCondCounter;

            rInterfacePart.AddConditions(condition_ids);
        } else {
            if (SimplestGeometry == false) {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);

                CreateNewCondition(pThisProperties, rEdgeGeometry, rConditionId, r_condition);
                condition_ids[0] = rConditionId;
                ++rCondCounter;

                rInterfacePart.AddConditions(condition_ids);
            } else {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(2);

                Line2D2< NodeType > lin_1(rEdgeGeometry(0), rEdgeGeometry(1));
                CreateNewCondition(pThisProperties, lin_1, rConditionId, r_condition);
                condition_ids[0] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Line2D2< NodeType > lin_2(rEdgeGeometry(1), rEdgeGeometry(2));
                CreateNewCondition(pThisProperties, lin_2, rConditionId, r_condition);
                condition_ids[1] = rConditionId;
                ++rCondCounter;

                rInterfacePart.AddConditions(condition_ids);
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

inline void InterfacePreprocessCondition::GenerateFaceCondition(
    ModelPart& rInterfacePart,
    Properties::Pointer pThisProperties,
    const GeometryType& rFaceGeometry,
    const bool SimplestGeometry,
    IndexType& rCondCounter,
    IndexType& rConditionId
    )
{
    KRATOS_TRY;

    // First we check if the nodes belong to the interface modelpart
    const auto& r_nodes_array = rInterfacePart.Nodes();
    for (auto& r_node : rFaceGeometry) {
        auto it_found = r_nodes_array.find(r_node.Id());
        if(it_found == r_nodes_array.end()) { // Node does not exist in the top model part
            return void();
        }
    }

    IndexType count = 0;
    const IndexType number_of_points = rFaceGeometry.PointsNumber();
    for (IndexType it_node = 0; it_node < number_of_points; ++it_node) {
        if (rFaceGeometry[it_node].IsDefined(INTERFACE))
            if (rFaceGeometry[it_node].Is(INTERFACE)) ++count;
    }

    const std::string condition_name = (number_of_points == 3 || SimplestGeometry) ? "SurfaceCondition3D3N" : (number_of_points == 4) ? "SurfaceCondition3D4N" : (number_of_points == 6) ? "SurfaceCondition3D6N" : (number_of_points == 8) ? "SurfaceCondition3D8N" : "SurfaceCondition3D9N";

    Condition const& r_condition =  KratosComponents<Condition>::Get(condition_name);

    if (count == number_of_points) {
        ++rConditionId;
        if (number_of_points == 3) {
            // We initialize a vector for the IDs
            std::vector<std::size_t> condition_ids(1);

            CreateNewCondition(pThisProperties, rFaceGeometry, rConditionId, r_condition);
            condition_ids[0] = rConditionId;
            ++rCondCounter;

            rInterfacePart.AddConditions(condition_ids);
        } else if (number_of_points == 4) {
            if (SimplestGeometry == false) {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);

                CreateNewCondition(pThisProperties, rFaceGeometry, rConditionId, r_condition);
                condition_ids[0] = rConditionId;
                ++rCondCounter;

                rInterfacePart.AddConditions(condition_ids);
            } else {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(2);

                Triangle3D3< NodeType > tri_1(rFaceGeometry(0), rFaceGeometry(1), rFaceGeometry(2));
                CreateNewCondition(pThisProperties, tri_1, rConditionId, r_condition);
                condition_ids[0] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_2(rFaceGeometry(2), rFaceGeometry(3), rFaceGeometry(0));
                CreateNewCondition(pThisProperties, tri_2, rConditionId, r_condition);
                condition_ids[1] = rConditionId;
                ++rCondCounter;

                rInterfacePart.AddConditions(condition_ids);
            }
        } else if (number_of_points == 6) {
            if (SimplestGeometry == false) {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);

                CreateNewCondition(pThisProperties, rFaceGeometry, rConditionId, r_condition);
                condition_ids[0] = rConditionId;
                ++rCondCounter;

                rInterfacePart.AddConditions(condition_ids);
            } else {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(4);

                Triangle3D3< NodeType > tri_1(rFaceGeometry(0), rFaceGeometry(1), rFaceGeometry(5));
                CreateNewCondition(pThisProperties, tri_1, rConditionId, r_condition);
                condition_ids[0] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_2(rFaceGeometry(1), rFaceGeometry(2), rFaceGeometry(3));
                CreateNewCondition(pThisProperties, tri_2, rConditionId, r_condition);
                condition_ids[1] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_3(rFaceGeometry(1), rFaceGeometry(3), rFaceGeometry(5));
                CreateNewCondition(pThisProperties, tri_3, rConditionId, r_condition);
                condition_ids[2] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_4(rFaceGeometry(3), rFaceGeometry(4), rFaceGeometry(5));
                CreateNewCondition(pThisProperties, tri_4, rConditionId, r_condition);
                condition_ids[3] = rConditionId;
                ++rCondCounter;

                rInterfacePart.AddConditions(condition_ids);
            }
        } else if (number_of_points == 8) {
            if (SimplestGeometry == false) {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);

                CreateNewCondition(pThisProperties, rFaceGeometry, rConditionId, r_condition);
                condition_ids[0] = rConditionId;
                ++rCondCounter;

                rInterfacePart.AddConditions(condition_ids);
            } else {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(6);

                Triangle3D3< NodeType > tri_1(rFaceGeometry(0), rFaceGeometry(1), rFaceGeometry(7));
                CreateNewCondition(pThisProperties, tri_1, rConditionId, r_condition);
                condition_ids[0] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_2(rFaceGeometry(1), rFaceGeometry(5), rFaceGeometry(7));
                CreateNewCondition(pThisProperties, tri_2, rConditionId, r_condition);
                condition_ids[1] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_3(rFaceGeometry(1), rFaceGeometry(3), rFaceGeometry(5));
                CreateNewCondition(pThisProperties, tri_3, rConditionId, r_condition);
                condition_ids[2] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_4(rFaceGeometry(1), rFaceGeometry(2), rFaceGeometry(3));
                CreateNewCondition(pThisProperties, tri_4, rConditionId, r_condition);
                condition_ids[3] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_5(rFaceGeometry(3), rFaceGeometry(4), rFaceGeometry(5));
                CreateNewCondition(pThisProperties, tri_5, rConditionId, r_condition);
                condition_ids[4] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_6(rFaceGeometry(5), rFaceGeometry(6), rFaceGeometry(7));
                CreateNewCondition(pThisProperties, tri_6, rConditionId, r_condition);
                condition_ids[5] = rConditionId;
                ++rCondCounter;

                rInterfacePart.AddConditions(condition_ids);
            }
        } else { // Assuming it will not be a very weird geometry
            if (SimplestGeometry == false) {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);

                CreateNewCondition(pThisProperties, rFaceGeometry, rConditionId, r_condition);
                condition_ids[0] = rConditionId;
                ++rCondCounter;

                rInterfacePart.AddConditions(condition_ids);
            } else {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(8);

                Triangle3D3< NodeType > tri_1(rFaceGeometry(0), rFaceGeometry(1), rFaceGeometry(8));
                CreateNewCondition(pThisProperties, tri_1, rConditionId, r_condition);
                condition_ids[0] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_2(rFaceGeometry(1), rFaceGeometry(2), rFaceGeometry(3));
                CreateNewCondition(pThisProperties, tri_2, rConditionId, r_condition);
                condition_ids[1] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_3(rFaceGeometry(1), rFaceGeometry(3), rFaceGeometry(8));
                CreateNewCondition(pThisProperties, tri_3, rConditionId, r_condition);
                condition_ids[2] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_4(rFaceGeometry(8), rFaceGeometry(3), rFaceGeometry(4));
                CreateNewCondition(pThisProperties, tri_4, rConditionId, r_condition);
                condition_ids[3] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_5(rFaceGeometry(8), rFaceGeometry(4), rFaceGeometry(5));
                CreateNewCondition(pThisProperties, tri_5, rConditionId, r_condition);
                condition_ids[4] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_6(rFaceGeometry(5), rFaceGeometry(6), rFaceGeometry(7));
                CreateNewCondition(pThisProperties, tri_6, rConditionId, r_condition);
                condition_ids[5] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_7(rFaceGeometry(5), rFaceGeometry(7), rFaceGeometry(8));
                CreateNewCondition(pThisProperties, tri_7, rConditionId, r_condition);
                condition_ids[6] = rConditionId;
                ++rCondCounter;
                ++rConditionId;
                Triangle3D3< NodeType > tri_8(rFaceGeometry(0), rFaceGeometry(8), rFaceGeometry(7));
                CreateNewCondition(pThisProperties, tri_8, rConditionId, r_condition);
                condition_ids[7] = rConditionId;
                ++rCondCounter;

                rInterfacePart.AddConditions(condition_ids);
            }
        }
    }

    KRATOS_CATCH("");
}

}  // namespace Kratos.

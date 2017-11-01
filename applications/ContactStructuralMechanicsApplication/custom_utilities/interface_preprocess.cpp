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
        "condition_name"                       : "", 
        "final_string"                         : "", 
        "simplify_geometry"                    : false
    })" );
    
    ThisParameters.ValidateAndAssignDefaults(default_parameters);
    
    const std::string& condition_name = ThisParameters["condition_name"].GetString();
    const std::string& final_string = ThisParameters["final_string"].GetString();
    const bool simplest_geometry = ThisParameters["simplify_geometry"].GetBool();
    
    unsigned int cond_counter = 0;
    
    // We reorder the conditions
    unsigned int cond_id = ReorderConditions();
    
    // Check if it is a mortar condition
    const bool is_mortar = (condition_name.find("Mortar") != std::string::npos) ? true : false;
    
    // Generate Conditions from original the edges that can be considered interface
    for (ModelPart::ElementsContainerType::const_iterator it_elem = rOriginPart.ElementsBegin(); it_elem != rOriginPart.ElementsEnd(); ++it_elem)
    {
        GeometryType& this_geometry = it_elem->GetGeometry();
        
        for (unsigned int i_edge = 0; i_edge < this_geometry.EdgesNumber(); ++i_edge)
        {
            GenerateEdgeCondition(rInterfacePart, *(it_elem.base()), this_geometry.Edges()[i_edge], condition_name, final_string, simplest_geometry, cond_counter, cond_id, is_mortar);
        }
    }
    
    
    // NOTE: Reorder ID if parallellization
    
    const unsigned int num_nodes = static_cast<int>(mrMainModelPart.Nodes().size());
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
        "condition_name"                       : "", 
        "final_string"                         : "", 
        "simplify_geometry"                    : false
    })" );
    
    ThisParameters.ValidateAndAssignDefaults(default_parameters);
    
    const std::string& condition_name = ThisParameters["condition_name"].GetString();
    const std::string& final_string = ThisParameters["final_string"].GetString();
    const bool simplest_geometry = ThisParameters["simplify_geometry"].GetBool();
    
    unsigned int cond_counter = 0;
    
    // We reorder the conditions
    unsigned int cond_id = ReorderConditions();

    // Check if it is a mortar condition
    const bool is_mortar = (condition_name.find("Mortar") != std::string::npos) ? true : false;
    
    // Generate Conditions from original the faces that can be considered interface
    for (ModelPart::ElementsContainerType::const_iterator it_elem = rOriginPart.ElementsBegin(); it_elem != rOriginPart.ElementsEnd(); ++it_elem)
    {          
        GeometryType& this_geometry = it_elem->GetGeometry();
        
        if (this_geometry.LocalSpaceDimension() == 3)
        {
            for (unsigned int i_face = 0; i_face < this_geometry.FacesNumber(); ++i_face)
            {
                GenerateFaceCondition(rInterfacePart, *(it_elem.base()), this_geometry.Faces()[i_face], condition_name, final_string, simplest_geometry, cond_counter, cond_id, is_mortar);
            }
        }
        else
        {
            GenerateFaceCondition(rInterfacePart, *(it_elem.base()), this_geometry, condition_name, final_string, simplest_geometry, cond_counter, cond_id, is_mortar);
        }
    }
    
    // NOTE: Reorder ID if parallellization
    
    const unsigned int num_nodes = static_cast<int>(mrMainModelPart.Nodes().size());
    PrintNodesAndConditions(num_nodes, cond_counter);
    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void InterfacePreprocessCondition::CreateNewCondition(
        Element::Pointer rpElem,
        Geometry<Node<3> > & rGeometry,
        const unsigned int CondId,
        const std::string& ConditionName,
        const bool IsMortar
        )
{
    KRATOS_TRY;

    Condition const& rCondition = KratosComponents<Condition>::Get(ConditionName); 
    Condition::Pointer p_cond = Condition::Pointer(rCondition.Create(CondId, rGeometry, rpElem->pGetProperties()));
    GeometryType& this_geometry = p_cond->GetGeometry();
    mrMainModelPart.AddCondition(p_cond);
    if (IsMortar == true)
    {
        p_cond->SetValue(ELEMENT_POINTER, rpElem);
        
        // We set the condition as master or slave (master by default)
        bool is_slave = true;
        for (unsigned int it_node = 0; it_node < this_geometry.size(); ++it_node)
        {
            if (this_geometry[it_node].Is(SLAVE) == false) is_slave = false;
        }
        if (is_slave == true)  p_cond->Set(SLAVE, true);
        else  p_cond->Set(MASTER, true);
    }
    
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void InterfacePreprocessCondition::PrintNodesAndConditions(
        const int& NodesCounter,
        const int& CondCounter
        )
{
    std::cout << "\t" << NodesCounter << " nodes ";
    std::cout << "and " << CondCounter <<  " conditions found." << std::endl;

    // Check that we actually found something
    if( NodesCounter == 0)
    {
        KRATOS_ERROR << "No interface nodes found. Please check that nodes on both sides of the interface have been assigned Is(INTERFACE) = true." << std::endl;
    }
    if( CondCounter == 0)
    {
        KRATOS_ERROR << "No interface conditions found. Please check that nodes on both sides of the interface have been assigned Is(INTERFACE) = true and that the contact surfaces have been assigned conditions." << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

unsigned int InterfacePreprocessCondition::ReorderConditions()
{
    // We reorder the conditions
    ConditionsArrayType& conditions_array = mrMainModelPart.Conditions();
    const unsigned int num_conditions = static_cast<int>(conditions_array.size());
    
    for(unsigned int i = 0; i < num_conditions; i++) 
    {
        auto it_condition = conditions_array.begin() + i;
        it_condition->SetId(i + 1);
    }
    
    return num_conditions;
}

/***********************************************************************************/
/***********************************************************************************/

inline void InterfacePreprocessCondition::GenerateEdgeCondition(
    ModelPart& rInterfacePart,
    Element::Pointer rpElem,
    GeometryType& EdgeGeometry,
    const std::string& ConditionName,
    const std::string& FinalString,
    const bool& SimplestGeometry,
    unsigned int& CondCounter,
    unsigned int& CondId,
    const bool IsMortar
    )
{
    unsigned int count = 0;
    const unsigned int number_of_points = EdgeGeometry.PointsNumber();
    for (unsigned int it_node = 0; it_node < number_of_points; ++it_node)
    {
        if (EdgeGeometry[it_node].IsDefined(INTERFACE) == true)  
        {
            if (EdgeGeometry[it_node].Is(INTERFACE) == true) ++count;
        }
    }
        
    if (count == number_of_points)
    {
        std::string EdgeConditionName = ConditionName;
        ++CondId; // NOTE: To paralellize be careful with this ID
        if (number_of_points == 2)
        {
            // We initialize a vector for the IDs
            std::vector<std::size_t> condition_ids(1);
    
            EdgeConditionName.append("Condition2D2N");
            EdgeConditionName.append(FinalString);
            
            CreateNewCondition(rpElem, EdgeGeometry, CondId, EdgeConditionName, IsMortar);
            condition_ids[0] = CondId;
            ++CondCounter;
            
            rInterfacePart.AddConditions(condition_ids);
        }
        else
        {                            
            if (SimplestGeometry == false)
            {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);
            
                EdgeConditionName.append("Condition2D3N"); 
                EdgeConditionName.append(FinalString); 
                
                CreateNewCondition(rpElem, EdgeGeometry, CondId, EdgeConditionName, IsMortar);
                condition_ids[0] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
            else
            {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(2);
            
                EdgeConditionName.append("Condition2D2N"); 
                EdgeConditionName.append(FinalString); 

                Line2D2< Node<3> > lin_1(EdgeGeometry(0), EdgeGeometry(1));
                CreateNewCondition(rpElem, lin_1, CondId, EdgeConditionName, IsMortar);
                condition_ids[0] = CondId;
                ++CondCounter;
                ++CondId;
                Line2D2< Node<3> > lin_2(EdgeGeometry(1), EdgeGeometry(2));
                CreateNewCondition(rpElem, lin_2, CondId, EdgeConditionName, IsMortar);
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
    Element::Pointer rpElem,
    GeometryType& FaceGeometry,
    const std::string& ConditionName,
    const std::string& FinalString,
    const bool& SimplestGeometry,
    unsigned int& CondCounter,
    unsigned int& CondId,
    const bool IsMortar
    )
{
    unsigned int count = 0;
    const unsigned int number_of_points = FaceGeometry.PointsNumber();
    for (unsigned int it_node = 0; it_node < number_of_points; ++it_node)
    {
        if (FaceGeometry[it_node].IsDefined(INTERFACE) == true)  
        {
            if (FaceGeometry[it_node].Is(INTERFACE) == true) ++count;
        }
    }
    
    if (count == number_of_points)
    {
        std::string face_condition_name = ConditionName;
        ++CondId;
        if (number_of_points == 3)
        {
            // We initialize a vector for the IDs
            std::vector<std::size_t> condition_ids(1);
                
            face_condition_name.append("Condition3D3N");
            face_condition_name.append(FinalString);
            
            CreateNewCondition(rpElem, FaceGeometry, CondId, face_condition_name, IsMortar);
            condition_ids[0] = CondId;
            ++CondCounter;
            
            rInterfacePart.AddConditions(condition_ids);
        }
        else if (number_of_points == 4)
        {
            if (SimplestGeometry == false)
            {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);
                
                face_condition_name.append("Condition3D4N");
                face_condition_name.append(FinalString);
            
                CreateNewCondition(rpElem, FaceGeometry, CondId, face_condition_name, IsMortar);
                condition_ids[0] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
            else
            {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(2);
                
                face_condition_name.append("Condition3D3N");
                face_condition_name.append(FinalString);
                
                Triangle3D3< Node<3> > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(2));
                CreateNewCondition(rpElem, tri_1, CondId, face_condition_name, IsMortar);
                condition_ids[0] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_2(FaceGeometry(2), FaceGeometry(3), FaceGeometry(0));
                CreateNewCondition(rpElem, tri_2, CondId, face_condition_name, IsMortar);
                condition_ids[1] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
        }
        else if (number_of_points == 6)
        {
            if (SimplestGeometry == false)
            {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);
                
                face_condition_name.append("Condition3D6N");
                face_condition_name.append(FinalString);
                
                CreateNewCondition(rpElem, FaceGeometry, CondId, face_condition_name, IsMortar);
                condition_ids[0] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
            else
            {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(4);
                
                face_condition_name.append("Condition3D3N");
                face_condition_name.append(FinalString);
                
                Triangle3D3< Node<3> > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(5));
                CreateNewCondition(rpElem, tri_1, CondId, face_condition_name, IsMortar);
                condition_ids[0] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_2(FaceGeometry(1), FaceGeometry(2), FaceGeometry(3));
                CreateNewCondition(rpElem, tri_2, CondId, face_condition_name, IsMortar);
                condition_ids[1] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_3(FaceGeometry(1), FaceGeometry(3), FaceGeometry(5));
                CreateNewCondition(rpElem, tri_3, CondId, face_condition_name, IsMortar);
                condition_ids[2] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_4(FaceGeometry(3), FaceGeometry(4), FaceGeometry(5));
                CreateNewCondition(rpElem, tri_4, CondId, face_condition_name, IsMortar);
                condition_ids[3] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
        }
        else if (number_of_points == 8)
        {
            if (SimplestGeometry == false)
            {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);
                
                face_condition_name.append("Condition3D8N");
                face_condition_name.append(FinalString);
                
                CreateNewCondition(rpElem, FaceGeometry, CondId, face_condition_name, IsMortar);
                condition_ids[0] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
            else
            {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(6);
                
                face_condition_name.append("Condition3D3N");
                face_condition_name.append(FinalString);

                Triangle3D3< Node<3> > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(7));
                CreateNewCondition(rpElem, tri_1, CondId, face_condition_name, IsMortar);
                condition_ids[0] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_2(FaceGeometry(1), FaceGeometry(5), FaceGeometry(7));
                CreateNewCondition(rpElem, tri_2, CondId, face_condition_name, IsMortar);
                condition_ids[1] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_3(FaceGeometry(1), FaceGeometry(3), FaceGeometry(5));
                CreateNewCondition(rpElem, tri_3, CondId, face_condition_name, IsMortar);
                condition_ids[2] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_4(FaceGeometry(1), FaceGeometry(2), FaceGeometry(3));
                CreateNewCondition(rpElem, tri_4, CondId, face_condition_name, IsMortar);
                condition_ids[3] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_5(FaceGeometry(3), FaceGeometry(4), FaceGeometry(5));
                CreateNewCondition(rpElem, tri_5, CondId, face_condition_name, IsMortar);
                condition_ids[4] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_6(FaceGeometry(5), FaceGeometry(6), FaceGeometry(7));
                CreateNewCondition(rpElem, tri_6, CondId, face_condition_name, IsMortar);
                condition_ids[5] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
        }
        else // Assuming it will not be a very weird geometry
        {
            if (SimplestGeometry == false)
            {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(1);
                
                face_condition_name.append("Condition3D4N");
                face_condition_name.append(FinalString);
                
                CreateNewCondition(rpElem, FaceGeometry, CondId, face_condition_name, IsMortar);
                condition_ids[0] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
            else
            {
                // We initialize a vector for the IDs
                std::vector<std::size_t> condition_ids(8);
                
                face_condition_name.append("Condition3D3N");
                face_condition_name.append(FinalString);

                Triangle3D3< Node<3> > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(8));
                CreateNewCondition(rpElem, tri_1, CondId, face_condition_name, IsMortar);
                condition_ids[0] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_2(FaceGeometry(1), FaceGeometry(2), FaceGeometry(3));
                CreateNewCondition(rpElem, tri_2, CondId, face_condition_name, IsMortar);
                condition_ids[1] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_3(FaceGeometry(1), FaceGeometry(3), FaceGeometry(8));
                CreateNewCondition(rpElem, tri_3, CondId, face_condition_name, IsMortar);
                condition_ids[2] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_4(FaceGeometry(8), FaceGeometry(3), FaceGeometry(4));
                CreateNewCondition(rpElem, tri_4, CondId, face_condition_name, IsMortar);
                condition_ids[3] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_5(FaceGeometry(8), FaceGeometry(4), FaceGeometry(5));
                CreateNewCondition(rpElem, tri_5, CondId, face_condition_name, IsMortar);
                condition_ids[4] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_6(FaceGeometry(5), FaceGeometry(6), FaceGeometry(7));
                CreateNewCondition(rpElem, tri_6, CondId, face_condition_name, IsMortar);
                condition_ids[5] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_7(FaceGeometry(5), FaceGeometry(7), FaceGeometry(8));
                CreateNewCondition(rpElem, tri_7, CondId, face_condition_name, IsMortar);
                condition_ids[6] = CondId;
                ++CondCounter;
                ++CondId;
                Triangle3D3< Node<3> > tri_8(FaceGeometry(0), FaceGeometry(8), FaceGeometry(7));
                CreateNewCondition(rpElem, tri_8, CondId, face_condition_name, IsMortar);
                condition_ids[7] = CondId;
                ++CondCounter;
                
                rInterfacePart.AddConditions(condition_ids);
            }
        }
    }
}
}  // namespace Kratos.

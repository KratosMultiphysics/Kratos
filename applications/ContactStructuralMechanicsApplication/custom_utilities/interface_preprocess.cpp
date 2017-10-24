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
        
        const std::string condition_name = ThisParameters["condition_name"].GetString();
        const std::string final_string = ThisParameters["final_string"].GetString();
        const bool simplest_geometry = ThisParameters["simplify_geometry"].GetBool();
        
        NodesArrayType& nodes_array = mrMainModelPart.Nodes();
        const unsigned int num_nodes = static_cast<int>(nodes_array.size());
        
        unsigned int cond_counter = 0;
        
        // We reorder the conditions
        unsigned int cond_id = ReorderConditions();
        
        // Generate Conditions from original the edges that can be considered interface
        for (ModelPart::ElementsContainerType::const_iterator it_elem = rOriginPart.ElementsBegin(); it_elem != rOriginPart.ElementsEnd(); ++it_elem)
        {
            GeometryType& this_geometry = it_elem->GetGeometry();
            
            for (unsigned int i_edge = 0; i_edge < this_geometry.EdgesNumber(); i_edge++)
            {
                unsigned int count = 0;
                const unsigned int number_of_points = this_geometry.Edges()[i_edge].PointsNumber();
                for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
                {
                    if (this_geometry.Edges()[i_edge][i_node].IsDefined(INTERFACE) == true)  
                    {
                        if (this_geometry.Edges()[i_edge][i_node].Is(INTERFACE) == true)  
                        {
                            count++;
                        }
                    }
                }
                
                if (count == number_of_points)
                {
                    std::string Edgecondition_name = condition_name;
                    ++cond_id; // NOTE: To paralellize be careful with this ID
                    if (number_of_points == 2)
                    {
                        Edgecondition_name.append("Condition2D2N");
                        Edgecondition_name.append(final_string);
                        
                        CreateNewCondition(rInterfacePart, *(it_elem.base()), this_geometry.Edges()[i_edge], cond_id, Edgecondition_name);
                        ++cond_counter;
                    }
                    else
                    {                            
                        if (simplest_geometry == false)
                        {
                            Edgecondition_name.append("Condition2D3N"); 
                            Edgecondition_name.append(final_string); 
                            
                            CreateNewCondition(rInterfacePart, *(it_elem.base()), this_geometry.Edges()[i_edge], cond_id, Edgecondition_name);
                            ++cond_counter;
                        }
                        else
                        {
                            Edgecondition_name.append("Condition2D2N"); 
                            Edgecondition_name.append(final_string); 

                            Line2D2< Node<3> > lin_1(this_geometry.Edges()[i_edge](0), this_geometry.Edges()[i_edge](1));
                            CreateNewCondition(rInterfacePart, *(it_elem.base()), lin_1, cond_id, Edgecondition_name);
                            ++cond_counter;
                            ++cond_id;
                            Line2D2< Node<3> > lin_2(this_geometry.Edges()[i_edge](1), this_geometry.Edges()[i_edge](2));
                            CreateNewCondition(rInterfacePart, *(it_elem.base()), lin_2, cond_id, Edgecondition_name);
                            ++cond_counter;
                        }
                    }
                }
            }
        }
       
      
        // NOTE: Reorder ID if parallellization
      
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
        const bool& simplest_geometry = ThisParameters["simplify_geometry"].GetBool();
        
        unsigned int cond_counter = 0;
        
        // We reorder the conditions
        unsigned int cond_id = ReorderConditions();

        // Generate Conditions from original the faces that can be considered interface
        for (ModelPart::ElementsContainerType::const_iterator it_elem = rOriginPart.ElementsBegin(); it_elem != rOriginPart.ElementsEnd(); ++it_elem)
        {          
            GeometryType& this_geometry = it_elem->GetGeometry();
            
            if (this_geometry.LocalSpaceDimension() == 3)
            {
                for (unsigned int i_face = 0; i_face < this_geometry.FacesNumber(); i_face++)
                {
                    GenerateFaceCondition(rInterfacePart, *(it_elem.base()), this_geometry.Faces()[i_face], condition_name, final_string, simplest_geometry, cond_counter, cond_id);
                }
            }
            else
            {
                GenerateFaceCondition(rInterfacePart, *(it_elem.base()), this_geometry, condition_name, final_string, simplest_geometry, cond_counter, cond_id);
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
            ModelPart& rInterfacePart,
            Element::Pointer rpElem,
            Geometry<Node<3> > & rGeometry,
            const unsigned int CondId,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;
    
        Condition const & rCondition = KratosComponents<Condition>::Get(ConditionName); 
        Condition::Pointer p_cond = Condition::Pointer(rCondition.Create(CondId, rGeometry, rpElem->pGetProperties()));
        rInterfacePart.AddCondition(p_cond);
        GeometryType& this_geometry = p_cond->GetGeometry();
        if (ConditionName.find("Mortar") != std::string::npos)
        {
            p_cond->GetValue(ELEMENT_POINTER) = rpElem;
            
            // We set the condition as master or slave (master by default)
            if (ConditionName.find("ALM") != std::string::npos || ConditionName.find("MeshTying") != std::string::npos)
            {
                bool is_slave = true;
                for (unsigned int i_node = 0; i_node < this_geometry.size(); ++i_node)
                {
                    if (this_geometry[i_node].Is(SLAVE) == false) is_slave = false;
                }
                if (is_slave == true)  p_cond->Set(SLAVE, true);
                else  p_cond->Set(MASTER, true);
            }
        }
        
        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    void InterfacePreprocessCondition::GenerateLine2NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;

        const unsigned int& dimension = mrMainModelPart.GetProcessInfo()[DOMAIN_SIZE];
                    
        // Store pointers to all interface nodes
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); ++i_node)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                ++nodes_counter;
            }
        }

        // Generate linear Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); ++it_cond)
        {
            if (
                ((*it_cond).GetGeometry()[0].Is(INTERFACE) == true) &&
                ((*it_cond).GetGeometry()[1].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                ++cond_counter;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

        if (dimension == 2) // By default, but someone can be interested in project values to a BEAM for example
        {
            GenerateLine2D2NConditions(aux, InterfacePart.Conditions(), ConditionName);
        }
        else
        {
            GenerateLine3D2NConditions(aux, InterfacePart.Conditions(), ConditionName);
        }

        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    void InterfacePreprocessCondition::GenerateLine3NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;
        
        const unsigned int& dimension = mrMainModelPart.GetProcessInfo()[DOMAIN_SIZE];
        
        // Store pointers to all interface nodes
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); ++i_node)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                ++nodes_counter;
            }
        }
        
        // Generate linear Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); ++it_cond)
        {
            GeometryType& this_geometry = (*it_cond).GetGeometry();
            if (
                (this_geometry[0].Is(INTERFACE) == true) &&
                (this_geometry[1].Is(INTERFACE) == true) &&
                (this_geometry[2].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                ++cond_counter;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

        if (dimension == 2) // By default, but someone can be interested in project values to a BEAM for example
        {
            GenerateLine2D3NConditions(aux, InterfacePart.Conditions(), ConditionName);
        }
        else
        {
            GenerateLine3D3NConditions(aux, InterfacePart.Conditions(), ConditionName);
        }

        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    void InterfacePreprocessCondition::GenerateTriangle3NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;
        
        // Store pointers to all interface nodes
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); ++i_node)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                ++nodes_counter;
            }
        }

        // Generate triangular Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); ++it_cond)
        {
            GeometryType& this_geometry = (*it_cond).GetGeometry();
            if (
                (this_geometry[0].Is(INTERFACE) == true) &&
                (this_geometry[1].Is(INTERFACE) == true) &&
                (this_geometry[2].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                ++cond_counter;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

        GenerateTriangular3D3NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void InterfacePreprocessCondition::GenerateTriangle6NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;

        // Store pointers to all interface nodes
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); ++i_node)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                ++nodes_counter;
            }
        }

        // Generate triangular Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); ++it_cond)
        {
            GeometryType& this_geometry = (*it_cond).GetGeometry();
            if (
                (this_geometry[0].Is(INTERFACE) == true) &&
                (this_geometry[1].Is(INTERFACE) == true) &&
                (this_geometry[2].Is(INTERFACE) == true) &&
                (this_geometry[3].Is(INTERFACE) == true) &&
                (this_geometry[4].Is(INTERFACE) == true) &&
                (this_geometry[5].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                ++cond_counter;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

        GenerateTriangular3D6NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    void InterfacePreprocessCondition::GenerateQuadrilateral4NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;

        // Store pointers to all interface nodes
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); ++i_node)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                ++nodes_counter;
            }
        }

        // Generate quadrilateral Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); ++it_cond)
        {
            GeometryType& this_geometry = (*it_cond).GetGeometry();
            if (
                (this_geometry[0].Is(INTERFACE) == true) &&
                (this_geometry[1].Is(INTERFACE) == true) &&
                (this_geometry[2].Is(INTERFACE) == true) &&
                (this_geometry[3].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                ++cond_counter;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

        GenerateQuadrilateral3D4NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    void InterfacePreprocessCondition::GenerateQuadrilateral8NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;

        // Store pointers to all interface nodes
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); ++i_node)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                ++nodes_counter;
            }
        }

        // Generate quadrilateral Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); ++it_cond)
        {
            GeometryType& this_geometry = (*it_cond).GetGeometry();
            if (
                (this_geometry[0].Is(INTERFACE) == true) &&
                (this_geometry[1].Is(INTERFACE) == true) &&
                (this_geometry[2].Is(INTERFACE) == true) &&
                (this_geometry[3].Is(INTERFACE) == true) &&
                (this_geometry[4].Is(INTERFACE) == true) &&
                (this_geometry[5].Is(INTERFACE) == true) &&
                (this_geometry[6].Is(INTERFACE) == true) &&
                (this_geometry[7].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                ++cond_counter;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

        GenerateQuadrilateral3D8NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    void InterfacePreprocessCondition::GenerateQuadrilateral9NInterfacePart(
            const ModelPart& rOriginPart,
            ModelPart& InterfacePart,
            const std::string ConditionName
            )
    {
        KRATOS_TRY;

        // Store pointers to all interface nodes
        unsigned int nodes_counter = 0;
        for (ModelPart::NodesContainerType::const_iterator i_node = rOriginPart.NodesBegin(); i_node != rOriginPart.NodesEnd(); ++i_node)
        {
            if (i_node->Is(INTERFACE) == true)
            {
                InterfacePart.Nodes().push_back( *(i_node.base()) );
                ++nodes_counter;
            }
        }

        // Generate quadrilateral Conditions from original interface conditions
        ModelPart::ConditionsContainerType aux;
        unsigned int cond_counter = 0;

        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginPart.ConditionsBegin(); it_cond != rOriginPart.ConditionsEnd(); ++it_cond)
        {
            GeometryType& this_geometry = (*it_cond).GetGeometry();
            if (
                (this_geometry[0].Is(INTERFACE) == true) &&
                (this_geometry[1].Is(INTERFACE) == true) &&
                (this_geometry[2].Is(INTERFACE) == true) &&
                (this_geometry[3].Is(INTERFACE) == true) &&
                (this_geometry[4].Is(INTERFACE) == true) &&
                (this_geometry[5].Is(INTERFACE) == true) &&
                (this_geometry[6].Is(INTERFACE) == true) &&
                (this_geometry[7].Is(INTERFACE) == true) &&
                (this_geometry[8].Is(INTERFACE) == true))
            {
                aux.push_back( *(it_cond.base()) );
                ++cond_counter;
            }
        }

        PrintNodesAndConditions(nodes_counter, cond_counter);

        GenerateQuadrilateral3D9NConditions(aux, InterfacePart.Conditions(), ConditionName);

        KRATOS_CATCH("");
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    void InterfacePreprocessCondition::GenerateLine2D2NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition2D2N"); 

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType lin_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginConds.begin(); it_cond != rOriginConds.end(); ++it_cond)
        {
            GeometryType& this_geometry = it_cond->GetGeometry();
            boost::shared_ptr<Kratos::Properties> this_properties = it_cond->pGetProperties();
            if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D2)
            {
                rLinConds.push_back( rCondition.Create(++lin_id, this_geometry, this_properties));
            }
            else if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D3)
            {
                Line2D2< Node<3> > lin_1(this_geometry(0), this_geometry(1));
                Line2D2< Node<3> > lin_2(this_geometry(1), this_geometry(2));

                rLinConds.push_back(rCondition.Create(++lin_id, lin_1, this_properties));
                rLinConds.push_back(rCondition.Create(++lin_id, lin_2, this_properties));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using linear lines " << std::endl;
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    void InterfacePreprocessCondition::GenerateLine3D2NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D2N");

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType lin_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginConds.begin(); it_cond != rOriginConds.end(); ++it_cond)
        {
            GeometryType& this_geometry = it_cond->GetGeometry();
            boost::shared_ptr<Kratos::Properties> this_properties = it_cond->pGetProperties();
            if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line3D2)
            {
                rLinConds.push_back( rCondition.Create(++lin_id, this_geometry, this_properties));
            }
            else if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line3D3)
            {
                Line3D2< Node<3> > lin_1(this_geometry(0), this_geometry(1));
                Line3D2< Node<3> > lin_2(this_geometry(1), this_geometry(2));

                rLinConds.push_back(rCondition.Create(++lin_id, lin_1, this_properties));
                rLinConds.push_back(rCondition.Create(++lin_id, lin_2, this_properties));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using linear lines " << std::endl;
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void InterfacePreprocessCondition::GenerateLine2D3NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition2D3N");

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType lin_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginConds.begin(); it_cond != rOriginConds.end(); ++it_cond)
        {
            GeometryType& this_geometry = it_cond->GetGeometry();
            boost::shared_ptr<Kratos::Properties> this_properties = it_cond->pGetProperties();
            if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D3)
            {
                rLinConds.push_back( rCondition.Create(++lin_id, this_geometry, this_properties));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using linear lines " << std::endl;
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void InterfacePreprocessCondition::GenerateLine3D3NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rLinConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D3N");

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType lin_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginConds.begin(); it_cond != rOriginConds.end(); ++it_cond)
        {
            GeometryType& this_geometry = it_cond->GetGeometry();
            boost::shared_ptr<Kratos::Properties> this_properties = it_cond->pGetProperties();
            if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line3D3)
            {
                rLinConds.push_back( rCondition.Create(++lin_id, this_geometry, this_properties));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using linear lines " << std::endl;
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    void InterfacePreprocessCondition::GenerateTriangular3D3NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rTriConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D"); 

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType tri_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginConds.begin(); it_cond != rOriginConds.end(); ++it_cond)
        {
            GeometryType& this_geometry = it_cond->GetGeometry();
            boost::shared_ptr<Kratos::Properties> this_properties = it_cond->pGetProperties();
            if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D3)
            {
                rTriConds.push_back( rCondition.Create(++tri_id, this_geometry, this_properties));
            }
            else if (this_geometry.PointsNumber() == 4)
            {
                Triangle3D3< Node<3> > tri_1(this_geometry(0), this_geometry(1), this_geometry(2));
                Triangle3D3< Node<3> > tri_2(this_geometry(2), this_geometry(3), this_geometry(0));

                rTriConds.push_back(rCondition.Create(++tri_id, tri_1, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_2, this_properties));
            }
            else if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D6)
            {
                Triangle3D3< Node<3> > tri_1(this_geometry(0), this_geometry(1), this_geometry(5));
                Triangle3D3< Node<3> > tri_2(this_geometry(1), this_geometry(2), this_geometry(3));
                Triangle3D3< Node<3> > tri_3(this_geometry(1), this_geometry(3), this_geometry(5));
                Triangle3D3< Node<3> > tri_4(this_geometry(3), this_geometry(4), this_geometry(5));

                rTriConds.push_back(rCondition.Create(++tri_id, tri_1, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_2, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_3, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_4, this_properties));
            }
            else if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8)
            {
                Triangle3D3< Node<3> > tri_1(this_geometry(0), this_geometry(1), this_geometry(7));
                Triangle3D3< Node<3> > tri_2(this_geometry(1), this_geometry(5), this_geometry(7));
                Triangle3D3< Node<3> > tri_3(this_geometry(1), this_geometry(3), this_geometry(5));
                Triangle3D3< Node<3> > tri_4(this_geometry(1), this_geometry(2), this_geometry(3));
                Triangle3D3< Node<3> > tri_5(this_geometry(3), this_geometry(4), this_geometry(5));
                Triangle3D3< Node<3> > tri_6(this_geometry(5), this_geometry(6), this_geometry(7));

                rTriConds.push_back(rCondition.Create(++tri_id, tri_1, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_2, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_3, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_4, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_5, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_6, this_properties));
            }
            else if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9)
            {
                Triangle3D3< Node<3> > tri_1(this_geometry(0), this_geometry(1), this_geometry(8));
                Triangle3D3< Node<3> > tri_2(this_geometry(1), this_geometry(2), this_geometry(3));
                Triangle3D3< Node<3> > tri_3(this_geometry(1), this_geometry(3), this_geometry(8));
                Triangle3D3< Node<3> > tri_4(this_geometry(8), this_geometry(3), this_geometry(4));
                Triangle3D3< Node<3> > tri_5(this_geometry(8), this_geometry(4), this_geometry(5));
                Triangle3D3< Node<3> > tri_6(this_geometry(5), this_geometry(6), this_geometry(7));
                Triangle3D3< Node<3> > tri_7(this_geometry(5), this_geometry(7), this_geometry(8));
                Triangle3D3< Node<3> > tri_8(this_geometry(0), this_geometry(8), this_geometry(7));

                rTriConds.push_back(rCondition.Create(++tri_id, tri_1, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_2, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_3, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_4, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_5, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_6, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_7, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_8, this_properties));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using linear triangles " << std::endl;
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    void InterfacePreprocessCondition::GenerateTriangular3D6NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rTriConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D6N"); 

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType tri_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginConds.begin(); it_cond != rOriginConds.end(); ++it_cond)
        {
            GeometryType& this_geometry = it_cond->GetGeometry();
            boost::shared_ptr<Kratos::Properties> this_properties = it_cond->pGetProperties();
            if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle3D6)
            {
                rTriConds.push_back( rCondition.Create(++tri_id, this_geometry, this_properties));
            }
            else if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9)
            {
                Triangle3D6< Node<3> > tri_1(this_geometry(0), this_geometry(1), this_geometry(3), this_geometry(4), this_geometry(8), this_geometry(7));
                Triangle3D6< Node<3> > tri_2(this_geometry(2), this_geometry(3), this_geometry(1), this_geometry(6), this_geometry(8), this_geometry(5));

                rTriConds.push_back(rCondition.Create(++tri_id, tri_1, this_properties));
                rTriConds.push_back(rCondition.Create(++tri_id, tri_2, this_properties));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using quadratic triangles " << std::endl;
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void InterfacePreprocessCondition::GenerateQuadrilateral3D4NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rQuadConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D4N"); 

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType quad_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginConds.begin(); it_cond != rOriginConds.end(); ++it_cond)
        {
            GeometryType& this_geometry = it_cond->GetGeometry();
            boost::shared_ptr<Kratos::Properties> this_properties = it_cond->pGetProperties();
            if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4)
            {
                rQuadConds.push_back( rCondition.Create(++quad_id, this_geometry, this_properties));
            }
            else if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9)
            {
                Quadrilateral3D4< Node<3> > quad_1(this_geometry(0), this_geometry(4), this_geometry(8), this_geometry(7));
                Quadrilateral3D4< Node<3> > quad_2(this_geometry(4), this_geometry(1), this_geometry(5), this_geometry(8));
                Quadrilateral3D4< Node<3> > quad_3(this_geometry(8), this_geometry(5), this_geometry(2), this_geometry(6));
                Quadrilateral3D4< Node<3> > quad_4(this_geometry(7), this_geometry(8), this_geometry(6), this_geometry(3));

                rQuadConds.push_back(rCondition.Create(++quad_id, quad_1, this_properties));
                rQuadConds.push_back(rCondition.Create(++quad_id, quad_2, this_properties));
                rQuadConds.push_back(rCondition.Create(++quad_id, quad_3, this_properties));
                rQuadConds.push_back(rCondition.Create(++quad_id, quad_4, this_properties));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using bilinear quadrilaterals " << std::endl;
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    void InterfacePreprocessCondition::GenerateQuadrilateral3D8NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rQuadConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D8N"); 

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType quad_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it_cond = rOriginConds.begin(); it_cond != rOriginConds.end(); ++it_cond)
        {
            if (it_cond->GetGeometry().GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8)
            {
                rQuadConds.push_back( rCondition.Create(++quad_id, it_cond->GetGeometry(), it_cond->pGetProperties()));
            }
            else
	    {
                KRATOS_ERROR << "The geometry can not be divided using quadratic (8 nodes) quadrilaterals " << std::endl;
	    }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    void InterfacePreprocessCondition::GenerateQuadrilateral3D9NConditions(
            const ModelPart::ConditionsContainerType& rOriginConds,
            ModelPart::ConditionsContainerType& rQuadConds,
            const std::string ConditionName
            )
    {
        // Define a condition to use as reference for all new triangle conditions
        const Condition& rCondition = KratosComponents<Condition>::Get(ConditionName); // The custom condition will be considered
//         const Condition& rCondition = KratosComponents<Condition>::Get("Condition3D9N"); 

        // Required information for new conditions: Id, geometry and properties
        Condition::IndexType quad_id = 1; // Id

        // Loop over origin conditions and create a set of triangular ones
        for (ModelPart::ConditionsContainerType::const_iterator it = rOriginConds.begin(); it != rOriginConds.end(); ++it)
        {
            if (it->GetGeometry().PointsNumber() == 9)
            {
                rQuadConds.push_back( rCondition.Create(++quad_id, it->GetGeometry(), it->pGetProperties()));
            }
            else
            {
                KRATOS_ERROR << "The geometry can not be divided using quadratic (9 nodes) quadrilaterals " << std::endl;
            }
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    void InterfacePreprocessCondition::PrintNodesAndConditions(
            const int NodesCounter,
            const int CondCounter
            )
    {
        std::cout << "    " << NodesCounter << " nodes ";
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

    inline void InterfacePreprocessCondition::GenerateFaceCondition(
        ModelPart& rInterfacePart,
        Element::Pointer rpElem,
        GeometryType& FaceGeometry,
        const std::string& ConditionName,
        const std::string& FinalString,
        const bool& SimplestGeometry,
        unsigned int& CondCounter,
        unsigned int& CondId
        )
    {
        unsigned int count = 0;
        const unsigned int number_of_points = FaceGeometry.PointsNumber();
        for (unsigned int i_node = 0; i_node < number_of_points; ++i_node)
        {
            if (FaceGeometry[i_node].IsDefined(INTERFACE) == true)  
            {
                if (FaceGeometry[i_node].Is(INTERFACE) == true) ++count;
            }
        }
        
        if (count == number_of_points)
        {
            std::string face_condition_name = ConditionName;
            ++CondId;
            if (number_of_points == 3)
            {
                face_condition_name.append("Condition3D3N");
                face_condition_name.append(FinalString);
                
                CreateNewCondition(rInterfacePart, rpElem, FaceGeometry, CondId, face_condition_name);
                ++CondCounter;
            }
            else if (number_of_points == 4)
            {
                if (SimplestGeometry == false)
                {
                    face_condition_name.append("Condition3D4N");
                    face_condition_name.append(FinalString);
                
                    CreateNewCondition(rInterfacePart, rpElem, FaceGeometry, CondId, face_condition_name);
                    ++CondCounter;
                }
                else
                {
                    face_condition_name.append("Condition3D3N");
                    face_condition_name.append(FinalString);
                    
                    Triangle3D3< Node<3> > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(2));
                    CreateNewCondition(rInterfacePart, rpElem, tri_1, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_2(FaceGeometry(2), FaceGeometry(3), FaceGeometry(0));
                    CreateNewCondition(rInterfacePart, rpElem, tri_2, CondId, face_condition_name);
                    ++CondCounter;
                }
            }
            else if (number_of_points == 6)
            {
                if (SimplestGeometry == false)
                {
                    face_condition_name.append("Condition3D6N");
                    face_condition_name.append(FinalString);
                    
                    CreateNewCondition(rInterfacePart, rpElem, FaceGeometry, CondId, face_condition_name);
                    ++CondCounter;
                }
                else
                {
                    face_condition_name.append("Condition3D3N");
                    face_condition_name.append(FinalString);
                    
                    Triangle3D3< Node<3> > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(5));
                    CreateNewCondition(rInterfacePart, rpElem, tri_1, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_2(FaceGeometry(1), FaceGeometry(2), FaceGeometry(3));
                    CreateNewCondition(rInterfacePart, rpElem, tri_2, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_3(FaceGeometry(1), FaceGeometry(3), FaceGeometry(5));
                    CreateNewCondition(rInterfacePart, rpElem, tri_3, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_4(FaceGeometry(3), FaceGeometry(4), FaceGeometry(5));
                    CreateNewCondition(rInterfacePart, rpElem, tri_4, CondId, face_condition_name);
                    ++CondCounter;
                }
            }
            else if (number_of_points == 8)
            {
                if (SimplestGeometry == false)
                {
                    face_condition_name.append("Condition3D8N");
                    face_condition_name.append(FinalString);
                    
                    CreateNewCondition(rInterfacePart, rpElem, FaceGeometry, CondId, face_condition_name);
                    ++CondCounter;
                }
                else
                {
                    face_condition_name.append("Condition3D3N");
                    face_condition_name.append(FinalString);

                    Triangle3D3< Node<3> > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(7));
                    CreateNewCondition(rInterfacePart, rpElem, tri_1, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_2(FaceGeometry(1), FaceGeometry(5), FaceGeometry(7));
                    CreateNewCondition(rInterfacePart, rpElem, tri_2, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_3(FaceGeometry(1), FaceGeometry(3), FaceGeometry(5));
                    CreateNewCondition(rInterfacePart, rpElem, tri_3, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_4(FaceGeometry(1), FaceGeometry(2), FaceGeometry(3));
                    CreateNewCondition(rInterfacePart, rpElem, tri_4, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_5(FaceGeometry(3), FaceGeometry(4), FaceGeometry(5));
                    CreateNewCondition(rInterfacePart, rpElem, tri_5, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_6(FaceGeometry(5), FaceGeometry(6), FaceGeometry(7));
                    CreateNewCondition(rInterfacePart, rpElem, tri_6, CondId, face_condition_name);
                    ++CondCounter;
                }
            }
            else // Assuming it will not be a very weird geometry
            {
                if (SimplestGeometry == false)
                {
                    face_condition_name.append("Condition3D4N");
                    face_condition_name.append(FinalString);
                    
                    CreateNewCondition(rInterfacePart, rpElem, FaceGeometry, CondId, face_condition_name);
                    ++CondCounter;
                }
                else
                {
                    face_condition_name.append("Condition3D3N");
                    face_condition_name.append(FinalString);

                    Triangle3D3< Node<3> > tri_1(FaceGeometry(0), FaceGeometry(1), FaceGeometry(8));
                    CreateNewCondition(rInterfacePart, rpElem, tri_1, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_2(FaceGeometry(1), FaceGeometry(2), FaceGeometry(3));
                    CreateNewCondition(rInterfacePart, rpElem, tri_2, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_3(FaceGeometry(1), FaceGeometry(3), FaceGeometry(8));
                    CreateNewCondition(rInterfacePart, rpElem, tri_3, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_4(FaceGeometry(8), FaceGeometry(3), FaceGeometry(4));
                    CreateNewCondition(rInterfacePart, rpElem, tri_4, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_5(FaceGeometry(8), FaceGeometry(4), FaceGeometry(5));
                    CreateNewCondition(rInterfacePart, rpElem, tri_5, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_6(FaceGeometry(5), FaceGeometry(6), FaceGeometry(7));
                    CreateNewCondition(rInterfacePart, rpElem, tri_6, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_7(FaceGeometry(5), FaceGeometry(7), FaceGeometry(8));
                    CreateNewCondition(rInterfacePart, rpElem, tri_7, CondId, face_condition_name);
                    ++CondCounter;
                    ++CondId;
                    Triangle3D3< Node<3> > tri_8(FaceGeometry(0), FaceGeometry(8), FaceGeometry(7));
                    CreateNewCondition(rInterfacePart, rpElem, tri_8, CondId, face_condition_name);
                    ++CondCounter;
                }
            }
        }
    }
}

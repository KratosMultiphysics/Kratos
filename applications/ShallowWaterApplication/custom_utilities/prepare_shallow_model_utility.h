//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#if !defined(KRATOS_PREPARE_SHALLOW_MODEL_UTILITY_H_INCLUDED)
#define KRATOS_PREPARE_SHALLOW_MODEL_UTILITY_H_INCLUDED


// System includes
#include <iostream>


// External includes


// Project includes
#include "includes/model_part.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"
#include "shallow_water_application.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Utilities for the prepare model process
/** This classes implement some utilities for the prepare model process
*/
class PrepareShallowModelUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// The type used for the node
    typedef Node<3> NodeType;

    // General containers type definitions
    typedef ModelPart::NodesContainerType                                     NodesArrayType;
    typedef ModelPart::ConditionsContainerType                           ConditionsArrayType;
    typedef ModelPart::ElementsContainerType                               ElementsArrayType;
    typedef ModelPart::MasterSlaveConstraintContainerType     MasterSlaveConstraintArrayType;

    // General containers iterators type definitions
    typedef NodesArrayType::iterator                                  IteratorNodesArrayType;
    typedef ConditionsArrayType::iterator                        IteratorConditionsArrayType;
    typedef ElementsArrayType::iterator                            IteratorElementsArrayType;
    typedef MasterSlaveConstraintArrayType::iterator IteratorMasterSlaveConstraintsArrayType;

    /// The type used to identify the size
    typedef std::size_t IndexType;

    /// The Types used by the AssignUniqueModelPartCollectionTag
    typedef std::unordered_map<IndexType, IndexType> IndexIndexMapType;
    typedef std::unordered_map<IndexType, std::vector<std::string>> IndexStringMapType;
    typedef std::unordered_map<IndexType, std::vector<IndexType>> IndexVectorMapType;

    /// Pointer definition of PrepareShallowModelUtility
    KRATOS_CLASS_POINTER_DEFINITION(PrepareShallowModelUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PrepareShallowModelUtility(ModelPart& rTopographicModelPart, ModelPart& rDestinationModelPart)
      : mrTopographicModelPart(rTopographicModelPart)
      , mrDestinationModelPart(rDestinationModelPart)
    {}

    /// Destructor.
    virtual ~PrepareShallowModelUtility(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void CopyTopographicModelPart()
    {
        // Get the map of the original model part
        IndexStringMapType collections;
        IndexIndexMapType nodes_tags, elems_tags, conds_tags;
        AssignUniqueModelPartCollectionTagUtility model_part_collection(mrTopographicModelPart);
        model_part_collection.ComputeTags(nodes_tags, conds_tags, elems_tags, collections);

        // Clone the nodes
        VariablesList* p_variables_list = &mrDestinationModelPart.GetNodalSolutionStepVariablesList();
        const IndexType num_nodes = mrTopographicModelPart.Nodes().size();
        const auto it_node_begin = mrTopographicModelPart.NodesBegin();
        for (IndexType i = 0; i < num_nodes; i++)
        {
            auto it_node = it_node_begin + i;
            NodeType::Pointer p_new_node = it_node->Clone();
            mrDestinationModelPart.AddNode(p_new_node);
            p_new_node->SetSolutionStepVariablesList(p_variables_list);
        }

        // Clone the elements
        const IndexType num_elements = mrTopographicModelPart.Elements().size();
        const auto it_elem_begin = mrTopographicModelPart.ElementsBegin();
        for (IndexType i = 0; i < num_elements; i++)
        {
            auto it_elem = it_elem_begin + i;
            // Build the geometry with the new nodes
            PointerVector<NodeType> nodes_array;
            for (NodeType& node : it_elem->GetGeometry().Points())
            {
                nodes_array.push_back(mrDestinationModelPart.pGetNode(node.Id()));
            }
            Element::Pointer p_new_element = it_elem->Clone(it_elem->Id(), nodes_array);
            mrDestinationModelPart.AddElement(p_new_element);
        }

        // Clone the conditions
        const IndexType num_conditions = mrTopographicModelPart.Conditions().size();
        const auto it_cond_begin = mrTopographicModelPart.ConditionsBegin();
        for (IndexType i = 0; i < num_conditions; i++)
        {
            auto it_cond = it_cond_begin + i;
            // Build the geometry with the new nodes
            PointerVector<NodeType> nodes_array;
            for (NodeType& node : it_cond->GetGeometry().Points())
            {
                nodes_array.push_back(mrDestinationModelPart.pGetNode(node.Id()));
            }
            Condition::Pointer p_new_condition = it_cond->Clone(it_cond->Id(), nodes_array);
            mrDestinationModelPart.AddCondition(p_new_condition);
        }

        // Invert the tags
        IndexVectorMapType nodes_collection, elems_collection, conds_collection;
        for (auto node_tag : nodes_tags)
        {
            nodes_collection[node_tag.second].push_back(node_tag.first);
        }
        for (auto elem_tag : elems_tags)
        {
            elems_collection[elem_tag.second].push_back(elem_tag.first);
        }
        for (auto cond_tag : conds_tags)
        {
            conds_collection[cond_tag.second].push_back(cond_tag.first);
        }

        // And populate the sub model parts
        // Nodes
        for (auto collection : nodes_collection)
        {
            IndexType tag = collection.first;
            if (tag != 0)
            {
                for (auto name : collections[tag])
                {
                    ModelPart& part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(mrDestinationModelPart, name);
                    part.AddNodes(collection.second);
                }
            }
        }

        // Elements
        for (auto collection : elems_collection)
        {
            IndexType tag = collection.first;
            if (tag != 0)
            {
                for (auto name : collections[tag])
                {
                    ModelPart& part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(mrDestinationModelPart, name);
                    part.AddElements(collection.second);
                }
            }
        }

        // Elements
        for (auto collection : conds_collection)
        {
            IndexType tag = collection.first;
            if (tag != 0)
            {
                for (auto name : collections[tag])
                {
                    ModelPart& part = AssignUniqueModelPartCollectionTagUtility::GetRecursiveSubModelPart(mrDestinationModelPart, name);
                    part.AddConditions(collection.second);
                }
            }
        }


        // Finally, we need to increase the ids and duplicate the properties, in order to have different meshes in the post-process
        // TODO: We should renumber the computing model part id, but the move particles utility will renumber the ids again. In the future the move particles utility should not renumber ids
        // NOTE: this is very dirty but it works in GiD
        Properties::Pointer topographic_property = mrTopographicModelPart.CreateNewProperties(100);
        IndexType node_id = num_nodes;
        IndexType elem_id = num_elements;
        IndexType cond_id = num_conditions;
        // Nodes
        const IndexType num_computing_nodes = mrTopographicModelPart.Nodes().size();
        const auto it_computing_node_begin = mrTopographicModelPart.NodesBegin();
        for (IndexType i = 0; i < num_computing_nodes; i++)
        {
            auto it_node = it_computing_node_begin + i;
            it_node->SetId(++node_id);
        }
        // Elements
        const IndexType num_computing_elements = mrTopographicModelPart.Elements().size();
        const auto it_computing_element_begin = mrTopographicModelPart.ElementsBegin();
        for (IndexType i = 0; i < num_computing_elements; i++)
        {
            auto it_elem = it_computing_element_begin + i;
            it_elem->SetId(++elem_id);
            it_elem->SetProperties(topographic_property);
        }
        // Conditions
        const IndexType num_computing_conditions = mrTopographicModelPart.Conditions().size();
        const auto it_computing_condition_begin = mrTopographicModelPart.ConditionsBegin();
        for (IndexType i = 0; i < num_computing_conditions; i++)
        {
            auto it_cond = it_computing_condition_begin + i;
            it_cond->SetId(++cond_id);
            it_cond->SetProperties(topographic_property);
        }

    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    // virtual std::string Info() const;

    // /// Print information about this object.
    // virtual void PrintInfo(std::ostream& rOStream) const;

    // /// Print object's data.
    // virtual void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrTopographicModelPart;
    ModelPart& mrDestinationModelPart;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    // PrepareShallowModelUtility& operator=(PrepareShallowModelUtility const& rOther);

    /// Copy constructor.
    PrepareShallowModelUtility(PrepareShallowModelUtility const& rOther);


    ///@}

}; // Class PrepareShallowModelUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                 PrepareShallowModelUtility& rThis);

// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                 const PrepareShallowModelUtility& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PREPARE_SHALLOW_MODEL_UTILITY_H_INCLUDED defined

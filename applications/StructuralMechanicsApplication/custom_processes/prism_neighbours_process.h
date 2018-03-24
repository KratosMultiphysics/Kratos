// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_PRISM_NEIGHBOURS_PROCESS_H_INCLUDED )
#define  KRATOS_PRISM_NEIGHBOURS_PROCESS_H_INCLUDED

// System includes

// External includes
#include <unordered_map>

// Project includes
#include "processes/process.h"
#include "includes/key_hash.h"
#include "includes/model_part.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

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

/**
 * @class PrismNeighboursProcess
 * @ingroup StructuralMechanicsApplication
 * @brief An algorithm that looks for neighbour nodes and elements in a mesh of prismatic elements
 * @details For that pourpose if builds an unordered map of the surrounding elements and nodes and performs different checks
 * @author Vicente Mataix Ferrandiz
 * @todo Remove the dependence of the boost::unordered_map
*/
class PrismNeighboursProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PrismNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(PrismNeighboursProcess);

    // General geometry type definitions
    typedef Node<3>                                          NodeType;
    typedef Geometry<NodeType>                           GeometryType;

    // Containers definition
    typedef ModelPart::NodesContainerType              NodesArrayType;
    typedef ModelPart::ConditionsContainerType    ConditionsArrayType;
    typedef ModelPart::ElementsContainerType        ElementsArrayType;

    // Containers iterators definition
    typedef NodesArrayType::iterator                NodesIterarorType;
    typedef ConditionsArrayType::iterator      ConditionsIteratorType;
    typedef ElementsArrayType::iterator          ElementsIteratorType;

    // Weak pointers vectors types
    typedef WeakPointerVector<NodeType> NodePointerVector;
    typedef WeakPointerVector<Element> ElementPointerVector;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    /// Definition of the vector indexes considered
    typedef vector<IndexType> VectorIndexType;

    /// Definition of the hasher considered
    typedef VectorIndexHasher<VectorIndexType> VectorIndexHasherType;

    /// Definition of the key comparor considered
    typedef VectorIndexComparor<VectorIndexType> VectorIndexComparorType;

    /// Define the map considered for indexes
    typedef std::unordered_map<VectorIndexType, IndexType, VectorIndexHasherType, VectorIndexComparorType > HashMapVectorIntIntType;

    /// Define the HashMapVectorIntIntType iterator type
    typedef HashMapVectorIntIntType::iterator HashMapVectorIntIntIteratorType;

    /// Define the map considered for elemento pointers
    typedef std::unordered_map<VectorIndexType, Element::Pointer, VectorIndexHasherType, VectorIndexComparorType > HashMapVectorIntElementPointerType;

    /// Define the HashMapVectorIntElementPointerType iterator type
    typedef HashMapVectorIntElementPointerType::iterator HashMapVectorIntElementPointerIteratorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModelPart The model part where the search of neighbours is performed
     */
    PrismNeighboursProcess(ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    virtual ~PrismNeighboursProcess() {}

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method executes the algorithm that looks for neighbour nodes and elements in a  mesh of prismatic elements
     */
    void Execute() override
    {
        KRATOS_TRY;

        for(NodesIterarorType it_node = mrModelPart.NodesBegin(); it_node!=mrModelPart.NodesEnd(); it_node++) {
            (it_node->GetValue(NEIGHBOUR_NODES)).reserve(6); // Just it_node-plane neighbours
            NodePointerVector& r_neighbour_nodes = it_node->GetValue(NEIGHBOUR_NODES);
            r_neighbour_nodes.erase(r_neighbour_nodes.begin(),r_neighbour_nodes.end() );

            (it_node->GetValue(NEIGHBOUR_ELEMENTS)).reserve(3); // Just it_node-plane neighbours
            ElementPointerVector& r_neighbour_elements = it_node->GetValue(NEIGHBOUR_ELEMENTS);
            r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end() );
        }

        for (ElementsIteratorType it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++) {
            (it_elem->GetValue(NEIGHBOUR_NODES)).reserve(6); // Just in-plane neighbours
            NodePointerVector& r_neighbour_nodes = it_elem->GetValue(NEIGHBOUR_NODES);
            r_neighbour_nodes.erase(r_neighbour_nodes.begin(),r_neighbour_nodes.end() );

            (it_elem->GetValue(NEIGHBOUR_ELEMENTS)).reserve(3); // Just in-plane neighbours
            ElementPointerVector& r_neighbour_elements = it_elem->GetValue(NEIGHBOUR_ELEMENTS);
            r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end() );
        }

        /* NEIGHBOUR ELEMENTS */
        // Create the HashMapVectorIntElementPointerType
        HashMapVectorIntElementPointerType face_map;
        
        for (ElementsIteratorType it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++) {
            GeometryType& geom = it_elem->GetGeometry();
             
            // Insert a pointer to the element identified by the hash value ids if it doesn't exist
            ElementPointerVector& r_neighbour_elements = it_elem->GetValue(NEIGHBOUR_ELEMENTS);

            /* IN-PLANE FACES */
            VectorIndexType ids_4(4);
            
            /* FACE 1 */
            ids_4[0] = geom[0].Id();
            ids_4[1] = geom[1].Id();
            ids_4[2] = geom[3].Id();
            ids_4[3] = geom[4].Id();
            
            /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
            std::sort(ids_4.begin(), ids_4.end());
            // Check if the elements already exist in the HashMapVectorIntElementPointerType
            HashMapVectorIntElementPointerIteratorType it_check = face_map.find(ids_4);

            if(it_check != face_map.end() ) {
                // If it exists the node is added as a neighbour, reciprocally
                r_neighbour_elements.push_back(it_check->second);
                ElementPointerVector& aux_3 = (it_check->second)->GetValue(NEIGHBOUR_ELEMENTS);
                aux_3.push_back(*it_elem.base());
            } else {
                // If it doesn't exist it is added to the database
                face_map.insert( HashMapVectorIntElementPointerType::value_type(ids_4, *it_elem.base()) );
            }

            /* FACE 2 */
            ids_4[0] = geom[1].Id();
            ids_4[1] = geom[2].Id();
            ids_4[2] = geom[4].Id();
            ids_4[3] = geom[5].Id();

            /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
            std::sort(ids_4.begin(), ids_4.end());
            // Check if the elements already exist in the HashMapVectorIntElementPointerType
            it_check = face_map.find(ids_4);

            if(it_check != face_map.end() ) {
                // If it exists the node is added as a neighbour, reciprocally
                r_neighbour_elements.push_back(it_check->second);
                ElementPointerVector& aux_el_2 = (it_check->second)->GetValue(NEIGHBOUR_ELEMENTS);
                aux_el_2.push_back(*it_elem.base());
            } else {
                // If it doesn't exist it is added to the database
                face_map.insert( HashMapVectorIntElementPointerType::value_type(ids_4, *it_elem.base()) );
            }

            /* FACE 3 */
            ids_4[0] = geom[0].Id();
            ids_4[1] = geom[2].Id();
            ids_4[2] = geom[3].Id();
            ids_4[3] = geom[5].Id();

            /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
            std::sort(ids_4.begin(), ids_4.end());
            // Check if the elements already exist in the HashMapVectorIntElementPointerType
            it_check = face_map.find(ids_4);

            if(it_check != face_map.end() ) {
                // If it exists the node is added as a neighbour, reciprocally
                r_neighbour_elements.push_back(it_check->second);
                ElementPointerVector& aux_el_3 = (it_check->second)->GetValue(NEIGHBOUR_ELEMENTS);
                aux_el_3.push_back(*it_elem.base());
            } else {
                // If it doesn't exist it is added to the database
                face_map.insert( HashMapVectorIntElementPointerType::value_type(ids_4, *it_elem.base()) );
            }
        }

        /* NEIGHBOURS NODES */

        // Create the ids and aux vectors
        VectorIndexType ids(2),  aux_1(2), aux_2(2), aux_3(2);

        // Search the neighbour nodes(for elements)
        for (ElementsIteratorType it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++) {
            GeometryType& geom = it_elem->GetGeometry();
            NodePointerVector& neighb_nodes = it_elem->GetValue(NEIGHBOUR_NODES);
            neighb_nodes.resize(6);

            for (IndexType fill = 0; fill < 6; fill++) {
                neighb_nodes(fill) = NodeType::WeakPointer(geom(fill));
            }

            // Just upper nodes, the others are +3 IDs
            ElementPointerVector& r_neighbour_elements = it_elem->GetValue(NEIGHBOUR_ELEMENTS);
            for (ElementPointerVector::iterator it_neigh_elem = r_neighbour_elements.begin(); it_neigh_elem != r_neighbour_elements.end(); it_neigh_elem++) {
                HashMapVectorIntIntType node_map;
                GeometryType& geom_neig = it_neigh_elem->GetGeometry();

                // Edge 1
                ids[0] = geom_neig[0].Id();
                ids[1] = geom_neig[1].Id();
                std::sort(ids.begin(), ids.end());
                node_map.insert( HashMapVectorIntIntType::value_type(ids, 2) );

                // Edge 2
                ids[0] = geom_neig[1].Id();
                ids[1] = geom_neig[2].Id();
                std::sort(ids.begin(), ids.end());
                node_map.insert( HashMapVectorIntIntType::value_type(ids, 0) );

                // Edge 3
                ids[0] = geom_neig[2].Id();
                ids[1] = geom_neig[0].Id();
                std::sort(ids.begin(), ids.end());
                node_map.insert( HashMapVectorIntIntType::value_type(ids, 1) );

                aux_1[0] = geom[1].Id();
                aux_1[1] = geom[2].Id();

                aux_2[0] = geom[2].Id();
                aux_2[1] = geom[0].Id();

                aux_3[0] = geom[0].Id();
                aux_3[1] = geom[1].Id();

                std::sort(aux_1.begin(), aux_1.end());
                HashMapVectorIntIntIteratorType it_1 = node_map.find(aux_1);
                std::sort(aux_2.begin(), aux_2.end());
                HashMapVectorIntIntIteratorType it_2 = node_map.find(aux_2);
                std::sort(aux_3.begin(), aux_3.end());
                HashMapVectorIntIntIteratorType it_3 = node_map.find(aux_3);

                if(it_1 != node_map.end() ) {
                    neighb_nodes(0) = NodeType ::WeakPointer(geom_neig(it_1->second));
                    neighb_nodes(3) = NodeType ::WeakPointer(geom_neig(it_1->second + 3));
                } else if(it_2 != node_map.end() ) {
                    neighb_nodes(1) = NodeType ::WeakPointer(geom_neig(it_2->second));
                    neighb_nodes(4) = NodeType ::WeakPointer(geom_neig(it_2->second + 3));
                } else if(it_3 != node_map.end() ) {
                    neighb_nodes(2) = NodeType ::WeakPointer(geom_neig(it_3->second));
                    neighb_nodes(5) = NodeType ::WeakPointer(geom_neig(it_3->second + 3));
                }
            }
        }

        // Add the neighbour elements to all the nodes in the mesh
        for (ElementsIteratorType it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++) {
            GeometryType& geom = it_elem->GetGeometry();
            for(IndexType i = 0; i < geom.size(); i++) {
                (geom[i].GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( *(it_elem.base()) ) );
            }
        }

        // Adding the neighbouring nodes (in the same face)
        for (ElementsIteratorType it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++) {
            GeometryType& geom = it_elem->GetGeometry();

            NodePointerVector& neighb_nodes = it_elem->GetValue(NEIGHBOUR_NODES);

            for (IndexType i = 0; i < 3; i++) {
                NodeType::WeakPointer temp;

                // Adding nodes from the element
                IndexType aux_index1, aux_index2;

                if (i == 0) {
                    aux_index1 = 2;
                    aux_index2 = 1;
                } else if (i == 1) {
                    aux_index1 = 0;
                    aux_index2 = 2;
                } else {
                    aux_index1 = 1;
                    aux_index2 = 0;
                }

                // Lower face
                temp = geom(aux_index1);
                AddUniqueWeakPointer< NodeType >(geom[i].GetValue(NEIGHBOUR_NODES), temp);
                temp = geom(aux_index2);
                AddUniqueWeakPointer< NodeType >(geom[i].GetValue(NEIGHBOUR_NODES), temp);

                // Upper face
                temp = geom(aux_index1 + 3);
                AddUniqueWeakPointer< NodeType >(geom[i + 3].GetValue(NEIGHBOUR_NODES), temp);
                temp = geom(aux_index2 + 3);
                AddUniqueWeakPointer< NodeType >(geom[i + 3].GetValue(NEIGHBOUR_NODES), temp);

                // Adding neighbour elements
                if (neighb_nodes[aux_index1].Id() != geom[i].Id()) {
                    // Lower face
                    temp = neighb_nodes(aux_index1);
                    AddUniqueWeakPointer< NodeType >(geom[i].GetValue(NEIGHBOUR_NODES), temp);
                    // Upper face
                    temp = neighb_nodes(aux_index1 + 3);
                    AddUniqueWeakPointer< NodeType >(geom[i + 3].GetValue(NEIGHBOUR_NODES), temp);
                }
                if (neighb_nodes[aux_index2].Id() != geom[i].Id()) {
                    // Lower face
                    temp = neighb_nodes(aux_index2);
                    AddUniqueWeakPointer< NodeType >(geom[i].GetValue(NEIGHBOUR_NODES), temp);
                    // Upper face
                    temp = neighb_nodes(aux_index2 + 3);
                    AddUniqueWeakPointer< NodeType >(geom[i + 3].GetValue(NEIGHBOUR_NODES), temp);
                }
            }
        }
        
        KRATOS_CATCH("");
    }

    /**
     * @brief This method should be called in case that the current list of neighbour must be drop
     */
    void ClearNeighbours()
    {
        for(NodesIterarorType it_node = mrModelPart.NodesBegin(); it_node!=mrModelPart.NodesEnd(); it_node++) {
            NodePointerVector& r_neighbour_nodes = it_node->GetValue(NEIGHBOUR_NODES);
            r_neighbour_nodes.erase(r_neighbour_nodes.begin(),r_neighbour_nodes.end() );

            ElementPointerVector& r_neighbour_elements = it_node->GetValue(NEIGHBOUR_ELEMENTS);
            r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end());
        }

        for (ElementsIteratorType it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); it_elem++) {
            ElementPointerVector& r_neighbour_elements = it_elem->GetValue(NEIGHBOUR_ELEMENTS);
            r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end());
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
    virtual std::string Info() const override
    {
        return "PrismNeighboursProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PrismNeighboursProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }

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

    ModelPart& mrModelPart; /// The main model part

    ///@}
    ///@name Private Operators
    ///@{
    
    template< class TDataType > void  AddUniqueWeakPointer
    (WeakPointerVector< TDataType >& v, const typename TDataType::WeakPointer candidate)
    {
        typename WeakPointerVector< TDataType >::iterator i     = v.begin();
        typename WeakPointerVector< TDataType >::iterator endit = v.end();
        while ( i != endit && (i)->Id() != (candidate.lock())->Id()) {
            i++;
        }
        if( i == endit ) {
            v.push_back(candidate);
        }

    }

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
    PrismNeighboursProcess& operator=(PrismNeighboursProcess const& rOther);

    ///@}

}; // Class PrismNeighboursProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  PrismNeighboursProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PrismNeighboursProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_PRISM_NEIGHBOURS_PROCESS_H_INCLUDED  defined

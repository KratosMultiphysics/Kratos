// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_SPRISM_NEIGHBOURS_H_INCLUDED )
#define  KRATOS_SPRISM_NEIGHBOURS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "utilities/math_utils.h"

#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp> 
#include <utility>

#include "structural_mechanics_application_variables.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
typedef  ModelPart::NodesContainerType NodesContainerType;
typedef  ModelPart::ElementsContainerType ElementsContainerType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{
#if !defined(KEY_COMPAROR)
#define KEY_COMPAROR
struct KeyComparor
{
    bool operator()(const vector<int>& lhs, const vector<int>& rhs) const
    {
        if(lhs.size() != rhs.size())
            return false;

        for(unsigned int i = 0; i < lhs.size(); i++)
        {
            if(lhs[i] != rhs[i]) return false;
        }

        return true;
    }
};
#endif
    
#if !defined(KEY_HASHER)
#define KEY_HASHER
struct KeyHasher
{
    std::size_t operator()(const vector<int>& k) const
    {
        return boost::hash_range(k.begin(), k.end());
    }
};
#endif

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class SprismNeighbours
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SprismNeighbours(ModelPart& model_part)
        : mr_model_part(model_part)
    {
    }

    /// Destructor.
    virtual ~SprismNeighbours() {}

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

    virtual void Execute()
    {
        KRATOS_TRY;

        NodesContainerType& rNodes = mr_model_part.Nodes();
        ElementsContainerType& rElems = mr_model_part.Elements();

        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            (in->GetValue(NEIGHBOUR_NODES)).reserve(6); // Just in-plane neighbours
            WeakPointerVector<Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);
            rN.erase(rN.begin(),rN.end() );

            (in->GetValue(NEIGHBOUR_ELEMENTS)).reserve(3); // Just in-plane neighbours
            WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(),rE.end() );
        }

        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
        {
            (ie->GetValue(NEIGHBOUR_NODES)).reserve(6); // Just in-plane neighbours
            WeakPointerVector<Node<3> >& rN = ie->GetValue(NEIGHBOUR_NODES);
            rN.erase(rN.begin(),rN.end() );

            (ie->GetValue(NEIGHBOUR_ELEMENTS)).reserve(3); // Just in-plane neighbours
            WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(),rE.end() );
        }

        /* NEIGHBOUR ELEMENTS */
        // Create the hashmap_el
        typedef boost::unordered_map<vector<int>, Element::Pointer, KeyHasher, KeyComparor > hashmap_el;
        hashmap_el face_map;
        
        for (ModelPart::ElementIterator itElem = mr_model_part.ElementsBegin(); itElem != mr_model_part.ElementsEnd(); itElem++)
        {
            Geometry< Node<3> >& geom = itElem->GetGeometry();
             
            // Insert a pointer to the element identified by the hash value ids if it doesn't exist
            WeakPointerVector< Element >& neighb_elems = itElem->GetValue(NEIGHBOUR_ELEMENTS);

            /* IN-PLANE FACES */
            vector<int> ids_4(4);
            
            /* FACE 1 */
            ids_4[0] = geom[0].Id();
            ids_4[1] = geom[1].Id();
            ids_4[2] = geom[3].Id();
            ids_4[3] = geom[4].Id();
            
            /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
            std::sort(ids_4.begin(), ids_4.end());
            // Check if the elements already exist in the hashmap_el
            hashmap_el::iterator it_check = face_map.find(ids_4);

            if(it_check != face_map.end() )
            {
                // If it exists the node is added as a neighbour, reciprocally
                neighb_elems.push_back(it_check->second);
                WeakPointerVector< Element >& aux_3 = (it_check->second)->GetValue(NEIGHBOUR_ELEMENTS);
                aux_3.push_back(*itElem.base());
            }
            else
            {
                // If it doesn't exist it is added to the database
                face_map.insert( hashmap_el::value_type(ids_4, *itElem.base()) );
            }

            /* FACE 2 */
            ids_4[0] = geom[1].Id();
            ids_4[1] = geom[2].Id();
            ids_4[2] = geom[4].Id();
            ids_4[3] = geom[5].Id();

            /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
            std::sort(ids_4.begin(), ids_4.end());
            // Check if the elements already exist in the hashmap_el
            it_check = face_map.find(ids_4);

            if(it_check != face_map.end() )
            {
                // If it exists the node is added as a neighbour, reciprocally
                neighb_elems.push_back(it_check->second);
                WeakPointerVector< Element >& aux_el_2 = (it_check->second)->GetValue(NEIGHBOUR_ELEMENTS);
                aux_el_2.push_back(*itElem.base());
            }
            else
            {
                // If it doesn't exist it is added to the database
                face_map.insert( hashmap_el::value_type(ids_4, *itElem.base()) );
            }

            /* FACE 3 */
            ids_4[0] = geom[0].Id();
            ids_4[1] = geom[2].Id();
            ids_4[2] = geom[3].Id();
            ids_4[3] = geom[5].Id();

            /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
            std::sort(ids_4.begin(), ids_4.end());
            // Check if the elements already exist in the hashmap_el
            it_check = face_map.find(ids_4);

            if(it_check != face_map.end() )
            {
                // If it exists the node is added as a neighbour, reciprocally
                neighb_elems.push_back(it_check->second);
                WeakPointerVector< Element >& aux_el_3 = (it_check->second)->GetValue(NEIGHBOUR_ELEMENTS);
                aux_el_3.push_back(*itElem.base());
            }
            else
            {
                // If it doesn't exist it is added to the database
                face_map.insert( hashmap_el::value_type(ids_4, *itElem.base()) );
            }
        }

//         /* Print the result */
//         for (ModelPart::ElementIterator itElem = mr_model_part.ElementsBegin(); itElem != mr_model_part.ElementsEnd(); itElem++)
//         {
//            KRATOS_WATCH(itElem->Id());
//            WeakPointerVector<Element >& rE = itElem->GetValue(NEIGHBOUR_ELEMENTS);
//            for (WeakPointerVector<Element>::iterator ttt = rE.begin(); ttt != rE.end(); ttt++)
//            {
//                KRATOS_WATCH(ttt->Id());
//            }
//         }

        /* NEIGHBOURS NODES */
        // Create the hashmap_no
        typedef boost::unordered_map<vector<int>, int, KeyHasher, KeyComparor > hashmap_no;

        // Create the ids and aux vectors
        vector<int> ids(2);
        vector<int> aux_1(2);
        vector<int> aux_2(2);
        vector<int> aux_3(2);

        // Search the neighbour nodes(for elements)
        for (ModelPart::ElementIterator itElem = mr_model_part.ElementsBegin(); itElem != mr_model_part.ElementsEnd(); itElem++)
        {
            Geometry< Node<3> >& geom = itElem->GetGeometry();
            WeakPointerVector< Node < 3 > >& neighb_nodes = itElem->GetValue(NEIGHBOUR_NODES);
            neighb_nodes.resize(6);

            for (unsigned int fill = 0; fill < 6; fill++)
            {
                neighb_nodes(fill) = Node < 3 > ::WeakPointer(geom(fill));
            }

            // Just upper nodes, the others are +3 IDs
            WeakPointerVector<Element >& rE = itElem->GetValue(NEIGHBOUR_ELEMENTS);
            for (WeakPointerVector<Element>::iterator nel = rE.begin(); nel != rE.end(); nel++)
            {
                hashmap_no node_map;
                Geometry< Node<3> >& geom_neig = nel->GetGeometry();

                // Edge 1
                ids[0] = geom_neig[0].Id();
                ids[1] = geom_neig[1].Id();
                std::sort(ids.begin(), ids.end());
                node_map.insert( hashmap_no::value_type(ids, 2) );

                // Edge 2
                ids[0] = geom_neig[1].Id();
                ids[1] = geom_neig[2].Id();
                std::sort(ids.begin(), ids.end());
                node_map.insert( hashmap_no::value_type(ids, 0) );

                // Edge 3
                ids[0] = geom_neig[2].Id();
                ids[1] = geom_neig[0].Id();
                std::sort(ids.begin(), ids.end());
                node_map.insert( hashmap_no::value_type(ids, 1) );

                aux_1[0] = geom[1].Id();
                aux_1[1] = geom[2].Id();

                aux_2[0] = geom[2].Id();
                aux_2[1] = geom[0].Id();

                aux_3[0] = geom[0].Id();
                aux_3[1] = geom[1].Id();

                std::sort(aux_1.begin(), aux_1.end());
                hashmap_no::iterator it_1 = node_map.find(aux_1);
                std::sort(aux_2.begin(), aux_2.end());
                hashmap_no::iterator it_2 = node_map.find(aux_2);
                std::sort(aux_3.begin(), aux_3.end());
                hashmap_no::iterator it_3 = node_map.find(aux_3);

                if(it_1 != node_map.end() )
                {
                    neighb_nodes(0) = Node < 3 > ::WeakPointer(geom_neig(it_1->second));
                    neighb_nodes(3) = Node < 3 > ::WeakPointer(geom_neig(it_1->second + 3));
                }
                else if(it_2 != node_map.end() )
                {
                    neighb_nodes(1) = Node < 3 > ::WeakPointer(geom_neig(it_2->second));
                    neighb_nodes(4) = Node < 3 > ::WeakPointer(geom_neig(it_2->second + 3));
                }
                else if(it_3 != node_map.end() )
                {
                    neighb_nodes(2) = Node < 3 > ::WeakPointer(geom_neig(it_3->second));
                    neighb_nodes(5) = Node < 3 > ::WeakPointer(geom_neig(it_3->second + 3));
                }
            }
        }

//         /* Print the result */
//         for (ModelPart::ElementIterator itElem = mr_model_part.ElementsBegin(); itElem != mr_model_part.ElementsEnd(); itElem++)
//         {
//            KRATOS_WATCH(itElem->Id());
//            Geometry< Node<3> >& geom = itElem->GetGeometry();
//            WeakPointerVector< Node < 3 > >& neighb_nodes = itElem->GetValue(NEIGHBOUR_NODES);
//            for (unsigned int i = 0; i < 6; i++)
//            {
//                std::cout << "Node " << i << " " << geom[i].Id() << " " << neighb_nodes[i].Id() <<std::endl;
//            }
//         }

        // Add the neighbour elements to all the nodes in the mesh
        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
        {
            Element::GeometryType& pGeom = ie->GetGeometry();
            for(unsigned int i = 0; i < pGeom.size(); i++)
            {
                (pGeom[i].GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( *(ie.base()) ) );
            }
        }

//        // Adding the neighbouring nodes
//        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
//        {
//            WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);

//            for(unsigned int ie = 0; ie < rE.size(); ie++)
//            {
//                Element::GeometryType& pGeom = rE[ie].GetGeometry();
//                for(unsigned int i = 0; i < pGeom.size(); i++)
//                {
//                    if(pGeom[i].Id() != in->Id() )
//                    {
//                        Element::NodeType::WeakPointer temp = pGeom(i);
//                        AddUniqueWeakPointer< Node<3> >(in->GetValue(NEIGHBOUR_NODES), temp);
//                    }
//                }
//            }
//        }

        // Adding the neighbouring nodes (in the same face)
        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
        {
            Element::GeometryType& pGeom = ie->GetGeometry();

            WeakPointerVector< Node < 3 > >& neighb_nodes = ie->GetValue(NEIGHBOUR_NODES);

            for (unsigned int i = 0; i < 3; i++)
            {
                Element::NodeType::WeakPointer temp;

                // Adding nodes from the element
                int aux_index1;
                int aux_index2;

                if (i == 0)
                {
                    aux_index1 = 2;
                    aux_index2 = 1;
                }
                else if (i == 1)
                {
                    aux_index1 = 0;
                    aux_index2 = 2;
                }
                else
                {
                    aux_index1 = 1;
                    aux_index2 = 0;
                }

                // Lower face
                temp = pGeom(aux_index1);
                AddUniqueWeakPointer< Node<3> >(pGeom[i].GetValue(NEIGHBOUR_NODES), temp);
                temp = pGeom(aux_index2);
                AddUniqueWeakPointer< Node<3> >(pGeom[i].GetValue(NEIGHBOUR_NODES), temp);

                // Upper face
                temp = pGeom(aux_index1 + 3);
                AddUniqueWeakPointer< Node<3> >(pGeom[i + 3].GetValue(NEIGHBOUR_NODES), temp);
                temp = pGeom(aux_index2 + 3);
                AddUniqueWeakPointer< Node<3> >(pGeom[i + 3].GetValue(NEIGHBOUR_NODES), temp);

                // Adding neighbour elements
                if (neighb_nodes[aux_index1].Id() != pGeom[i].Id())
                {
                    // Lower face
                    temp = neighb_nodes(aux_index1);
                    AddUniqueWeakPointer< Node<3> >(pGeom[i].GetValue(NEIGHBOUR_NODES), temp);
                    // Upper face
                    temp = neighb_nodes(aux_index1 + 3);
                    AddUniqueWeakPointer< Node<3> >(pGeom[i + 3].GetValue(NEIGHBOUR_NODES), temp);
                }
                if (neighb_nodes[aux_index2].Id() != pGeom[i].Id())
                {
                    // Lower face
                    temp = neighb_nodes(aux_index2);
                    AddUniqueWeakPointer< Node<3> >(pGeom[i].GetValue(NEIGHBOUR_NODES), temp);
                    // Upper face
                    temp = neighb_nodes(aux_index2 + 3);
                    AddUniqueWeakPointer< Node<3> >(pGeom[i + 3].GetValue(NEIGHBOUR_NODES), temp);
                }
            }

        }

//        /* Print the result */
//        for (ModelPart::NodeIterator itNode = mr_model_part.NodesBegin(); itNode != mr_model_part.NodesEnd(); itNode++)
//        {
//            WeakPointerVector< Node < 3 > >& aux_neighb = itNode->GetValue(NEIGHBOUR_NODES);
//            unsigned int aux_index = aux_neighb.size();

//            std::cout << "Node " << itNode->Id() << std::endl;
//            for (unsigned int i = 0; i < aux_index; i++)
//            {
//                std::cout << "Neighbour " << i << " " << aux_neighb[i].Id() << std::endl;
//            }
//        }
        
        KRATOS_CATCH("");
    }

    void ClearNeighbours()
    {
        NodesContainerType& rNodes = mr_model_part.Nodes();
        for(NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        {
            WeakPointerVector<Element >& rE = in->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(),rE.end());

            WeakPointerVector<Node<3> >& rN = in->GetValue(NEIGHBOUR_NODES);
            rN.erase(rN.begin(),rN.end() );
        }

        ElementsContainerType& rElems = mr_model_part.Elements();
        for(ElementsContainerType::iterator ie = rElems.begin(); ie!=rElems.end(); ie++)
        {
            WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);
            rE.erase(rE.begin(),rE.end());
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
    virtual std::string Info() const
    {
        return "SprismNeighbours";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SprismNeighbours";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
    ModelPart& mr_model_part;

    ///@}
    ///@name Private Operators
    ///@{
    
    template< class TDataType > void  AddUniqueWeakPointer
    (WeakPointerVector< TDataType >& v, const typename TDataType::WeakPointer candidate)
    {
        typename WeakPointerVector< TDataType >::iterator i     = v.begin();
        typename WeakPointerVector< TDataType >::iterator endit = v.end();
        while ( i != endit && (i)->Id() != (candidate.lock())->Id())
        {
            i++;
        }
        if( i == endit )
        {
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
    SprismNeighbours& operator=(SprismNeighbours const& rOther);

    ///@}

}; // Class SprismNeighbours

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SprismNeighbours& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SprismNeighbours& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SPRISM_NEIGHBOURS_H_INCLUDED  defined

/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_EDGE_SWAPPING_2D_MODELER_H_INCLUDED )
#define  KRATOS_EDGE_SWAPPING_2D_MODELER_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "modeler/modeler.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "utilities/timer.h"
#include "processes/find_elements_neighbours_process.h"
#include "processes/node_erase_process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Optimizes a 2D mesh by swapping the edges between elements.
/** This class can be used to optimized a 2D mesh modifying the connectivities of the elements.
    This modeler also collapses the nodes which are marked to be deleted. Nodes must be marked
  by setting the TO_ERASE in data value container (using Node::SetValue) to true.
  The algorithm consist of doing iteration over following steps:
	- Marking the pair of elements to be swapped by checking the Delaunay criteria
	- Swapping all marked pair of elements
	- Collapsing the mark nodes for erasing with TO_ERASE set to true.

  @see Remesh
*/
class EdgeSwapping2DModeler : public Modeler
{


    struct SwappingData
    {
    public:
        SwappingData()
        {
            Reset();
        }
        void Reset()
        {
            for(int i = 0 ; i < 3 ; i++)
            {
                NeighbourElements[i] = -1;
                OppositeNodes[i] = -1;
                OppositeEdge[i] = -1;
                IsInside[i] = false;
            }
            IsSwapCandidate = false;
            IsElementErased = false;
            SwapWith = -1;
        }
        array_1d<int, 3> NeighbourElements;
        array_1d<int, 3> OppositeNodes;
        array_1d<int, 3> OppositeEdge;
        array_1d<bool, 3> IsInside;
        bool IsSwapCandidate;
        bool IsElementErased;
        int SwapWith;
        int SwapEdge;
    };

    struct CollapsingData
    {
    public:
        CollapsingData()
        {
            Reset();
        }
        void Reset()
        {
            for(int i = 0 ; i < 2 ; i++)
            {
                CollapseElements[i] = -1;
            }
            MinimumDistance = -1;
            NearestNodeId = -1;
        }
        double MinimumDistance;
        int NearestNodeId;
        array_1d<int, 2> CollapseElements;
    };



public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EdgeSwapping2DModeler
    KRATOS_CLASS_POINTER_DEFINITION(EdgeSwapping2DModeler);

    /** Definition of the Modeler class as BaseType of EdgeSwapping2DModeler
    */
    typedef Modeler BaseType;

    typedef Point<3> PointType;

    /** Defining 3D Node class as nodes of the mesh
    */
    typedef Node<3> NodeType;

    /** Defining Geometry of 3D Node class as geometry type
    */
    typedef Geometry<NodeType> GeometryType;


    /** A PointerVector of 3D Node class as array of nodes.
    */
    typedef PointerVector<NodeType> NodesVectorType;

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Empty default constructor.
    EdgeSwapping2DModeler()
    {
    }

    /// Empty destructor.
    virtual ~EdgeSwapping2DModeler() {}



    ///@}
    ///@name Operations
    ///@{



    /** Remesh is the main method of this class and perform all the swapping and collapsing process.
        It accepts a ModelPart as its inputs and performs the remeshing process over its current mesh.
      The iterative process of finding pair of swapping elements and swap them is limited to
      maximum_repeat_number which is set to 10.

      @param rThisModelPart The model part containing the mesh to be remeshed

    */
    void Remesh(ModelPart& rThisModelPart)
    {
        Timer::Start("Edge Swapping");
        const int maximum_repeat_number = 10;

        unsigned int number_of_bad_quality_elements = 1 ;//MarkBadQualityElements(rThisModelPart);
        //unsigned int number_of_bad_quality_elements_old = number_of_bad_quality_elements;


        //  	    FindElementalNeighboursProcess find_element_neighbours_process(rThisModelPart, 2);


        // 	    find_element_neighbours_process.Execute();


        //  WriteMesh(rThisModelPart, "/home/pooyan/kratos/tests/before_remesh.post.msh");
        for(int repeat = 0 ; (repeat < maximum_repeat_number) & (number_of_bad_quality_elements != 0) ; repeat++)
        {
            ModelPart::NodesContainerType::ContainerType& nodes_array = rThisModelPart.NodesArray();
            ModelPart::ElementsContainerType::ContainerType& elements_array = rThisModelPart.ElementsArray();

            const int number_of_nodes = rThisModelPart.NumberOfNodes();
            const int number_of_elements = rThisModelPart.NumberOfElements();

            #pragma omp parallel for
            for(int i = 0 ; i < number_of_nodes ; i++)
                nodes_array[i]->SetId(i+1);

            #pragma omp parallel for
            for(int i = 0 ; i < number_of_elements ; i++)
                elements_array[i]->SetId(i+1);

            std::cout << "Edge swapping iteration #" << repeat; // << " : No. of bad quality elements = " << number_of_bad_quality_elements;
            SetSwappingData(rThisModelPart);


            unsigned int swap_counter = 0;
            for(int i = 0 ; i < number_of_elements ; i++)
            {
                if(mSwappingData[i].IsSwapCandidate == false)
                    //if(mBadQuality[i])
                {
                    for(int j = 0 ; j < 3 ; j++)
                    {
                        const int neighbour_index = mSwappingData[i].NeighbourElements[j] - 1;
                        if(neighbour_index >= 0)
                        {
                            if(mSwappingData[neighbour_index].IsSwapCandidate == false)
                            {
                                // std::cout << "check element #" << elements_array[i]->Id() << " with nodes: " << (*(elements_array[i])).GetGeometry()[0].Id()  << "," << (*(elements_array[i])).GetGeometry()[1].Id() << ","  << (*(elements_array[i])).GetGeometry()[2].Id() << std::endl;
                                if(mSwappingData[i].IsInside[j])
                                    // if(IsInElementSphere(*(elements_array[i]), *(nodes_array[mSwappingData[i].OppositeNodes[j]-1])))
                                {
                                    swap_counter++;
                                    mSwappingData[neighbour_index].IsSwapCandidate = true;
                                    mSwappingData[neighbour_index].SwapEdge = mSwappingData[i].OppositeEdge[j];
                                    mSwappingData[i].IsSwapCandidate = true;
                                    mSwappingData[i].SwapWith = neighbour_index;
                                    mSwappingData[i].SwapEdge = j;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            std::cout << " number of swaps = " << swap_counter << std::endl;
            if(swap_counter == 0)
                break;

            // 		for(int i = 0 ; i < number_of_elements ; i++)
            // 		{
            // 			if(mSwappingData[i].SwapWith != -1)
            // 			{
            // 				std::cout << "Element #" << i + 1 << " : " ;
            // 				std::cout << mSwappingData[i].NeighbourElements << ",";
            // 				std::cout << mSwappingData[i].OppositeNodes << ",";
            // 				std::cout << mSwappingData[i].SwapWith;
            // 				std::cout << std::endl;
            // 			}
            // 		}


            for(int i = 0 ; i < number_of_elements ; i++)
            {
                if(mSwappingData[i].SwapWith != -1)
                {
                    //if((elements_array[i]->Id() > 3100 ) && (elements_array[i]->Id() < 3200 ))
                    // std::cout << "swapping element #" << elements_array[i]->Id() << " with element #" << elements_array[mSwappingData[i].SwapWith]->Id() << std::endl;
                    Swap(*(elements_array[i]),*(elements_array[mSwappingData[i].SwapWith]), mSwappingData[i].SwapEdge, mSwappingData[mSwappingData[i].SwapWith].SwapEdge);
                }
            }
            //KRATOS_WATCH("Swapping done!!");
            //number_of_bad_quality_elements = MarkBadQualityElements(rThisModelPart);
            //if(number_of_bad_quality_elements == number_of_bad_quality_elements_old)
            //{
            // break;
            //}
            //else
            //{
            // number_of_bad_quality_elements_old=number_of_bad_quality_elements;
            //}
            //WriteMesh(rThisModelPart, "/home/pooyan/kratos/tests/before_collapse.post.msh");

            CollapseNodes(rThisModelPart);
            //WriteMesh(rThisModelPart, "/home/pooyan/kratos/tests/after_collapse.post.msh");

        }
        //WriteMesh(rThisModelPart, "/home/pooyan/kratos/tests/after_remesh.post.msh");

        //set the coordinates to the original value
        /*  for( ModelPart::NodeIterator it =  rThisModelPart.NodesBegin(); it!= rThisModelPart.NodesEnd(); it++)
        {
        	const array_1d<double,3>& disp = (it)->FastGetSolutionStepValue(DISPLACEMENT);
        	(it)->X0() = (it)->X() - disp[0];
        	(it)->Y0() = (it)->Y() - disp[1];
        	(it)->Z0() = 0.0;
        }
        */
        Timer::Stop("Edge Swapping");
    }

    /** An auxiliary method for writing the mesh for GiD for debugging purpose
    */
    void WriteMesh(ModelPart& rThisModelPart, std::string Filename)
    {
        std::ofstream temp_file(Filename.c_str());

        temp_file << "MESH \"Kratos_Triangle2D3_Mesh\" dimension 2 ElemType Triangle Nnode 3" << std::endl;

        temp_file << "Coordinates" << std::endl;

        for(ModelPart::NodeIterator i_node = rThisModelPart.NodesBegin() ; i_node != rThisModelPart.NodesEnd() ; i_node++)
            temp_file << i_node->Id() << " " << i_node->X() << " " << i_node->Y() << " " << i_node->Z() << std::endl;

        temp_file << "End Coordinates" << std::endl;

        temp_file << "Elements" << std::endl;

        for(ModelPart::ElementIterator i_element = rThisModelPart.ElementsBegin() ; i_element != rThisModelPart.ElementsEnd() ; i_element++)
            temp_file << i_element->Id() << " " << i_element->GetGeometry()[0].Id() << " " << i_element->GetGeometry()[1].Id() << " " << i_element->GetGeometry()[2].Id() << " 0"<< std::endl;


        temp_file << "End Elements" << std::endl;



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
        return "EdgeSwapping2DModeler";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }




    ///@}

private:


    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    std::vector<int> mBadQuality;
    std::vector<std::vector<int> > mNodalNeighbourElements;
    std::vector<CollapsingData> mCollapsingData;
    std::vector<SwappingData > mSwappingData;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void CollapseNodes(ModelPart& rThisModelPart)
    {
        FindNodalNeighbours(rThisModelPart);
        SetCollapsingData(rThisModelPart);
        //const double threshold = 0.2;

        ModelPart::NodesContainerType::ContainerType& nodes_array = rThisModelPart.NodesArray();
        ModelPart::ElementsContainerType::ContainerType& elements_array = rThisModelPart.ElementsArray();

        const int number_of_nodes = rThisModelPart.NumberOfNodes();
        const int number_of_elements = rThisModelPart.NumberOfElements();

        for(int i = 0 ; i < number_of_nodes ; i++)
        {
            NodeType& r_node = *(nodes_array[i]);

            //			  std::cout << "node #" << r_node.Id() << " TO_ERASE : " << r_node.Is(TO_ERASE) << std::endl;

            if(r_node.Is(TO_ERASE))
            {
                const int node_index = r_node.Id() - 1;
                const int nearest_node_id = mCollapsingData[i].NearestNodeId;

                if(nearest_node_id != -1) // Exist some node to collapsed to
                {
                    NodeType::Pointer p_node = nodes_array[nearest_node_id - 1];
                    //look for the elements around node
                    for(std::vector<int>::iterator i = mNodalNeighbourElements[node_index].begin() ; i != mNodalNeighbourElements[node_index].end() ; i++)
                    {
                        //look for the this node in neighbour elements
                        Geometry<Node<3> >& r_neighbour_element_geometry = elements_array[*i-1]->GetGeometry();
                        for( unsigned int node_i = 0 ; node_i < r_neighbour_element_geometry.size(); node_i++)
                        {
                            int other_node_id = r_neighbour_element_geometry[node_i].Id();
                            if(other_node_id == static_cast<int>(r_node.Id()))
                            {
                                r_neighbour_element_geometry(node_i) = p_node;
                            }
                            else if(other_node_id == nearest_node_id)
                            {
                                mSwappingData[*i-1].IsElementErased = true;
                            }

                        }
                    }
                }
                else
                {
                    r_node.Set(TO_ERASE, false);
                }
            }
        }

        NodeEraseProcess node_erase_process(rThisModelPart);

        node_erase_process.Execute();


        ModelPart::ElementsContainerType temp_elements_container;
        temp_elements_container.reserve(number_of_elements);

        temp_elements_container.swap(rThisModelPart.Elements());



        for(ModelPart::ElementsContainerType::iterator i_element = temp_elements_container.begin() ; i_element != temp_elements_container.end() ; i_element++)
        {
            if( static_cast<bool>(mSwappingData[i_element->Id()-1].IsElementErased) == false)
                (rThisModelPart.Elements()).push_back(*(i_element.base()));
        }

    }



    void SetCollapsingData(ModelPart& rThisModelPart)
    {
        const int number_of_nodes = static_cast<int>(rThisModelPart.NumberOfNodes());
        ModelPart::NodesContainerType::ContainerType& nodes_array = rThisModelPart.NodesArray();

        if(static_cast<int>(mCollapsingData.size()) != number_of_nodes)
            mCollapsingData.resize(number_of_nodes);

        for(int i = 0 ; i < number_of_nodes ; i++)
            mCollapsingData[i].Reset();

        for(int i = 0 ; i < number_of_nodes ; i++)
        {
            NodeType& r_node = *(nodes_array[i]);

            if(r_node.Is(TO_ERASE))
            {
                SetNodalCollapsingData(r_node, rThisModelPart);
            }
        }
    }

    void SetNodalCollapsingData(NodeType& rThisNode, ModelPart& rThisModelPart)
    {
        ModelPart::ElementsContainerType::ContainerType& elements_array = rThisModelPart.ElementsArray();
        //const unsigned int number_of_elements = rThisModelPart.NumberOfElements();

        const int node_index = rThisNode.Id() - 1;

        int next[3] = {1,2,0};
        int previous[3] = {2,0,1};

        //look for the elements around node
        for(std::vector<int>::iterator i = mNodalNeighbourElements[node_index].begin() ; i != mNodalNeighbourElements[node_index].end() ; i++)
        {
            //look for the nodes of the neighbour faces
            Geometry<Node<3> >& r_neighbour_element_geometry = elements_array[*i-1]->GetGeometry();
            for( unsigned int i_node = 0 ; i_node < r_neighbour_element_geometry.size(); i_node++)
            {
                if(r_neighbour_element_geometry[i_node].Is(TO_ERASE) == 0) // can be used for collapse and is not the same node!
                {
                    int other_node_id = r_neighbour_element_geometry[i_node].Id();
                    if(other_node_id == mCollapsingData[node_index].NearestNodeId)
                    {
                        mCollapsingData[node_index].CollapseElements[1] = *i;
                    }
                    else if(r_neighbour_element_geometry[next[i_node]].Id() == rThisNode.Id())
                    {
                        double d=Distance2(r_neighbour_element_geometry[i_node], r_neighbour_element_geometry[next[i_node]]);

                        if((d < mCollapsingData[node_index].MinimumDistance) || (mCollapsingData[node_index].MinimumDistance < 0.00))
                        {
                            mCollapsingData[node_index].MinimumDistance = d;
                            mCollapsingData[node_index].NearestNodeId = other_node_id;
                            mCollapsingData[node_index].CollapseElements[0] = *i;
                        }
                    }
                    else if(r_neighbour_element_geometry[previous[i_node]].Id() == rThisNode.Id())
                    {
                        double d=Distance2(r_neighbour_element_geometry[i_node], r_neighbour_element_geometry[previous[i_node]]);

                        if((d < mCollapsingData[node_index].MinimumDistance) || (mCollapsingData[node_index].MinimumDistance < 0.00))
                        {
                            mCollapsingData[node_index].MinimumDistance = d;
                            mCollapsingData[node_index].NearestNodeId = other_node_id;
                            mCollapsingData[node_index].CollapseElements[0] = *i;
                        }
                    }

                }
            }

        }
    }



    double Distance2(NodeType& rNode1, NodeType& rNode2)
    {
        double dx = rNode1.X() - rNode2.X();
        double dy = rNode1.Y() - rNode2.Y();

        return dx*dx + dy*dy;
    }

    void Swap(Element& rElement1, Element& rElement2, int Edge1, int Edge2)
    {

        int next2[3] = {2,0,1};

        /* 		  std::cout << "Element1 #" << rElement1.Id() << " : " << rElement1.GetGeometry()[0].Id()  << " : " << rElement1.GetGeometry()[1].Id() << " : " << rElement1.GetGeometry()[2].Id() << std::endl; */
        /* 		  std::cout << "Element2 #" << rElement2.Id() << " : " << rElement2.GetGeometry()[0].Id()  << " : " << rElement2.GetGeometry()[1].Id() << " : " << rElement2.GetGeometry()[2].Id() << std::endl; */
        //KRATOS_WATCH(Edge1);
        //KRATOS_WATCH(Edge2);

        //KRATOS_WATCH(rElement1.GetGeometry().Area())
        //KRATOS_WATCH(rElement2.GetGeometry().Area())

        rElement1.GetGeometry()(next2[Edge1]) = rElement2.GetGeometry()(Edge2);
        rElement2.GetGeometry()(next2[Edge2]) = rElement1.GetGeometry()(Edge1);

        //KRATOS_WATCH(rElement1.GetGeometry().Area())
        //KRATOS_WATCH(rElement2.GetGeometry().Area())


    }

    unsigned int MarkBadQualityElements(ModelPart& rThisModelPart)
    {
// 		  unsigned int counter = 0;
        unsigned int marked = 0;

        const double threshold = 0.25;

        if(mBadQuality.size() != rThisModelPart.NumberOfElements())
            mBadQuality.resize(rThisModelPart.NumberOfElements());

        ModelPart::ElementsContainerType::ContainerType& elements_array = rThisModelPart.ElementsArray();
        const int number_of_elements = rThisModelPart.NumberOfElements();

        #pragma omp parallel for
        for(int i = 0 ; i < number_of_elements ; i++)
        {
            if(ElementQuality(*elements_array[i]) < threshold)
            {
                marked++;
                mBadQuality[i] = true;
            }
        }

// To be parallelized.
// 		  for(ModelPart::ElementIterator i_element = rThisModelPart.ElementsBegin() ; i_element != rThisModelPart.ElementsEnd() ; i_element++)
// 		  {
// 			  if(ElementQuality(*i_element) < threshold)
// 			  {
// 				  marked++;
// 				  mBadQuality[counter] = true;
// 			  }
// 			  counter++;
// 		  }

        return marked;
    }

    double ElementQuality(Element& rThisElement)
    {
        double h_min;
        double h_max;
        double area = rThisElement.GetGeometry().Area();

        GeometryUtils::SideLenghts2D(rThisElement.GetGeometry(), h_min, h_max);

        return (area / (h_max*h_max));
    }

    void FindNodalNeighbours(ModelPart& rThisModelPart)
    {
        ModelPart::ElementsContainerType::ContainerType& elements_array = rThisModelPart.ElementsArray();
        const int number_of_nodes = static_cast<int>(rThisModelPart.NumberOfNodes());
        const int number_of_elements = static_cast<int>(rThisModelPart.NumberOfElements());

        if(static_cast<int>(mNodalNeighbourElements.size()) != number_of_nodes)
            mNodalNeighbourElements.resize(number_of_nodes);

        for(int i = 0 ; i < number_of_nodes ; i++)
            mNodalNeighbourElements[i].clear();

        for(int i = 0 ; i < number_of_elements ; i++)
        {
            Element::GeometryType& r_geometry = elements_array[i]->GetGeometry();
            for(unsigned int j = 0; j < r_geometry.size(); j++)
            {
                mNodalNeighbourElements[r_geometry[j].Id()-1].push_back(elements_array[i]->Id());
            }
        }

        for(int i = 0 ; i < number_of_nodes ; i++)
        {
            std::sort(mNodalNeighbourElements[i].begin(),mNodalNeighbourElements[i].end());
            std::vector<int>::iterator new_end = std::unique(mNodalNeighbourElements[i].begin(),mNodalNeighbourElements[i].end());
            mNodalNeighbourElements[i].erase(new_end, mNodalNeighbourElements[i].end());
        }
    }

    void SetSwappingData(ModelPart& rThisModelPart)
    {
        ModelPart::NodesContainerType::ContainerType& nodes_array = rThisModelPart.NodesArray();
        ModelPart::ElementsContainerType::ContainerType& elements_array = rThisModelPart.ElementsArray();
        const int number_of_elements = static_cast<int>(rThisModelPart.NumberOfElements());

        FindNodalNeighbours(rThisModelPart);

        if(static_cast<int>(mSwappingData.size()) != number_of_elements)
            mSwappingData.resize(number_of_elements);

        #pragma omp parallel for
        for(int i = 0 ; i < number_of_elements ; i++)
            mSwappingData[i].Reset();

        #pragma omp parallel for
        for(int i = 0 ; i < number_of_elements ; i++)
        {
            Element::GeometryType& r_geometry = elements_array[i]->GetGeometry();
            int id = elements_array[i]->Id();
            FindNeighbourElement(r_geometry[1].Id(), r_geometry[2].Id(), id, elements_array, mSwappingData[i], 0);
            FindNeighbourElement(r_geometry[2].Id(), r_geometry[0].Id(), id, elements_array, mSwappingData[i], 1);
            FindNeighbourElement(r_geometry[0].Id(), r_geometry[1].Id(), id, elements_array, mSwappingData[i], 2);
            for(unsigned int j = 0; j < r_geometry.size(); j++)
            {
                mSwappingData[i].IsInside[j] = IsInElementSphere(*(elements_array[i]), *(nodes_array[mSwappingData[i].OppositeNodes[j]-1]));
            }
        }


    }


    void FindNeighbourElement(unsigned int NodeId1, unsigned int NodeId2, int ElementId, ModelPart::ElementsContainerType::ContainerType const& ElementsArray, SwappingData& rSwappingData, int EdgeIndex)
    {
        const int node_index_1 = NodeId1 - 1;

        rSwappingData.NeighbourElements[EdgeIndex] = -1;


        //look for the elements around node NodeId1
        for(std::vector<int>::iterator i = mNodalNeighbourElements[node_index_1].begin() ; i != mNodalNeighbourElements[node_index_1].end() ; i++)
        {
            //look for the nodes of the neighbour faces
            Geometry<Node<3> >& r_neighbour_element_geometry = ElementsArray[*i-1]->GetGeometry();
            rSwappingData.OppositeNodes[EdgeIndex] = -1;
            for( int node_i = 0 ; node_i < static_cast<int>(r_neighbour_element_geometry.size()); node_i++)
            {
                int other_node_id = r_neighbour_element_geometry[node_i].Id();
                if ((other_node_id !=  static_cast<int>(NodeId1)) && (other_node_id !=  static_cast<int>(NodeId2)))
                {
                    rSwappingData.OppositeNodes[EdgeIndex] = other_node_id;
                    rSwappingData.OppositeEdge[EdgeIndex] = node_i;
                }

                if (r_neighbour_element_geometry[node_i].Id() == NodeId2)
                {
                    if(*i != ElementId)
                    {
                        rSwappingData.NeighbourElements[EdgeIndex] =  *i;
                    }
                }
            }
            if(rSwappingData.NeighbourElements[EdgeIndex] != -1)
                return;
        }
        return;
    }

    void PrepareForSwapping(ModelPart& rThisModelPart)
    {
        for(ModelPart::ElementIterator i_element = rThisModelPart.ElementsBegin() ; i_element != rThisModelPart.ElementsEnd() ; i_element++)
        {
        }
    }

    bool IsInElementSphere(Element& rThisElement, NodeType& rThisNode)
    {
        Element::GeometryType& r_geometry = rThisElement.GetGeometry();


        double ax = r_geometry[0].X();
        double ay = r_geometry[0].Y();

        double bx = r_geometry[1].X();
        double by = r_geometry[1].Y();

        double cx = r_geometry[2].X();
        double cy = r_geometry[2].Y();

        double a2 = ax*ax + ay*ay;
        double b2 = bx*bx + by*by;
        double c2 = cx*cx + cy*cy;

        double vx = rThisNode.X();
        double vy = rThisNode.Y();

        double v2 = vx*vx + vy*vy;


        double a11 = r_geometry[0].X()-rThisNode.X();
        double a12 = r_geometry[0].Y()-rThisNode.Y();
        double a13 = a2 - v2;

        double a21 = r_geometry[1].X()-rThisNode.X();
        double a22 = r_geometry[1].Y()-rThisNode.Y();
        double a23 = b2 - v2;

        double a31 = r_geometry[2].X()-rThisNode.X();
        double a32 = r_geometry[2].Y()-rThisNode.Y();
        double a33 = c2 - v2;

        return ( ( a11*(a22*a33-a23*a32) + a12*(a23*a31-a21*a33) + a13*(a21*a32-a22*a31) ) > 0.0 );
    }


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
    EdgeSwapping2DModeler& operator=(EdgeSwapping2DModeler const& rOther);

    /// Copy constructor.
    EdgeSwapping2DModeler(EdgeSwapping2DModeler const& rOther);


    ///@}

}; // Class EdgeSwapping2DModeler

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  EdgeSwapping2DModeler& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const EdgeSwapping2DModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_EDGE_SWAPPING_2D_MODELER_H_INCLUDED  defined 



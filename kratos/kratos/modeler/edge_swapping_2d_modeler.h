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
#include "processes/find_elements_neighbours_process.h"

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
  
  /// Short class definition.
  /** Detail class definition.
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
				}
				IsSwapCandidate = false;
				SwapWith = -1;
			}
			array_1d<int, 3> NeighbourElements;
			array_1d<int, 3> OppositeNodes;
			array_1d<int, 3> OppositeEdge;
			bool IsSwapCandidate;
			int SwapWith;
			int SwapEdge;
		};

	
	
	public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of EdgeSwapping2DModeler
      KRATOS_CLASS_POINTER_DEFINITION(EdgeSwapping2DModeler);

	  typedef Modeler BaseType;

	  typedef Point<3> PointType;
	
	  typedef Node<3> NodeType;

	  typedef Geometry<NodeType> GeometryType;

	  typedef PointerVector<NodeType> NodesVectorType;

	  typedef std::size_t SizeType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// constructor.
	  EdgeSwapping2DModeler()
	  {
	  }

      /// Destructor.
	  virtual ~EdgeSwapping2DModeler(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      


	  void Remesh(ModelPart& rThisModelPart)
	  {
	    Timer::Start("Edge Swapping");
  	  
		ModelPart::NodesContainerType::ContainerType& nodes_array = rThisModelPart.NodesArray();
		ModelPart::ElementsContainerType::ContainerType& elements_array = rThisModelPart.ElementsArray();

		const int number_of_elements = rThisModelPart.NumberOfElements(); 

	    unsigned int number_of_bad_quality_elements = MarkBadQualityElements(rThisModelPart);

	    KRATOS_WATCH(number_of_bad_quality_elements);

//  	    FindElementalNeighboursProcess find_element_neighbours_process(rThisModelPart, 2); 

		SetSwappingData(rThisModelPart);


// 	    find_element_neighbours_process.Execute(); 


		for(int i = 0 ; i < number_of_elements ; i++)
		{
			if(mBadQuality[i])
			{
				for(int j = 0 ; j < 3 ; j++)
				{
					const int neighbour_index = mSwappingData[i].NeighbourElements[j] - 1;
					if(neighbour_index >= 0)
					{
						if(mSwappingData[neighbour_index].IsSwapCandidate == false)
						{
							if(IsInElementSphere(*(elements_array[i]), *(nodes_array[mSwappingData[i].OppositeNodes[j]])))
							{
								mSwappingData[neighbour_index].IsSwapCandidate = true;
								mSwappingData[neighbour_index].SwapEdge = mSwappingData[i].OppositeEdge[j];
								mSwappingData[i].SwapWith = neighbour_index;
								mSwappingData[i].SwapEdge = j;
							}
						}
					}
				}
			}
		}

		for(int i = 0 ; i < number_of_elements ; i++)
		{
			if(mSwappingData[i].SwapWith != -1)
			{
				std::cout << "Element #" << i + 1 << " : " ;
				std::cout << mSwappingData[i].NeighbourElements << ",";
				std::cout << mSwappingData[i].OppositeNodes << ",";
				std::cout << mSwappingData[i].SwapWith;
				std::cout << std::endl;
			}
		}


		for(int i = 0 ; i < number_of_elements ; i++)
		{
			if(mSwappingData[i].SwapWith != -1)
			{
				Swap(*(elements_array[i]),*(elements_array[mSwappingData[i].SwapWith]), mSwappingData[i].SwapEdge, mSwappingData[mSwappingData[i].SwapWith].SwapEdge);
			}
		}

		Timer::Stop("Edge Swapping");
	  }
      
	  void Swap(Element& rElement1, Element& rElement2, int Edge1, int Edge2)
	  {

		  int next2[3] = {2,0,1};

		  rElement1.GetGeometry()(next2[Edge1]) = rElement2.GetGeometry()(Edge2);
		  rElement2.GetGeometry()(next2[Edge2]) = rElement1.GetGeometry()(Edge1);

	  }

	  unsigned int MarkBadQualityElements(ModelPart& rThisModelPart)
	  {
	    unsigned int counter = 0;
	    unsigned int marked = 0;

	    const double threshold = 0.1;

	    if(mBadQuality.size() != rThisModelPart.NumberOfElements())
	      mBadQuality.resize(rThisModelPart.NumberOfElements());
	    
	    // To be parallelized.
	    for(ModelPart::ElementIterator i_element = rThisModelPart.ElementsBegin() ; i_element != rThisModelPart.ElementsEnd() ; i_element++)
	      {
		if(ElementQuality(*i_element) < threshold)
		  {
		    marked++;
		    mBadQuality[counter] = true;
		  }
		counter++;
	      }

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
		  const int number_of_nodes = rThisModelPart.NumberOfNodes(); 
		  const int number_of_elements = rThisModelPart.NumberOfElements(); 

		  if(mNodalNeighbourElements.size() != number_of_nodes)
			  mNodalNeighbourElements.resize(number_of_nodes);
		  else
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
		  ModelPart::ElementsContainerType::ContainerType& elements_array = rThisModelPart.ElementsArray();
		  const int number_of_nodes = rThisModelPart.NumberOfNodes(); 
		  const int number_of_elements = rThisModelPart.NumberOfElements(); 

		  FindNodalNeighbours(rThisModelPart);

		  if(mSwappingData.size() != number_of_elements)
			  mSwappingData.resize(number_of_elements);
		  else
			  for(int i = 0 ; i < number_of_elements ; i++)
				  mSwappingData[i].Reset();

		  for(int i = 0 ; i < number_of_elements ; i++)
		  {
			  Element::GeometryType& r_geometry = elements_array[i]->GetGeometry();
			  int id = elements_array[i]->Id();
			  for(unsigned int j = 0; j < r_geometry.size(); j++)
			  {
				FindNeighbourElement(r_geometry[1].Id(), r_geometry[2].Id(), id, elements_array, mSwappingData[i], 0);
				FindNeighbourElement(r_geometry[2].Id(), r_geometry[0].Id(), id, elements_array, mSwappingData[i], 1);
				FindNeighbourElement(r_geometry[0].Id(), r_geometry[1].Id(), id, elements_array, mSwappingData[i], 2);
			  }
		  }


	  }


	  void FindNeighbourElement(unsigned int NodeId1, unsigned int NodeId2, int ElementId, ModelPart::ElementsContainerType::ContainerType const& ElementsArray, SwappingData& rSwappingData, int EdgeIndex)
	  {
		  const int node_index_1 = NodeId1 - 1;

		  rSwappingData.NeighbourElements[EdgeIndex] = -1;


		  //look for the elements around node NodeId1
 		  for(std::vector<int>::iterator i = mNodalNeighbourElements[node_index_1].begin() ; i != mNodalNeighbourElements[node_index_1].end() ; i++)
		  {	//look for the nodes of the neighbour faces
			  Geometry<Node<3> >& r_neighbour_element_geometry = ElementsArray[*i-1]->GetGeometry();
			  rSwappingData.OppositeNodes[EdgeIndex] = -1;
			  for( unsigned int node_i = 0 ; node_i < r_neighbour_element_geometry.size(); node_i++) 
			  {	
				  int other_node_id = r_neighbour_element_geometry[node_i].Id();
				  if ((other_node_id != NodeId1) && (other_node_id != NodeId2))
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

	    double a11 = r_geometry[0].X()-rThisNode.X();
	    double a12 = r_geometry[0].Y()-rThisNode.Y();
	    double a13 = 1.00;
	    
	    double a21 = r_geometry[1].X()-rThisNode.X();
	    double a22 = r_geometry[1].Y()-rThisNode.Y();
	    double a23 = 1.00;
	    
	    double a31 = r_geometry[2].X()-rThisNode.X();
	    double a32 = r_geometry[2].Y()-rThisNode.Y();
	    double a33 = 1.00;
	    
	    return ( ( a11*(a22*a33-a23*a32) + a12*(a23*a31-a21*a33) + a13*(a21*a32-a22*a31) ) > 0.0 );
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
        
      std::vector<int> mBadQuality;
      std::vector<std::vector<int> > mNodalNeighbourElements;
      std::vector<SwappingData > mSwappingData;
        
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



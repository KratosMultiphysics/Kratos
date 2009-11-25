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


#if !defined(KRATOS_CARTESIAN_MESH_GENERATOR_H_INCLUDED )
#define  KRATOS_CARTESIAN_MESH_GENERATOR_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "modeler/modeler.h"
#include "spatial_containers/spatial_containers.h"

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
	class CartesianMeshGeneratorModeler : public Modeler
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of CartesianMeshGeneratorModeler
      KRATOS_CLASS_POINTER_DEFINITION(CartesianMeshGeneratorModeler);

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
	  CartesianMeshGeneratorModeler(ModelPart& rSourceModelPart, double ElementSize) :
		mrModelPart(rSourceModelPart), mElementSize(ElementSize)
	  {
	  }

      /// Destructor.
	  virtual ~CartesianMeshGeneratorModeler(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
	  void GenerateMesh(ModelPart& rThisModelPart, Element const& rReferenceElement)
	  {
	    const unsigned int dimension = rReferenceElement.GetGeometry().Dimension();

	    KRATOS_WATCH(dimension);

	  	  Timer::Start("Generating Mesh");

		  CalculateBoundingBox(mrModelPart, mMinPoint, mMaxPoint);
		  CalculateDivisionNumbers();

		  unsigned int start_node_id = mrModelPart.NumberOfNodes() + 1;
		  unsigned int start_element_id = 1;
 		  //unsigned int segment_number_1 =  mSegmentsNumber[0] + 1;
 		  //unsigned int segment_number_2 =  mSegmentsNumber[1] + 1;
 		  //unsigned int segment_number_3 =  mSegmentsNumber[2] + 1;



		  const unsigned int number_of_nodes =  mDivisionsNumber[0] * mDivisionsNumber[1] * mDivisionsNumber[2];
		  //const unsigned int number_of_nodes =  segment_number_1 * segment_number_2 * segment_number_3;

		  const unsigned int number_of_elements =  mSegmentsNumber[0] * mSegmentsNumber[1] * mSegmentsNumber[2];

		  KRATOS_WATCH(number_of_nodes);

		  KRATOS_WATCH(number_of_elements);

		  CalculateBoundaryIntersections(mrModelPart);

		  CalculateIsInside(mrModelPart);

		  CalculateNormals();

	  	  Timer::Start("Generating Nodes");

		  double x0 =  mMinPoint.X();
		  double y0 =  mMinPoint.Y();
		  double z0 =  mMinPoint.Z();

		  ModelPart::NodesContainerType::ContainerType& nodes_array = rThisModelPart.NodesArray();
		  ModelPart::NodesContainerType::ContainerType temp_nodes_array(number_of_nodes);
				  
		  ModelPart::ElementsContainerType::ContainerType& elements_array = rThisModelPart.ElementsArray();
		  ModelPart::ElementsContainerType::ContainerType temp_elements_array(number_of_elements);
				  
		  for(unsigned int k = 0 ; k < mDivisionsNumber[2] ; k++)
			  for(unsigned int j = 0 ; j < mDivisionsNumber[1] ; j++)
				  for(unsigned int i = 0 ; i < mDivisionsNumber[0] ; i++)
					  temp_nodes_array[NodeIndex(i,j,k)] = NodeType::Pointer(new NodeType(start_node_id++, x0 + i * mElementSize, y0 + j * mElementSize, z0 + k * mElementSize));

		  unsigned int number_of_active_nodes = 0;

		  for(unsigned int i = 0 ; i < number_of_nodes ; i++)
			  if(mIsInside[i])
				  number_of_active_nodes++;

		  nodes_array.resize(number_of_active_nodes);

		  unsigned int index = 0;

		  for(unsigned int i = 0 ; i < number_of_nodes ; i++)
			  if(mIsInside[i])
				nodes_array[index++]=temp_nodes_array[i];

	  	  Timer::Stop("Generating Nodes");

	  	  Timer::Start("Generating Elements");

		  Element::NodesArrayType element_nodes(4);

 		  if(dimension == 2)
		  {
			  unsigned int number_of_active_elements = 0;
			  for(unsigned int j = 0 ; j < mSegmentsNumber[1] ; j++)
				  for(unsigned int i = 0 ; i < mSegmentsNumber[0] ; i++)
				  {
					  if(mIsInside[NodeIndex(i,j,0)] & mIsInside[NodeIndex(i+1,j,0)] & mIsInside[NodeIndex(i+1,j+1,0)] & mIsInside[NodeIndex(i,j+1,0)])
						  number_of_active_elements++;
				  }

			  elements_array.resize(number_of_active_elements);

			  unsigned int counter = 0;

			  for(unsigned int j = 0 ; j < mSegmentsNumber[1] ; j++)
				  for(unsigned int i = 0 ; i < mSegmentsNumber[0] ; i++)
				  {
					  if(mIsInside[NodeIndex(i,j,0)] & mIsInside[NodeIndex(i+1,j,0)] & mIsInside[NodeIndex(i+1,j+1,0)] & mIsInside[NodeIndex(i,j+1,0)])
					  {
						  element_nodes(0) = temp_nodes_array[NodeIndex(i,j,0)];
						  element_nodes(1) = temp_nodes_array[NodeIndex(i+1,j,0)];
						  element_nodes(2) = temp_nodes_array[NodeIndex(i+1,j+1,0)];
						  element_nodes(3) = temp_nodes_array[NodeIndex(i,j+1,0)];

						  elements_array[counter++] = rReferenceElement.Create(start_element_id++, element_nodes, rReferenceElement.pGetProperties());
						  //elements_array[ElementIndex(i,j,0)] = rReferenceElement.Create(start_element_id++, element_nodes, rReferenceElement.pGetProperties());
					  }
				  }
				  for(ModelPart::ElementIterator i_element = mrModelPart.ElementsBegin() ; i_element != mrModelPart.ElementsEnd() ; i_element++)
				  {
					  Element::GeometryType& r_geometry = i_element->GetGeometry();

					  PointType point1 = r_geometry[0] + mNormals[r_geometry[0].Id()-1];
					  PointType point2 = r_geometry[1] + mNormals[r_geometry[1].Id()-1];

					  unsigned int index1 = FindNearestNodeIndex(point1,mNormals[r_geometry[0].Id()-1]);
					  unsigned int index2 = FindNearestNodeIndex(point2,mNormals[r_geometry[1].Id()-1]);

					  element_nodes(0) = r_geometry(0);
					  element_nodes(1) = r_geometry(1);
					  element_nodes(2) = temp_nodes_array[index2];
					  element_nodes(3) = temp_nodes_array[index1];

					  elements_array.push_back(rReferenceElement.Create(start_element_id++, element_nodes, rReferenceElement.pGetProperties()));
				  }
		  }

		  // TODO: assigning nodes and elements to the modelpart
		  // TODO: adding source boundary mesh nods to modelpart

	  	  Timer::Stop("Generating Elements");

		  

	  	  Timer::Stop("Generating Mesh");
	  }


	  void CalculateNormals()
	  {
		  array_1d<double,3> zero = ZeroVector(3);

		  if(mNormals.size() != mrModelPart.NumberOfNodes())
			  mNormals.resize(mrModelPart.NumberOfNodes(), zero);
		  else
			  std::fill(mNormals.begin(), mNormals.end(), zero);

		  double coefficient = mElementSize / 2.00;

		  for(ModelPart::ElementIterator i_element = mrModelPart.ElementsBegin() ; i_element != mrModelPart.ElementsEnd() ; i_element++)
		  {
			  Element::GeometryType& r_geometry = i_element->GetGeometry();
			  array_1d<double,3> normal;

			  normal[0] =    r_geometry[1].Y() - r_geometry[0].Y();
			  normal[1] = - (r_geometry[1].X() - r_geometry[0].X());

			  normal *= coefficient / r_geometry.Length();


			  mNormals[r_geometry[0].Id()-1] += normal ;
			  mNormals[r_geometry[1].Id()-1] += normal ;
		  }
		  for(ModelPart::NodeIterator i_node = mrModelPart.NodesBegin() ; i_node != mrModelPart.NodesEnd() ; i_node++)
			noalias(i_node->FastGetSolutionStepValue(NORMAL)) = mNormals[i_node->Id()-1] ;
	  }

	  unsigned int FindNearestNodeIndex(PointType& rThisPoint, array_1d<double,3>& rNormal)
	  {
		  double x = (rThisPoint.X() - mMinPoint.X()) / mElementSize;
		  double y = (rThisPoint.Y() - mMinPoint.Y()) / mElementSize;

		  unsigned int i = static_cast<unsigned int>(x);
		  unsigned int j = static_cast<unsigned int>(y);

		  if(mIsInside[NodeIndex(i,j,0)])
			return NodeIndex(i,j,0);


		  if(rNormal[0] >= 0.00)
		  {
			  if(rNormal[1] >= 0.00)
			  {
				  if(rNormal[0] > rNormal[1])
					  i++;
				  else
					  j++;
			  }
			  else
			  {
				  if(rNormal[0] > -rNormal[1])
					  i++;
				  else
					  j--;
			  }
		  }
		  else
		  {
			  if(rNormal[1] >= 0.00)
			  {
				  if(-rNormal[0] > rNormal[1])
					  i--;
				  else
					  j++;
			  }
			  else
			  {
				  if(-rNormal[0] > -rNormal[1])
					  i--;
				  else
					  j--;
			  }
		  }


		  return NodeIndex(i,j,0);

	  }

	  void CalculateIsInside(ModelPart& rThisModelPart)
	  {
		  const unsigned int number_of_nodes =  mDivisionsNumber[0] * mDivisionsNumber[1] * mDivisionsNumber[2];

		  if(mIsInside.size() != number_of_nodes)
		  {
			  mIsInside.resize(number_of_nodes, 0);
		  }
		  else
			  std::fill(mIsInside.begin(),mIsInside.end(), 0);

		  int size = mSegmentsNumber[1] + 1;

		  for(int j = 0 ; j < size ; j++)
		  {
			  bool is_inside = false;
			  std::vector<double>& j_intersections = mIntersections[j];

			  for(std::vector<double>::iterator j_x = j_intersections.begin() ; j_x != j_intersections.end() ; j_x++)
			  {
				  std::vector<double>::iterator start = j_x++;

				unsigned int i_start = static_cast<unsigned int>(*start / mElementSize);
				unsigned int i_end = static_cast<unsigned int>(*j_x / mElementSize);

				for(unsigned int i = i_start + 1 ; i < i_end ; i++)
					mIsInside[NodeIndex(i,j,0)]=true;
			  }
		  }
	  }

	  void CalculateBoundaryIntersections(ModelPart& rThisModelPart)
	  {
		  std::vector<int> number_of_intersections(mSegmentsNumber[1] + 1, 0); // number of intersections per row for mSegmentsNumber[1] + 1 rows in x direction
		  if(mIntersections.size() != mSegmentsNumber[1] + 1)
			  mIntersections.resize(mSegmentsNumber[1] + 1); // intersection points x coordinate per row for mSegmentsNumber[1] + 1 rows 

		  for(ModelPart::ElementIterator i_element = rThisModelPart.ElementsBegin() ; i_element != rThisModelPart.ElementsEnd() ; i_element++)
		  {
			  Element::GeometryType& r_geometry = i_element->GetGeometry();

			  double x1;
			  double x2;
			  double y1;
			  double y2;
			  if(r_geometry[0].Y() < r_geometry[1].Y())
			  {
				  x1 = (r_geometry[0].X() - mMinPoint.X());
				  x2 = (r_geometry[1].X() - mMinPoint.X());
				  y1 = (r_geometry[0].Y() - mMinPoint.Y());
				  y2 = (r_geometry[1].Y() - mMinPoint.Y());
			  }
			  else
			  {
				  x1 = (r_geometry[1].X() - mMinPoint.X());
				  x2 = (r_geometry[0].X() - mMinPoint.X());
				  y1 = (r_geometry[1].Y() - mMinPoint.Y());
				  y2 = (r_geometry[0].Y() - mMinPoint.Y());
			  }

			  unsigned int i_start = static_cast<unsigned int>(y1 / mElementSize);
			  unsigned int i_end = static_cast<unsigned int>(y2 / mElementSize);
				  
			  if (i_start*mElementSize < y1)
					  i_start++;

			  double m = 0.00;
			  if(y1 != y2)
				m = (x2 - x1) / (y2 - y1);

			  double delta_x = m*mElementSize;
			  double x = x1 + (i_start * mElementSize - y1) * m ;
			  for(unsigned int i = i_start  ; i <= i_end ; i++)
			  {
				  number_of_intersections[i]++;
				  mIntersections[i].push_back(x);
				  x += delta_x;
			  }
		  }

		  for(int i = 0 ; i < mIntersections.size() ; i++)
			  std::sort(mIntersections[i].begin(), mIntersections[i].end());
	  }

	  void CalculateBoundingBox(ModelPart& rThisModelPart, Point<3>& rMinPoint, Point<3>& rMaxPoint)
	  {
		  if(rThisModelPart.NumberOfElements() == 0)
		  {
			  rMinPoint = PointType();
			  rMaxPoint = PointType();
			  return;
		  }

		  if(rThisModelPart.ElementsBegin()->GetGeometry().empty())
		  {
			  rMinPoint = PointType();
			  rMaxPoint = PointType();
			  return;
		  }

		  
		  rMinPoint = rThisModelPart.ElementsBegin()->GetGeometry()[0];
		  rMaxPoint = rMinPoint;

		  for(ModelPart::ElementIterator i_element = rThisModelPart.ElementsBegin() ; 
			  i_element != rThisModelPart.ElementsEnd() ; i_element++)
			  for(Element::GeometryType::iterator i_point = i_element->GetGeometry().begin() ; i_point != i_element->GetGeometry().end() ; i_point++)
				  for(unsigned int i = 0 ; i < PointType::Dimension() ; i++)
				  {
					  if(rMinPoint[i] > (*i_point)[i])
						  rMinPoint[i] = (*i_point)[i];

					  if(rMaxPoint[i] < (*i_point)[i])
						  rMaxPoint[i] = (*i_point)[i];
				  }
	  }


	  void CalculateDivisionNumbers()
	  {
		  if(mElementSize == 0.00)
			  return;

		  for(unsigned int i = 0 ; i < PointType::Dimension() ; i++)
		  {
			  double delta = mMaxPoint[i] - mMinPoint[i];
			  int segments_number = static_cast<int>(delta / mElementSize);

			  if (((segments_number * mElementSize) < delta))
				  segments_number++;
				
			  mSegmentsNumber[i] = segments_number;
			  KRATOS_WATCH(mSegmentsNumber[i]);

			  mDivisionsNumber[i] = segments_number + 1;

			  if(mSegmentsNumber[i] == 0)
				  mSegmentsNumber[i]++;
		  }

	  }

	  virtual void GenerateNodes(ModelPart& ThisModelPart)
	  {
		  //std::vector<PointType> 
	  }
      

		unsigned int ElementIndex(unsigned int i, unsigned int j, unsigned int k)
		{
			return i + mDivisionsNumber[0] * j + mDivisionsNumber[0] * mDivisionsNumber[1] * k;
		}

		unsigned int NodeIndex(unsigned int i, unsigned int j, unsigned int k)
		{
			return i + mSegmentsNumber[0] * j + mSegmentsNumber[0] * mSegmentsNumber[1] * k;
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
		  return "CartesianMeshGeneratorModeler";
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
        
	  ModelPart& mrModelPart;

	  double mElementSize;
	  PointType mMinPoint;
	  PointType mMaxPoint;

	  unsigned int mSegmentsNumber[3];

	  unsigned int mDivisionsNumber[3];

	  std::vector<std::vector<double> > mIntersections;
	  std::vector<int> mIsInside;

	  std::vector<array_1d<double,3> > mNormals;
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        

	  void GenerateNodes(ModelPart& ThisModelPart, GeometryType& rGeometry, SizeType NumberOfSegments, SizeType StartNodeId)
	  {
		  double x1 = rGeometry[0][0];
		  double y1 = rGeometry[0][1];
		  double z1 = rGeometry[0][2];
		  double x2 = rGeometry[1][0];
		  double y2 = rGeometry[1][1];
		  double z2 = rGeometry[1][2];

		  double dx = (x2 - x1) / NumberOfSegments;
		  double dy = (y2 - y1) / NumberOfSegments;
		  double dz = (z2 - z1) / NumberOfSegments;

		  for(SizeType i = 1 ; i < NumberOfSegments - 1 ; i++)
		  {
			  ThisModelPart.CreateNewNode(StartNodeId++, x1 + i * dx, y1 + i * dy, z1 + i * dz);
		  }
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
      CartesianMeshGeneratorModeler& operator=(CartesianMeshGeneratorModeler const& rOther);

      /// Copy constructor.
      CartesianMeshGeneratorModeler(CartesianMeshGeneratorModeler const& rOther);

        
      ///@}    
        
    }; // Class CartesianMeshGeneratorModeler 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    CartesianMeshGeneratorModeler& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const CartesianMeshGeneratorModeler& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_CARTESIAN_MESH_GENERATOR_H_INCLUDED  defined 



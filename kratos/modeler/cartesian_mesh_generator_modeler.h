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

		  unsigned int start_node_id = 1;
		  unsigned int start_element_id = 1;
 		  unsigned int segment_number_1 =  mSegmentsNumber[0] + 1;
 		  unsigned int segment_number_2 =  mSegmentsNumber[1] + 1;
 		  unsigned int segment_number_3 =  mSegmentsNumber[2] + 1;

		  for(int i = 0 ; i < 3 ; i++)
		    if(mSegmentsNumber[i] == 0)
		      mSegmentsNumber[i]++;


		  const unsigned int number_of_nodes =  segment_number_1 * segment_number_2 * segment_number_3;

		  const unsigned int number_of_elements =  mSegmentsNumber[0] * mSegmentsNumber[1] * mSegmentsNumber[2];

		  KRATOS_WATCH(number_of_nodes);

		  KRATOS_WATCH(number_of_elements);

	  	  Timer::Start("Generating Nodes");

		  double x0 =  mMinPoint.X();
		  double y0 =  mMinPoint.Y();
		  double z0 =  mMinPoint.Z();

		  ModelPart::NodesContainerType::ContainerType& nodes_array = rThisModelPart.NodesArray();
		  nodes_array.resize(number_of_nodes);
				  
		  ModelPart::ElementsContainerType::ContainerType& elements_array = rThisModelPart.ElementsArray();
		  elements_array.resize(number_of_elements);
				  
		  for(unsigned int i = 0 ; i < segment_number_1 ; i++)
		    for(unsigned int j = 0 ; j < segment_number_2 ; j++)
		      for(unsigned int k = 0 ; k < segment_number_3 ; k++)
			nodes_array[i*segment_number_2*segment_number_3+j*segment_number_3+k] = NodeType::Pointer(new NodeType(start_node_id++, x0 + i * mElementSize, y0 + j * mElementSize, z0 + k * mElementSize));

	  	  Timer::Stop("Generating Nodes");

	  	  Timer::Start("Generating Elements");

		  Element::NodesArrayType element_nodes(4);

 		  if(dimension == 2)
 		    for(unsigned int i = 0 ; i < mSegmentsNumber[0] ; i++)
 		      for(unsigned int j = 0 ; j < mSegmentsNumber[1] ; j++)
			{
			  element_nodes(0) = nodes_array[i* segment_number_2 + j];
			  element_nodes(1) = nodes_array[(i+1)*(mSegmentsNumber[1] + 1) + j];
			  element_nodes(2) = nodes_array[(i+1)*(mSegmentsNumber[1] + 1) + j+1];
			  element_nodes(3) = nodes_array[(i)*(mSegmentsNumber[1] + 1) + j+1];

			  elements_array[i* mSegmentsNumber[1] + j] = rReferenceElement.Create(start_element_id++, element_nodes, rReferenceElement.pGetProperties());
			}

	  	  Timer::Stop("Generating Elements");

		  

	  	  Timer::Stop("Generating Mesh");
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
			  for(GeometryType::iterator i_point = i_element->GetGeometry().begin() ; i_point != i_element->GetGeometry().end() ; i_point++)
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
		  }
	  }

	  virtual void GenerateNodes(ModelPart& ThisModelPart)
	  {
		  //std::vector<PointType> 
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



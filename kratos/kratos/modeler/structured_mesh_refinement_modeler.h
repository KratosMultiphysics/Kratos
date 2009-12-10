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


#if !defined(KRATOS_STRUCTURED_MESH_REFINEMENT_H_INCLUDED )
#define  KRATOS_STRUCTURED_MESH_REFINEMENT_H_INCLUDED



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
	class StructuredMeshRefinementModeler : public Modeler
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of StructuredMeshRefinementModeler
      KRATOS_CLASS_POINTER_DEFINITION(StructuredMeshRefinementModeler);

	  typedef Modeler BaseType;

	  typedef Node<3> NodeType;

	  typedef PointerVector<NodeType> NodesVectorType;

	  typedef std::size_t SizeType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
	  StructuredMeshRefinementModeler(){}

      /// Destructor.
	  virtual ~StructuredMeshRefinementModeler(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
	  void GenerateMesh(ModelPart& rThisModelPart, Element const& rReferenceElement)
	  {
		ModelPart::ElementsContainerType old_elements;
		old_elements.swap(rThisModelPart.Elements());

		ModelPart::ConditionsContainerType old_conditions;
		old_conditions.swap(rThisModelPart.Conditions());

		const double tolerance = 1e-9;

		// Creating a bins for searching the nodes to be collapsed
		//typedef Bucket<3, NodeType, ModelPart::NodesContainerType::ContainerType> bucket_type;
		typedef BinsDynamic<3, NodeType, ModelPart::NodesContainerType::ContainerType> bins_type;

		bins_type nodes_bins(rThisModelPart.Nodes().ptr_begin(), rThisModelPart.Nodes().ptr_end()); 
		ModelPart::NodesContainerType::ContainerType founded_nodes(1);


		  for(ModelPart::ElementIterator i_element = old_elements.begin() ; i_element != old_elements.end() ; i_element++)
		  {
			  // Getting the number of divisions;
			  vector<int>& number_of_divisions = i_element->GetValue(STRUCTURED_MESH_DIVISIONS);

			  //int node_id = rThisModelPart.Nodes().size() + 1;
			  int node_id = rThisModelPart.Nodes().back().Id() + 1;
			  int start_element_id = 1;
			  if(rThisModelPart.Elements().size() != 0)
				  start_element_id = rThisModelPart.Elements().back().Id() + 1;


			  // Calcualting the size of the nodes after refining this element
			  SizeType nodes_vector_size = number_of_divisions[0] + 1;
			  for(SizeType i = 1 ; i < number_of_divisions.size() ; i++)
				  nodes_vector_size *= number_of_divisions[i] + 1;

			  if(nodes_vector_size > i_element->GetGeometry().size()) // if we need to create more nodes than actual we have to refine
			  {
				  NodesVectorType nodes_vector;
				  NodesVectorType surface_nodes_vector;
				  GenerateNodes(*i_element,  nodes_vector, surface_nodes_vector, number_of_divisions);
				  for(SizeType i = 0 ; i < nodes_vector.size() ; i++)
				  {
					  NodeType::Pointer p_node = nodes_vector(i);
					  if(nodes_bins.SearchInRadius(*p_node, tolerance, founded_nodes.begin(), 1) > 0) // It founds a node in the radius of tolerance so there was a node there
					  {
						  nodes_vector(i) = founded_nodes[0];
					  }
					  else
					  {
						  p_node->SetId(node_id++);
					  	  rThisModelPart.AddNode(p_node);
						  nodes_bins.AddPoint(p_node);
					  }
				  }
				  GenerateElements(rThisModelPart, *i_element,  nodes_vector, surface_nodes_vector, number_of_divisions, start_element_id);

			  }
			  else
			  {
				  i_element->SetId(start_element_id);
				  rThisModelPart.Elements().push_back(*(i_element.base()));
			  }

		  }

		  for(ModelPart::ConditionIterator i_condition = old_conditions.begin() ; i_condition != old_conditions.end() ; i_condition++)
		  {
			  // Getting the number of divisions;
			  vector<int>& number_of_divisions = i_condition->GetValue(STRUCTURED_MESH_DIVISIONS);

			  int node_id = rThisModelPart.Nodes().size() + 1;
			  int start_condition_id = rThisModelPart.Conditions().size() + 1;
			  // Calcualting the size of the nodes after refining this condition
			  SizeType nodes_vector_size = number_of_divisions[0] + 1;
			  for(SizeType i = 1 ; i < i_condition->GetGeometry().LocalSpaceDimension() ; i++)
				  nodes_vector_size *= number_of_divisions[i] + 1;

			  if(nodes_vector_size > i_condition->GetGeometry().size()) // if we need to create more nodes than actual we have to refine
			  {
				  NodesVectorType nodes_vector;
				  NodesVectorType surface_nodes_vector;
				  GenerateNodes(*i_condition,  nodes_vector, surface_nodes_vector, number_of_divisions);
				  for(SizeType i = 0 ; i < nodes_vector.size() ; i++)
				  {
					  NodeType::Pointer p_node = nodes_vector(i);
					  if(nodes_bins.SearchInRadius(*p_node, tolerance, founded_nodes.begin(), 1) > 0) // It founds a node in the radius of tolerance so there was a node there
					  {
						  nodes_vector(i) = founded_nodes[0];
					  }
					  else
					  {
						  p_node->SetId(node_id++);
					  	  rThisModelPart.AddNode(p_node);
						  nodes_bins.AddPoint(p_node);
					  }
				  }
				  GenerateConditions(rThisModelPart, *i_condition,  nodes_vector, surface_nodes_vector, number_of_divisions, start_condition_id);

			  }
			  else
			  {
				  i_condition->SetId(start_condition_id);
				  rThisModelPart.Conditions().push_back(*(i_condition.base()));
			  }
		  }
	  }

	  template<class ComponentType>
	  void GenerateNodes(ComponentType& rThisComponent, NodesVectorType& VolumeNodesVector, NodesVectorType& SurfaceNodesVector, vector<int>& number_of_divisions)
	  {
		  SizeType local_dimension = rThisComponent.GetGeometry().LocalSpaceDimension();
		  Element::GeometryType::CoordinatesArrayType local_coordinates(local_dimension);
		  Element::GeometryType::CoordinatesArrayType global_coordinates;
		  
		  if(local_dimension == 2)
		  {
			  for(int i = 0 ; i <= number_of_divisions[0] ; i++)
			  {
				  local_coordinates[0] = double(-1.00) + double(2.00 * i) / number_of_divisions[0];
				  for(int j = 0 ; j <= number_of_divisions[1] ; j++)
				  {
					  local_coordinates[1] = double(-1.00) + double(2.00 * j) / number_of_divisions[1];
					  GenerateNode(rThisComponent, VolumeNodesVector, SurfaceNodesVector, local_coordinates);
				  }
			  }
		  }
		  if(local_dimension == 3)
		  {
			  for(int i = 0 ; i <= number_of_divisions[0] ; i++)
			  {
				  local_coordinates[0] = double(-1.00) + double(2.00 * i) / number_of_divisions[0];
				  for(int j = 0 ; j <= number_of_divisions[1] ; j++)
				  {
					  local_coordinates[1] = double(-1.00) + double(2.00 * j) / number_of_divisions[1];
					  for(int k = 0 ; k <= number_of_divisions[2] ; k++)
					  {
						  local_coordinates[2] = double(-1.00) + double(2.00 * k) / number_of_divisions[2];
						  GenerateNode(rThisComponent, VolumeNodesVector, SurfaceNodesVector, local_coordinates);
					  }
				  }
			  }
		  }
	  }

	  template<class ComponentType>
	  void GenerateNode(ComponentType& rThisComponent, NodesVectorType& VolumeNodesVector, NodesVectorType& SurfaceNodesVector, Element::GeometryType::CoordinatesArrayType& rLocalCoordinates)
	  {
		  SizeType local_dimension = rThisComponent.GetGeometry().LocalSpaceDimension();
		  typename ComponentType::GeometryType::CoordinatesArrayType global_coordinates;

		  typename ComponentType::GeometryType& r_geometry = rThisComponent.GetGeometry();

		  const SizeType components_nodes_number = r_geometry.size();

		  Vector shape_functions_values(components_nodes_number);

		  // Getting the shape function value for given local coordinates
		  for(SizeType h = 0 ; h < components_nodes_number ; h++)
			  shape_functions_values[h] = r_geometry.ShapeFunctionValue(h, rLocalCoordinates);


		  // Calculating the GlobalCoordinates
		  r_geometry.GlobalCoordinates(global_coordinates, rLocalCoordinates);


		  // Interpolating the Nodal data
		  SizeType ndoal_data_size = r_geometry[0].SolutionStepData().TotalDataSize() / sizeof(double);

		  typedef VariablesListDataValueContainer::BlockType block_type;

		  block_type* nodal_data = new block_type[ndoal_data_size];

		  // Initializing to zero
		  for(SizeType i = 0 ; i < ndoal_data_size ; i++)
				nodal_data[i] = block_type();
			  
		  for(SizeType i = 0 ; i < components_nodes_number ; i++)
		  {
			  block_type* p_data = r_geometry[i].SolutionStepData().Data();
			  for(SizeType j = 0 ; j < ndoal_data_size ; j++)
			  {
				  nodal_data[j] += shape_functions_values[i] * p_data[j];
			  }
		  }

		  // Creating the new node
		  NodeType::Pointer p_node(new NodeType(0, global_coordinates[0], global_coordinates[1], global_coordinates[2], r_geometry[0].pGetVariablesList(), nodal_data, r_geometry[0].GetBufferSize())); 

		  SizeType number_of_dofs = r_geometry[0].GetDofs().size();

		  for(NodeType::DofsContainerType::iterator i_dof = r_geometry[0].GetDofs().begin() ; i_dof != r_geometry[0].GetDofs().end() ; i_dof++)
		  {
			  VariableData const& r_dof_variable = i_dof->GetVariable();
			  double dof_is_fixed = double();

			  for(SizeType i = 0 ; i < components_nodes_number ; i++)
			  {
				  dof_is_fixed += shape_functions_values[i] * static_cast<double>(r_geometry[i].IsFixed(r_dof_variable));				  
			  }
			  if(dof_is_fixed > 0.999999)
			  {
				  p_node->pAddDof(*i_dof)->FixDof();
			  }
		  }

			

		  VolumeNodesVector.push_back(p_node);
	  }

	  void GenerateElements(ModelPart& rThisModelPart, Element& rThisElement, NodesVectorType& VolumeNodesVector, NodesVectorType& SurfaceNodesVector, vector<int>& number_of_divisions, SizeType StartElementId)
	  {
		  SizeType local_dimension = rThisElement.GetGeometry().LocalSpaceDimension();
		  Element::NodesArrayType element_nodes(rThisElement.GetGeometry().size());
		  SizeType nodes_number_i = number_of_divisions[0] + 1; // Number of nodes in i direction
		  SizeType nodes_number_j = number_of_divisions[1] + 1; // Number of nodes in j direction
		  SizeType nodes_number_k = number_of_divisions[2] + 1; // Number of nodes in k direction

		  if(local_dimension == 3)
		  {
			  for(int i = 0 ; i < number_of_divisions[0] ; i++)
			  {
				  for(int j = 0 ; j < number_of_divisions[1] ; j++)
				  {
					  for(int k = 0 ; k < number_of_divisions[2] ; k++)
					  {
						  // This is done only for Hexahedra
						  SizeType local_id = i * nodes_number_k * nodes_number_j + j * nodes_number_k + k;
						  element_nodes(0) = VolumeNodesVector(local_id);
						  element_nodes(1) = VolumeNodesVector(local_id + nodes_number_j * nodes_number_k);
						  element_nodes(2) = VolumeNodesVector(local_id + nodes_number_j * nodes_number_k + nodes_number_k);
						  element_nodes(3) = VolumeNodesVector(local_id + nodes_number_k);
						  element_nodes(4) = VolumeNodesVector(local_id + 1);
						  element_nodes(5) = VolumeNodesVector(local_id + nodes_number_j * nodes_number_k + 1);
						  element_nodes(6) = VolumeNodesVector(local_id + nodes_number_j * nodes_number_k + nodes_number_k + 1);
						  element_nodes(7) = VolumeNodesVector(local_id + nodes_number_k + 1);

						  rThisModelPart.Elements().push_back(rThisElement.Create(local_id + StartElementId, element_nodes, rThisElement.pGetProperties()));
					
						  
					  }
				  }
			  }
		  }
	  }
      
	  void GenerateConditions(ModelPart& rThisModelPart, Condition& rThisCondition, NodesVectorType& VolumeNodesVector, NodesVectorType& SurfaceNodesVector, vector<int>& number_of_divisions, SizeType StartConditionId)
	  {
		  SizeType local_dimension = rThisCondition.GetGeometry().LocalSpaceDimension();
		  Condition::NodesArrayType condition_nodes(rThisCondition.GetGeometry().size());
		  SizeType nodes_number_i = number_of_divisions[0] + 1; // Number of nodes in i direction
		  SizeType nodes_number_j = number_of_divisions[1] + 1; // Number of nodes in j direction

		  if(local_dimension == 2)
		  {
			  for(int i = 0 ; i < number_of_divisions[0] ; i++)
			  {
				  for(int j = 0 ; j < number_of_divisions[1] ; j++)
				  {
						  // This is done only for quadrilateral
						  SizeType local_id = i * nodes_number_j + j;
						  condition_nodes(0) = VolumeNodesVector(local_id);
						  condition_nodes(1) = VolumeNodesVector(local_id + nodes_number_j);
						  condition_nodes(2) = VolumeNodesVector(local_id + nodes_number_j + 1);
						  condition_nodes(3) = VolumeNodesVector(local_id + 1);

						  rThisModelPart.Conditions().push_back(rThisCondition.Create(local_id + StartConditionId, condition_nodes, rThisCondition.pGetProperties()));
					
				  }
			  }
		  }
	  }
      
	  void GenerateNodes(ModelPart& ThisModelPart, SizeType NumberOfSegments)
	  {
	  }

/*
		void Interpolate(Element& rThisElement 
				ModelPart::ElementsContainerType::iterator el_it, 
				const array_1d<double,3>& N, 
				int step_data_size,
      				Node<3>::Pointer pnode)
		{
			//Geometry element of the rOrigin_ModelPart
			Geometry< Node<3> >& geom = el_it->GetGeometry();
			
			unsigned int buffer_size = pnode->GetBufferSize();
			
			for(unsigned int step = 0; step<buffer_size; step++)
			{
				//getting the data of the solution step
				double* step_data = (pnode)->SolutionStepData().Data(step);
				
				double* node0_data = geom[0].SolutionStepData().Data(step);
				double* node1_data = geom[1].SolutionStepData().Data(step);
				double* node2_data = geom[2].SolutionStepData().Data(step);
					
				//copying this data in the position of the vector we are interested in
				for(int j= 0; j< step_data_size; j++)
				{
					step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j];
				}						
			}				
		}

*/









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
		  return "StructuredMeshRefinementModeler";
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
      StructuredMeshRefinementModeler& operator=(StructuredMeshRefinementModeler const& rOther);

      /// Copy constructor.
      StructuredMeshRefinementModeler(StructuredMeshRefinementModeler const& rOther);

        
      ///@}    
        
    }; // Class StructuredMeshRefinementModeler 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    StructuredMeshRefinementModeler& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const StructuredMeshRefinementModeler& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_STRUCTURED_MESH_REFINEMENT_H_INCLUDED  defined 



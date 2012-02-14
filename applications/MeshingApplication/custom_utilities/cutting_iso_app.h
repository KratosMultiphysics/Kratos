//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: jirazabal $
//   Date:                $Date: 2011-12-05 11:51:35 $
//   Revision:            $Revision: 1.5 $
//
//

#if !defined(KRATOS_CUTTING_ISOSURFACE_APPLICATION)
#define  KRATOS_CUTTING_ISOSURFACE_APPLICATION


#ifdef _OPENMP
#include <omp.h> 
#endif

#include "boost/smart_ptr.hpp"
#include <boost/timer.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>


// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>
#include <cmath>
#include <algorithm>


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/dof.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "processes/find_nodal_neighbours_process.h"
#include "processes/find_elements_neighbours_process.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"
#include "utilities/split_tetrahedra.h"
#include "utilities/split_triangle.c"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_3d_3.h"
#include "processes/node_erase_process.h" 
#include "spatial_containers/spatial_containers.h"

namespace Kratos
{

	/// CUTTING ISOSURFACE APPLICATION. 
	/** The cutting app is used to create cutting isosurfaces in a 3D domain. 
	     Despite this can be done in the postprocessing (using Kratos), on large domains the output data would be really large and include info that is not useful for the user.
	     This app creates nodes and conditions (in a new, empty model part)intersecting the domain with as many surfaces as the user want. It is also possible to use different layers.
	     If conditions are triangles, they can be added too to the new model part.
	     After each time step UpdateCutData must be called to save the new info into the new model part, interpolating from the tetraeda elements.
	     After each time step DeleteCutData must also be called to have a clean model part for the next step.
	     NOTE: it only work with tetrahedra elements, generating triangles (1 or 2) from the tetrahedra - plane intersection.
	 */
	class Cutting_Isosurface_Application
	{
	public:
    
    ///@name Type Definitions
    ///@{
	typedef ModelPart::NodesContainerType NodesArrayType;           
	typedef ModelPart::ElementsContainerType ElementsArrayType;
	typedef ModelPart::ConditionsContainerType ConditionsArrayType;
	typedef boost::numeric::ublas::vector<Matrix> Matrix_Order_Tensor;
	typedef boost::numeric::ublas::vector<Vector> Vector_Order_Tensor;
	typedef boost::numeric::ublas::vector<Vector_Order_Tensor> Node_Vector_Order_Tensor;
	typedef Node < 3 > PointType;
	typedef Node < 3 > ::Pointer PointPointerType;
	typedef std::vector<PointType::Pointer> PointVector;
     
    ///@}
    ///@name Life Cycle 
    ///@{ 

    /// Default constructor.
	Cutting_Isosurface_Application(){}  

    /// Destructor.
	virtual ~Cutting_Isosurface_Application(){}
		
		
    ///This function Creates cutting isosurfaces by creating nodes and conditions (to define the conectivities) in a different model part. (new_model_part)
    /** each time it is called a new cutting isosurface is created and therefore new nodes and conditions are added to the new model part
        WARNING: the cutting isosurface MUST cut the domain in at least one triangle, otherwise a segmentiation fault might appear
	* @param mr_model_part . original model part
	* @param mr__new_model_part . destinantion model part
	* @param variable. variable to define the isosurface 
	* @param isovalue. value of the variable to define the isosurface
	* @param isosurface_number . layer to add the conditions (integer)
	* @param tolerance factor .
	**/
	
	template <class TDataType>    
    
	void GenerateVariableCut(ModelPart& mr_model_part, ModelPart& mr_new_model_part, Variable<TDataType>& variable, double isovalue, int isosurface_number, float tolerance)
        {
		KRATOS_WATCH("Generating Isosurface with the following data:")
		KRATOS_WATCH(variable);
		KRATOS_WATCH(isovalue);
		KRATOS_TRY
		                      
		boost::numeric::ublas::vector<array_1d<int, 2 > > Position_Node;
		boost::numeric::ublas::vector<int> List_New_Nodes;
		        
		compressed_matrix<int> Coord;
		ModelPart& this_model_part = mr_model_part;
		ModelPart& new_model_part = mr_new_model_part;
		
		vector<double> Value_Nodes(this_model_part.Elements().size());                 //value of the variable in each node
		vector<int>  Elems_In_Isosurface( this_model_part.Elements().size());          //our (int) vector, where we write 1 when the element is cut by the isosurface, when it is 2, it means we have 2 triangles (4 cutting points)
		int number_of_triangles = 0;
		
		//finished creating the variables and vector/matrix needed. now we have to call the subroutines.
		   
		CSR_Row_Matrix_Mod(this_model_part, Coord);     
		
		FirstLoop(this_model_part, Coord, variable, isovalue, number_of_triangles, Elems_In_Isosurface, tolerance);
		
		Create_List_Of_New_Nodes_Mod(this_model_part,  new_model_part, Coord, List_New_Nodes, Position_Node);
		   
		Calculate_Coordinate_And_Insert_New_Nodes_Mod(this_model_part, new_model_part, Position_Node, List_New_Nodes, variable, isovalue, tolerance);
		    
		GenerateElements (this_model_part, new_model_part, Position_Node, Elems_In_Isosurface, Coord, variable, isosurface_number);
		
		KRATOS_WATCH("Finished generating isosurface")
		KRATOS_CATCH("")
	}
	
	template <class TAdaptatorType>
	
	void GenerateVectorialComponentVariableCut(ModelPart& mr_model_part, ModelPart& mr_new_model_part, VariableComponent < VectorComponentAdaptor < array_1d<double,3> > >& variable, double isovalue, int isosurface_number, float tolerance)
        {
		KRATOS_WATCH("Generating Isosurface with the following data:")
		KRATOS_WATCH(variable);
		KRATOS_WATCH(isovalue);
		KRATOS_TRY
		                      
		boost::numeric::ublas::vector<array_1d<int, 2 > > Position_Node;
		boost::numeric::ublas::vector<int> List_New_Nodes;
		        
		compressed_matrix<int> Coord;
		ModelPart& this_model_part = mr_model_part;
		ModelPart& new_model_part = mr_new_model_part;
		
		vector<double> Value_Nodes(this_model_part.Elements().size());                 //value of the variable in each node
		vector<int>  Elems_In_Isosurface( this_model_part.Elements().size());          //our (int) vector, where we write 1 when the element is cut by the isosurface, when it is 2, it means we have 2 triangles (4 cutting points)
		int number_of_triangles = 0;
		
		//finished creating the variables and vector/matrix needed. now we have to call the subroutines.
		   
		CSR_Row_Matrix_Mod(this_model_part, Coord);     
		
		FirstLoop(this_model_part, Coord, variable, isovalue, number_of_triangles, Elems_In_Isosurface, tolerance);
		
		Create_List_Of_New_Nodes_Mod(this_model_part,  new_model_part, Coord, List_New_Nodes, Position_Node);
		     
		Calculate_Coordinate_And_Insert_New_Nodes_Mod(this_model_part, new_model_part, Position_Node, List_New_Nodes, variable, isovalue, tolerance);
		     
		GenerateElements (this_model_part, new_model_part, Position_Node, Elems_In_Isosurface, Coord, variable, isosurface_number);
		
		KRATOS_WATCH("Finished generating isosurface")
		KRATOS_CATCH("")
	}
	        
        ///************************************************************************************************
        ///************************************************************************************************
	       
        /// ADDMODELPARTELEMENTS: THIS FUNCTION ADDS TO THE NEW MODEL PART THE DATA BELONGING TO THE OLD MODEL PART.
		/** this function adds the elements to the new_model_part. 
		    WARNING: They have to be tetrahedras and it CAN'T be empty, otherwise a segmentation fault will appear
		 * @param mr_model_part . original model part
		 * @param mr_new_model_part . destinantion model part
		 * @param number . layer to add the conditions (integer)
		 **/
	void AddModelPartElements(ModelPart& mr_model_part, ModelPart& mr_new_model_part, int number )
	{
		ModelPart& this_model_part = mr_model_part;
		ModelPart& new_model_part = mr_new_model_part;
		
		KRATOS_WATCH("Adding Skin Conditions to the new model part, added in layer:")
		KRATOS_WATCH(number)
		KRATOS_TRY
		
		vector<int>  Element_Nodes( this_model_part.Nodes().size());          //our (int) vector, where we write -1 when the node is part of the condition faces
		int number_of_tetrahedras = 0;
		int number_of_nodes = 0; //same as above
		int number_of_previous_nodes = 0; // nodes from the previous conditions (planes) created
		
		ElementsArrayType& rElements = this_model_part.Elements();
		
		if (new_model_part.Nodes().size()!=0) //it means this is not the first plane that has to be created
		{
			NodesArrayType& rNodes_new = new_model_part.Nodes();        //i need the model part just to check the id of the new nodes.
			NodesArrayType::iterator it_end_node_new = rNodes_new.ptr_end();
			NodesArrayType::iterator it_begin_node_new = rNodes_new.ptr_begin();
			number_of_previous_nodes=(it_end_node_new-it_begin_node_new);
			KRATOS_WATCH(number_of_previous_nodes)
		}
		else //nothing. number of triangles and nodes already set to 0
		{ 
			KRATOS_WATCH("First Cutting Plane"); 
		}
		
		//we loop all the nodes in the old model part and copy the ones used by conditions to the new model part:
		for (unsigned int index=0 ; index != this_model_part.Nodes().size() ; ++index)
		{
				++number_of_nodes; //one new node!
				Element_Nodes[index] = number_of_nodes + number_of_previous_nodes; //we give this node consecutives ids. now we create the new node
				ModelPart::NodesContainerType::iterator it_node = this_model_part.Nodes().begin()+index;
				Node < 3 > ::Pointer pnode = new_model_part.CreateNewNode(number_of_nodes+number_of_previous_nodes, it_node->X(), it_node->Y(), it_node->Z());  //recordar que es el nueevo model part!!
				pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());
				pnode->GetValue(FATHER_NODES).resize(0);
				pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node.base() ) );       // we keep the same size despite we only need one. to have everyhing with the same size
				pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node.base() ) );
				pnode-> GetValue(WEIGHT_FATHER_NODES) = 1.0;  //lo anterior no anduvo... creo q es así entonces. o hace falta el GetValue?
				
				pnode->X0() = it_node->X0();    
				pnode->Y0() = it_node->Y0();    
				pnode->Z0() = it_node->Z0();
		}//finished the list of nodes to be added. 
		
		//now to the conditions!
		vector<int>  tetrahedra_nodes(4); //here we'll save the nodes' ids with the new node names
		Element const& rReferenceElement = KratosComponents<Element>::Get("Fluid3D");         //element type
		Properties::Pointer properties = this_model_part.GetMesh().pGetProperties(number);        //this will allow us later to turn this layer on/off in GID
		
		for(ModelPart::ElementsContainerType::iterator i_element = rElements.begin() ; i_element != rElements.end() ; i_element++) //looping all the elements
		{
			Geometry<Node<3> >&geom = i_element->GetGeometry(); //current element(nodes, etc)
			for(unsigned int i = 0; i < i_element->GetGeometry().size() ; i++)         //looping the nodes
			{
				int position = (geom[i].Id()) - 1; //the id of the node minus one is the position in the Condition_Node array (position0 = node1)
				tetrahedra_nodes[i]=Element_Nodes[position]; // saving the i nodeId
			} //nodes id saved. now we have to create the element.
			
			Tetrahedra3D4<Node<3> > geometry(
						new_model_part.Nodes()(tetrahedra_nodes[0]),  //condition to be added
						new_model_part.Nodes()(tetrahedra_nodes[1]),
						new_model_part.Nodes()(tetrahedra_nodes[2]),
						new_model_part.Nodes()(tetrahedra_nodes[3]) 
						);

			Element::Pointer p_element = rReferenceElement.Create(number_of_tetrahedras+1, geometry, properties);
			new_model_part.Elements().push_back(p_element);
			++number_of_tetrahedras;
		}
		KRATOS_WATCH(number_of_tetrahedras)
		KRATOS_WATCH("Finished copying elements")
		KRATOS_CATCH("")
	}

	
	///************************************************************************************************
	///************************************************************************************************
	       
	/// ADDSKINCONDITIONS: THIS FUNCTION ADDS TO THE NEW MODEL PART THE DATA OF THE CONDITIONS BELONGING TO THE OLD MODEL PART, BASICALLY THE SAME AS THE PREVIOUS FUNCTION BUT THERE'S NO NEED TO INTERPOLATE SINCE THE NODES COORDINATES ALREADY EXIST. WE ONLY NEED TO COPY THEM TO THE NEW MODEL PART.
		/** this function adds the skin condition. 
		    WARNING: They have to be triangles and it CAN'T be empty, otherwise a segmentation fault will appear
		 * @param mr_model_part . original model part
		 * @param mr_new_model_part . destinantion model part
		 * @param number . layer to add the conditions (integer)
		 **/
	void AddSkinConditions(ModelPart& mr_model_part, ModelPart& mr_new_model_part, int number )
	{
		ModelPart& this_model_part = mr_model_part;
		ModelPart& new_model_part = mr_new_model_part;
		
		KRATOS_WATCH("Adding Skin Conditions to the new model part, added in layer:")
		KRATOS_WATCH(number)
		KRATOS_TRY
		vector<int>  Condition_Nodes( this_model_part.Nodes().size());          //our (int) vector, where we write -1 when the node is part of the condition faces
		int number_of_triangles = 0; //we set it to zero to start
		int number_of_nodes = 0; //same as above
		int number_of_previous_nodes = 0; // nodes from the previous conditions (planes) created

		new_model_part.GetNodalSolutionStepVariablesList() = this_model_part.GetNodalSolutionStepVariablesList();
		new_model_part.AddNodalSolutionStepVariable(FATHER_NODES);
		new_model_part.AddNodalSolutionStepVariable(WEIGHT_FATHER_NODES);
		
		ConditionsArrayType& rConditions = this_model_part.Conditions();
		//ConditionsArrayType::iterator cond_it_begin = rConditions.ptr_begin();
		//ConditionsArrayType::iterator cond_it_end = rConditions.ptr_end();
		            
		ConditionsArrayType& rConditions_new = new_model_part.Conditions();
		ConditionsArrayType::iterator cond_it_end_new = rConditions_new.ptr_end();
		ConditionsArrayType::iterator cond_it_begin_new = rConditions_new.ptr_begin();
		
		if (new_model_part.Nodes().size()!=0) //it means this is not the first plane that has to be created
		{
			NodesArrayType& rNodes_new = new_model_part.Nodes();        //i need the model part just to check the id of the new nodes.
			NodesArrayType::iterator it_end_node_new = rNodes_new.ptr_end();
			NodesArrayType::iterator it_begin_node_new = rNodes_new.ptr_begin();
			number_of_previous_nodes=(it_end_node_new-it_begin_node_new);
			number_of_triangles=cond_it_end_new - cond_it_begin_new ;
			KRATOS_WATCH(number_of_triangles)
			KRATOS_WATCH(number_of_previous_nodes)
			//KRATOS_CATCH("")
		}
		else //nothing. number of triangles and nodes already set to 0
		{ 
			KRATOS_WATCH("First Cutting Plane"); 
		}
		
		for (unsigned int index=0 ; index != this_model_part.Nodes().size() ; ++index) Condition_Nodes[index]=0; //initializing in zero the whole vector (meaning no useful nodes for the condition layer)
		
		for(ModelPart::ConditionsContainerType::iterator i_condition = rConditions.begin() ; i_condition != rConditions.end() ; i_condition++)
		{
			Geometry<Node<3> >&geom = i_condition->GetGeometry();
			for(unsigned int i = 0; i < i_condition->GetGeometry().size() ; i++)         
			{
				int position = (geom[i].Id()) - 1; //the id of the node minus one is the position in the Condition_Node array (position0 = node1)
				Condition_Nodes[position] =  -1 ; //a -1 means we need this node!
			}
		}//done. now we know all the nodes that will have to be added to the new model part.
		
		//we loop all the nodes in the old model part and copy the ones used by conditions to the new model part:
		for (unsigned int index=0 ; index != this_model_part.Nodes().size() ; ++index)
		{
			if (Condition_Nodes[index]==-1)
			{
				++number_of_nodes; //one new node!
				Condition_Nodes[index] = number_of_nodes + number_of_previous_nodes; //we give this node consecutives ids. now we create the new node
				ModelPart::NodesContainerType::iterator it_node = this_model_part.Nodes().begin()+index;
				Node < 3 > ::Pointer pnode = new_model_part.CreateNewNode(number_of_nodes+number_of_previous_nodes, it_node->X(), it_node->Y(), it_node->Z());  //recordar que es el nueevo model part!!
				pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());
				pnode->GetValue(FATHER_NODES).resize(0);
				pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node.base() ) );       // we keep the same size despite we only need one. to have everyhing with the same size
				pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node.base() ) );
				pnode-> GetValue(WEIGHT_FATHER_NODES) = 1.0;  //lo anterior no anduvo... creo q es así entonces. o hace falta el GetValue?
				
				pnode->X0() = it_node->X0();    
				pnode->Y0() = it_node->Y0();    
				pnode->Z0() = it_node->Z0(); 
			}
		}//finished the list of nodes to be added. 
		 
		 //now to the conditions!
		vector<int>  triangle_nodes(3); //here we'll save the nodes' ids with the new node names
		Condition const& rReferenceCondition = KratosComponents<Condition>::Get("Condition3D");         //condition type
		Properties::Pointer properties = this_model_part.GetMesh().pGetProperties(number); 		//this will allow us later to turn this layer on/off in GID
		 
		for(ModelPart::ConditionsContainerType::iterator i_condition = rConditions.begin() ; i_condition != rConditions.end() ; i_condition++) //looping all the conditions
		{
			Geometry<Node<3> >&geom = i_condition->GetGeometry(); //current condition(nodes, etc)
			for(unsigned int i = 0; i < i_condition->GetGeometry().size() ; i++)         //looping the nodes
			{
				int position = (geom[i].Id()) - 1; //the id of the node minus one is the position in the Condition_Node array (position0 = node1)
				triangle_nodes[i]=Condition_Nodes[position]; // saving the i nodeId
			} //nodes id saved. now we have to create the element.
			Triangle3D3<Node<3> > geometry(
						new_model_part.Nodes()(triangle_nodes[0]),  //condition to be added
						new_model_part.Nodes()(triangle_nodes[1]),
						new_model_part.Nodes()(triangle_nodes[2])	 
						);
			Condition::Pointer p_condition = rReferenceCondition.Create(number_of_triangles+1, geometry, properties); //está bien? acá la verdad ni idea. sobre todo number_of_triangles (la posición). o debe ser un puntero de elemento en vez de un entero?
			new_model_part.Conditions().push_back(p_condition);
			++number_of_triangles;
		}
		KRATOS_WATCH("Finished copying conditions surfaces")
		KRATOS_CATCH("")
	}
	
        ///************************************************************************************************
        ///************************************************************************************************
	
    ///LIST OF SUBROUTINES

	void CSR_Row_Matrix_Mod(ModelPart& this_model_part, compressed_matrix<int>& Coord)
	{
		NodesArrayType& pNodes = this_model_part.Nodes();
		Coord.resize(pNodes.size(), pNodes.size());
		NodesArrayType::iterator i_begin = pNodes.ptr_begin();
		NodesArrayType::iterator i_end = pNodes.ptr_end();
		
		std::vector<unsigned int> aux(10000);
		
		for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i)
		{
			int index_i = i->Id() - 1;
			WeakPointerVector< Node < 3 > >& neighb_nodes = i->GetValue(NEIGHBOUR_NODES);
			Coord.push_back(index_i, index_i, -1);        //only modification added, now the diagonal is filled with -1 too.
			                
			unsigned int active = 0;
			for (WeakPointerVector< Node < 3 > >::iterator inode = neighb_nodes.begin(); inode != neighb_nodes.end(); inode++)
			{
				int index_j = inode->Id() - 1;
				if (index_j > index_i)
				{
					aux[active] = index_j;
					active++;
				}
			}
			std::sort(aux.begin(), aux.begin() + active);
			
			for (unsigned int k = 0; k < active; k++)
			{
				Coord.push_back(index_i, aux[k], -1);
			}
		}
	}
	
    ///************************************************************************************************      
	        
	void FirstLoop(ModelPart& this_model_part, compressed_matrix<int>& Coord, Variable<double>& variable, double isovalue, int number_of_triangles, vector<int>&  Elems_In_Isosurface, float tolerance)//
	{                                           //varible is the variable that defines the isosurface
						    //isovalue is the value of the variable that defines the isosurface

		ElementsArrayType& rElements = this_model_part.Elements();

		ElementsArrayType::iterator it_begin = rElements.ptr_begin();
		ElementsArrayType::iterator it_end = rElements.ptr_end(); 
            
            
		array_1d<double, 4 > node_value;                     // value of the variable in the node
		array_1d<double, 4 > neigh_value;                    // value of the variable in thr other node of the edge
		array_1d<double, 4 > diff_node_value;                // difference between the imposed value of the variable and its value in the node
		array_1d<double, 4 > diff_neigh_value;               // difference beween the variable value on other node of the edge (neighbour) and the imposed value                      
		array_1d<unsigned int, 4 > list_matching_nodes;      //used to save the new nodes that match exactly old nodes (very unlikely, but might be 4 for very plane elements)
		unsigned int exact_nodes=0;
            
		number_of_triangles = 0;
		int current_element= 0; //current element. it's a position. NOT ID!
            
		int number_of_cuts= 0; //this is the counter explained in the following lines
            
		for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) //looping all the elements
		{
			++current_element;
			number_of_cuts = 0;
			exact_nodes = 0 ;       
			Geometry<Node<3> >&geom = it->GetGeometry(); //geometry of the element
			for(unsigned int i = 0; i < it->GetGeometry().size() ; i++)          //size = 4 ; nodes per element. NOTICE WE'LL BE LOOPING THE EDGES TWICE. THIS IS A WASTE OF TIME BUT MAKES IT EASIER TO IDENTITY ELEMENTS. LOOK BELOW.
			//when we have a triangle inside a thetraedra, its edges (or nodes) must be cut 3 times by the plane. if we loop all 2 times we can have a counter. when it's = 6 then we have a triangle. when tetrahedras are cutted 8 times then we have 2 triangles (or a cuatrilateral, the same)
			{
				node_value[i] = geom[i].FastGetSolutionStepValue(variable);
				
				diff_node_value[i] = isovalue - node_value[i];
				
				for(unsigned int j = 0; j < it->GetGeometry().size() ; j++)          //  looping on the neighbours
				{
					if (i != j)  //(cant link node with itself)
					{
						neigh_value[j] = geom[j].FastGetSolutionStepValue(variable);
						
						diff_neigh_value[j] = isovalue - neigh_value[j];
						
						bool isovernode=false;   // if true, then it can't be between the nodes
												
						if (fabs(diff_node_value[i]) < (tolerance)) //then our node is part of the plane (this should have been done before the loop on neighbours, but this way it is easier to read .
						{
							int index_i = geom[i].Id() -1 ; // i node id
							Coord(index_i, index_i) = -2;	              //saving a -2 in the diagonal
							isovernode=true;
							number_of_cuts += 2; //since its neighbour wont take this case as a cut, we must save 2 cuts instead of one. (to reach number_of_cuts=6), 
							++exact_nodes;
							list_matching_nodes[i]= geom[i].Id();
							break;
						}
						//check this last condition, used to avoid talking points that belong to other node. might cause some problems when the plane is almost paralel to an edge. to be improved! (it seems to be working correcly even when the edge is part of the plane.)
						if ((diff_node_value[i]*diff_neigh_value[j]) < 0.0 && isovernode==false && fabs(diff_neigh_value[j]) > tolerance) // this means one has a bigger value than the imposed values and the other smaller, no need to do more checks, between there is a point that should match the value!
						{							
							int index_i = geom[i].Id() - 1;     //i node id
							int index_j = geom[j].Id() - 1;     //j node id
							if (index_j > index_i)  Coord(index_i, index_j) = -2;    //saving a -2 in the upper side of the matrix
							else Coord(index_j, index_i) = -2;
							number_of_cuts += 1;
						}
					} //closing the i!=j if
				} //closing the neighbour loop
			} //closing the nodes loop
			if (tolerance > (diff_node_value[0]) && tolerance > (diff_node_value[1]) && tolerance > (diff_node_value[2]) && tolerance > (diff_node_value[3])) //if all the nodes have the same the same value as the isovalue we do not take into account the triangles
			{
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++)
				{
					number_of_cuts -= 2;
					--exact_nodes;
				}  
			} 
			//now we have to save the data. we should get a list with the elements that will genereate triangles and the total number of triangles
			Elems_In_Isosurface[current_element-1] = 0 ; //we initialize as 0
			
			if (number_of_cuts == 6)     //it can be 8, in that case we have 2 triangles (the cut generates a square)
			{
				number_of_triangles +=1;
				Elems_In_Isosurface[current_element-1] = 1 ; //i still don't know the number of the node so i'll have to do another loop later to assign to define node id's of each triangular element  
			}
			else if (number_of_cuts == 8 ) // 2 triangles in the element!
			{
				number_of_triangles +=2;
				Elems_In_Isosurface[ current_element-1] = 2 ;
			}
				
		} //closing the elem loop
		KRATOS_WATCH(number_of_triangles);
	} //closing "FirstLoop"
	        
	void FirstLoop(ModelPart& this_model_part, compressed_matrix<int>& Coord, VariableComponent < VectorComponentAdaptor < array_1d<double,3> > >& variable, double isovalue, int number_of_triangles, vector<int>&  Elems_In_Isosurface, float tolerance)//
	{                                           //varible is the variable that defines the isosurface
						    //isovalue is the value of the variable that defines the isosurface

		ElementsArrayType& rElements = this_model_part.Elements();

		ElementsArrayType::iterator it_begin = rElements.ptr_begin();
		ElementsArrayType::iterator it_end = rElements.ptr_end(); 
            
            
		array_1d<double, 4 > node_value;                     // value of the variable in the node
		array_1d<double, 4 > neigh_value;                    // value of the variable in thr other node of the edge
		array_1d<double, 4 > diff_node_value;                // difference between the imposed value of the variable and its value in the node
		array_1d<double, 4 > diff_neigh_value;               // difference beween the variable value on other node of the edge (neighbour) and the imposed value                      
		array_1d<unsigned int, 4 > list_matching_nodes;      //used to save the new nodes that match exactly old nodes (very unlikely, but might be 4 for very plane elements)
		unsigned int exact_nodes=0;
            
		number_of_triangles = 0;
		int current_element= 0; //current element. it's a position. NOT ID!
            
		int number_of_cuts= 0; //this is the counter explained in the following lines
            
		for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) //looping all the elements
		{
			++current_element;
			number_of_cuts = 0;
			
			exact_nodes = 0 ;       
			Geometry<Node<3> >&geom = it->GetGeometry(); //geometry of the element
			for(unsigned int i = 0; i < it->GetGeometry().size() ; i++)          //size = 4 ; nodes per element. NOTICE WE'LL BE LOOPING THE EDGES TWICE. THIS IS A WASTE OF TIME BUT MAKES IT EASIER TO IDENTITY ELEMENTS. LOOK BELOW.
			//when we have a triangle inside a thetraedra, its edges (or nodes) must be cut 3 times by the plane. if we loop all 2 times we can have a counter. when it's = 6 then we have a triangle. when tetrahedras are cutted 8 times then we have 2 triangles (or a cuatrilateral, the same)
			{
				node_value[i] = geom[i].FastGetSolutionStepValue(variable);
				
				diff_node_value[i] = isovalue - node_value[i];
				
				for(unsigned int j = 0; j < it->GetGeometry().size() ; j++)          //  looping on the neighbours
				{
					if (i != j)  //(cant link node with itself)
					{
						neigh_value[j] = geom[j].FastGetSolutionStepValue(variable);
						
						diff_neigh_value[j] = isovalue - neigh_value[j];
						
						bool isovernode=false;   // if true, then it can't be between the nodes
												
						if (fabs(diff_node_value[i]) < (tolerance)) //then our node is part of the plane (this should have been done before the loop on neighbours, but this way it is easier to read .
						{
							int index_i = geom[i].Id() -1 ; // i node id
							Coord(index_i, index_i) = -2;	              //saving a -2 in the diagonal
							isovernode=true;
							number_of_cuts += 2; //since its neighbour wont take this case as a cut, we must save 2 cuts instead of one. (to reach number_of_cuts=6), 
							++exact_nodes;
							list_matching_nodes[i]= geom[i].Id();
							break;
						}
						//check this last condition, used to avoid talking points that belong to other node. might cause some problems when the plane is almost paralel to an edge. to be improved! (it seems to be working correcly even when the edge is part of the plane.)
						if ((diff_node_value[i]*diff_neigh_value[j]) < 0.0 && isovernode==false && fabs(diff_neigh_value[j]) > tolerance) // this means one has a bigger value than the imposed values and the other smaller, no need to do more checks, between there is a point that should match the value!
						{							
							int index_i = geom[i].Id() - 1;     //i node id
							int index_j = geom[j].Id() - 1;     //j node id
							if (index_j > index_i)  Coord(index_i, index_j) = -2;    //saving a -2 in the upper side of the matrix
							else Coord(index_j, index_i) = -2;
							number_of_cuts += 1;
						}
					} //closing the i!=j if
				} //closing the neighbour loop
			} //closing the nodes loop
			if (tolerance > (diff_node_value[0]) && tolerance > (diff_node_value[1]) && tolerance > (diff_node_value[2]) && tolerance > (diff_node_value[3])) //if all the nodes have the same the same value as the isovalue we do not take into account the triangles
			{
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++)
				{
					number_of_cuts -= 2;
					--exact_nodes;
				}  
			} 
			//now we have to save the data. we should get a list with the elements that will genereate triangles and the total number of triangles
			Elems_In_Isosurface[current_element-1] = 0 ; //we initialize as 0
			
			if (number_of_cuts == 6)     //it can be 8, in that case we have 2 triangles (the cut generates a square)
			{
				number_of_triangles +=1;
				Elems_In_Isosurface[current_element-1] = 1 ; //i still don't know the number of the node so i'll have to do another loop later to assign to define node id's of each triangular element  
			}
			else if (number_of_cuts == 8 ) // 2 triangles in the element!
			{
				number_of_triangles +=2;
				Elems_In_Isosurface[ current_element-1] = 2 ;
			}
				
		} //closing the elem loop
		KRATOS_WATCH(number_of_triangles);
	} //closing "FirstLoop"
	
	void FirstLoop(ModelPart& this_model_part, compressed_matrix<int>& Coord, Variable < array_1d<double,3> >& variable, double isovalue, int number_of_triangles, vector<int>&  Elems_In_Isosurface, float tolerance)//
	{                                           //varible is the variable that defines the isosurface
						    //isovalue is the value of the variable that defines the isosurface

		ElementsArrayType& rElements = this_model_part.Elements();

		ElementsArrayType::iterator it_begin = rElements.ptr_begin();
		ElementsArrayType::iterator it_end = rElements.ptr_end(); 
            
		array_1d< array_1d<double,3>, 4 > vector_node_value;  // value of the variable in the node
		array_1d< array_1d<double,3>, 4 > vector_neigh_value; // value of the variable in the other node of the edge
		array_1d<long double, 4 > squared_node_value;
		array_1d<long double, 4 > squared_neigh_value;
		double squared_isovalue = isovalue * isovalue;
		double aux_value;
		double a;
		double b;
		double c;
		double weight1;
		double weight2;
		array_1d<double, 4 > node_value;                      // module of the variable in the node
		array_1d<double, 4 > neigh_value;                     // module of the variable in the other node of the edge
		array_1d<double, 4 > diff_node_value;                 // difference between the imposed value of the variable and its value in the node
		array_1d<double, 4 > diff_neigh_value;                // difference beween the variable value on other node of the edge (neighbour) and the imposed value                     
		array_1d<unsigned int, 4 > list_matching_nodes;      //used to save the new nodes that match exactly old nodes (very unlikely, but might be 4 for very plane elements)
		unsigned int exact_nodes=0;
            
		number_of_triangles = 0;
		int current_element= 0; //current element. it's a position. NOT ID!
            
		int number_of_cuts= 0; //this is the counter explained in the following lines
            
		for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) //looping all the elements
		{
			++current_element;
			number_of_cuts = 0;
			
			exact_nodes = 0 ;       
			Geometry<Node<3> >&geom = it->GetGeometry(); //geometry of the element
			for(unsigned int i = 0; i < it->GetGeometry().size() ; i++)          //size = 4 ; nodes per element. NOTICE WE'LL BE LOOPING THE EDGES TWICE. THIS IS A WASTE OF TIME BUT MAKES IT EASIER TO IDENTITY ELEMENTS. LOOK BELOW.
			//when we have a triangle inside a thetraedra, its edges (or nodes) must be cut 3 times by the plane. if we loop all 2 times we can have a counter. when it's = 6 then we have a triangle. when tetrahedras are cutted 8 times then we have 2 triangles (or a cuatrilateral, the same)
			{
				vector_node_value[i] = geom[i].FastGetSolutionStepValue(variable);
				
				squared_node_value[i] = (vector_node_value[i][0] * vector_node_value[i][0]) + (vector_node_value[i][1] * vector_node_value[i][1]) + (vector_node_value[i][2] * vector_node_value[i][2]);
				node_value[i] = sqrt(squared_node_value[i]);
				
				diff_node_value[i] = isovalue - node_value[i];
				
				for(unsigned int j = 0; j < it->GetGeometry().size() ; j++)          //  looping on the neighbours
				{
					if (i != j)  //(cant link node with itself)
					{
						vector_neigh_value[j] = geom[j].FastGetSolutionStepValue(variable);
						
						squared_neigh_value[j] = (vector_neigh_value[j][0] * vector_neigh_value[j][0]) + (vector_neigh_value[j][1] * vector_neigh_value[j][1]) + (vector_neigh_value[j][2] * vector_neigh_value[j][2]);
						neigh_value[j] = sqrt(squared_neigh_value[j]);
						
						diff_neigh_value[j] = isovalue - neigh_value[j];
						
						aux_value = (vector_node_value[j][0] * vector_neigh_value[j][0]) + (vector_node_value[j][1] * vector_neigh_value[j][1]) + (vector_node_value[j][2] * vector_neigh_value[j][2]);
						
						a = squared_node_value[i] + squared_neigh_value[j] - (2 * aux_value);
						b = 2 * (aux_value - squared_node_value[i]);
						c = squared_node_value[i] - squared_isovalue;
						
						float tol = 0.0000000000000000000000000000000000001;
						
						if (fabs(a) < tol)
						{  
							if (a < 0) a = tol;
							if (a < 0) a = -tol;
						}
						if (fabs(b) < tol) b = tol;
						{  
							if (b < 0) b = tol;
							if (b < 0) b = -tol;
						}
						if (fabs(c) < tol) c = tol;
						{  
							if (c < 0) c = tol;
							if (c < 0) c = -tol;
						}
						
						weight1 = (-b + sqrt((b * b) - (4 * a * c))) / (2 * a);
						weight2 = (-b - sqrt((b * b) - (4 * a * c))) / (2 * a);
						
						bool isovernode=false;   // if true, then it can't be between the nodes
							
						if (fabs(diff_node_value[i]) <= (tolerance)) //then our node is part of the plane (this should have been done before the loop on neighbours, but this way it is easier to read .
						{
							int index_i = geom[i].Id() -1 ; // i node id
							Coord(index_i, index_i) = -2;	              //saving a -2 in the diagonal
							isovernode=true;
							number_of_cuts += 2; //since its neighbour wont take this case as a cut, we must save 2 cuts instead of one. (to reach number_of_cuts=6), 
							++exact_nodes;
							list_matching_nodes[i]= geom[i].Id();
							break;
						}
						//check this last condition, used to avoid talking points that belong to other node. might cause some problems when the plane is almost paralel to an edge. to be improved! (it seems to be working correcly even when the edge is part of the plane.)
						if (weight1 > 0.0 && weight1 < 1.0 && fabs(diff_neigh_value[j]) > tolerance  && isovernode==false) // this means one has a bigger value than the imposed values and the other smaller, no need to do more checks, between there is a point that should match the value!
						{
							int index_i = geom[i].Id() - 1;     //i node id
							int index_j = geom[j].Id() - 1;     //j node id
							if (index_j > index_i)  Coord(index_i, index_j) = -2;    //saving a -2 in the upper side of the matrix
							else Coord(index_j, index_i) = -2;
							isovernode=true;
							number_of_cuts += 1;
							break;
						}
						if (weight2 > 0.0 && weight2 < 1.0 && fabs(diff_neigh_value[j]) > tolerance  && isovernode==false)
						{
							int index_i = geom[i].Id() - 1;     //i node id
							int index_j = geom[j].Id() - 1;     //j node id
							if (index_j > index_i)  Coord(index_i, index_j) = -2;    //saving a -2 in the upper side of the matrix
							else Coord(index_j, index_i) = -2;
							number_of_cuts += 1;
						}  
					} //closing the i!=j if
				} //closing the neighbour loop
			} //closing the nodes loop
			if (tolerance > (diff_node_value[0]) && tolerance > (diff_node_value[1]) && tolerance > (diff_node_value[2]) && tolerance > (diff_node_value[3])) //if all the nodes have the same the same value as the isovalue we do not take into account the triangles
			{
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++)
				{
					number_of_cuts -= 2;
					--exact_nodes;
				}  
			} 
			//now we have to save the data. we should get a list with the elements that will genereate triangles and the total number of triangles
			Elems_In_Isosurface[current_element-1] = 0 ; //we initialize as 0
			
			if (number_of_cuts == 6)     //it can be 8, in that case we have 2 triangles (the cut generates a square)
			{
				number_of_triangles +=1;
				Elems_In_Isosurface[current_element-1] = 1 ; //i still don't know the number of the node so i'll have to do another loop later to assign to define node id's of each triangular element   
			}
			else if (number_of_cuts == 8 ) // 2 triangles in the element!
			{
				number_of_triangles +=2;
				Elems_In_Isosurface[ current_element-1] = 2 ;
			}
		} //closing the elem loop
		KRATOS_WATCH(number_of_triangles);
	} //closing "FirstLoop"
	
    ///************************************************************************************************
	
	void Create_List_Of_New_Nodes_Mod(ModelPart& this_model_part, ModelPart& new_model_part, compressed_matrix<int>& Coord, boost::numeric::ublas::vector<int> &List_New_Nodes, boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node) //surface number = 1 -- inf (but the first one should be always one!
	{
		unsigned int number_of_new_nodes = 0;
		//NodesArrayType& pNodes = this_model_part.Nodes();
		typedef compressed_matrix<int>::iterator1 i1_t;
		typedef compressed_matrix<int>::iterator2 i2_t;
	
	    ///*WARNING
		for (i1_t i1 = Coord.begin1(); i1 != Coord.end1(); ++i1)
		{
			for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
			{
				if (Coord(i2.index1(), i2.index2()) == -2)
				{
					number_of_new_nodes++;                          //this should work without any change
				}
			}
		}
		if (number_of_new_nodes != 0)
		{  
		    ///* New Id de los Nodos
			List_New_Nodes.resize(number_of_new_nodes);
			//int total_node = pNodes.size();
			
			if (new_model_part.Nodes().size()!=0) //it means this is not the first surface that has to be created
			{
				NodesArrayType& rNodes_new = new_model_part.Nodes();        //i need the model part just to check the id of the new nodes.
				NodesArrayType::iterator it_end_node_new = rNodes_new.ptr_end();
				NodesArrayType::iterator it_begin_node_new = rNodes_new.ptr_begin();
				List_New_Nodes[0]=(it_end_node_new-it_begin_node_new)+1;
				KRATOS_WATCH("New Cutting Surface")
				first_cutting_surface = false;
			}
			else
			{
				List_New_Nodes[0]=1; 
				KRATOS_WATCH("First Cutting Surface"); 
				first_cutting_surface = true;
			}
			
			for (unsigned int i = 1; i < number_of_new_nodes; i++)
			{
				List_New_Nodes[i] = List_New_Nodes[0] + i;      //just a list. necessary because other cutting surfaces might have been generated before	      
			}                                  
			
		    ///* setting edges -2 to the new id of the new node
			Position_Node.resize(number_of_new_nodes);
			unsigned int index = 0;
			for (i1_t i1 = Coord.begin1(); i1 != Coord.end1(); ++i1)
			{
				for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
				{
					if (Coord(i2.index1(), i2.index2()) == -2)
					{
						Coord(i2.index1(), i2.index2()) = List_New_Nodes[index];
						Position_Node[index][0] = i2.index1() + 1;
						Position_Node[index][1] = i2.index2() + 1;                     //in the diagonal terms both indexes should point to the same node but still should work i guess
						index++;
					}
				}
			}
		}
	}

    ///************************************************************************************************ 
	        
	void Calculate_Coordinate_And_Insert_New_Nodes_Mod(ModelPart& this_model_part, ModelPart& new_model_part, const boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node, const boost::numeric::ublas::vector<int> &List_New_Nodes, Variable<double>& variable, double isovalue, float tolerance)//,
	{
		array_1d<double, 3 > Coord_Node_1;
		array_1d<double, 3 > Coord_Node_2;
		array_1d<double, 3 > temp_dist;
		array_1d<double, 3 > vector_distance;
		double node_value;
		double neigh_value;
		double diff_node_value;
		double diff_neigh_value;
		double diff_neigh_node_value;
		double weight;
		boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node;
		Coordinate_New_Node.resize(Position_Node.size());
		
		new_model_part.GetNodalSolutionStepVariablesList() = this_model_part.GetNodalSolutionStepVariablesList();
		new_model_part.AddNodalSolutionStepVariable(FATHER_NODES);
		new_model_part.AddNodalSolutionStepVariable(WEIGHT_FATHER_NODES);
		
		for (unsigned int i = 0; i < Position_Node.size(); i++) //looping the new nodes
		{
			/// calculating the coordinate of the new nodes
			const int& node_i = Position_Node[i][0];
			const int& node_j = Position_Node[i][1];
			ModelPart::NodesContainerType::iterator it_node1 = this_model_part.Nodes().find(node_i);
			noalias(Coord_Node_1) = it_node1->Coordinates();
			ModelPart::NodesContainerType::iterator it_node2 = this_model_part.Nodes().find(node_j);
			noalias(Coord_Node_2) = it_node2->Coordinates();
			//ok, now we have both coordinates. now we must define a weight coefficient based on the distance.
			//this coeff will be =node_plane_distance/node_neigh_distance (linear interpolation)
			node_value = it_node1->FastGetSolutionStepValue(variable);
			neigh_value = it_node2->FastGetSolutionStepValue(variable);
			
			diff_node_value = isovalue - node_value;
			diff_neigh_value = isovalue - neigh_value;
			diff_neigh_node_value = neigh_value - node_value;
			
			vector_distance = Coord_Node_2 - Coord_Node_1;
			
			if (fabs(diff_node_value) <  tolerance) weight = 0.0;
			else if (fabs(diff_neigh_value) <  tolerance) weight = 1.0;
			else weight = diff_node_value / diff_neigh_node_value;
			
			if (weight > 1.00005) KRATOS_WATCH("**** something's wrong! weight higher than 1! ****");
			
			for (unsigned int index=0; index!=3; ++index) //we loop the 3 coordinates)
			{
				if (Position_Node[i][0]!=Position_Node[i][1])  Coordinate_New_Node[i][index] = Coord_Node_1[index] + vector_distance[index] * weight; 
				else   Coordinate_New_Node[i][index] = Coord_Node_1[index]; //when both nodes are the same it doesnt make any sense to interpolate
				
				/// inserting the new node in the model part
				Node < 3 > ::Pointer pnode = new_model_part.CreateNewNode(List_New_Nodes[i], Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2]);  //remember that is the new model part!!
			
				pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());
				
				pnode->GetValue(FATHER_NODES).resize(0);
				pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node1.base() ) );       //saving data about fathers in the model part
				pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node2.base() ) );
				pnode-> GetValue(WEIGHT_FATHER_NODES) = weight;  
				
				pnode->X0() = (1.0 - weight) * (it_node1->X0())  +  weight * it_node2->X0();   
				pnode->Y0() = (1.0 - weight) * (it_node1->Y0())  +  weight * it_node2->Y0();    
				pnode->Z0() = (1.0 - weight) * (it_node1->Z0())  +  weight * it_node2->Z0();
			}
		}
	} 
	        
	void Calculate_Coordinate_And_Insert_New_Nodes_Mod(ModelPart& this_model_part, ModelPart& new_model_part, const boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node, const boost::numeric::ublas::vector<int> &List_New_Nodes, VariableComponent < VectorComponentAdaptor < array_1d<double,3> > >& variable, double isovalue, float tolerance)//,
	{
		array_1d<double, 3 > Coord_Node_1;
		array_1d<double, 3 > Coord_Node_2;
		array_1d<double, 3 > temp_dist;
		array_1d<double, 3 > vector_distance;
		double node_value;
		double neigh_value;
		double diff_node_value;
		double diff_neigh_value;
		double diff_neigh_node_value;
		double weight;
		boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node;
		Coordinate_New_Node.resize(Position_Node.size());
		
		new_model_part.GetNodalSolutionStepVariablesList() = this_model_part.GetNodalSolutionStepVariablesList();
		new_model_part.AddNodalSolutionStepVariable(FATHER_NODES);
		new_model_part.AddNodalSolutionStepVariable(WEIGHT_FATHER_NODES);
		
		for (unsigned int i = 0; i < Position_Node.size(); i++) //looping the new nodes
		{
			/// calculating the coordinate of the new nodes
			const int& node_i = Position_Node[i][0];
			const int& node_j = Position_Node[i][1];
			ModelPart::NodesContainerType::iterator it_node1 = this_model_part.Nodes().find(node_i);
			noalias(Coord_Node_1) = it_node1->Coordinates();
			ModelPart::NodesContainerType::iterator it_node2 = this_model_part.Nodes().find(node_j);
			noalias(Coord_Node_2) = it_node2->Coordinates();
			//ok, now we have both coordinates. now we must define a weight coefficient based on the distance.
			//this coeff will be =node_plane_distance/node_neigh_distance (linear interpolation)
			node_value = it_node1->FastGetSolutionStepValue(variable);
			neigh_value = it_node2->FastGetSolutionStepValue(variable);
			
			diff_node_value = isovalue - node_value;
			diff_neigh_value = isovalue - neigh_value;
			diff_neigh_node_value = neigh_value - node_value;
			
			vector_distance = Coord_Node_2 - Coord_Node_1;
			
			if (fabs(diff_node_value) <  tolerance) weight = 0.0;
			else if (fabs(diff_neigh_value) <  tolerance) weight = 1.0;
			else weight = diff_node_value / diff_neigh_node_value;
			
			if (weight > 1.00005) KRATOS_WATCH("**** something's wrong! weight higher than 1! ****");
			
			for (unsigned int index=0; index!=3; ++index) //we loop the 3 coordinates)
			{
				if (Position_Node[i][0]!=Position_Node[i][1])  Coordinate_New_Node[i][index] = Coord_Node_1[index] + vector_distance[index] * weight; 
				else   Coordinate_New_Node[i][index] = Coord_Node_1[index]; //when both nodes are the same it doesnt make any sense to interpolate
				
				/// inserting the new node in the model part
				Node < 3 > ::Pointer pnode = new_model_part.CreateNewNode(List_New_Nodes[i], Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2]);  //remember that is the new model part!!
			
				pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());
				
				pnode->GetValue(FATHER_NODES).resize(0);
				pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node1.base() ) );       //saving data about fathers in the model part
				pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node2.base() ) );
				pnode-> GetValue(WEIGHT_FATHER_NODES) = weight;  
				
				pnode->X0() = (1.0 - weight) * (it_node1->X0())  +  weight * it_node2->X0();   
				pnode->Y0() = (1.0 - weight) * (it_node1->Y0())  +  weight * it_node2->Y0();    
				pnode->Z0() = (1.0 - weight) * (it_node1->Z0())  +  weight * it_node2->Z0();
			}
		}
	}
	
	void Calculate_Coordinate_And_Insert_New_Nodes_Mod(ModelPart& this_model_part, ModelPart& new_model_part, const boost::numeric::ublas::vector<array_1d<int, 2 > >& Position_Node, const boost::numeric::ublas::vector<int> &List_New_Nodes, Variable< array_1d< double, 3> >& variable, double isovalue, float tolerance)//,
	{
		array_1d<double, 3 > Coord_Node_1;
		array_1d<double, 3 > Coord_Node_2;
		array_1d<double, 3 > temp_dist;
		array_1d<double, 3 > vector_distance;
		array_1d<double, 3 > vector_node_value;
		array_1d<double, 3 > vector_neigh_value;
		double squared_node_value;
		double squared_neigh_value;
		double node_value;
		double neigh_value;
		double aux_value;
		double squared_isovalue;
		double a;
		double b;
		double c;
		double diff_node_value;
		double diff_neigh_value;
		double weight;
		double weight1;
		double weight2;
		boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node;
		Coordinate_New_Node.resize(Position_Node.size());
		
		new_model_part.GetNodalSolutionStepVariablesList() = this_model_part.GetNodalSolutionStepVariablesList();
		new_model_part.AddNodalSolutionStepVariable(FATHER_NODES);
		new_model_part.AddNodalSolutionStepVariable(WEIGHT_FATHER_NODES);
		
		for (unsigned int i = 0; i < Position_Node.size(); i++) //looping the new nodes
		{
			/// calculating the coordinate of the new nodes
			const int& node_i = Position_Node[i][0];
			const int& node_j = Position_Node[i][1];
			ModelPart::NodesContainerType::iterator it_node1 = this_model_part.Nodes().find(node_i);
			noalias(Coord_Node_1) = it_node1->Coordinates();
			ModelPart::NodesContainerType::iterator it_node2 = this_model_part.Nodes().find(node_j);
			noalias(Coord_Node_2) = it_node2->Coordinates();
			//ok, now we have both coordinates. now we must define a weight coefficient based on the distance.
			//this coeff will be =node_plane_distance/node_neigh_distance (linear interpolation)
			
			vector_node_value = it_node1->FastGetSolutionStepValue(variable);
			squared_node_value = (vector_node_value[0] * vector_node_value[0]) + (vector_node_value[1] * vector_node_value[1]) + (vector_node_value[2] * vector_node_value[2]);
			node_value = sqrt(squared_node_value);
			
			vector_neigh_value = it_node2->FastGetSolutionStepValue(variable);
			squared_neigh_value = (vector_neigh_value[0] * vector_neigh_value[0]) + (vector_neigh_value[1] * vector_neigh_value[1]) + (vector_neigh_value[2] * vector_neigh_value[2]);
			neigh_value = sqrt(squared_neigh_value);
			
			aux_value = (vector_node_value[0] * vector_neigh_value[0]) + (vector_node_value[1] * vector_neigh_value[1]) + (vector_node_value[2] * vector_neigh_value[2]);
			
			squared_isovalue = isovalue * isovalue;
			
			a = squared_node_value + squared_neigh_value - (2 * aux_value);
			b = 2 * (aux_value - squared_node_value);
			c = squared_node_value - squared_isovalue;
			
			float tol = 0.0000000000000000000000000000000000001;
			
			if (fabs(a) < tol)
			{  
				if (a < 0) a = tol;
				if (a < 0) a = -tol;
			}
			if (fabs(b) < tol) b = tol;
			{  
				if (b < 0) b = tol;
				if (b < 0) b = -tol;
			}
			if (fabs(c) < tol) c = tol;
			{  
				if (c < 0) c = tol;
				if (c < 0) c = -tol;
			}
			
			diff_node_value = isovalue - node_value;
			diff_neigh_value = isovalue - neigh_value;
			
			vector_distance = Coord_Node_2 - Coord_Node_1; 
			
			if (fabs(diff_node_value) <  tolerance) weight = 0.0;
			else if (fabs(diff_neigh_value) <  tolerance) weight = 1.0;
			else
			{
				weight1 = (-b + sqrt((b * b) - (4 * a * c))) / (2 * a);
				weight2 = (-b - sqrt((b * b) - (4 * a * c))) / (2 * a);
				if (weight1 < 1.0 && weight1 > 0.0) weight = weight1;
				else if (weight2 < 1.0 && weight2 > 0.0) weight = weight2;
			}
			if (weight > 1.00005) KRATOS_WATCH("**** something's wrong! weight higher than 1! ****");
			
			for (unsigned int index=0; index!=3; ++index) //we loop the 3 coordinates)
			{
				if (Position_Node[i][0]!=Position_Node[i][1])  Coordinate_New_Node[i][index] = Coord_Node_1[index] + vector_distance[index] * weight; 
				else   Coordinate_New_Node[i][index] = Coord_Node_1[index]; //when both nodes are the same it doesnt make any sense to interpolate
				
				/// inserting the new node in the model part
				Node < 3 > ::Pointer pnode = new_model_part.CreateNewNode(List_New_Nodes[i], Coordinate_New_Node[i][0], Coordinate_New_Node[i][1], Coordinate_New_Node[i][2]);  //remember that is the new model part!!
			
				pnode->SetBufferSize(this_model_part.NodesBegin()->GetBufferSize());
				
				pnode->GetValue(FATHER_NODES).resize(0);
				pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node1.base() ) );       //saving data about fathers in the model part
				pnode->GetValue(FATHER_NODES).push_back( Node<3>::WeakPointer( *it_node2.base() ) );
				pnode-> GetValue(WEIGHT_FATHER_NODES) = weight;  
				
				pnode->X0() = (1.0 - weight) * (it_node1->X0())  +  weight * it_node2->X0();   
				pnode->Y0() = (1.0 - weight) * (it_node1->Y0())  +  weight * it_node2->Y0();    
				pnode->Z0() = (1.0 - weight) * (it_node1->Z0())  +  weight * it_node2->Z0();
			}
			weight = 0;
			weight1 = 0;
			weight2 = 0;
		}
	}
	
    ///************************************************************************************************ 
	
	void GenerateElements (ModelPart& this_model_part, ModelPart& new_model_part, const boost::numeric::ublas::vector<array_1d<int, 4 > >& Position_Node, vector<int> Elems_In_Isosurface, compressed_matrix<int>& Coord, Variable<double>& variable, int surface_number)
	{
		array_1d<double, 3 > temp_vector1; 
		array_1d<double, 3 > temp_vector2;
		array_1d<double, 3 > temp_vector3;
		array_1d<double, 3 > temp_vector4;
		array_1d<double, 3 > temp_vector5;
		array_1d<double, 3 > Coord_Node;
		array_1d<double, 3 > Coord_Node_1;
		array_1d<double, 3 > Coord_Node_2;
		array_1d<double, 3 > Coord_Node_3;
		array_1d<double, 3 > Coord_Node_4;
		array_1d<double, 3 > vector_max;
		array_1d<double, 3 > vector_min;
		array_1d<double, 3 > vector_normal;
		double dist2; //to be used when there are 2 triangles in the tetrahedra
		double dist3;
		double control;
		array_1d<int, 6 > nodes_for_2triang; //to be used when there are 2 triangles
		array_1d<double, 4 > node_value;
		unsigned int temp_int;
		unsigned int first_element=0;
		
		ElementsArrayType& rElements_old = this_model_part.Elements();
		ElementsArrayType::iterator it_begin_old = rElements_old.ptr_begin();
		ElementsArrayType::iterator it_end_old = rElements_old.ptr_end(); 
		            
		ElementsArrayType& rElements_new = new_model_part.Elements();        
		ElementsArrayType::iterator it_begin_new = rElements_new.ptr_begin(); 
		ElementsArrayType::iterator it_end_new = rElements_new.ptr_end();
		
		boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node;
		Coordinate_New_Node.resize(Position_Node.size());

		if (first_cutting_surface == false) //it means this is not the first surface that has to be created
		{
			first_element=(it_end_new-it_begin_new);
			KRATOS_WATCH("New Cutting Surface")
			KRATOS_WATCH(first_element)
		}
		else
		{
			first_element=0;                      //or we're creating the first elements of the model part
			KRATOS_WATCH("First Cutting Surface"); 
		}
		                       
		Condition const& rReferenceCondition = KratosComponents<Condition>::Get("Condition3D");         
		Properties::Pointer properties = this_model_part.GetMesh().pGetProperties(surface_number); 
		            
		int number_of_triangles =  0; 
		int current_element= 0; 
		unsigned int triangle_nodes= 0;  // number of nodes already saved (of the current element)   
		bool new_node = false ; //used to check whether the current node has been saved or not      
		          
		array_1d<int, 4 > TriangleNodesArray; //nodes of the element to be generated. 4 in case they're 2 triangles
		for (unsigned int k=0; k!=4; ++k) { TriangleNodesArray[k] = -1;} //initializing in -1, meaning we have no nodes yet
		            
		///we enter the element loop
		for (ElementsArrayType::iterator it = it_begin_old; it != it_end_old; ++it)
		{
			++current_element;
			triangle_nodes = 0; //starting, no nodes yet
			///we enter in the if for only one triangle in the tetrahedra
			if (Elems_In_Isosurface[current_element-1] == 1 ) //do not forget than can be both 1 or 2 triangles per tetrahedra. this is the simplest case. no need to check anything, we just create an element with the 3 nodes
			{  
				//checking element conectivities
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++) 
				{               
					Geometry<Node<3> >&geom = it->GetGeometry(); //i node of the element
					for(unsigned int j = 0; j < it->GetGeometry().size() ; j++) //j node of the element
					{
						new_node= true; //by default it's a new node
						int index_i = geom[i].Id() - 1;     //i node id
						int index_j = geom[j].Id() - 1; 
						for (unsigned int l=0; l!=3; ++l)
						{
							if(TriangleNodesArray[l]==Coord(index_i, index_j) || Coord(index_i, index_j) <1 ) new_node=false; //if we have already saved this node or it has not been cutted, then we have no new node to add (coord(i,j)=-1)
						}
					        //if it's a new node and the indexes are correct:
					        if (new_node && index_i<=index_j)
						{
							TriangleNodesArray[triangle_nodes]=Coord(index_i, index_j) ;
							triangle_nodes++;
						}
						if (triangle_nodes ==3) break;               //if we have already found 3 nodes then we can exit
					} //closing j node loop
					if (triangle_nodes ==3) break;               //egal
				} //closing i node loop	  
				//now we have to check that the normal of the element is in the direction of the gradient of the scalar field. To do this the value of the variable in all the nodes of the element of the old model part is evaluated. With the biggest and the smallest values a vector between those coordinates is calculated. The inner product of that vector and the normal vector of the element should be positive.
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++)          //size = 4 ; nodes per element. NOTICE WE'LL BE LOOPING THE EDGES TWICE. THIS IS A WASTE OF TIME BUT MAKES IT EASIER TO IDENTITY ELEMENTS. LOOK BELOW.
				{
				Geometry<Node<3> >&geom = it->GetGeometry();
				node_value[i] = geom[i].FastGetSolutionStepValue(variable);
				
				Coord_Node[0] = geom[i].X();
				Coord_Node[1] = geom[i].Y();
				Coord_Node[2] = geom[i].Z();
				
				if(i == 0) { Coord_Node_1 = Coord_Node; }
				if(i == 1) { Coord_Node_2 = Coord_Node; }
				if(i == 2) { Coord_Node_3 = Coord_Node; }
				if(i == 3) { Coord_Node_4 = Coord_Node; }
				}
                        
				if (node_value[0] >= node_value[1] && node_value[0] >= node_value[2] && node_value[0] >= node_value[3]) vector_max=Coord_Node_1;
				if (node_value[1] > node_value[0] && node_value[1] > node_value[2] && node_value[1] > node_value[3]) vector_max=Coord_Node_2;
				if (node_value[2] > node_value[0] && node_value[2] > node_value[1] && node_value[2] > node_value[3]) vector_max=Coord_Node_3;
				if (node_value[3] > node_value[0] && node_value[3] > node_value[1] && node_value[3] > node_value[2]) vector_max=Coord_Node_4;
				
				if (node_value[0] < node_value[1] && node_value[0] < node_value[2] && node_value[0] < node_value[3]) vector_min=Coord_Node_1;
				if (node_value[1] < node_value[0] && node_value[1] < node_value[2] && node_value[1] < node_value[3]) vector_min=Coord_Node_2;
				if (node_value[2] < node_value[0] && node_value[2] < node_value[1] && node_value[2] < node_value[3]) vector_min=Coord_Node_3;
				if (node_value[3] <= node_value[0] && node_value[3] <= node_value[1] && node_value[3] <= node_value[2]) vector_min=Coord_Node_4;
				
				vector_normal = vector_max - vector_min;
				
				ModelPart::NodesContainerType::iterator new_it_node1 =  new_model_part.Nodes().find(TriangleNodesArray[0]);
				noalias(temp_vector1) = new_it_node1->Coordinates(); //node 1
				
				ModelPart::NodesContainerType::iterator new_it_node2 = new_model_part.Nodes().find(TriangleNodesArray[1]);
				noalias(temp_vector2) = new_it_node2->Coordinates(); //node 2
				                        
				ModelPart::NodesContainerType::iterator new_it_node3 = new_model_part.Nodes().find(TriangleNodesArray[2]);
				noalias(temp_vector3) = new_it_node3->Coordinates(); //nodo 3
				
				temp_vector3 -=temp_vector1; //first edge
				temp_vector2 -=temp_vector1; //second edge
				
				MathUtils<double>::CrossProduct(temp_vector4, temp_vector2 , temp_vector3) ; //multiplying the 2 edges gives us a normal vector to the element
				
				if (inner_prod(temp_vector4,vector_normal)<0.0) //if the signs do not match then they have opposite directions
				{
					temp_int= TriangleNodesArray[2];
					TriangleNodesArray[2] =  TriangleNodesArray[1];
					TriangleNodesArray[1] =  temp_int;              //we switch 2 nodes and ready
				}
				        
				//generate new Elements
				Triangle3D3<Node<3> > geom(
						new_model_part.Nodes()(TriangleNodesArray[0]),  
						new_model_part.Nodes()(TriangleNodesArray[1]),
						new_model_part.Nodes()(TriangleNodesArray[2])	 
						);
				
				Condition::Pointer p_condition = rReferenceCondition.Create(number_of_triangles+1+first_element, geom, properties); //creating the element using the reference element. notice we are using the first element to avoid overriting nodes created by other cutting planes
				new_model_part.Conditions().push_back(p_condition);
				++number_of_triangles; 
				          
				for (int counter=0; counter!=4; ++counter) TriangleNodesArray[counter]=0;     
			}//closing if elems_in_plane==6 (1 triangle)
			
			///entering now the if for 2 triangles inside the tetrahedra.
			if (Elems_In_Isosurface[current_element-1] == 2)  //now we have 2 elements. we cant just create 2 elements with a random node order because they might overlap and not cover the whole area defined by the trapezoid
			{   //to fix this we'll first create a plane. see below
			    //checking conectivities to find nodes
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++) //node i 
				{
					Geometry<Node<3> >&geom = it->GetGeometry();
					for(unsigned int j = 0; j < it->GetGeometry().size() ; j++) //node j
					{
						new_node= true;
						int index_i = geom[i].Id() - 1;    
						int index_j = geom[j].Id() - 1; 
						for (unsigned int l=0; l!=3; ++l)
						{
							if(TriangleNodesArray[l]==Coord(index_i, index_j) || Coord(index_i, index_j) < 1 )  //same as the part with only one triangle (look above)
							{       
								new_node=false; 
							}  
						}
						if (new_node && index_i<index_j)
						{
							TriangleNodesArray[triangle_nodes]=Coord(index_i, index_j) ;
							triangle_nodes++;
						}
						if (triangle_nodes ==4) break;               //once we've found the 4 nodes we can exit
					} //closing i loop
					if (triangle_nodes ==4) break;        
				} 
				
				//now we have to start checking the angles. the easiest way (i think) is creating a new plane using the original plane and a segment created by 2 nodes
				// using crossproduct we get a perpendicular plane. (either point can be used as origin).
				//since the cuadrilateral is created by the cut of teatraedra, 
				//none of its internal angles can exceed 180 degrees and hence our new plane divides the cuadrilateral into 2 triangles if the distances to the other points have different signs (one on top and the other on the bottom of this new plane)
				//otherwise this edge is just an edge of the cuadrilateral and we have to look for another.
				//so let's begin! (we'll keep an origin node and we'll loop different nodes as the end of the segment till we find one that satisfies our criteria)
				
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++)          //size = 4 ; nodes per element. NOTICE WE'LL BE LOOPING THE EDGES TWICE. THIS IS A WASTE OF TIME BUT MAKES IT EASIER TO IDENTITY ELEMENTS. LOOK BELOW.
				{
					Geometry<Node<3> >&geom = it->GetGeometry();
					node_value[i] = geom[i].FastGetSolutionStepValue(variable);
					
					Coord_Node[0] = geom[i].X();
					Coord_Node[1] = geom[i].Y();
					Coord_Node[2] = geom[i].Z();
					
					if(i == 0) Coord_Node_1 = Coord_Node;
					if(i == 1) Coord_Node_2 = Coord_Node;
					if(i == 2) Coord_Node_3 = Coord_Node;
					if(i == 3) Coord_Node_4 = Coord_Node;
				}
				
				if (node_value[0] >= node_value[1] && node_value[0] >= node_value[2] && node_value[0] >= node_value[3]) vector_max=Coord_Node_1;
				if (node_value[1] > node_value[0] && node_value[1] > node_value[2] && node_value[1] > node_value[3]) vector_max=Coord_Node_2;
				if (node_value[2] > node_value[0] && node_value[2] > node_value[1] && node_value[2] > node_value[3]) vector_max=Coord_Node_3;
				if (node_value[3] > node_value[0] && node_value[3] > node_value[1] && node_value[3] > node_value[2]) vector_max=Coord_Node_4;
				
				if (node_value[0] < node_value[1] && node_value[0] < node_value[2] && node_value[0] < node_value[3]) vector_min=Coord_Node_1;
				if (node_value[1] < node_value[0] && node_value[1] < node_value[2] && node_value[1] < node_value[3]) vector_min=Coord_Node_2;
				if (node_value[2] < node_value[0] && node_value[2] < node_value[1] && node_value[2] < node_value[3]) vector_min=Coord_Node_3;
				if (node_value[3] <= node_value[0] && node_value[3] <= node_value[1] && node_value[3] <= node_value[2]) vector_min=Coord_Node_4;
				
				vector_normal = vector_max - vector_min;
				ModelPart::NodesContainerType::iterator it_node1 =  new_model_part.Nodes().find(TriangleNodesArray[0]);
				noalias(temp_vector1) = it_node1->Coordinates(); //nodo 1 (origin)
				int jjj;
				int kkk;
				
				for (int iii=1; iii!=4; ++iii)    //end node of the segment that will be used to create the plane (will be contained too)
				{
					it_node1=new_model_part.Nodes().find(TriangleNodesArray[iii]); //i node. we always keep node 0 as origin
					noalias(temp_vector2) = it_node1->Coordinates(); //node2 (end)
					noalias(temp_vector3) = temp_vector2 - temp_vector1; //segment 1-2
					//now i have to create the new plane
					MathUtils<double>::CrossProduct(temp_vector4, vector_normal , temp_vector3); //done. now temp_vector4 is the (normal to the) new plane, perpendicular to the one containing the triangles
					//the origin of the plane is temp_vector1 (temp_vector2 could also be used)
					//now we need to check distances to the other nodes (i+2 (let's call them jjj and i+3=kkk since we can't go futher than i=3)
					
					if (iii==1) { jjj=2; kkk=3;}
					else if (iii==2) { jjj=3; kkk=1;}
					else { jjj=2; kkk=1;} 
				
					it_node1=new_model_part.Nodes().find(TriangleNodesArray[jjj]);
					noalias(temp_vector2) = it_node1->Coordinates(); //one of the remaining nodes;

					it_node1=new_model_part.Nodes().find(TriangleNodesArray[kkk]);				
					noalias(temp_vector3) = it_node1->Coordinates(); //the other remaining node;
				
					noalias(temp_vector2) -=temp_vector1;   // minus origin point of the plane
					noalias(temp_vector3) -=temp_vector1;

					dist2= inner_prod(temp_vector2,temp_vector4);  // dot product
					dist3= inner_prod(temp_vector3,temp_vector4);
					control=dist2*dist3;   
					//and that's it. we now have to check if the distance have different signs. to do so we multiply :
					if (control<0.0) //we have the right one! one node on each side of the plane generated by nodes 0 and iii
					{
						nodes_for_2triang[0] = TriangleNodesArray[0] ;
						nodes_for_2triang[1] = TriangleNodesArray[jjj];
						nodes_for_2triang[2] = TriangleNodesArray[iii]; //finish first triangle
						nodes_for_2triang[3] = TriangleNodesArray[iii];
						nodes_for_2triang[4] = TriangleNodesArray[kkk];
						nodes_for_2triang[5] = TriangleNodesArray[0]; //finish 2nd triangle
						break; //no need to keep looking, we can exit the loop	
					} //closing the if
				}//by the time this finishes i should already have TriangleNodesArray
        
				//checking if the normal to our element is oriented correctly, just as we did when we had only 1 triangle (not commented here)
				for (int index=0 ; index !=2 ; ++index) //for triangle 1 and 2
				{
					ModelPart::NodesContainerType::iterator new_it_node1 =  new_model_part.Nodes().find(nodes_for_2triang[index*3+0]);
					noalias(temp_vector1) = new_it_node1->Coordinates(); //node 1
					
					ModelPart::NodesContainerType::iterator new_it_node2 = new_model_part.Nodes().find(nodes_for_2triang[index*3+1]);
					noalias(temp_vector2) = new_it_node2->Coordinates(); //node 2
					
					ModelPart::NodesContainerType::iterator new_it_node3 = new_model_part.Nodes().find(nodes_for_2triang[index*3+2]);
					noalias(temp_vector3) = new_it_node3->Coordinates(); //nodo 3
					
					temp_vector3 -=temp_vector1;
					temp_vector2 -=temp_vector1;
					MathUtils<double>::CrossProduct(temp_vector4, temp_vector2 , temp_vector3) ; 
					
					if (inner_prod(temp_vector4,vector_normal)<0.0)
					{
						temp_int= nodes_for_2triang[index*3+2];
						nodes_for_2triang[index*3+2] =  nodes_for_2triang[index*3+1];
						nodes_for_2triang[index*3+1] =  temp_int;              
					}
					
					Triangle3D3<Node<3> > geom(
							new_model_part.Nodes()(nodes_for_2triang[index*3+0]), 
							new_model_part.Nodes()(nodes_for_2triang[index*3+1]),
							new_model_part.Nodes()(nodes_for_2triang[index*3+2])	 
							);
					
					Condition::Pointer p_condition = rReferenceCondition.Create(number_of_triangles+1+first_element, geom, properties); 
					
					new_model_part.Conditions().push_back(p_condition);
					++number_of_triangles; 
					for (int counter_1=0; counter_1 !=4; ++counter_1) TriangleNodesArray[counter_1]=0;//resetting, just in case
					//for (int counter_2=0; counter_2 !=6; ++counter_2) nodes_for_2triang[counter_2]=0;//resetting, just in case
				}//cierro el index
			}//closing if elems_in_surface=2
	  }//closing element loops

	}
	        
	void GenerateElements (ModelPart& this_model_part, ModelPart& new_model_part, const boost::numeric::ublas::vector<array_1d<int, 4 > >& Position_Node, vector<int> Elems_In_Isosurface, compressed_matrix<int>& Coord, VariableComponent < VectorComponentAdaptor < array_1d<double,3> > >& variable, int surface_number)
	{
		array_1d<double, 3 > temp_vector1; 
		array_1d<double, 3 > temp_vector2;
		array_1d<double, 3 > temp_vector3;
		array_1d<double, 3 > temp_vector4;
		array_1d<double, 3 > temp_vector5;
		array_1d<double, 3 > Coord_Node;
		array_1d<double, 3 > Coord_Node_1;
		array_1d<double, 3 > Coord_Node_2;
		array_1d<double, 3 > Coord_Node_3;
		array_1d<double, 3 > Coord_Node_4;
		array_1d<double, 3 > vector_max;
		array_1d<double, 3 > vector_min;
		array_1d<double, 3 > vector_normal;
		double dist2; //to be used when there are 2 triangles in the tetrahedra
		double dist3;
		double control;
		array_1d<int, 6 > nodes_for_2triang; //to be used when there are 2 triangles
		array_1d< array_1d<double,3> , 4 > vector_node_value;
		array_1d<double, 3 > node_value;
		unsigned int temp_int;
		unsigned int first_element=0;
		
		ElementsArrayType& rElements_old = this_model_part.Elements();
		ElementsArrayType::iterator it_begin_old = rElements_old.ptr_begin();
		ElementsArrayType::iterator it_end_old = rElements_old.ptr_end(); 
		            
		ElementsArrayType& rElements_new = new_model_part.Elements();        
		ElementsArrayType::iterator it_begin_new = rElements_new.ptr_begin(); 
		ElementsArrayType::iterator it_end_new = rElements_new.ptr_end();
		
		boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node;
		Coordinate_New_Node.resize(Position_Node.size());
		
		if (first_cutting_surface == false)   //it means this is not the first surface that has to be created
		{
			first_element=(it_end_new-it_begin_new);
			KRATOS_WATCH("New Cutting Surface")
			KRATOS_WATCH(first_element)
		}
		
		else
		{
			first_element=0;                      //or we're creating the first elements of the model part
			KRATOS_WATCH("First Cutting Surface"); 
		}
		                       
		Condition const& rReferenceCondition = KratosComponents<Condition>::Get("Condition3D");         
		Properties::Pointer properties = this_model_part.GetMesh().pGetProperties(surface_number); 
		            
		int number_of_triangles =  0; 
		int current_element= 0; 
		unsigned int triangle_nodes= 0;  // number of nodes already saved (of the current element)   
		bool new_node = false ; //used to check whether the current node has been saved or not      
		          
		array_1d<int, 4 > TriangleNodesArray; //nodes of the element to be generated. 4 in case they're 2 triangles
		for (unsigned int k=0; k!=4; ++k) { TriangleNodesArray[k] = -1;} //initializing in -1, meaning we have no nodes yet
		            
	     ///we enter the element loop
		for (ElementsArrayType::iterator it = it_begin_old; it != it_end_old; ++it)
		{
			++current_element;
			triangle_nodes = 0; //starting, no nodes yet
		     ///we enter in the if for only one triangle in the tetrahedra
			if (Elems_In_Isosurface[current_element-1] == 1 ) //do not forget than can be both 1 or 2 triangles per tetrahedra. this is the simplest case. no need to check anything, we just create an element with the 3 nodes
			{  
				//checking element conectivities
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++) 
				{ 
					Geometry<Node<3> >&geom = it->GetGeometry(); //i node of the element
					for(unsigned int j = 0; j < it->GetGeometry().size() ; j++) //j node of the element
					{
						new_node= true; //by default it's a new node
						int index_i = geom[i].Id() - 1;     //i node id
						int index_j = geom[j].Id() - 1; 
						for (unsigned int l=0; l!=3; ++l)
						{
							if(TriangleNodesArray[l]==Coord(index_i, index_j) || Coord(index_i, index_j) < 1 ) new_node=false;  //if we have already saved this node or it has not been cutted, then we have no new node to add (coord(i,j)=-1)
						}
					        //if it's a new node and the indexes are correct:
					        if (new_node && index_i<=index_j)
						{
							TriangleNodesArray[triangle_nodes]=Coord(index_i, index_j) ;
							triangle_nodes++;
						}
						if (triangle_nodes ==3) break;               //if we have already found 3 nodes then we can exit
					} //closing j node loop
					if (triangle_nodes ==3) break;               //egal
				} //closing i node loop	  
				//now we have to check that the normal of the element is in the direction of the gradient of the scalar field. To do this the value of the variable in all the nodes of the element of the old model part is evaluated. With the biggest and the smallest values a vector between those coordinates is calculated. The inner product of that vector and the normal vector of the element should be positive.
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++)          //size = 4 ; nodes per element. NOTICE WE'LL BE LOOPING THE EDGES TWICE. THIS IS A WASTE OF TIME BUT MAKES IT EASIER TO IDENTITY ELEMENTS. LOOK BELOW.
				{
					Geometry<Node<3> >&geom = it->GetGeometry();
					node_value[i] = geom[i].FastGetSolutionStepValue(variable);
					
					Coord_Node[0] = geom[i].X();
					Coord_Node[1] = geom[i].Y();
					Coord_Node[2] = geom[i].Z();
					
					if(i == 0) { Coord_Node_1 = Coord_Node; }
					if(i == 1) { Coord_Node_2 = Coord_Node; }
					if(i == 2) { Coord_Node_3 = Coord_Node; }
					if(i == 3) { Coord_Node_4 = Coord_Node; }
				}
				                    
				if (node_value[0] >= node_value[1] && node_value[0] >= node_value[2] && node_value[0] >= node_value[3]) vector_max=Coord_Node_1;
				if (node_value[1] > node_value[0] && node_value[1] > node_value[2] && node_value[1] > node_value[3]) vector_max=Coord_Node_2;
				if (node_value[2] > node_value[0] && node_value[2] > node_value[1] && node_value[2] > node_value[3]) vector_max=Coord_Node_3;
				if (node_value[3] > node_value[0] && node_value[3] > node_value[1] && node_value[3] > node_value[2]) vector_max=Coord_Node_4;
				
				if (node_value[0] < node_value[1] && node_value[0] < node_value[2] && node_value[0] < node_value[3]) vector_min=Coord_Node_1;
				if (node_value[1] < node_value[0] && node_value[1] < node_value[2] && node_value[1] < node_value[3]) vector_min=Coord_Node_2;
				if (node_value[2] < node_value[0] && node_value[2] < node_value[1] && node_value[2] < node_value[3]) vector_min=Coord_Node_3;
				if (node_value[3] <= node_value[0] && node_value[3] <= node_value[1] && node_value[3] <= node_value[2]) vector_min=Coord_Node_4;
				
				vector_normal = vector_max - vector_min;
				
				ModelPart::NodesContainerType::iterator new_it_node1 =  new_model_part.Nodes().find(TriangleNodesArray[0]);
				noalias(temp_vector1) = new_it_node1->Coordinates(); //node 1
				                        
				ModelPart::NodesContainerType::iterator new_it_node2 = new_model_part.Nodes().find(TriangleNodesArray[1]);
				noalias(temp_vector2) = new_it_node2->Coordinates(); //node 2
				        
				ModelPart::NodesContainerType::iterator new_it_node3 = new_model_part.Nodes().find(TriangleNodesArray[2]);
				noalias(temp_vector3) = new_it_node3->Coordinates(); //nodo 3
				
				temp_vector3 -=temp_vector1; //first edge
				temp_vector2 -=temp_vector1; //second edge
				
				MathUtils<double>::CrossProduct(temp_vector4, temp_vector2 , temp_vector3) ; //multiplying the 2 edges gives us a normal vector to the element
				
				if (inner_prod(temp_vector4,vector_normal)<0.0) //if the signs do not match then they have opposite directions
				{
					temp_int= TriangleNodesArray[2];
					TriangleNodesArray[2] =  TriangleNodesArray[1];
					TriangleNodesArray[1] =  temp_int;              //we switch 2 nodes and ready
				}
				        
				//generate new Elements
				Triangle3D3<Node<3> > geom(
							new_model_part.Nodes()(TriangleNodesArray[0]),  
							new_model_part.Nodes()(TriangleNodesArray[1]),
							new_model_part.Nodes()(TriangleNodesArray[2])	 
							);
				
				Condition::Pointer p_condition = rReferenceCondition.Create(number_of_triangles+1+first_element, geom, properties); //creating the element using the reference element. notice we are using the first element to avoid overriting nodes created by other cutting planes
				new_model_part.Conditions().push_back(p_condition);
				++number_of_triangles; 
				
				for (int counter=0; counter!=4; ++counter) TriangleNodesArray[counter]=0;
			}//closing if elems_in_plane==6 (1 triangle)
		
			///entering now the if for 2 triangles inside the tetrahedra.
			if (Elems_In_Isosurface[current_element-1] == 2)  //now we have 2 elements. we cant just create 2 elements with a random node order because they might overlap and not cover the whole area defined by the trapezoid
			{   //to fix this we'll first create a plane. see below
			    //checking conectivities to find nodes
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++) //node i 
				{
					Geometry<Node<3> >&geom = it->GetGeometry();
					for(unsigned int j = 0; j < it->GetGeometry().size() ; j++) //node j
					{
						new_node= true;
						int index_i = geom[i].Id() - 1;    
						int index_j = geom[j].Id() - 1; 
						for (unsigned int l=0; l!=3; ++l)
						{
							if(TriangleNodesArray[l]==Coord(index_i, index_j) || Coord(index_i, index_j) < 1 )  //same as the part with only one triangle (look above)
							{       
								new_node=false; 
							}  
						}
						if (new_node && index_i<index_j)
						{
							TriangleNodesArray[triangle_nodes]=Coord(index_i, index_j) ;
							triangle_nodes++;
						}
						if (triangle_nodes ==4) break;               //once we've found the 4 nodes we can exit
					} //closing j loop
					if (triangle_nodes ==4) break;        
				} //closing i loop
				
				//now we have to start checking the angles. the easiest way (i think) is creating a new plane using the original plane and a segment created by 2 nodes
				// using crossproduct we get a perpendicular plane. (either point can be used as origin).
				//since the cuadrilateral is created by the cut of teatraedra, 
				//none of its internal angles can exceed 180 degrees and hence our new plane divides the cuadrilateral into 2 triangles if the distances to the other points have different signs (one on top and the other on the bottom of this new plane)
				//otherwise this edge is just an edge of the cuadrilateral and we have to look for another.
				//so let's begin! (we'll keep an origin node and we'll loop different nodes as the end of the segment till we find one that satisfies our criteria)
				
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++)          //size = 4 ; nodes per element. NOTICE WE'LL BE LOOPING THE EDGES TWICE. THIS IS A WASTE OF TIME BUT MAKES IT EASIER TO IDENTITY ELEMENTS. LOOK BELOW.
				{
					Geometry<Node<3> >&geom = it->GetGeometry();
					node_value[i] = geom[i].FastGetSolutionStepValue(variable);
					
					Coord_Node[0] = geom[i].X();
					Coord_Node[1] = geom[i].Y();
					Coord_Node[2] = geom[i].Z();
					
					if(i == 0) Coord_Node_1 = Coord_Node;
					if(i == 1) Coord_Node_2 = Coord_Node;
					if(i == 2) Coord_Node_3 = Coord_Node;
					if(i == 3) Coord_Node_4 = Coord_Node;
				}
				
				if (node_value[0] >= node_value[1] && node_value[0] >= node_value[2] && node_value[0] >= node_value[3]) vector_max=Coord_Node_1;
				if (node_value[1] > node_value[0] && node_value[1] > node_value[2] && node_value[1] > node_value[3]) vector_max=Coord_Node_2;
				if (node_value[2] > node_value[0] && node_value[2] > node_value[1] && node_value[2] > node_value[3]) vector_max=Coord_Node_3;
				if (node_value[3] > node_value[0] && node_value[3] > node_value[1] && node_value[3] > node_value[2]) vector_max=Coord_Node_4;
				
				if (node_value[0] < node_value[1] && node_value[0] < node_value[2] && node_value[0] < node_value[3]) vector_min=Coord_Node_1;
				if (node_value[1] < node_value[0] && node_value[1] < node_value[2] && node_value[1] < node_value[3]) vector_min=Coord_Node_2;
				if (node_value[2] < node_value[0] && node_value[2] < node_value[1] && node_value[2] < node_value[3]) vector_min=Coord_Node_3;
				if (node_value[3] <= node_value[0] && node_value[3] <= node_value[1] && node_value[3] <= node_value[2]) vector_min=Coord_Node_4;
				
				vector_normal = vector_max - vector_min;
				ModelPart::NodesContainerType::iterator it_node1 =  new_model_part.Nodes().find(TriangleNodesArray[0]);
				noalias(temp_vector1) = it_node1->Coordinates(); //nodo 1 (origin)
				int jjj;
				int kkk;
				
				for (int iii=1; iii!=4; ++iii)    //end node of the segment that will be used to create the plane (will be contained too)
				{
					it_node1=new_model_part.Nodes().find(TriangleNodesArray[iii]); //i node. we always keep node 0 as origin
					noalias(temp_vector2) = it_node1->Coordinates(); //node2 (end)
					noalias(temp_vector3) = temp_vector2 - temp_vector1; //segment 1-2
					//now i have to create the new plane
					MathUtils<double>::CrossProduct(temp_vector4, vector_normal , temp_vector3); //done. now temp_vector4 is the (normal to the) new plane, perpendicular to the one containing the triangles
					//the origin of the plane is temp_vector1 (temp_vector2 could also be used)
					//now we need to check distances to the other nodes (i+2 (let's call them jjj and i+3=kkk since we can't go futher than i=3)
					
					if (iii==1) { jjj=2; kkk=3;}
					else if (iii==2) { jjj=3; kkk=1;}
					else { jjj=2; kkk=1;} 
				
					it_node1=new_model_part.Nodes().find(TriangleNodesArray[jjj]);
					noalias(temp_vector2) = it_node1->Coordinates(); //one of the remaining nodes;

					it_node1=new_model_part.Nodes().find(TriangleNodesArray[kkk]);				
					noalias(temp_vector3) = it_node1->Coordinates(); //the other remaining node;
				
					noalias(temp_vector2) -=temp_vector1;   // minus origin point of the plane
					noalias(temp_vector3) -=temp_vector1;

					dist2= inner_prod(temp_vector2,temp_vector4);  // dot product
					dist3= inner_prod(temp_vector3,temp_vector4);
					control=dist2*dist3;   
					//and that's it. we now have to check if the distance have different signs. to do so we multiply :
					if (control<0.0) //we have the right one! one node on each side of the plane generated by nodes 0 and iii
					{
						nodes_for_2triang[0] = TriangleNodesArray[0] ;
						nodes_for_2triang[1] = TriangleNodesArray[jjj];
						nodes_for_2triang[2] = TriangleNodesArray[iii]; //finish first triangle
						nodes_for_2triang[3] = TriangleNodesArray[iii];
						nodes_for_2triang[4] = TriangleNodesArray[kkk];
						nodes_for_2triang[5] = TriangleNodesArray[0]; //finish 2nd triangle
						break; //no need to keep looking, we can exit the loop	
					} //closing the if
				}//by the time this finishes i should already have TriangleNodesArray
        
				//checking if the normal to our element is oriented correctly, just as we did when we had only 1 triangle (not commented here)
				for (int index=0 ; index !=2 ; ++index) //for triangle 1 and 2
				{
					ModelPart::NodesContainerType::iterator new_it_node1 =  new_model_part.Nodes().find(nodes_for_2triang[index*3+0]);
					noalias(temp_vector1) = new_it_node1->Coordinates(); //node 1
					
					ModelPart::NodesContainerType::iterator new_it_node2 = new_model_part.Nodes().find(nodes_for_2triang[index*3+1]);
					noalias(temp_vector2) = new_it_node2->Coordinates(); //node 2
					
					ModelPart::NodesContainerType::iterator new_it_node3 = new_model_part.Nodes().find(nodes_for_2triang[index*3+2]);
					noalias(temp_vector3) = new_it_node3->Coordinates(); //nodo 3
					
					temp_vector3 -=temp_vector1;
					temp_vector2 -=temp_vector1;
					MathUtils<double>::CrossProduct(temp_vector4, temp_vector2 , temp_vector3) ; 
					
					if (inner_prod(temp_vector4,vector_normal)<0.0)
					{
						temp_int= nodes_for_2triang[index*3+2];
						nodes_for_2triang[index*3+2] =  nodes_for_2triang[index*3+1];
						nodes_for_2triang[index*3+1] =  temp_int;              
					}
					
					Triangle3D3<Node<3> > geom(
							new_model_part.Nodes()(nodes_for_2triang[index*3+0]), 
							new_model_part.Nodes()(nodes_for_2triang[index*3+1]),
							new_model_part.Nodes()(nodes_for_2triang[index*3+2])	 
							);
					
					Condition::Pointer p_condition = rReferenceCondition.Create(number_of_triangles+1+first_element, geom, properties); 
					
					new_model_part.Conditions().push_back(p_condition);
					++number_of_triangles; 
					for (int counter=0; counter!=4; ++counter) TriangleNodesArray[counter]=0; //resetting, just in case
				}//cierro el index
			}//closing if elems_in_surface=2
		}//closing element loops
	}
	
	void GenerateElements (ModelPart& this_model_part, ModelPart& new_model_part, const boost::numeric::ublas::vector<array_1d<int, 4 > >& Position_Node, vector<int> Elems_In_Isosurface, compressed_matrix<int>& Coord, Variable<array_1d<double, 3 > >& variable, int surface_number)
	{
		array_1d<double, 3 > temp_vector1; 
		array_1d<double, 3 > temp_vector2;
		array_1d<double, 3 > temp_vector3;
		array_1d<double, 3 > temp_vector4;
		array_1d<double, 3 > temp_vector5;
		array_1d<double, 3 > Coord_Node;
		array_1d<double, 3 > Coord_Node_1;
		array_1d<double, 3 > Coord_Node_2;
		array_1d<double, 3 > Coord_Node_3;
		array_1d<double, 3 > Coord_Node_4;
		array_1d<double, 3 > vector_max;
		array_1d<double, 3 > vector_min;
		array_1d<double, 3 > vector_normal;
		array_1d<int, 6 > nodes_for_2triang; //to be used when there are 2 triangles
		array_1d< array_1d <double, 3 >, 4 > vector_node_value;
		array_1d<double, 4 > node_value;
		unsigned int temp_int;
		unsigned int first_element=0;
		
		ElementsArrayType& rElements_old = this_model_part.Elements();
		ElementsArrayType::iterator it_begin_old = rElements_old.ptr_begin();
		ElementsArrayType::iterator it_end_old = rElements_old.ptr_end(); 
		            
		ElementsArrayType& rElements_new = new_model_part.Elements();        
		ElementsArrayType::iterator it_begin_new = rElements_new.ptr_begin(); 
		ElementsArrayType::iterator it_end_new = rElements_new.ptr_end();
		
		boost::numeric::ublas::vector< array_1d<double, 3 > > Coordinate_New_Node;
		Coordinate_New_Node.resize(Position_Node.size());
		
		if (first_cutting_surface == false)  //it means this is not the first surface that has to be created
		{
			first_element=(it_end_new-it_begin_new);
			KRATOS_WATCH("New Cutting Surface")
			KRATOS_WATCH(first_element)
		}
		else
		{
			first_element=0;                      //or we're creating the first elements of the model part
			KRATOS_WATCH("First Cutting Surface"); 
		}
		                       
		Condition const& rReferenceCondition = KratosComponents<Condition>::Get("Condition3D");         
		Properties::Pointer properties = this_model_part.GetMesh().pGetProperties(surface_number); 
		            
		int number_of_triangles =  0; 
		int current_element= 0; 
		unsigned int triangle_nodes= 0;  // number of nodes already saved (of the current element)   
		bool new_node = false ; //used to check whether the current node has been saved or not      
		
		array_1d<int, 4 > TriangleNodesArray; //nodes of the element to be generated. 4 in case they're 2 triangles
		for (unsigned int k=0; k!=4; ++k) TriangleNodesArray[k] = -1; //initializing in -1, meaning we have no nodes yet
		            
		///we enter the element loop
		for (ElementsArrayType::iterator it = it_begin_old; it != it_end_old; ++it)
		{
			++current_element;
			triangle_nodes = 0; //starting, no nodes yet
			
			///we enter in the if for only one triangle in the tetrahedra
			if (Elems_In_Isosurface[current_element-1] == 1 ) //do not forget than can be both 1 or 2 triangles per tetrahedra. this is the simplest case. no need to check anything, we just create an element with the 3 nodes
			{                                                 //checking element conectivities
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++) 
				{               
					Geometry<Node<3> >&geom = it->GetGeometry(); //i node of the element
					for(unsigned int j = 0; j < it->GetGeometry().size() ; j++) //j node of the element
					{
						new_node= true; //by default it's a new node
						int index_i = geom[i].Id() - 1;     //i node id
						int index_j = geom[j].Id() - 1; 
						for (unsigned int l=0; l!=3; ++l)
						{
							if(TriangleNodesArray[l]==Coord(index_i, index_j) || Coord(index_i, index_j) <1 ) new_node=false; //if we have already saved this node or it has not been cutted, then we have no new node to add (coord(i,j)=-1)      
						}
						if (new_node && index_i<=index_j)  //if it's a new node and the indexes are correct:
						{
							TriangleNodesArray[triangle_nodes]=Coord(index_i, index_j) ;
							triangle_nodes++;
						}
						if (triangle_nodes ==3) break;               //if we have already found 3 nodes then we can exit
					} //closing j node loop
					if (triangle_nodes ==3) break;               //egal
				} //closing i node loop	  
				//now we have to check that the normal of the element is in the direction of the gradient of the scalar field. To do this the value of the variable in all the nodes of the element of the old model part is evaluated. With the biggest and the smallest values a vector between those coordinates is calculated. The inner product of that vector and the normal vector of the element should be positive.
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++)          //size = 4 ; nodes per element. NOTICE WE'LL BE LOOPING THE EDGES TWICE. THIS IS A WASTE OF TIME BUT MAKES IT EASIER TO IDENTITY ELEMENTS. LOOK BELOW.
				{
					Geometry<Node<3> >&geom = it->GetGeometry();
					vector_node_value[i] = geom[i].FastGetSolutionStepValue(variable);
					
					node_value[i] = sqrt((vector_node_value[i][0] * vector_node_value[i][0]) + (vector_node_value[i][1] * vector_node_value[i][1]) + (vector_node_value[i][2] * vector_node_value[i][2]));
					
					Coord_Node[0] = geom[i].X();
					Coord_Node[1] = geom[i].Y();
					Coord_Node[2] = geom[i].Z();
					
					if(i == 0) { Coord_Node_1 = Coord_Node; }
					if(i == 1) { Coord_Node_2 = Coord_Node; }
					if(i == 2) { Coord_Node_3 = Coord_Node; }
					if(i == 3) { Coord_Node_4 = Coord_Node; }
				}
				                        
				if (node_value[0] >= node_value[1] && node_value[0] >= node_value[2] && node_value[0] >= node_value[3]) vector_max=Coord_Node_1;
				if (node_value[1] > node_value[0] && node_value[1] > node_value[2] && node_value[1] > node_value[3]) vector_max=Coord_Node_2;
				if (node_value[2] > node_value[0] && node_value[2] > node_value[1] && node_value[2] > node_value[3]) vector_max=Coord_Node_3;
				if (node_value[3] > node_value[0] && node_value[3] > node_value[1] && node_value[3] > node_value[2]) vector_max=Coord_Node_4;
				
				if (node_value[0] < node_value[1] && node_value[0] < node_value[2] && node_value[0] < node_value[3]) vector_min=Coord_Node_1;
				if (node_value[1] < node_value[0] && node_value[1] < node_value[2] && node_value[1] < node_value[3]) vector_min=Coord_Node_2;
				if (node_value[2] < node_value[0] && node_value[2] < node_value[1] && node_value[2] < node_value[3]) vector_min=Coord_Node_3;
				if (node_value[3] <= node_value[0] && node_value[3] <= node_value[1] && node_value[3] <= node_value[2]) vector_min=Coord_Node_4;
				
				vector_normal = vector_max - vector_min;
				
				ModelPart::NodesContainerType::iterator new_it_node1 =  new_model_part.Nodes().find(TriangleNodesArray[0]);
				noalias(temp_vector1) = new_it_node1->Coordinates(); //node 1
				                       
				ModelPart::NodesContainerType::iterator new_it_node2 = new_model_part.Nodes().find(TriangleNodesArray[1]);
				noalias(temp_vector2) = new_it_node2->Coordinates(); //node 2
				                       
				ModelPart::NodesContainerType::iterator new_it_node3 = new_model_part.Nodes().find(TriangleNodesArray[2]);
				noalias(temp_vector3) = new_it_node3->Coordinates(); //nodo 3
				
				temp_vector3 -=temp_vector1; //first edge
				temp_vector2 -=temp_vector1; //second edge
				
				MathUtils<double>::CrossProduct(temp_vector4, temp_vector2 , temp_vector3) ; //multiplying the 2 edges gives us a normal vector to the element
				
				if (inner_prod(temp_vector4,vector_normal)<0.0) //if the signs do not match then they have opposite directions
				{
					temp_int= TriangleNodesArray[2];
					TriangleNodesArray[2] =  TriangleNodesArray[1];
					TriangleNodesArray[1] =  temp_int;              //we switch 2 nodes and ready
				}
				        
					//generate new Elements
					Triangle3D3<Node<3> > geom(
							new_model_part.Nodes()(TriangleNodesArray[0]),  
							new_model_part.Nodes()(TriangleNodesArray[1]),
							new_model_part.Nodes()(TriangleNodesArray[2])	 
							);
					
					Condition::Pointer p_condition = rReferenceCondition.Create(number_of_triangles+1+first_element, geom, properties); //creating the element using the reference element. notice we are using the first element to avoid overriting nodes created by other cutting planes
					new_model_part.Conditions().push_back(p_condition);
					++number_of_triangles; 
				          
					for (int counter=0; counter!=4; ++counter) TriangleNodesArray[counter]=0;     
			}//closing if elems_in_plane==6 (1 triangle)
			
			///entering now the if for 2 triangles inside the tetrahedra.
			if (Elems_In_Isosurface[current_element-1] == 2)  //now we have 2 elements. we cant just create 2 elements with a random node order because they might overlap and not cover the whole area defined by the trapezoid
			{                                                 //to fix this we'll first create a plane. see below
				//checking conectivities to find nodes
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++) //node i 
				{
					Geometry<Node<3> >&geom = it->GetGeometry();
					for(unsigned int j = 0; j < it->GetGeometry().size() ; j++) //node j
					{
						new_node= true;
						int index_i = geom[i].Id() - 1;    
						int index_j = geom[j].Id() - 1; 
						for (unsigned int l=0; l!=3; ++l)
						{
							if(TriangleNodesArray[l]==Coord(index_i, index_j) || Coord(index_i, index_j) < 1 ) new_node=false; //same as the part with only one triangle (look above)
						}
						if (new_node && index_i<index_j)
						{
							TriangleNodesArray[triangle_nodes]=Coord(index_i, index_j) ;
							triangle_nodes++;
						}
						if (triangle_nodes ==4) break;               //once we've found the 4 nodes we can exit
					} //closing i loop
					if (triangle_nodes ==4) break;        
				} 
				//now we have to start checking the angles. the easiest way (i think) is creating a new plane using the original plane and a segment created by 2 nodes
				// using crossproduct we get a perpendicular plane. (either point can be used as origin).
				//since the cuadrilateral is created by the cut of teatraedra, 
				//none of its internal angles can exceed 180 degrees and hence our new plane divides the cuadrilateral into 2 triangles if the distances to the other points have different signs (one on top and the other on the bottom of this new plane)
				//otherwise this edge is just an edge of the cuadrilateral and we have to look for another.
				//so let's begin! (we'll keep an origin node and we'll loop different nodes as the end of the segment till we find one that satisfies our criteria)
				int jjj;
				int kkk;
				for (int iii=1; iii!=4; ++iii)    //end node of the segment that will be used to create the plane (will be contained too)
				{
					if (iii==1) { jjj=2; kkk=3;}
					else if (iii==2) { jjj=3; kkk=1;}
					else { jjj=2; kkk=1;} 
					
					nodes_for_2triang[0] = TriangleNodesArray[0] ;
					nodes_for_2triang[1] = TriangleNodesArray[jjj];
					nodes_for_2triang[2] = TriangleNodesArray[iii]; //finish first triangle
					nodes_for_2triang[3] = TriangleNodesArray[iii];
					nodes_for_2triang[4] = TriangleNodesArray[kkk];
					nodes_for_2triang[5] = TriangleNodesArray[0]; //finish 2nd triangle
				}//by the time this finishes i should already have TriangleNodesArray
				for(unsigned int i = 0; i < it->GetGeometry().size() ; i++)          //size = 4 ; nodes per element. NOTICE WE'LL BE LOOPING THE EDGES TWICE. THIS IS A WASTE OF TIME BUT MAKES IT EASIER TO IDENTITY ELEMENTS. LOOK BELOW.
				{
					Geometry<Node<3> >&geom = it->GetGeometry();
					vector_node_value[i] = geom[i].FastGetSolutionStepValue(variable);
					
					node_value[i] = sqrt((vector_node_value[i][0] * vector_node_value[i][0]) + (vector_node_value[i][1] * vector_node_value[i][1]) + (vector_node_value[i][2] * vector_node_value[i][2]));
					
					Coord_Node[0] = geom[i].X();
					Coord_Node[1] = geom[i].Y();
					Coord_Node[2] = geom[i].Z();
					
					if(i == 0) Coord_Node_1 = Coord_Node;
					if(i == 1) Coord_Node_2 = Coord_Node;
					if(i == 2) Coord_Node_3 = Coord_Node;
					if(i == 3) Coord_Node_4 = Coord_Node;
				}
				
				if (node_value[0] >= node_value[1] && node_value[0] >= node_value[2] && node_value[0] >= node_value[3]) vector_max=Coord_Node_1;
				if (node_value[1] > node_value[0] && node_value[1] > node_value[2] && node_value[1] > node_value[3]) vector_max=Coord_Node_2;
				if (node_value[2] > node_value[0] && node_value[2] > node_value[1] && node_value[2] > node_value[3]) vector_max=Coord_Node_3;
				if (node_value[3] > node_value[0] && node_value[3] > node_value[1] && node_value[3] > node_value[2]) vector_max=Coord_Node_4;
				
				if (node_value[0] < node_value[1] && node_value[0] < node_value[2] && node_value[0] < node_value[3]) vector_min=Coord_Node_1;
				if (node_value[1] < node_value[0] && node_value[1] < node_value[2] && node_value[1] < node_value[3]) vector_min=Coord_Node_2;
				if (node_value[2] < node_value[0] && node_value[2] < node_value[1] && node_value[2] < node_value[3]) vector_min=Coord_Node_3;
				if (node_value[3] <= node_value[0] && node_value[3] <= node_value[1] && node_value[3] <= node_value[2]) vector_min=Coord_Node_4;
				
				vector_normal = vector_max - vector_min;
				
				//checking if the normal to our element is oriented correctly, just as we did when we had only 1 triangle (not commented here)
				
				for (int index=0 ; index !=2 ; ++index) //for triangle 1 and 2
				{
					ModelPart::NodesContainerType::iterator new_it_node1 =  new_model_part.Nodes().find(nodes_for_2triang[index*3+0]);
					noalias(temp_vector1) = new_it_node1->Coordinates(); //node 1
					
					ModelPart::NodesContainerType::iterator new_it_node2 = new_model_part.Nodes().find(nodes_for_2triang[index*3+1]);
					noalias(temp_vector2) = new_it_node2->Coordinates(); //node 2
					
					ModelPart::NodesContainerType::iterator new_it_node3 = new_model_part.Nodes().find(nodes_for_2triang[index*3+2]);
					noalias(temp_vector3) = new_it_node3->Coordinates(); //nodo 3
				        
					temp_vector3 -=temp_vector1;
					temp_vector2 -=temp_vector1;
					MathUtils<double>::CrossProduct(temp_vector4, temp_vector2 , temp_vector3) ; 
					
					if (inner_prod(temp_vector4,vector_normal)<0.0)
					{
						temp_int= nodes_for_2triang[index*3+2];
						nodes_for_2triang[index*3+2] =  nodes_for_2triang[index*3+1];
						nodes_for_2triang[index*3+1] =  temp_int;              
					}
					
					Triangle3D3<Node<3> > geom(
							new_model_part.Nodes()(nodes_for_2triang[index*3+0]), 
							new_model_part.Nodes()(nodes_for_2triang[index*3+1]),
							new_model_part.Nodes()(nodes_for_2triang[index*3+2]) 
							);
					
					Condition::Pointer p_condition = rReferenceCondition.Create(number_of_triangles+1+first_element, geom, properties); 
					
					new_model_part.Conditions().push_back(p_condition);
					++number_of_triangles; 
					for (int counter_1=0; counter_1!=4; ++counter_1) TriangleNodesArray[counter_1]=0;//resetting, just in case
					//for (int counter_2=0; counter_2!=6; ++counter_2) nodes_for_2triang[counter_2]=0;//resetting, just in case
				}//cierro el index
			}//closing if elems_in_surface=2
		}//closing element loops
	}
	
     ///************************************************************************************************
     ///************************************************************************************************       

     ///UPDATECUTDATA: THIS FUNCTION UPDATES THE DATA OF THE NEW MODEL PART READING FROM THE DATA OF THE FATHER NODES (and weight factor)   
		/** interpolates data form the origin model part into the new model part (containing the cutting planes)
		 * @param new_model_part . destination model part
		 * @param old_model_part . origin model part
		 **/
     void UpdateCutData( ModelPart& new_model_part, ModelPart& old_model_part) 
     {
		KRATOS_WATCH("Updating Cut Data");
		int step_data_size = old_model_part.GetNodalSolutionStepDataSize();
		
		//looping the nodes, no data is assigned to elements
		for (ModelPart::NodesContainerType::iterator it = new_model_part.NodesBegin(); it != new_model_part.NodesEnd(); it++)
		{	
			double* node0_data = it->GetValue(FATHER_NODES)[0].SolutionStepData().Data(0); //current step only, (since we'll call this every timestep
			double* node1_data = it->GetValue(FATHER_NODES)[1].SolutionStepData().Data(0);
			double    weight   = it->GetValue(WEIGHT_FATHER_NODES);
			double* step_data = (it)->SolutionStepData().Data(0);
				
			//now we only have to copy the information from node_data to step_data
			for(int j= 0; j< step_data_size; j++)  //looping all the variables and interpolating using weight
			{
				step_data[j] = (1.0-weight)*node0_data[j] + weight*node1_data[j]; 
			}
		}//closing node loop
	}//closing subroutine

     ///************************************************************************************************
     ///************************************************************************************************
     
     ///DELETECUTDATA: THIS FUNCTION DELETES THE MESH FROM THE PREVIOUS TIME STEP   
		/** deletes data
		 * @param new_model_part . cut model part
		 **/
		
	void DeleteCutData(ModelPart& new_model_part) 
	{
		 KRATOS_WATCH("Deleting Cut Data");
		 new_model_part.Nodes().clear();
		 new_model_part.Conditions().clear();
		 new_model_part.Elements().clear();
	}
	
     ///************************************************************************************************
     ///************************************************************************************************  
	private:
	bool first_cutting_surface; // to avoid checking again if we're working on the new cutting surface or if some have already been created
	
	};
}

#endif
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Pavel Ryzhakov



#if !defined(KRATOS_ADD_PERIODIC_CONDITIONS_PROCESS_INCLUDED )
#define  KRATOS_ADD_PERIODIC_CONDITIONS_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
//#include "custom_elements/updated_lagrangian_fluid.h"
//#include "custom_elements/updated_lagrangian_fluid3D.h"
//#include "custom_elements/updated_lagrangian_fluid_inc.h"
//#include "custom_elements/updated_lagrangian_fluid3D_inc.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/model_part.h"
#include "includes/condition.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"

#include "ULF_application.h"
#include "custom_conditions/tangent_velocity_periodic_condition_2D.h"
#include "custom_conditions/tangent_velocity_periodic_condition_3D.h"

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
	Update the PRESSURE_FORCE on the nodes


*/

class AddPeriodicConditionsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(AddPeriodicConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AddPeriodicConditionsProcess(ModelPart& model_part): mr_model_part(model_part) 
    //ModelPart& fluid_model_part, ModelPart& structure_model_part, ModelPart& combined_model_part)
    //: mr_fluid_model_part(fluid_model_part), mr_structure_model_part(structure_model_part), mr_combined_model_part(combined_model_part)
    {
	//Check();

	KRATOS_CHECK_VARIABLE_KEY( NODE_PAIR_X_COMPONENT )
        KRATOS_CHECK_VARIABLE_KEY( NODE_PAIR_Y_COMPONENT )
        KRATOS_CHECK_VARIABLE_KEY( NODE_PAIR_PRESSURE)
        KRATOS_CHECK_VARIABLE_KEY( NODE_PAIR_X_COMPONENT_ANTIPERIODIC)
        KRATOS_CHECK_VARIABLE_KEY( NODE_PAIR_Y_COMPONENT_ANTIPERIODIC)
	//KRATOS_WATCH(" INSIDE ADD WALL NODES CONSTRUCTOR") 
    }

    /// Destructor.
    ~AddPeriodicConditionsProcess() //override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    //	void operator()()
    //	{
    //		MergeParts();
    //	}


    ///@}
    ///@name Operations
    ///@{

 
    //works on a cube. Note that we dont prescribe periodics on the corner nodes!
    void AddTangentConditions3D(const double x_low ,const double x_high, const double y_low, const double y_high, const double z_low, const double z_high  ) 
		{
			KRATOS_TRY
			
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing LAGRANGE_MULTIPLIER_TANGENT_VELOCITY variable on solution step data","");

			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(LAGRANGE_MULTIPLIER_SECOND_TANGENT_VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing LAGRANGE_MULTIPLIER_TANGENT_VELOCITY variable on solution step data","");

			
			unsigned int condition_number=1;
			if(mr_model_part.Conditions().size()!=0)
			{
				ModelPart::ConditionsContainerType::iterator lastcondition = (mr_model_part.ConditionsEnd()-1);
				condition_number = lastcondition->Id();
			}
			

				
			//now we have to loop the elements and create the conditions! no easy task!
			//LEFT AND RIGHT WALLS PAIR
			std::vector<double> left_wall_y_coord; //distance from the 0,0 point to the closest point of the line segment
			std::vector<double> right_wall_y_coord; //distance from the 0,0 point to the furthest point of the line segment
			std::vector<double> left_wall_z_coord; //distance from the 0,0 point to the closest point of the line segment
			std::vector<double> right_wall_z_coord; //distance from the 0,0 point to the furthest point of the line segment
			std::vector<Node<3>::Pointer>left_wall_pointers; //just pointers to nodes
			std::vector<Node<3>::Pointer>right_wall_pointers; //just pointers to nodes
			
			//UP AND LOW WALLS PAIR
			std::vector<double> low_wall_x_coord; //distance from the 0,0 point to the closest point of the line segment
			std::vector<double> up_wall_x_coord; //distance from the 0,0 point to the furthest point of the line segment
			std::vector<double> low_wall_z_coord; //distance from the 0,0 point to the closest point of the line segment
			std::vector<double> up_wall_z_coord; //distance from the 0,0 point to the furthest point of the line segment
			std::vector<Node<3>::Pointer>low_wall_pointers; //just pointers to nodes
			std::vector<Node<3>::Pointer>up_wall_pointers; //just pointers to nodes

			//UP AND LOW WALLS PAIR
			std::vector<double> back_wall_x_coord; //distance from the 0,0 point to the closest point of the line segment
			std::vector<double> front_wall_x_coord; //distance from the 0,0 point to the furthest point of the line segment
			std::vector<double> back_wall_y_coord; //distance from the 0,0 point to the closest point of the line segment
			std::vector<double> front_wall_y_coord; //distance from the 0,0 point to the furthest point of the line segment
			std::vector<Node<3>::Pointer>back_wall_pointers; //just pointers to nodes
			std::vector<Node<3>::Pointer>front_wall_pointers; //just pointers to nodes

			std::vector<bool> taken_vector_right; //we have to make sure that all the elements have found their matches
			std::vector<bool> taken_vector_up; //we have to make sure that all the elements have found their matches
			std::vector<bool> taken_vector_front; //we have to make sure that all the elements have found their matches

			//first loop, creating list of nodes in each of the four walls (should be modified to work also in 3d, adding 2 more walls.
			//first loop, creating list of nodes in each of the 6 walls
			for(unsigned int ii=0; ii<mr_model_part.Nodes().size(); ii++)
			{
				ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin() + ii;
				if(inode->X()<=x_low) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					left_wall_y_coord.push_back(inode->Y());
					left_wall_z_coord.push_back(inode->Z());
					left_wall_pointers.push_back(*inode.base());
					//taken_right.push_back(false);
				}
				else if(inode->X()>=x_high) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					right_wall_y_coord.push_back(inode->Y());
					right_wall_z_coord.push_back(inode->Z());
					right_wall_pointers.push_back(*inode.base());
					taken_vector_right.push_back(false);
				}
				if(inode->Y()<=y_low) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					low_wall_x_coord.push_back(inode->X());
					low_wall_z_coord.push_back(inode->Z());
					low_wall_pointers.push_back(*inode.base());
					//taken_right.push_back(false);
				}
				else if(inode->Y()>=y_high) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					up_wall_x_coord.push_back(inode->X());
					up_wall_z_coord.push_back(inode->Z());
					up_wall_pointers.push_back(*inode.base());
					taken_vector_up.push_back(false);
				}

				if(inode->Z()<=z_low) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					back_wall_x_coord.push_back(inode->X());
					back_wall_y_coord.push_back(inode->Y());
					back_wall_pointers.push_back(*inode.base());
					//taken_right.push_back(false);
				}
				else if(inode->Z()>=z_high) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					front_wall_x_coord.push_back(inode->X());
					front_wall_y_coord.push_back(inode->Y());
					front_wall_pointers.push_back(*inode.base());
					taken_vector_front.push_back(false);
				}
				
			}
			KRATOS_WATCH("Left and right walls")
			KRATOS_WATCH(left_wall_y_coord)
			KRATOS_WATCH(right_wall_y_coord)
			KRATOS_WATCH(left_wall_z_coord)
			KRATOS_WATCH(right_wall_z_coord)

			KRATOS_WATCH("Up and down walls")
			KRATOS_WATCH(up_wall_x_coord)
			KRATOS_WATCH(low_wall_x_coord)
			KRATOS_WATCH(up_wall_z_coord)			
			KRATOS_WATCH(low_wall_z_coord)

			KRATOS_WATCH("Back and front walls")
			KRATOS_WATCH(back_wall_x_coord)
			KRATOS_WATCH(front_wall_x_coord)
			KRATOS_WATCH(back_wall_y_coord)
			KRATOS_WATCH(front_wall_y_coord)



			//after looping all the elements, we should have found all the segments that are needed to create the conditions.
			//they must match EXACTLY! otherwise this will return an error.				
			if (up_wall_x_coord.size()==low_wall_x_coord.size()     && up_wall_z_coord.size()==low_wall_z_coord.size() &&
			    right_wall_y_coord.size()==left_wall_y_coord.size() && right_wall_z_coord.size()==left_wall_z_coord.size() &&  
			    back_wall_x_coord.size()==front_wall_x_coord.size() && back_wall_y_coord.size()==front_wall_y_coord.size())
			{
				
				Condition const& rReferenceCondition = KratosComponents<Condition>::Get("TangentVelocityPeriodicCondition3D2N");         //conditiontype for conds with nodes on cube faces
				Condition const& rReferenceConditionEdge = KratosComponents<Condition>::Get("TangentVelocityPeriodicEdgeCondition3D2N"); //condition type for conds with nodes on cube edges
				Condition const& rReferenceConditionVertex = KratosComponents<Condition>::Get("TangentVelocityPeriodicVertexCondition3D2N"); //condition type for verices
				Condition const& rReferenceSecondConditionVertex = KratosComponents<Condition>::Get("SecondTangentVelocityPeriodicVertexCondition3D2N"); //condition type for conds with nodes on cube

				Condition const& rReferenceConditionDummy = KratosComponents<Condition>::Get("Condition3D2N");         //condition type

				Condition const& rReferenceConditionNormalEdge = KratosComponents<Condition>::Get("TangentVelocityPeriodicNormalToEdgeCondition3D2N");

				Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0); 		//this will allow us later to turn this layer on/off in GID

				//we loop on the LEFT boundary, trying to find the one that matches with the cooridantes of their nodes in the right boundary:
				for(unsigned int ii=0; ii<left_wall_y_coord.size(); ii++)
				{
					const double this_y_coordinate = left_wall_y_coord[ii];
					const double this_z_coordinate = left_wall_z_coord[ii];
					bool found_match=false;
					for(unsigned int jj=0; jj<right_wall_y_coord.size(); jj++)
					{
						if(taken_vector_right[jj]==false) //ok, this one is still free! it could be our match
						{
							if( fabs(this_y_coordinate-right_wall_y_coord[jj])<1.0e-6 &&  fabs(this_z_coordinate-right_wall_z_coord[jj])<1.0e-6  )
							{
								taken_vector_right[jj]=true;
								//adding the condition to the mdpa:
								//if (coordinate > y_low && coordinate < y_high) //then it means we are not corner nodes.
								//{
								Line3D2<Node<3> > geometry(
									left_wall_pointers[ii],  //condition to be added
									right_wall_pointers[jj]);
								
								left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =right_wall_pointers[jj]->Id();
								right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = left_wall_pointers[ii]->Id();
								left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_PRESSURE) =right_wall_pointers[jj]->Id();
								right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_PRESSURE) = left_wall_pointers[ii]->Id();

								//now we check if we can also do it for the tangential component: //NO CORNER NODES!
								if (this_y_coordinate > y_low && this_y_coordinate < y_high && this_z_coordinate > z_low && this_z_coordinate < z_high) 
								//WE INCLUDE THE inner ppart of left and right wall
								//if (this_z_coordinate > z_low && this_z_coordinate < z_high)
								{
									
									Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									//left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = right_wall_pointers[jj]->Id();
									//right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = left_wall_pointers[ii]->Id();
									

								}
								//Front and back vertical edges of the left and right walls: need to be treated in a specialy way.. as they have only one tangent velocity component: uy (other components are normal, thus periodic and are taken care inside the builder and solver. Note that the corner vertices are still exclused: there all velocity components are normal to either of the cube faces. thus its taken care by builder and solver
								else if ((this_y_coordinate > y_low && this_y_coordinate < y_high) || (this_z_coordinate > z_low && this_z_coordinate < z_high)) 
								//INNER PART OF THE BACK AND FRONT WALL EDGES
								//if (this_x_coordinate > x_low && this_x_coordinate < x_high) 								
								{
									
									//Y component									
									Condition::Pointer p_condition = rReferenceConditionEdge.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);	
									if (this_z_coordinate<0.0000001 || (this_z_coordinate > z_low && this_z_coordinate < z_high))// && this_y_coordinate>y_high) )	
										{//FOR Z COMPONENT
									Condition::Pointer p_condition1 = rReferenceConditionNormalEdge.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition1);	
										}						
									//low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =up_wall_pointers[jj]->Id();
									//up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = low_wall_pointers[ii]->Id();
									 

								}
								//VERTEX CONDITIONS
								//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
								//VERTEX CONDITIONS for Y: low edge (z-orienteation) of left wall
								else if (this_y_coordinate <= y_low && this_z_coordinate <= z_low) //in this node cond for Y only
								{
									// FOR Y VEL : 1 conditions will be created... they are located at the vertices lower edge of the left wall..
									Condition::Pointer p_condition = rReferenceConditionVertex.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
								}
								else if (this_y_coordinate <= y_low && this_z_coordinate >= z_high) //IN THIS NODE: conds for both Y and Z
								{
									// FOR Y VEL 1 conditions will be created... 
									Condition::Pointer p_condition = rReferenceConditionVertex.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);

									// FOR Z VEL : 1 conditions will be created... they are located at the vertices lower edge of the left wall..
									Condition::Pointer p_condition1 = rReferenceSecondConditionVertex.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition1);
								}

								//VERTEX CONDITIONS for Z: front vertical edge (y-orientation) of left wall
								else if (this_y_coordinate >= y_high && this_z_coordinate >= z_high) //in this node: con for z only
								{	

									// FOR Z VEL : 1 conditions will be created... they are located at the vertices lower edge of the left wall..
									Condition::Pointer p_condition = rReferenceSecondConditionVertex.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);

								}
								/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
								else
								{
									Condition::Pointer p_condition = rReferenceConditionDummy.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
								}
								

								found_match=true;							
								break;
							}	
						}
					}
					//if we got to the end of the list and still found no matching, we have a problem!
					if(found_match==false)
						KRATOS_THROW_ERROR(std::logic_error, "FINISHED LOOPING LEFT and RIGHT wall elements AND FOUND NO MATCHING. low_coordinate = ", this_y_coordinate );
				}
				//now for the upper and lower
				for(unsigned int ii=0; ii<low_wall_x_coord.size(); ii++)
				{
					const double this_x_coordinate = low_wall_x_coord[ii];
					const double this_z_coordinate = low_wall_z_coord[ii];
					
					bool found_match=false;
					for(unsigned int jj=0; jj<up_wall_x_coord.size(); jj++)
					{
						if(taken_vector_up[jj]==false) //ok, this one is still free! it could be our match
						{
							if( fabs(this_x_coordinate-up_wall_x_coord[jj])<1.0e-6 && fabs(this_z_coordinate-up_wall_z_coord[jj])<1.0e-6)
							{	
								taken_vector_up[jj]=true;
								//adding the condition to the mdpa:
								//if (coordinate > y_low && coordinate < y_high) //then it means we are not corner nodes.
								//{
								Line3D2<Node<3> > geometry(
									low_wall_pointers[ii],  //condition to be added
									up_wall_pointers[jj]);
								//Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
								//mr_model_part.Conditions().push_back(p_condition);
								
								low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) =up_wall_pointers[jj]->Id();
								up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = low_wall_pointers[ii]->Id();
								low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_PRESSURE) =up_wall_pointers[jj]->Id();
								up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_PRESSURE) = low_wall_pointers[ii]->Id();
								//now we check if we can also do it for the tangential component: NO CORNER NODES!
								if (this_x_coordinate > x_low && this_x_coordinate < x_high && this_z_coordinate > z_low && this_z_coordinate < z_high) 
								//WE INCLUDE THE MIDDLE PART OF THE LOWE AND UPPER WALLS INCLUDING THE EDGES
								//if (this_x_coordinate > x_low && this_x_coordinate < x_high)
								{
									
									Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									//low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =up_wall_pointers[jj]->Id();
									//up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = low_wall_pointers[ii]->Id();
									 

								}
								//Back and front horizonatl edges of the up and top walls: need to be treated in a specialy way.. as they have only one tangent velocity component: ux (other components are normal, thus periodic and are taken care inside the builder and solver. Note that the corner vertices are still exclused: there all velocity components are normal to either of the cube faces. thus its taken care by builder and solver
								else if ((this_x_coordinate > x_low && this_x_coordinate < x_high) || (this_z_coordinate > z_low && this_z_coordinate < z_high))						
								{
									
									Condition::Pointer p_condition = rReferenceConditionEdge.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									//low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =up_wall_pointers[jj]->Id();
									//up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = low_wall_pointers[ii]->Id();
									if (this_z_coordinate>0.999999 || (this_z_coordinate > z_low && this_z_coordinate < z_high && this_x_coordinate<0.00001 ))// && this_y_coordinate>y_high) )	
										{//FOR Z COMPONENT
									Condition::Pointer p_condition1 = rReferenceConditionNormalEdge.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition1);	
										}				
									 

								}
								//VERTEX CONDITIONS
								////////////////////////////////////////////////////////////////////////////////////////
								else if ((this_x_coordinate >= x_high && this_z_coordinate <= z_low) || (this_x_coordinate >= x_high && this_z_coordinate >= z_high) ) 
								{
									//2 conditions will be created for x-vel..they are located at the vertices on the rright side (Lag dof will be stored up (sec node))
									Condition::Pointer p_condition = rReferenceConditionVertex.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									//and only in one of thse nodes we add the condition for Z:
									if (this_x_coordinate >= x_high && this_z_coordinate <= z_low)	
										{
										Condition::Pointer p_condition1 = rReferenceSecondConditionVertex.Create(++condition_number, geometry, properties);
										mr_model_part.Conditions().push_back(p_condition1);	
										}
								}								
								////////////////////////////////////////////////////////////////////////////////////////////
								else
								{
									Condition::Pointer p_condition = rReferenceConditionDummy.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
								}
								
								found_match=true;							
								break;
							}	
						}
					}
					//if we got to the end of the list and still found no matching, we have a problem!
					if(found_match==false)
						KRATOS_THROW_ERROR(std::logic_error, "FINISHED LOOPING UP and DOWN wall elements AND FOUND NO MATCHING. low_coordinate = ", this_x_coordinate );
				}
				//now for the back and front
				
				for(unsigned int ii=0; ii<back_wall_x_coord.size(); ii++)
				{
					const double this_x_coordinate = back_wall_x_coord[ii];
					const double this_y_coordinate = back_wall_y_coord[ii];
					
					bool found_match=false;
					for(unsigned int jj=0; jj<front_wall_x_coord.size(); jj++)
					{
						if(taken_vector_front[jj]==false) //ok, this one is still free! it could be our match
						{
							if( fabs(this_x_coordinate-front_wall_x_coord[jj])<1.0e-6 && fabs(this_y_coordinate-front_wall_y_coord[jj])<1.0e-6)
							{	
								taken_vector_front[jj]=true;
								//adding the condition to the mdpa:
								//if (coordinate > y_low && coordinate < y_high) //then it means we are not corner nodes.
								//{
								Line3D2<Node<3> > geometry(
									back_wall_pointers[ii],  //condition to be added
									front_wall_pointers[jj]);
								//Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
								//mr_model_part.Conditions().push_back(p_condition);
								
								back_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Z_COMPONENT) =front_wall_pointers[jj]->Id();
								front_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Z_COMPONENT) = back_wall_pointers[ii]->Id();
								back_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_PRESSURE) =front_wall_pointers[jj]->Id();
								front_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_PRESSURE) = back_wall_pointers[ii]->Id();
								//now we check if we can also do it for the tangential component: NO CORNER NODES!
								if (this_x_coordinate > x_low && this_x_coordinate < x_high && this_y_coordinate > y_low && this_y_coordinate < y_high) 
								//INNER PART OF THE BACK AND FRONT WALL EDGES
								//if (this_x_coordinate > x_low && this_x_coordinate < x_high) 								
								{
									
									Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									//low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =up_wall_pointers[jj]->Id();
									//up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = low_wall_pointers[ii]->Id();
									 

								}
								//Upper and lower edges of the front and back walls: need to be treated in a specialy way.. as they have only one tangent velocity component: ux (other components are normal, thus periodic and are taken care inside the builder and solver. Note that the corner vertices are still exclused: there all velocity components are normal to either of the cube faces. thus its taken care by builder and solver
								else if ((this_x_coordinate > x_low && this_x_coordinate < x_high && this_y_coordinate>0.99999)  || (this_y_coordinate > y_low && this_y_coordinate < y_high && this_x_coordinate>0.999999) )
								//INNER PART OF THE BACK AND FRONT WALL EDGES
								//if (this_x_coordinate > x_low && this_x_coordinate < x_high) 								
								{
								//KRATOS_WATCH("NOT DOING ANYTHING-.... JUST TO CHECK!!")	
									Condition::Pointer p_condition = rReferenceConditionEdge.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);


									Condition::Pointer p_condition1 = rReferenceConditionNormalEdge.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition1);
									//low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =up_wall_pointers[jj]->Id();
									//up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = low_wall_pointers[ii]->Id();
									 

								}
								//VERTEX CONDITIONS: FIRST CASE: condition for Y (var: LAG), SECOND CASE: cond. for X (var: SECOND_LAG)
								else if (this_x_coordinate >= x_high && this_y_coordinate <= y_low) 
								{
									//1 cond for vy will be created...it  is located at righ down vertex (dof is stored at the front wall-in the second node of the cond)
									Condition::Pointer p_condition = rReferenceConditionVertex.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);	

									//1 cond for vx be created... it  is located at righ down vertex (dof is stored at the front wall-in the second node of the cond)
									Condition::Pointer p_condition1 = rReferenceSecondConditionVertex.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition1);
								}
								else
								{
									Condition::Pointer p_condition = rReferenceConditionDummy.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
								}
								
								found_match=true;							
								break;
							}	
						}
					}
					//if we got to the end of the list and still found no matching, we have a problem!
					if(found_match==false)
						KRATOS_THROW_ERROR(std::logic_error, "FINISHED LOOPING BACK and FRONT wall elements FOUND NO MATCHING. low_coordinate = ", this_x_coordinate);
				}
				
			}
			else //too bad! the number does not match, the geometry is not symmetric or maybe the tolerance was wrong?
			{
				KRATOS_WATCH(up_wall_x_coord.size())
				KRATOS_WATCH(low_wall_x_coord.size())
				KRATOS_WATCH(up_wall_z_coord.size())
				KRATOS_WATCH(low_wall_z_coord.size())

				KRATOS_WATCH(right_wall_y_coord.size())
				KRATOS_WATCH(left_wall_y_coord.size())
				KRATOS_WATCH(right_wall_z_coord.size())
				KRATOS_WATCH(left_wall_z_coord.size())

				KRATOS_WATCH(back_wall_x_coord.size())
				KRATOS_WATCH(front_wall_x_coord.size())
				KRATOS_WATCH(back_wall_y_coord.size())
				KRATOS_WATCH(front_wall_y_coord.size())
				KRATOS_THROW_ERROR(std::logic_error, "DIFFERENT NUMBER OF NODES ON BOUNDARIES!!", "");
			}
			
			
			std::cout << "Finished adding periodic tangent velocity condtions with a jump " << condition_number << std::endl;

			
			KRATOS_CATCH("")
		} 

	//////////////////////////////////////////////////////////////////////////////////////
	void AddTangentConditions2D(const double x_low ,const double x_high, const double y_low, const double y_high  ) 
		{
			KRATOS_TRY
			
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing LAGRANGE_MULTIPLIER_TANGENT_VELOCITY variable on solution step data","");

			
			unsigned int condition_number=1;
			if(mr_model_part.Conditions().size()!=0)
			{
				ModelPart::ConditionsContainerType::iterator lastcondition = (mr_model_part.ConditionsEnd()-1);
				condition_number = lastcondition->Id();
			}
			

				
			//now we have to loop the elements and create the conditions! no easy task!
			std::vector<double> left_wall_y_coord; //distance from the 0,0 point to the closest point of the line segment
			std::vector<double> right_wall_y_coord; //distance from the 0,0 point to the furthest point of the line segment
			std::vector<Node<3>::Pointer>left_wall_pointers; //just pointers to nodes
			std::vector<Node<3>::Pointer>right_wall_pointers; //just pointers to nodes
			std::vector<double> low_wall_x_coord; //distance from the 0,0 point to the closest point of the line segment
			std::vector<double> up_wall_x_coord; //distance from the 0,0 point to the furthest point of the line segment
			std::vector<Node<3>::Pointer>low_wall_pointers; //just pointers to nodes
			std::vector<Node<3>::Pointer>up_wall_pointers; //just pointers to nodes
			std::vector<bool> taken_vector_right; //we have to make sure that all the elements have found their matches
			std::vector<bool> taken_vector_up; //we have to make sure that all the elements have found their matches

			//first loop, creating list of nodes in each of the four walls (should be modified to work also in 3d, adding 2 more walls.
			for(unsigned int ii=0; ii<mr_model_part.Nodes().size(); ii++)
			{
				ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin() + ii;
				if(inode->X()<=x_low) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					left_wall_y_coord.push_back(inode->Y());
					left_wall_pointers.push_back(*inode.base());
					//taken_right.push_back(false);
				}
				else if(inode->X()>=x_high) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					right_wall_y_coord.push_back(inode->Y());
					right_wall_pointers.push_back(*inode.base());
					taken_vector_right.push_back(false);
				}
				if(inode->Y()<=y_low) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					low_wall_x_coord.push_back(inode->X());
					low_wall_pointers.push_back(*inode.base());
					//taken_right.push_back(false);
				}
				else if(inode->Y()>=y_high) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					up_wall_x_coord.push_back(inode->X());
					up_wall_pointers.push_back(*inode.base());
					taken_vector_up.push_back(false);
				}
				
			}
			
			
			//after looping all the elements, we should have found all the segments that are needed to create the conditions.
			//they must match EXACTLY! otherwise this will return an error.				
			if (up_wall_x_coord.size()==low_wall_x_coord.size() && right_wall_y_coord.size()==left_wall_y_coord.size()) //good! it maches!
			{
				
				Condition const& rReferenceCondition = KratosComponents<Condition>::Get("TangentVelocityPeriodicCondition2D2N");         //condition type
				Condition const& rReferenceConditionDummy = KratosComponents<Condition>::Get("Condition2D2N");         //condition type

				Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0); 		//this will allow us later to turn this layer on/off in GID

				//we loop on the left boundary, trying to find the one that matches with the cooridantes of their nodes in the right boundary:
				for(unsigned int ii=0; ii<left_wall_y_coord.size(); ii++)
				{
					const double coordinate = left_wall_y_coord[ii];
					bool found_match=false;
					for(unsigned int jj=0; jj<right_wall_y_coord.size(); jj++)
					{
						if(taken_vector_right[jj]==false) //ok, this one is still free! it could be our match
						{
							if( fabs(coordinate-right_wall_y_coord[jj])<1.0e-6 )
							{
								taken_vector_right[jj]=true;
								//adding the condition to the mdpa:
								//if (coordinate > y_low && coordinate < y_high) //then it means we are not corner nodes.
								//{
								Line2D2<Node<3> > geometry(
									left_wall_pointers[ii],  //condition to be added
									right_wall_pointers[jj]);
								
								left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =right_wall_pointers[jj]->Id();
								right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = left_wall_pointers[ii]->Id();
								left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_PRESSURE) =right_wall_pointers[jj]->Id();
								right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_PRESSURE) = left_wall_pointers[ii]->Id();

								//now we check if we can also do it for the tangential component:
								if (coordinate > y_low) // && coordinate < y_high) //then it means we are not corner nodes.
								{
									
									Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									//left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = right_wall_pointers[jj]->Id();
									//right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = left_wall_pointers[ii]->Id();
									
									/*
									Condition::Pointer p_condition = rReferenceConditionDummy.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC) = right_wall_pointers[jj]->Id();
									right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC) = left_wall_pointers[ii]->Id();
									*/ 
								}
								else
								{
									Condition::Pointer p_condition = rReferenceConditionDummy.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
								}
								

								found_match=true;							
								break;
							}	
						}
					}
					//if we got to the end of the list and still found no matching, we have a problem!
					if(found_match==false)
						KRATOS_THROW_ERROR(std::logic_error, "FINISHED LOOPING UPPER ELEMENTS AND FOUND NO MATCHING. low_coordinate = ", coordinate );
				}
				//now for the upper and lower
				for(unsigned int ii=0; ii<low_wall_x_coord.size(); ii++)
				{
					const double coordinate = low_wall_x_coord[ii];
					bool found_match=false;
					for(unsigned int jj=0; jj<up_wall_x_coord.size(); jj++)
					{
						if(taken_vector_up[jj]==false) //ok, this one is still free! it could be our match
						{
							if( fabs(coordinate-up_wall_x_coord[jj])<1.0e-6 )
							{
								taken_vector_up[jj]=true;
								//adding the condition to the mdpa:
								//if (coordinate > y_low && coordinate < y_high) //then it means we are not corner nodes.
								//{
								Line2D2<Node<3> > geometry(
									low_wall_pointers[ii],  //condition to be added
									up_wall_pointers[jj]);
								//Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
								//mr_model_part.Conditions().push_back(p_condition);
								
								low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) =up_wall_pointers[jj]->Id();
								up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = low_wall_pointers[ii]->Id();
								low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_PRESSURE) =up_wall_pointers[jj]->Id();
								up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_PRESSURE) = low_wall_pointers[ii]->Id();
								//now we check if we can also do it for the tangential component:
								if (coordinate > x_low ) //&& coordinate < x_high) //then it means we are not corner nodes.
								{
									
									Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									//low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =up_wall_pointers[jj]->Id();
									//up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = low_wall_pointers[ii]->Id();
									 
									/*
									Condition::Pointer p_condition = rReferenceConditionDummy.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT_ANTIPERIODIC) =up_wall_pointers[jj]->Id();
									up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT_ANTIPERIODIC) = low_wall_pointers[ii]->Id();
									*/ 
								}
								else
								{
									Condition::Pointer p_condition = rReferenceConditionDummy.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
								}
								
								found_match=true;							
								break;
							}	
						}
					}
					//if we got to the end of the list and still found no matching, we have a problem!
					if(found_match==false)
						KRATOS_THROW_ERROR(std::logic_error, "FINISHED LOOPING UPPER ELEMENTS AND FOUND NO MATCHING. low_coordinate = ", coordinate );
				}
			}
			else //too bad! the number does not match, the geometry is not symmetric or maybe the tolerance was wrong?
			{
				KRATOS_WATCH(up_wall_x_coord.size())
				KRATOS_WATCH(low_wall_x_coord.size())
				KRATOS_WATCH(right_wall_y_coord.size())
				KRATOS_WATCH(left_wall_y_coord.size())
				KRATOS_THROW_ERROR(std::logic_error, "DIFFERENT NUMBER OF NODES ON BOUNDARIES!!", "");
			}
			
			
			std::cout << "Finished adding periodic tangent velocity condtions with a jump " << condition_number << std::endl;

			
			KRATOS_CATCH("")
		} 
		/////
		void AddNormalAndTangentConditions2D(const double x_low ,const double x_high, const double y_low, const double y_high  ) 
		{
			KRATOS_TRY
			
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing LAGRANGE_MULTIPLIER_TANGENT_VELOCITY variable on solution step data","");

			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(LAGRANGE_MULTIPLIER_NORMAL_VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing LAGRANGE_MULTIPLIER_NORMAL_VELOCITY variable on solution step data","");


			
			unsigned int condition_number=1;
			if(mr_model_part.Conditions().size()!=0)
			{
				ModelPart::ConditionsContainerType::iterator lastcondition = (mr_model_part.ConditionsEnd()-1);
				condition_number = lastcondition->Id();
			}
			

				
			//now we have to loop the elements and create the conditions! no easy task!
			std::vector<double> left_wall_y_coord; //distance from the 0,0 point to the closest point of the line segment
			std::vector<double> right_wall_y_coord; //distance from the 0,0 point to the furthest point of the line segment
			std::vector<Node<3>::Pointer>left_wall_pointers; //just pointers to nodes
			std::vector<Node<3>::Pointer>right_wall_pointers; //just pointers to nodes
			std::vector<double> low_wall_x_coord; //distance from the 0,0 point to the closest point of the line segment
			std::vector<double> up_wall_x_coord; //distance from the 0,0 point to the furthest point of the line segment
			std::vector<Node<3>::Pointer>low_wall_pointers; //just pointers to nodes
			std::vector<Node<3>::Pointer>up_wall_pointers; //just pointers to nodes
			std::vector<bool> taken_vector_right; //we have to make sure that all the elements have found their matches
			std::vector<bool> taken_vector_up; //we have to make sure that all the elements have found their matches

			//first loop, creating list of nodes in each of the four walls (should be modified to work also in 3d, adding 2 more walls.
			for(unsigned int ii=0; ii<mr_model_part.Nodes().size(); ii++)
			{
				ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin() + ii;
				if(inode->X()<=x_low) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					left_wall_y_coord.push_back(inode->Y());
					left_wall_pointers.push_back(*inode.base());
					//taken_right.push_back(false);
				}
				else if(inode->X()>=x_high) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					right_wall_y_coord.push_back(inode->Y());
					right_wall_pointers.push_back(*inode.base());
					taken_vector_right.push_back(false);
				}
				if(inode->Y()<=y_low) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					low_wall_x_coord.push_back(inode->X());
					low_wall_pointers.push_back(*inode.base());
					//taken_right.push_back(false);
				}
				else if(inode->Y()>=y_high) //good!, we have found an element that is on the lower boundary of the slice (y=0)
				{
					up_wall_x_coord.push_back(inode->X());
					up_wall_pointers.push_back(*inode.base());
					taken_vector_up.push_back(false);
				}
				
			}
			
			
			//after looping all the elements, we should have found all the segments that are needed to create the conditions.
			//they must match EXACTLY! otherwise this will return an error.				
			if (up_wall_x_coord.size()==low_wall_x_coord.size() && right_wall_y_coord.size()==left_wall_y_coord.size()) //good! it maches!
			{
				
				Condition const& rReferenceCondition = KratosComponents<Condition>::Get("TangentVelocityPeriodicCondition2D2N");         //condition type

				Condition const& rReferenceNormalCondition = KratosComponents<Condition>::Get("NormalVelocityPeriodicCondition2D2N");         //conditiontype for conds with nodes on cube faces
				Condition const& rReferenceConditionDummy = KratosComponents<Condition>::Get("Condition2D2N");         //condition type

				Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0); 		//this will allow us later to turn this layer on/off in GID

				//we loop on the left boundary, trying to find the one that matches with the cooridantes of their nodes in the right boundary:
				for(unsigned int ii=0; ii<left_wall_y_coord.size(); ii++)
				{
					const double coordinate = left_wall_y_coord[ii];
					bool found_match=false;
					for(unsigned int jj=0; jj<right_wall_y_coord.size(); jj++)
					{
						if(taken_vector_right[jj]==false) //ok, this one is still free! it could be our match
						{
							if( fabs(coordinate-right_wall_y_coord[jj])<1.0e-6 )
							{
								taken_vector_right[jj]=true;
								//adding the condition to the mdpa:
								//if (coordinate > y_low && coordinate < y_high) //then it means we are not corner nodes.
								//{
								Line2D2<Node<3> > geometry(
									left_wall_pointers[ii],  //condition to be added
									right_wall_pointers[jj]);
								
								left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =right_wall_pointers[jj]->Id();
								right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = left_wall_pointers[ii]->Id();
								left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_PRESSURE) =right_wall_pointers[jj]->Id();
								right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_PRESSURE) = left_wall_pointers[ii]->Id();

								//now we check if we can also do it for the tangential component:
								if (coordinate > y_low) // && coordinate < y_high) //then it means we are not corner nodes.
								{
									
									//tangential velocity cond
									Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									//normal velocity cond
									Condition::Pointer p_condition_normal = rReferenceNormalCondition.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition_normal);


									//left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = right_wall_pointers[jj]->Id();
									//right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = left_wall_pointers[ii]->Id();
									
									/*
									Condition::Pointer p_condition = rReferenceConditionDummy.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC) = right_wall_pointers[jj]->Id();
									right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC) = left_wall_pointers[ii]->Id();
									*/ 
								}
								else
								{
									Condition::Pointer p_condition = rReferenceConditionDummy.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
								}
								

								found_match=true;							
								break;
							}	
						}
					}
					//if we got to the end of the list and still found no matching, we have a problem!
					if(found_match==false)
						KRATOS_THROW_ERROR(std::logic_error, "FINISHED LOOPING UPPER ELEMENTS AND FOUND NO MATCHING. low_coordinate = ", coordinate );
				}
				//now for the upper and lower
				for(unsigned int ii=0; ii<low_wall_x_coord.size(); ii++)
				{
					const double coordinate = low_wall_x_coord[ii];
					bool found_match=false;
					for(unsigned int jj=0; jj<up_wall_x_coord.size(); jj++)
					{
						if(taken_vector_up[jj]==false) //ok, this one is still free! it could be our match
						{
							if( fabs(coordinate-up_wall_x_coord[jj])<1.0e-6 )
							{
								taken_vector_up[jj]=true;
								//adding the condition to the mdpa:
								//if (coordinate > y_low && coordinate < y_high) //then it means we are not corner nodes.
								//{
								Line2D2<Node<3> > geometry(
									low_wall_pointers[ii],  //condition to be added
									up_wall_pointers[jj]);
								//Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
								//mr_model_part.Conditions().push_back(p_condition);
								
								low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) =up_wall_pointers[jj]->Id();
								up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = low_wall_pointers[ii]->Id();
								low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_PRESSURE) =up_wall_pointers[jj]->Id();
								up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_PRESSURE) = low_wall_pointers[ii]->Id();
								//now we check if we can also do it for the tangential component:
								if (coordinate > x_low ) //&& coordinate < x_high) //then it means we are not corner nodes.
								{
									//tangential velocity cond
									Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									//normal velocity cond
									Condition::Pointer p_condition_normal = rReferenceNormalCondition.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition_normal);

									//low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =up_wall_pointers[jj]->Id();
									//up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = low_wall_pointers[ii]->Id();
									 
									/*
									Condition::Pointer p_condition = rReferenceConditionDummy.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT_ANTIPERIODIC) =up_wall_pointers[jj]->Id();
									up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT_ANTIPERIODIC) = low_wall_pointers[ii]->Id();
									*/ 
								}
								else
								{
									Condition::Pointer p_condition = rReferenceConditionDummy.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
								}
								
								found_match=true;							
								break;
							}	
						}
					}
					//if we got to the end of the list and still found no matching, we have a problem!
					if(found_match==false)
						KRATOS_THROW_ERROR(std::logic_error, "FINISHED LOOPING UPPER ELEMENTS AND FOUND NO MATCHING. low_coordinate = ", coordinate );
				}
			}
			else //too bad! the number does not match, the geometry is not symmetric or maybe the tolerance was wrong?
			{
				KRATOS_WATCH(up_wall_x_coord.size())
				KRATOS_WATCH(low_wall_x_coord.size())
				KRATOS_WATCH(right_wall_y_coord.size())
				KRATOS_WATCH(left_wall_y_coord.size())
				KRATOS_THROW_ERROR(std::logic_error, "DIFFERENT NUMBER OF NODES ON BOUNDARIES!!", "");
			}
			
			
			std::cout << "Finished adding periodic tangent velocity condtions with a jump " << condition_number << std::endl;

			
			KRATOS_CATCH("")
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
    std::string Info() const override
    {
        return "AddPeriodicConditionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AddPeriodicConditionsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:
//THE CHECK IS NOW DONE IN THE CONSTRUCTOR
/*
		void Check()
		{
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODE_PAIR_X_COMPONENT) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODE_PAIR_X_COMPONENT variable on solution step data","");
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODE_PAIR_Y_COMPONENT) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODE_PAIR_Y_COMPONENT variable on solution step data","");
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODE_PAIR_PRESSURE) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODE_PAIR_PRESSURE variable on solution step data","");
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODE_PAIR_X_COMPONENT_ANTIPERIODIC) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODE_PAIR_X_COMPONENT_ANTIPERIODIC variable on solution step data","");
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODE_PAIR_Y_COMPONENT_ANTIPERIODIC variable on solution step data","");

		}

*/
	private:
		ModelPart& mr_model_part;

	};
	
	


    ///@}
    ///@name Member Variables
    ///@{
    //ModelPart& mr_fluid_model_part;
    //ModelPart& mr_structure_model_part;
    //ModelPart& mr_combined_model_part;

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
//		AddPeriodicConditionsProcess& operator=(AddPeriodicConditionsProcess const& rOther);

    /// Copy constructor.
//		AddPeriodicConditionsProcess(AddPeriodicConditionsProcess const& rOther);


    ///@}

 // Class AddPeriodicConditionsProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*
inline std::istream& operator >> (std::istream& rIStream,
                                  AddPeriodicConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AddPeriodicConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
*/
///@}


}  // namespace Kratos.

#endif // KRATOS_ADD_PERIODIC_CONDITIONS_PROCESS_INCLUDED  defined 



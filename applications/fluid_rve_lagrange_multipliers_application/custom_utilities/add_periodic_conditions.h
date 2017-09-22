

#if !defined(KRATOS_ADD_PERIODIC_CONDITIONS_NORMAL_ONLY_UTILITY_INCLUDED )
#define KRATOS_ADD_PERIODIC_CONDITIONS_NORMAL_ONLY_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// Project includes 
#include "includes/define.h"
#include "fluid_rve_lagrange_multipliers_application.h"
#include "custom_conditions/inverse_tangent_velocity_periodic_condition_2d.h"
#include "custom_conditions/inverse_normal_velocity_periodic_condition_2d.h"

#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/line_2d_2.h"



namespace Kratos
{
	//this class is to be modified by the user to customize the interpolation process
	//template< unsigned int TDim>
	class AddPeriodicConditionsNormalOnly2D
	{
	public:
	
		KRATOS_CLASS_POINTER_DEFINITION(AddMeanVelocityLagrangeMultiplierConditions2D);

		AddPeriodicConditionsNormalOnly2D(ModelPart& model_part)
			: mr_model_part(model_part) 
		{
			KRATOS_TRY	
			//std::cout << "Hello, I am the constructor of the Fixed Velocity 2d Utility" << std::endl;
			Check();		
			KRATOS_CATCH("")	
		}
		

		~AddPeriodicConditionsNormalOnly2D()
		{}

		
		void AddThem(const double x_low ,const double x_high, const double y_low, const double y_high  ) 
		{
			KRATOS_TRY
			
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
				Condition const& rReferenceCondition = KratosComponents<Condition>::Get("Condition2D2N");         //condition type
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
								Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
								mr_model_part.Conditions().push_back(p_condition);
								
								left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =right_wall_pointers[jj]->Id();
								right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = left_wall_pointers[ii]->Id();
								
								//now we check if we can also do it for the tangential component:
								if (coordinate > y_low && coordinate < y_high) //then it means we are not corner nodes.
								{
									//left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = right_wall_pointers[jj]->Id();
									//right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = left_wall_pointers[ii]->Id();
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
								Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
								mr_model_part.Conditions().push_back(p_condition);
								
								low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) =up_wall_pointers[jj]->Id();
								up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = low_wall_pointers[ii]->Id();
								
								//now we check if we can also do it for the tangential component:
								if (coordinate > x_low && coordinate < x_high) //then it means we are not corner nodes.
								{
									//low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =up_wall_pointers[jj]->Id();
									//up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = low_wall_pointers[ii]->Id();
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
			
			
			std::cout << "Finished adding periodic normal condtions" << condition_number << std::endl;

			
			KRATOS_CATCH("")
		} 
		
		void AddThemWithTangentInversePeriodicity(const double x_low ,const double x_high, const double y_low, const double y_high  ) 
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
				
				Condition const& rReferenceCondition = KratosComponents<Condition>::Get("InverseTangentVelocityPeriodicCondition2D2N");         //condition type
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
			
			
			std::cout << "Finished adding periodic normal condtions" << condition_number << std::endl;

			
			KRATOS_CATCH("")
		} 
		
		
		void AddThemWithNormalInversePeriodicity(const double x_low ,const double x_high, const double y_low, const double y_high  ) 
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
				
				Condition const& rReferenceCondition = KratosComponents<Condition>::Get("InverseNormalVelocityPeriodicCondition2D2N");         //condition type
				//Condition const& rReferenceConditionDummy = KratosComponents<Condition>::Get("Condition2D2N");         //condition type

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
								
								//left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =right_wall_pointers[jj]->Id();
								//right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = left_wall_pointers[ii]->Id();
								
								//now we check if we can also do it for the tangential component:
								if (coordinate > y_low && coordinate < y_high) //then it means we are not corner nodes.
								{
									
									Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = right_wall_pointers[jj]->Id();
									right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = left_wall_pointers[ii]->Id();
									
									/*
									Condition::Pointer p_condition = rReferenceConditionDummy.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									left_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC) = right_wall_pointers[jj]->Id();
									right_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC) = left_wall_pointers[ii]->Id();
									*/ 
								}
								else
								{
									Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
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
								
								//low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) =up_wall_pointers[jj]->Id();
								//up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_Y_COMPONENT) = low_wall_pointers[ii]->Id();
								
								//now we check if we can also do it for the tangential component:
								if (coordinate > x_low && coordinate < x_high) //then it means we are not corner nodes.
								{
									
									Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) =up_wall_pointers[jj]->Id();
									up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT) = low_wall_pointers[ii]->Id();
									 
									/*
									Condition::Pointer p_condition = rReferenceConditionDummy.Create(++condition_number, geometry, properties);
									mr_model_part.Conditions().push_back(p_condition);
									low_wall_pointers[ii]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT_ANTIPERIODIC) =up_wall_pointers[jj]->Id();
									up_wall_pointers[jj]->GetSolutionStepValue(NODE_PAIR_X_COMPONENT_ANTIPERIODIC) = low_wall_pointers[ii]->Id();
									*/ 
								}
								else
								{
									Condition::Pointer p_condition = rReferenceCondition.Create(++condition_number, geometry, properties);
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
			
			
			std::cout << "Finished adding periodic normal condtions" << condition_number << std::endl;

			
			KRATOS_CATCH("")
		} 
		
			
	protected:

		void Check()
		{
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODE_PAIR_X_COMPONENT) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODE_PAIR_X_COMPONENT variable on solution step data","");
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODE_PAIR_Y_COMPONENT) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODE_PAIR_Y_COMPONENT variable on solution step data","");
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODE_PAIR_PRESSURE) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODE_PAIR_PRESSURE variable on solution step data","");
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODE_PAIR_X_COMPONENT_ANTIPERIODIC) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODE_PAIR_X_COMPONENT_ANTIPERIODIC variable on solution step data","");
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODE_PAIR_Y_COMPONENT_ANTIPERIODIC variable on solution step data","");

		}

	private:
		ModelPart& mr_model_part;

	};
	
	
	
	

}  // namespace Kratos.

#endif // KRATOS_ADD_PERIODIC_CONDITIONS_NORMAL_ONLY_UTILITY_INCLUDED



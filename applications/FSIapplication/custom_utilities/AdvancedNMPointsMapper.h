////////////////////////////////////////////////////////// 
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:42 $
//   Revision:            $Revision: 1.1.1.1 $
//
#if !defined(KRATOS_ADVANCED_NM_POINTS_MAPPER_H_INCLUDED )
#define  KRATOS_ADVANCED_NM_POINTS_MAPPER_H_INCLUDED

//
#include <ANN/ANN.h>					// ANN declarations
#include <boost/timer.hpp> 
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "utilities/normal_calculation_utils.h" 
#include "fsi_application.h" 

//
namespace Kratos
{
//we want to enable the use of any kind of tree
//and transfer any kind of variable (e.g. pressure, displacement...)

//template<class TTree, class TVar>
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//template<class TDataType>
	class AdvancedNMPointsMapper
	{
    public:
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//typedef Variable<TDataType> Tvar
		typedef double TScalar;
		typedef array_1d<double,3> TVector;
		//cosntructor
		AdvancedNMPointsMapper(ModelPart& dest_model_part, ModelPart& origin_model_part);
		//!!!!!!!!!!!!!!!!!!!!!!
		~AdvancedNMPointsMapper();
		template<class TDataType>
		void ScalarMap(const Variable<TDataType>& rOriginVariable, const Variable<TDataType>& rDestinationVariable, int max_iter, double tol_iter)
		{
				boost::timer Mapini_timer;
				OriginScalarValue = new TScalar[n]; //if we have a template class later we should allocate
													//this memory in constructor and overwrite entries in each
													//timestep (not delete), because then we already will know if
													//we have to deal with scalar or vector values.
				for(ModelPart::NodeIterator i = DestinationModelPart.NodesBegin() ; 
				i != DestinationModelPart.NodesEnd() ; ++i)
					{
						if ((i->FastGetSolutionStepValue(IS_INTERFACE)==1.0))   
							{
								i->FastGetSolutionStepValue(NODAL_AREA)=0.0;
							}
					}
				double area;
				for(int i=0;i<n;i+=3)
				{
					area =DestArea[i];
					DestinationConds[i]->GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA)+=0.333333333333333*area;
					//we set destinationVarible to zero, remove this later
					DestinationConds[i]->GetGeometry()[0].FastGetSolutionStepValue(rDestinationVariable)=0.0;
					
					DestinationConds[i]->GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA)+=0.333333333333333*area;
					DestinationConds[i]->GetGeometry()[1].FastGetSolutionStepValue(rDestinationVariable)=0.0;

					DestinationConds[i]->GetGeometry()[2].FastGetSolutionStepValue(NODAL_AREA)+=0.333333333333333*area;
					DestinationConds[i]->GetGeometry()[2].FastGetSolutionStepValue(rDestinationVariable)=0.0;

				}
				for(int i=0;i<n;i++)
				{
					if(check[i]!=0)
					{
						OriginScalarValue[i] =	OriginConds[i]->GetGeometry()[0].FastGetSolutionStepValue(rOriginVariable)*OriginProjectionPoints[i][3]
											+ OriginConds[i]->GetGeometry()[1].FastGetSolutionStepValue(rOriginVariable)*OriginProjectionPoints[i][4]
											+ OriginConds[i]->GetGeometry()[2].FastGetSolutionStepValue(rOriginVariable)*(1.0-OriginProjectionPoints[i][3]-OriginProjectionPoints[i][4]);
						//					std::cout<<i<<"  "<<OriginScalarValue[i]<<std::endl;
					}
					else {OriginScalarValue[i] =0.0;}
				}
				
				//later, if we have a template class, make this outside of function, because we don´t want to do t
				//create the arrays in each timestep
				array_1d<double,3> loc_RHS;
				array_1d<double,3> x_vect;
				array_1d<double,3> p_orig;
				Matrix M_cons(3,3);
							M_cons(0,0)=2.0;	M_cons(0,1)=1.0;	M_cons(0,2)=1.0;
							M_cons(1,0)=1.0;	M_cons(1,1)=2.0;	M_cons(1,2)=1.0;
							M_cons(2,0)=1.0;	M_cons(2,1)=1.0;	M_cons(2,2)=2.0;

				double norm_dx=0.0;
				int inodes = 0;
				std::cout << "  Map Initialization TIME = " << Mapini_timer.elapsed() << std::endl;
				///////////////////////////////////////////////////////////////////////////////////////
				boost::timer Mapiter_timer;
				double norm_x = 0;
				//std::cout<<"Start Iterations"<<std::endl;
				int iter = 0;
				for(int k=0;k<max_iter;k++) //Iteration
				{
					iter += 1;
					for(int i=0;i<n;i+=3) //enough for all 3 steps
					{
							DestinationConds[i]->GetGeometry()[0].GetValue(AUX)=0.0;
							DestinationConds[i]->GetGeometry()[1].GetValue(AUX)=0.0;
							DestinationConds[i]->GetGeometry()[2].GetValue(AUX)=0.0;
					}			
					for(int i=0;i<n;i+=3)	
					{
							loc_RHS[0]=0.0; loc_RHS[1]=0.0; loc_RHS[2]=0.0;
						
							//gather
							x_vect[0] = DestinationConds[i]->GetGeometry()[0].FastGetSolutionStepValue(rDestinationVariable);
							x_vect[1] = DestinationConds[i]->GetGeometry()[1].FastGetSolutionStepValue(rDestinationVariable);
							x_vect[2] = DestinationConds[i]->GetGeometry()[2].FastGetSolutionStepValue(rDestinationVariable);
							
							//make this better with matrix operation
							//bei vererbten klassen, eine funktion zur berechnung von p_orig schreiben, so dass VectorMap gleich bleibt
							p_orig[0]=DestArea[i]/24.0*(6.0*OriginScalarValue[i]+1.0*OriginScalarValue[i+1]+1.0*OriginScalarValue[i+2]); 
							p_orig[1]=DestArea[i]/24.0*(1.0*OriginScalarValue[i]+6.0*OriginScalarValue[i+1]+1.0*OriginScalarValue[i+2]);
							p_orig[2]=DestArea[i]/24.0*(1.0*OriginScalarValue[i]+1.0*OriginScalarValue[i+1]+6.0*OriginScalarValue[i+2]);
							

							//M_cons*=DestArea[i]/(12.0);
							
							loc_RHS = p_orig - DestArea[i]/(12.0)*prod(M_cons, x_vect);

							DestinationConds[i]->GetGeometry()[0].GetValue(AUX)+=loc_RHS[0];
							DestinationConds[i]->GetGeometry()[1].GetValue(AUX)+=loc_RHS[1];
							DestinationConds[i]->GetGeometry()[2].GetValue(AUX)+=loc_RHS[2];
					}
					//////////////////////////////////////////////////////////////
					norm_dx=0.0;
					inodes=0;
					for(ModelPart::NodeIterator id = DestinationModelPart.NodesBegin() ; 
					id != DestinationModelPart.NodesEnd() ; ++id)
					{
						if ((id->FastGetSolutionStepValue(IS_INTERFACE)==1.0))   
						{
							const double& A_I = id->FastGetSolutionStepValue(NODAL_AREA);
							id->FastGetSolutionStepValue(rDestinationVariable)+=id->GetValue(AUX)/A_I;//p+=dp
							norm_dx+=(id->GetValue(AUX)/A_I)*(id->GetValue(AUX)/A_I);
							inodes++;
						}
					}

					
					if(iter == 1)
					{
						norm_x = 0.00;
						double temp;
						for(ModelPart::NodeIterator id = DestinationModelPart.NodesBegin() ; 
						id != DestinationModelPart.NodesEnd() ; ++id)
						{
							if ((id->FastGetSolutionStepValue(IS_INTERFACE)==1.0))   
							{
								temp = id->FastGetSolutionStepValue(rDestinationVariable);
								norm_x += temp*temp;
							}
						}
					}
					
					double ratio = 0.00;
					if (norm_x > 1e-15)
						ratio = sqrt(norm_dx/norm_x);
					std::cout << "ratio = " << ratio << std::endl;
					
					if ((norm_dx/inodes)<=pow(tol_iter*0.001,2) ||  ratio <= tol_iter) break;
						
				///////////////////////////////////////////////////////////////
				}//Iteration
				delete[] OriginScalarValue;

				std::cout << "  Map Iteration TIME for "<<iter<<" Iterations = " << Mapiter_timer.elapsed() << std::endl;
				
	


		}//ScalarMap
		
		template<class TDataType>
		void VectorMap(const Variable<TDataType>& rOriginVariable, const Variable<TDataType>& rDestinationVariable, int max_iter, double tol_iter)
		{
				boost::timer Mapini_timer;
				//in template class this will be defined in constuctor an overwritten in eacg timestep. Deletion is not necessary then
				OriginVectorValue = new TVector[n];
				TVector zero_vec; zero_vec[0]=0.0; zero_vec[1]=0.0; zero_vec[2]=0.0; 
				
				for(ModelPart::NodeIterator i = DestinationModelPart.NodesBegin() ; 
				i != DestinationModelPart.NodesEnd() ; ++i)
					{
						if (i->FastGetSolutionStepValue(IS_INTERFACE)==1.0)   
						{
							i->FastGetSolutionStepValue(NODAL_AREA)=0.0;
						}
					}
				////////////////////////////////////////////////////////////////
				//std::cout<<"Calculate Nodal Area for DestModelPart"<<std::endl;
				//double area;
				for(int i=0;i<n;i+=3)
					{
							//substitute with area!!!
							DestinationConds[i]->GetGeometry()[0].FastGetSolutionStepValue(NODAL_AREA)+=0.333333333333333*DestArea[i];
							DestinationConds[i]->GetGeometry()[0].FastGetSolutionStepValue(rDestinationVariable)=zero_vec;

							DestinationConds[i]->GetGeometry()[1].FastGetSolutionStepValue(NODAL_AREA)+=0.333333333333333*DestArea[i];
							DestinationConds[i]->GetGeometry()[1].FastGetSolutionStepValue(rDestinationVariable)=zero_vec;

							DestinationConds[i]->GetGeometry()[2].FastGetSolutionStepValue(NODAL_AREA)+=0.333333333333333*DestArea[i];
							DestinationConds[i]->GetGeometry()[2].FastGetSolutionStepValue(rDestinationVariable)=zero_vec;
					}
				for(int i=0;i<n;i++)
				{			
					//this is to filter out the Gauss point which do not have projectin point
					if(check[i]!=0)
					{	
							noalias( OriginVectorValue[i] ) = OriginConds[i]->GetGeometry()[0].FastGetSolutionStepValue(rOriginVariable)*OriginProjectionPoints[i][3]
											+ OriginConds[i]->GetGeometry()[1].FastGetSolutionStepValue(rOriginVariable)*OriginProjectionPoints[i][4]
											+ OriginConds[i]->GetGeometry()[2].FastGetSolutionStepValue(rOriginVariable)*(1.0-OriginProjectionPoints[i][3]-OriginProjectionPoints[i][4]);
					}
					else { noalias( OriginVectorValue[i] ) = zero_vec;}
				}
				//array_1d<double,3> loc_RHS_0; array_1d<double,3> loc_RHS_1; array_1d<double,3> loc_RHS_2;
				array_1d<double,3> x_vect_0; array_1d<double,3> x_vect_1; array_1d<double,3> x_vect_2;
				array_1d<double,3> p_orig_0; array_1d<double,3> p_orig_1; array_1d<double,3> p_orig_2;
				Matrix M_id(3,3);
							M_id(0,0)=1.0; M_id(0,1)=0.0; M_id(0,2)=0.0;
							M_id(1,0)=0.0; M_id(1,1)=1.0; M_id(1,2)=0.0; 
							M_id(2,0)=0.0; M_id(2,1)=0.0; M_id(2,2)=1.0; 

				double norm_dx=0.0;
				int inodes = 0;
				std::cout << "  Map Initialization TIME = " << Mapini_timer.elapsed() << std::endl;
				///////////////////////////////////////////////////////////////////////////////////////
				boost::timer Mapiter_timer;
				//std::cout<<"Start Iterations"<<std::endl;
				int iter = 0;
				for(int k=0;k<max_iter;k++) //Iteration
				{
					iter += 1;
					for(int i=0;i<n;i+=3)
						{
							//take zerovector from ublas
							noalias( DestinationConds[i]->GetGeometry()[0].GetValue(VAUX) ) = zero_vec; //(*ic).GetGeometry()[0].GetValue(AUX_Y)=0;(*ic).GetGeometry()[0].GetValue(AUXZ)=0;
							noalias( DestinationConds[i]->GetGeometry()[1].GetValue(VAUX) ) = zero_vec; //(*ic).GetGeometry()[1].GetValue(AUXY)=0;(*ic).GetGeometry()[1].GetValue(AUXZ)=0;
							noalias( DestinationConds[i]->GetGeometry()[2].GetValue(VAUX) ) = zero_vec; //(*ic).GetGeometry()[2].GetValue(AUXY)=0;(*ic).GetGeometry()[2].GetValue(AUXZ)=0;
						}
					for(int i=0;i<n;i+=3)	
						{
							//noalias( loc_RHS_0 ) = zero_vec;
							//noalias( loc_RHS_1 ) = zero_vec;
							//noalias( loc_RHS_2 ) = zero_vec;
							
							//gather
							noalias( x_vect_0 ) = DestinationConds[i]->GetGeometry()[0].FastGetSolutionStepValue(rDestinationVariable);
							noalias( x_vect_1 ) = DestinationConds[i]->GetGeometry()[1].FastGetSolutionStepValue(rDestinationVariable);
							noalias( x_vect_2 ) = DestinationConds[i]->GetGeometry()[2].FastGetSolutionStepValue(rDestinationVariable);
							

							p_orig_0=DestArea[i]/(24.0)*(prod(6.0*M_id, OriginVectorValue[i]) + prod(1.0*M_id, OriginVectorValue[i+1]) + prod(1.0*M_id, OriginVectorValue[i+2]) ); 
							p_orig_1=DestArea[i]/(24.0)*(prod(1.0*M_id, OriginVectorValue[i]) + prod(6.0*M_id, OriginVectorValue[i+1]) + prod(1.0*M_id, OriginVectorValue[i+2]) ); 
							p_orig_2=DestArea[i]/(24.0)*(prod(1.0*M_id, OriginVectorValue[i]) + prod(1.0*M_id, OriginVectorValue[i+1]) + prod(6.0*M_id, OriginVectorValue[i+2]) ); 

														
							//loc_RHS_0 = p_orig_0 - DestArea[i]/(12.0)*(prod(2.0*M_id, x_vect_0) + prod(1.0*M_id, x_vect_1) + prod(1.0*M_id, x_vect_2) ); 
							//loc_RHS_1 = p_orig_1 - DestArea[i]/(12.0)*(prod(1.0*M_id, x_vect_0) + prod(2.0*M_id, x_vect_1) + prod(1.0*M_id, x_vect_2) ); 
							//loc_RHS_2 = p_orig_2 - DestArea[i]/(12.0)*(prod(1.0*M_id, x_vect_0) + prod(1.0*M_id, x_vect_1) + prod(2.0*M_id, x_vect_2) );  

							
							noalias( DestinationConds[i]->GetGeometry()[0].GetValue(VAUX) ) += //loc_RHS_0;
								p_orig_0 - DestArea[i]/(12.0)*(prod(2.0*M_id, x_vect_0) + prod(1.0*M_id, x_vect_1) + prod(1.0*M_id, x_vect_2) ); 
							noalias( DestinationConds[i]->GetGeometry()[1].GetValue(VAUX) ) += //loc_RHS_1;
								p_orig_1 - DestArea[i]/(12.0)*(prod(1.0*M_id, x_vect_0) + prod(2.0*M_id, x_vect_1) + prod(1.0*M_id, x_vect_2) ); 
							noalias( DestinationConds[i]->GetGeometry()[2].GetValue(VAUX) ) += //loc_RHS_2;
								p_orig_2 - DestArea[i]/(12.0)*(prod(1.0*M_id, x_vect_0) + prod(1.0*M_id, x_vect_1) + prod(2.0*M_id, x_vect_2) );  
							
						}
					//////////////////////////////////////////////////////////////
					norm_dx=0.0;
					inodes=0;
					for(ModelPart::NodeIterator id = DestinationModelPart.NodesBegin() ; 
					id != DestinationModelPart.NodesEnd() ; ++id)
					{
						if ((id->FastGetSolutionStepValue(IS_INTERFACE)==1.0))   
						{
							const double& A_I = id->FastGetSolutionStepValue(NODAL_AREA);
							noalias( id->FastGetSolutionStepValue(rDestinationVariable) )+=id->GetValue(VAUX)/A_I;//p+=dp
							norm_dx+=pow((id->GetValue(VAUX_X)/A_I),2)+pow((id->GetValue(VAUX_Y)/A_I),2)+pow((id->GetValue(VAUX_Z)/A_I),2);
							inodes++;
											
						}
					}
					//std::cout<<norm_dx/inodes<<std::endl;
					if((norm_dx/inodes)<=pow(tol_iter,2))break;
				////////////////////////////////////////////////////////////////
				}//Iteration
				//different in template class
				delete[] OriginVectorValue;

				std::cout << "  Map Iteration TIME for "<<iter<<" Iterations = " << Mapiter_timer.elapsed() << std::endl;

		}//VectorMap
		
		//in template class template<TVar> makeMap and use Map<TScalar> or Map<TScalar> 
		//void Map(double, int, double);
		//void setDestGaussPoints();
			
		//it takes all neighbours of the baricenter of the origin element i, and finds all the dest. elements, whose
		//Gauss points lies within a radius alpha*r (r is the radius of the sphere, around the origin triangle i
		//with the center in the middle, and alpha is the tolerance parameter, alpha>1)
		//the search function returns an ARRAY of the indeces of the Gauss points of the obtained NEIGHBOURS (dest)
		//neighbours are Gauss points....
		//void FindNeighbours(ModelPart&);
		void FindNeighbours(double);
		
		//It takes NEIGHBOUR g and checks if the projection onto the respective origin element O from 
		//the g lies inside the origin element O.
		//if it lies inside then we check the entry of the CHECK array (if the Gauss point G already
		//has a projection point - if yes we do the distance comparison and take the shorter distance
		//if no - apply the fct. setOriginProjectionPoints (write the respective entries (projection
		//coordinates) in the array, and the shape functions values at the projection point (and also write 1 in the
        //CHECK array )
		//if it does not lie - do nothing
		
				
		void setOriginProjectionPoints(int, Condition&);
		//void setOriginProjectionPoints(int, Condition, ModelPart::ConditionsContainerType::iterator);
		
		//this is our old Map function for pressure, we can remove it now
		
		
		void checkBadGP();
		
		void setPPCoordinates(array_1d<double,3>, int,  Condition&);
		void setlocalPPCoordinates(array_1d<double,3>, int,  Condition&);
		void setlocalPPCoordinates(double, double, int,  Condition&);
		void setglobalPPCoordinates(double, double, int,  Condition&); //function overloading
		
		void calcLocalPPCoords(array_1d<double,3>&, int, Condition&);
		

		void GetMidPoint(Condition& cond, double*);
		double GetSquaredRadius(Condition& cond, double*);
		double GetSquaredDistance(double*, double*);
		double GetArea(Condition& cond);
		
		/*void operator()()
		{
			Execute();
		}


		virtual void Execute()
		{
			KRATOS_TRY;
		
			
			KRATOS_CATCH("")
		}*/


	private:
		//number of interface conditions
		int n;
		//array that stores the coordinates of GaussPoints of destination
		double** DestGaussPoints;
		//and the respective condition area
		double* DestArea;
		//this array is invented to check that a normal from the Gauss point of the 
		//destination conditions, crosses only one origin condition
		//when the entry of CHECK is more than one - it means an ambiguous case where 
		//distances have to be compared and only one has to be chosen for projection
		int* check;
		//this array is intended to store the coordinates of the origin points, and 
		//the values of the shape functions at these points which will be needed for 
		//projection (xi, eta, zeta)
		double** OriginProjectionPoints;
		
		double** DestCondsNormals;
		Condition** OriginConds;
		Condition** DestinationConds;
		
		//Maybe we can do this with class template like
		//template<class TVar>
		//create object MNPointsMapper<TScalar>(arguments...)
		//TVar* OriginVarValue;
		TScalar* OriginScalarValue;
		TVector* OriginVectorValue;

		//ModelPart* OriginModelPart;
		ModelPart& OriginModelPart;
		ModelPart& DestinationModelPart;
	};

}	 // namespace Kratos.

#endif // KRATOS_NM_POINTS_MAPPER_H_INCLUDED  defined 

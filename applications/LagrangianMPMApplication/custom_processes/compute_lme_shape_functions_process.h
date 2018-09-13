//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Zhiming Guo
//                   Riccardo Rossi
//





#if !defined(COMPUTE_LME_SHAPE_FUNCTIONS_PROCESS)
#define COMPUTE_LME_SHAPE_FUNCTIONS_PROCESS



/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/variables.h"


#include "lagrangian_mpm_application_variables.h"
#include "utilities/math_utils.h"
#include "custom_utilities/compute_mls_shape_functions_utility.h"
#include "custom_utilities/compute_lme_shape_functions_utility.h"
#include <cmath>
#include <algorithm>

#define PI 3.14159265

namespace Kratos
{



class ComputeLMEShapeFunctionsProcess : public Process
{

public:

    ComputeLMEShapeFunctionsProcess( ModelPart& model_part ) : Process(), mr_model_part( model_part )
    {
        /////////////

    }

    ~ComputeLMEShapeFunctionsProcess() {}



    virtual void Execute()

    {

        KRATOS_TRY;

        //const unsigned int dim = mr_model_part.GetProcessInfo()[DOMAIN_SIZE];
        const unsigned int dim = 2;
		
		
		//I am making this search to define the nodal spacing for each node
		
		typedef Node<3>                                         NodeType;
		typedef Node<3>::Pointer                                NodePointerType;
		typedef ModelPart::NodesContainerType::ContainerType    NodesContainerType;
		
		

		BinsDynamic<3, NodeType, NodesContainerType> bins(mr_model_part.NodesArray().begin(),mr_model_part.NodesArray().end());

		Node<3> aux_node(0,0.0,0.0,0.0);
		for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode!=mr_model_part.NodesEnd(); inode++)
        {
		
			double searching_radius = inode->FastGetSolutionStepValue(NODAL_RADIUS);
		
			unsigned int MaximumNumberOfResults = 1000;
            std::vector<NodePointerType>                            results(MaximumNumberOfResults);
            std::vector<double>                                     distances(MaximumNumberOfResults);
		
			aux_node.Coordinates() = inode->Coordinates();
			unsigned int nresults = bins.SearchInRadius(aux_node,searching_radius,results.begin(),distances.begin(),MaximumNumberOfResults);
			
			
			//if(inode->Id() == 817)
			//{
				//std::cout<<" nresults "<< nresults<<std::endl;
			//}
			int icount = 0;
			while(nresults<6 )//|| nresults>11)
			{
				//if(nresults<7)
				//{
				icount++;
				searching_radius *= 1.001; 
				nresults = bins.SearchInRadius(aux_node,searching_radius,results.begin(),distances.begin(),MaximumNumberOfResults);
				//}
				//else if(nresults>11)
				//{
				//searching_radius *= 0.9;
				
				//nresults = bins.SearchInRadius(aux_node,searching_radius,results.begin(),distances.begin(),MaximumNumberOfResults);
					
				//}
				if(icount>1000)
				{
					
					break;
				}
				
			}
			
			
			double distance = 0.0;
			unsigned int counter=0;
            
            for(NodesContainerType::iterator it = results.begin(); it != results.end(); ++it)
            {
                if(counter< nresults)
                {
                    //std::cout<<" test "<< (*(it.base()))->Id() <<std::endl;
                    //std::cout<<" test "<< (*(it.base()))->Coordinates() <<std::endl;
                    
                    double dist = 0.0;
                    for(unsigned int d = 0; d < dim; d++)
                    {
						dist += (inode->Coordinates()(d)-(*(it.base()))->Coordinates()(d)) * (inode->Coordinates()(d)-(*(it.base()))->Coordinates()(d));
					}
					distance += sqrt(dist);
					
                }
                 
                counter += 1;
               
            }
            //AVERAGE NODAL DISTANCE
            distance /= (nresults-1);
            
			
			inode->FastGetSolutionStepValue(NODAL_SPACING) = distance;
		
		}

        for(ModelPart::ElementsContainerType::iterator iel = mr_model_part.ElementsBegin();
                iel!=mr_model_part.ElementsEnd(); iel++)
        {
			
            //calculating shape functions values
            Geometry< Node<3> >& geom = iel->GetGeometry();
			
			
            Vector& Ng = iel->GetValue(SHAPE_FUNCTIONS);
            Matrix& DN_DX = iel->GetValue(SHAPE_FUNCTIONS_DERIVATIVES);
            const array_1d<double,3>& xg = iel->GetValue(GAUSS_POINT_COORDINATES);
            
            Vector h;
            h = zero_vector<double> (geom.size());
            //double h = iel->GetValue(EFFECTIVE_RADIUS);

            Ng = zero_vector<double> (geom.size());
            DN_DX = ZeroMatrix(geom.size(), dim);

            //KRATOS_WATCH(this->GetValue(TEMP_POS));
            Matrix node_coordinates(geom.size(),dim);
            for ( unsigned int i = 0; i < geom.size(); i++ )
            {
                for(unsigned int k=0; k<node_coordinates.size2(); k++)
                {
                    node_coordinates(i,k) = geom[i].Coordinates()[k];
                }
                h(i) = geom[i].FastGetSolutionStepValue(NODAL_SPACING);
            }

           
            int id_el = iel->Id();
            bool iflag = 0;
            
            
           
            
            LinearULFLME::ComputeLMEShapef(Ng, DN_DX, node_coordinates, xg, h,id_el,iflag);
            
            
			
			//IT IS POSSIBLE THAT ITERATIVE PROCEDURE DOES NOT WORK
			//IN THIS LOOP I START INCREASING THE SIZE OF CONNECTIVITY
            while(iflag == 1)
            {
				
				//std::cout<<"HERE "<<std::endl;
				unsigned int MaximumNumberOfResults = 1000;
				std::vector<NodePointerType>                            results(MaximumNumberOfResults);
				std::vector<double>                                     distances(MaximumNumberOfResults);

				//calculating shape functions values
				Geometry< Node<3> >::Pointer pgeom = iel->pGetGeometry();
				
				//clear the connectivity
				pgeom->clear();
				double search_radius = iel->GetValue(SEARCH_RADIUS);
				
				aux_node.Coordinates() =  xg;
				unsigned int nresults = bins.SearchInRadius(aux_node,search_radius,results.begin(),distances.begin(),MaximumNumberOfResults);
				
				 	
					
				while(nresults <= h.size())
					{
						search_radius *= 1.0001;
						//std::cout<<"search_radius old "<<search_radius<<std::endl;
						nresults = bins.SearchInRadius(aux_node,search_radius,results.begin(),distances.begin(),MaximumNumberOfResults);
						//std::cout<<"nresults new "<<nresults<<std::endl;
					}
					
					unsigned int counter=0;
            
					for(NodesContainerType::iterator it = results.begin(); it != results.end(); ++it)
					{
						if(counter< nresults)
						{
							pgeom->push_back( *(it.base()));
							
							
						}
					 
						counter += 1;
				   
					}
				
			
            iflag = 0;
            Geometry< Node<3> >& new_geom = iel->GetGeometry();
            Matrix new_node_coordinates(new_geom.size(),dim);
            
            
            h = zero_vector<double> (new_geom.size());
            

            Ng = zero_vector<double> (new_geom.size());
            DN_DX = ZeroMatrix(new_geom.size(), dim);
            
            
            for ( unsigned int i = 0; i < new_geom.size(); i++ )
            {
                for(unsigned int k=0; k<new_node_coordinates.size2(); k++)
                {
                    new_node_coordinates(i,k) = new_geom[i].Coordinates()[k];
                }
				h(i) = new_geom[i].FastGetSolutionStepValue(NODAL_SPACING);
            }
            
            
           
            //I compute the shape functions with a new connectivity
            LinearULFLME::ComputeLMEShapef(Ng, DN_DX, new_node_coordinates, xg, h,id_el,iflag);
            
            
            
            
            if(new_geom.size() > 30)
            {
				std::cout<<"CONVERGENCE NOT REACHED in 50 iterations"<<std::endl;
				break;
            }
//********************************************************************************            
            
            //If the iterative procedure does not work again I choose those nodes
            //of the connectivity which distance vectors do not create an angle 
            //less than 15 degree. In this way I choose only well located nodes.
            
            if(iflag != 0 && new_geom.size() > 5)
            {
				
				Vector v1 = zero_vector<double> (dim);
				Vector v2 = zero_vector<double> (dim);
				
				for(unsigned int i = 0; i<new_geom.size();i++)
				{
					
					for(unsigned int j = 0; j<new_geom.size(); j++)
					{
							
						if(i != j)
						{
							
							for(unsigned int k = 0; k<dim;k++)
							{
								v1[k] = new_geom[i].Coordinates()[k] - xg[k];
								v2[k] = new_geom[j].Coordinates()[k] - xg[k];
								
							}
							
							
					
							
							double Norm_v1 = MathUtils<double>::Norm(v1);
							double Norm_v2 = MathUtils<double>::Norm(v2);
							
							v1 /= Norm_v1;
							v2 /= Norm_v2;
							
							
							
							
							double cos_alfa = v1[0]*v2[0] + v1[1]*v2[1];
							
							if(acos(cos_alfa)*180/PI < 15)
							{
								if(Norm_v1 > Norm_v2)
								{
									new_geom[i].Set(MARKER, true);
									
								}
								else
								{
									new_geom[j].Set(MARKER, true);
								}
								
							} 
						}
					}
					
				}
				
				pgeom->clear();
				
				unsigned int counter=0;	
				
				for(NodesContainerType::iterator it = results.begin(); it != results.end(); ++it)
				{
					
					Node<3>::Pointer pnode = *(it.base());
					
					
					
					if(counter< nresults)
					{
						
					if(pnode->Is(MARKER) == 0) 
					{
						
						
						pgeom->push_back( pnode);
						
						
					}
					
					
					pnode->Set(MARKER, false);
					
					}
					counter += 1;
				}
				//std::cout<<"here check angle 4"<<std::endl;	
				iflag = 0;
				Geometry< Node<3> >& modified_geom = iel->GetGeometry();
				Matrix modified_node_coordinates(modified_geom.size(),dim);
				
				Vector h_modified = zero_vector<double> (modified_geom.size());
				

				Ng = zero_vector<double> (modified_geom.size());
				DN_DX = ZeroMatrix(modified_geom.size(), dim);
            
				
				for ( unsigned int i = 0; i < modified_geom.size(); i++ )
				{
					
					
					for(unsigned int k=0; k<modified_node_coordinates.size2(); k++)
					{
						modified_node_coordinates(i,k) = modified_geom[i].Coordinates()[k];
					}
					h_modified(i) = modified_geom[i].FastGetSolutionStepValue(NODAL_SPACING);
				}
            
				
				LinearULFLME::ComputeLMEShapef(Ng, DN_DX, modified_node_coordinates, xg, h_modified,id_el,iflag);
				
				//if(iflag == 0)
				//{
					//std::cout<<"the algorithm is verified"<<std::endl;
					//std::cout<<"ID ELEMENT "<<id_el<<std::endl;
				//}
				
				
			}
//*************************************************************************************************************            
            
           
            
            }
            
            
            
            
            double norm_N;
            
            norm_N = MathUtils<double>::Norm(Ng);
            
            
			//If the max-ent algorithm does not have any solution I eliminate the 
			// material point under considered
            if ( norm_N == 0.0 )//(std::isnan(norm_N))// &
            {			
				std::cout<<"ID ELEMENT"<<iel->Id()<<std::endl;
				std::cout<<"xg ELEMENT"<<xg<<std::endl;
				std::cout<<"node_coordinates"<<node_coordinates<<std::endl;
				//std::cout<<"new_node_coordinates"<<new_node_coordinates<<std::endl;
				std::cout<<"connectivity size"<<Ng.size()<<std::endl;
				
				std::cout<<"Ng"<<Ng<<std::endl;
				std::cout<<"iflag"<<iflag<<std::endl;
				//break;
				
				
				iel->Set(TO_ERASE, true);
      	
			}
			
            
        }
        mr_model_part.RemoveElements( TO_ERASE );
//********************code for lagrange multipliers*****************************************************************************************************

        //for(ModelPart::ConditionsContainerType::iterator icon = mr_model_part.ConditionsBegin();
                //icon!=mr_model_part.ConditionsEnd(); icon++)
        //{


            ////calculating shape functions values
            //Geometry< Node<3> >& geom = icon->GetGeometry();

            //Vector& Nc = icon->GetValue(SHAPE_FUNCTIONS);
            //Matrix& DN_DXc = icon->GetValue(SHAPE_FUNCTIONS_DERIVATIVES);

            //const double h = icon->GetValue(EFFECTIVE_RADIUS);




            //Nc = zero_vector<double> (geom.size());
            //DN_DXc = ZeroMatrix(geom.size(), dim);

            ////KRATOS_WATCH(this->GetValue(TEMP_POS));
            //Matrix node_coordinates(geom.size(),dim);
            //for ( unsigned int i = 0; i < geom.size(); i++ )
            //{
                //for(unsigned int k=0; k<node_coordinates.size2(); k++)
                //{
                    //node_coordinates(i,k) = geom[i].Coordinates()[k];
                    ////node_coordinates(i,k) = geom[i].GetInitialPosition()[k];
                //}
            //}
        ////KRATOS_WATCH("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb");
            ////const array_1d<double,dim> xc = row(node_coordinates,0);
            //const array_1d<double,dim> xc = icon->GetGeometry()[0].Coordinates();
        ////KRATOS_WATCH("rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr");
            ////Vector Ng( number_of_nodes);
            ////Matrix DN_DX( number_of_nodes, 2 );
            ////LinearMLSKernel::ComputeMLSKernel(Nc, DN_DXc, node_coordinates, xc, h);
            //LinearULFLME::ComputeLMEShapef(Nc, DN_DXc, node_coordinates, xc, h);
            ////KRATOS_WATCH(xc);

            ////KRATOS_WATCH(node_coordinates);

        //}
//*************************************************************************************************************
        KRATOS_CATCH("");
    }

    



private:

    ModelPart& mr_model_part;


};


}

#endif //COMPUTE_SHAPE_FUNCTIONS_PROCESS


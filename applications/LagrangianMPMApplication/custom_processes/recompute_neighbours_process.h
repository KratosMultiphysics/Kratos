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





#if !defined(RECOMPUTE_NEIGHBOURS_PROCESS)
#define RECOMPUTE_NEIGHBOURS_PROCESS



/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"

#include "lagrangian_mpm_application_variables.h"
#include "utilities/math_utils.h"
#include <cmath>
#include <algorithm>



namespace Kratos
{



class RecomputeNeighboursProcess : public Process
{

public:


    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;



    RecomputeNeighboursProcess( ModelPart& model_part ) : Process(), mr_model_part( model_part )
    {
        /////////////

    }

    ~RecomputeNeighboursProcess() {}

	virtual void EliminateNodes()
	
	{
		
		KRATOS_TRY;
		
		for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode!=mr_model_part.NodesEnd(); inode++)
        {
			//if(inode->Id() == 211)
			//{
				//std::cout<<" IS STRUCTURE "<< inode->FastGetSolutionStepValue(IS_STRUCTURE)<<std::endl;
				//std::cout<<" IS BOUNDARY "<< inode->FastGetSolutionStepValue(IS_BOUNDARY)<<std::endl;
			//}
			if( (inode->Y() < 0.0 && inode->GetSolutionStepValue(IS_BOUNDARY) != 1))//||(inode->Y() > 0.0 && inode->Y() < 0.017 && inode->X() > 2.117 && inode->GetSolutionStepValue(IS_BOUNDARY) != 1) || (inode->GetSolutionStepValue(IS_STRUCTURE) != 1 && inode->GetSolutionStepValue(IS_BOUNDARY) != 1))
			{
				std::cout<<" ID NODE "<< inode->Id()<<std::endl;
				
				
				inode->Set(TO_ERASE, true);
				
				
			}
			inode->FastGetSolutionStepValue(IS_STRUCTURE) = 0;
		}
		mr_model_part.RemoveNodesFromAllLevels(TO_ERASE);
		
		KRATOS_CATCH("");
		
	}

    virtual void Execute()

    {

        KRATOS_TRY;


        typedef Node<3>                                         NodeType;
        typedef Node<3>::Pointer                                NodePointerType;
        typedef ModelPart::NodesContainerType::ContainerType    NodesContainerType;
        ModelPart::ElementsContainerType                        temp_container;
        ModelPart::ConditionsContainerType                      temp_condition_container;



        BinsDynamic<3, NodeType, NodesContainerType> bins(mr_model_part.NodesArray().begin(),mr_model_part.NodesArray().end());

        //creating an auxiliary node for searching purposes
        Node<3> aux_node(0,0.0,0.0,0.0);

		//set INACTIVE all the nodes
		//for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                //inode!=mr_model_part.NodesEnd(); inode++)
        //{
			
			//inode->Reset(ACTIVE);
			
		//}



        for(ModelPart::ElementsContainerType::iterator iel = mr_model_part.ElementsBegin();
                iel!=mr_model_part.ElementsEnd(); iel++)
        {
            unsigned int MaximumNumberOfResults = 1000;
            std::vector<NodePointerType>                            results(MaximumNumberOfResults);
            std::vector<double>                                     distances(MaximumNumberOfResults);

            //calculating shape functions values
            Geometry< Node<3> >::Pointer pgeom = iel->pGetGeometry();
            pgeom->clear();

            double search_radius = iel->GetValue(SEARCH_RADIUS);//0.02;



            //get the coordinates of the gauss point as  the sum of N[i]*xnode
            const array_1d<double,3>& coordinates_of_gauss = iel->GetValue(GAUSS_POINT_COORDINATES);

            aux_node.Coordinates() =  coordinates_of_gauss;
            unsigned int nresults = bins.SearchInRadius(aux_node,search_radius,results.begin(),distances.begin(),MaximumNumberOfResults);
			


//If the moving least square method is used is better to comment this part
//**********************************************************************************************************************************
			int icount = 0;

			while(nresults != 3)
			{
				icount++;
				if(nresults<3)
				{
					search_radius *= 1.001;
					nresults = bins.SearchInRadius(aux_node,search_radius,results.begin(),distances.begin(),MaximumNumberOfResults);
					
				}
				else if(nresults>3)
				{
					search_radius *= 0.999;
					nresults = bins.SearchInRadius(aux_node,search_radius,results.begin(),distances.begin(),MaximumNumberOfResults);
				}
				if(icount>1000)
				{				
					//std::cout<<" ELEMENT ID "<< iel->Id()<<std::endl;
					//std::cout<<" icount "<<icount<<std::endl;
					break;
				}
				
				
			}
			
			
			
			iel->SetValue(SEARCH_RADIUS,search_radius);
			
//**********************************************************************************************************************************			
			
			
            unsigned int counter=0;
            
            for(NodesContainerType::iterator it = results.begin(); it != results.end(); ++it)
            {
                if(counter< nresults)
                {
                    pgeom->push_back( *(it.base()));
                    (*(it.base()))->FastGetSolutionStepValue(IS_STRUCTURE) = 1;
                    
                }
                 
                counter += 1;
               
            }
            

        }
        
        
        //I erase the nodes which do not belong to any clouds
        //for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                //inode!=mr_model_part.NodesEnd(); inode++)
        //{
			////if(inode->Id() == 211)
			////{
				////std::cout<<" IS STRUCTURE "<< inode->FastGetSolutionStepValue(IS_STRUCTURE)<<std::endl;
				////std::cout<<" IS BOUNDARY "<< inode->FastGetSolutionStepValue(IS_BOUNDARY)<<std::endl;
			////}
			//if( (inode->GetSolutionStepValue(IS_STRUCTURE) != 1 ))//(inode->Y() < 0.0 && inode->GetSolutionStepValue(IS_BOUNDARY) != 1)|| (inode->GetSolutionStepValue(IS_STRUCTURE) != 1 && inode->GetSolutionStepValue(IS_BOUNDARY) != 1)) //||(inode->Y() > 0.0 && inode->Y() < 0.016 && inode->X() > 2.117 && inode->GetSolutionStepValue(IS_BOUNDARY) != 1) 
			//{
				//std::cout<<" ID NODE "<< inode->Id()<<std::endl;
				
				
				//inode->Set(TO_ERASE, true);
				
				
			//}
			//inode->FastGetSolutionStepValue(IS_STRUCTURE) = 0;
		//}
		//mr_model_part.RemoveNodesFromAllLevels(TO_ERASE);
//********************code for lagrange multipliers*****************************************************************************************************

        ////creating an auxiliary node for searching purposes
        //Node<3> auxc_node(0,0.0,0.0,0.0);
        //for(ModelPart::ConditionsContainerType::iterator icon = mr_model_part.ConditionsBegin();
                //icon!=mr_model_part.ConditionsEnd(); icon++)
        //{
            //unsigned int MaximumNumberOfResults = 1000;
            //std::vector<NodePointerType>                            results(MaximumNumberOfResults);
            //std::vector<double>                                     distances(MaximumNumberOfResults);

            ////get the coordinates of the gauss point as  the sum of N[i]*xnode
            //const array_1d<double,3> xc = icon->GetGeometry()[0].Coordinates();
            //Geometry< Node<3> > rGeom = icon->GetGeometry();
                   ////KRATOS_WATCH(xc);
            //auxc_node.Coordinates() =  xc;
            //Geometry< Node<3> >::Pointer pcgeom = icon->pGetGeometry();
            //pcgeom->clear();

            ////Geometry< Node<3> >::Pointer pcgeom = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >());

            //for(Geometry< Node<3> >::iterator inode = rGeom.begin(); inode!=rGeom.end(); inode++)
            //{
                //if(inode->Id() == icon->GetValue(CENTER_ID))
                //{
                    //pcgeom->push_back(*(inode.base())) ;
                //}


            //}
			////if(icon->Id()==129)
					////{
						////for(unsigned int i = 0; i< icon->GetGeometry().size(); i++)
						////{
							////KRATOS_WATCH(icon->Id());
							////KRATOS_WATCH(icon->GetGeometry().size());
							
							////KRATOS_WATCH(icon->GetGeometry()[i].Id());
						////}
					////}

            ////Geometry< Node<3> >::Pointer pcgeom = icon->pGetGeometry();



            //double search_radius = icon->GetValue(SEARCH_RADIUS);//0.02;






        ////KRATOS_WATCH(search_radius);
        ////KRATOS_WATCH(auxc_node);

            //unsigned int nresults = bins.SearchInRadius(auxc_node,search_radius,results.begin(),distances.begin(),MaximumNumberOfResults);

            ////KRATOS_WATCH(nresults);
            //unsigned int counter=0;
            //for(NodesContainerType::iterator it = results.begin(); it != results.end(); ++it)
            //{

                //if(counter< nresults)
                //{
                    //if(icon->GetValue(CENTER_ID) != (*it)->Id() )
                    //{
                    
                    //pcgeom->push_back( *(it.base()));
					//}
                //}
                //counter += 1;
            //}

			//if(icon->Id()==129)
					//{
						//for(unsigned int i = 0; i< icon->GetGeometry().size(); i++)
						//{
							//KRATOS_WATCH(icon->Id());
							//KRATOS_WATCH(icon->GetGeometry().size());
							
							//KRATOS_WATCH(icon->GetGeometry()[i].Id());
						//}
					//}

        //}
//***********************************************************************************************************************************************

        KRATOS_CATCH("");
    }





private:

    ModelPart& mr_model_part;


};
}

#endif


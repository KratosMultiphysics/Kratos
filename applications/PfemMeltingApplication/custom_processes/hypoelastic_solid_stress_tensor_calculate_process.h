//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pavel Ryzhakov
//


#if !defined(KRATOS_HYPOELASTIC_CALCULATE_PROCESS_INCLUDED )
#define  KRATOS_HYPOELASTIC_CALCULATE_PROCESS_INCLUDED



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
#include "utilities/geometry_utilities.h"
#include "pfem_melting_application.h"
#include "custom_elements/hypo.h"


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

class HypoelasticStressCalculateProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of HypoelasticStressCalculateProcess
    KRATOS_CLASS_POINTER_DEFINITION(HypoelasticStressCalculateProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HypoelasticStressCalculateProcess(ModelPart& model_part, unsigned int domain_size)
        : mr_model_part(model_part),mdomain_size(domain_size)
    {
    }

    /// Destructor.
    ~HypoelasticStressCalculateProcess() override
    { 
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        KRATOS_TRY

        ProcessInfo& proc_info = mr_model_part.GetProcessInfo();
        
       
	if (mdomain_size==2)
	   {
           Matrix dummy=ZeroMatrix(2,2);
           double Area;
           array_1d<double, 3> N;
           
           BoundedMatrix<double, 3, 2> msDN_Dx;
           
           BoundedMatrix<double,3,3> Sigma=ZeroMatrix(2,2);

	   BoundedMatrix<double,2,2> CauchyStress=ZeroMatrix(2,2);

			

           
           //THIS SHOULD BE EXECUTED ONLY FOR THOSE ELEMENTS OF THE MONOLITHIC MODEL THAT ARE IDENTIFIED TO BELONG TO THE SOLID domain
	   //before the first step we initialize Cauchy stress to zero
	   if (proc_info[TIME]==0.0)
 	      {
              for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ;
                  im != mr_model_part.ElementsEnd() ; ++im)
                  {
    	          //IN A MONOLITHIC FLUID-SOLID MODEL WE WANT TO EXCEUTE THIS FUNCTION ONLY FOR THE SOLID ELEMENTS
                  //if(im->GetGeometry()[0].Is(STRUCTURE) && im->GetGeometry()[1].Is(STRUCTURE) && im->GetGeometry()[2].Is(STRUCTURE))
                      im->SetValue(CAUCHY_STRESS_TENSOR, dummy);
                      CauchyStress=im->GetValue(CAUCHY_STRESS_TENSOR);
                  }
                
  
                for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
            	{
                
                i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XX)=0.0;
                i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XY)=0.0;

                i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YX)=0.0;
                i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YY)=0.0;

  		 i->FastGetSolutionStepValue(DELTA_SIGMA_XX)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_XY)=0.0;

                i->FastGetSolutionStepValue(DELTA_SIGMA_YX)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_YY)=0.0;

          	}
          	    
                  
	      }
	   //and now we actually compute it
	   else
	     {     
	     
	     
	        for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
            	{
                i->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_XX)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_XY)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_YX)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_YY)=0.0;
          	}
          	
          	
             for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ; im != mr_model_part.ElementsEnd() ; ++im)
                 {
	         //IN A MONOLITHIC FLUID-SOLID MODEL WE WANT TO EXCEUTE THIS FUNCTION ONLY FOR THE SOLID ELEMENTS
                 //if(im->GetGeometry()[0].Is(STRUCTURE) && im->GetGeometry()[1].Is(STRUCTURE) && im->GetGeometry()[2].Is(STRUCTURE))
                 if(im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 )
                 
                    im->Calculate(CAUCHY_STRESS_TENSOR,dummy,proc_info);
                    
                    CauchyStress=im->GetValue(CAUCHY_STRESS_TENSOR);

		     //KRATOS_WATCH(CauchyStress)	
                    GeometryUtils::CalculateGeometryData(im->GetGeometry(), msDN_Dx, N, Area);
                    Geometry< Node < 3 > >& geom = im->GetGeometry();
                    
                    for (unsigned int k = 0; k < geom.size(); k++)
			{
                    		//TENSOR PROYECTED IN THE NODES!
                    		geom[k].SetLock();
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_XX) += N[k] * CauchyStress(0,0) * Area;
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_XY) += N[k] * CauchyStress(0,1) * Area;        
                                                              
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_YX) += N[k] * CauchyStress(1,0) * Area;
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_YY) += N[k] * CauchyStress(1,1) * Area;
                                                              
                               //DDDDDDDDDDDDDDDDDDD

                		// So it is safe to write in the node in OpenMP
                		if(Area==0){KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");}
                		double &area=geom[k].FastGetSolutionStepValue(NODAL_AREA); 
                		area +=Area * N[k];
				if (area==0) KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","")
                		geom[k].UnSetLock(); // Free the node for other threads
                               //DDDDDDDDDDDDDDDDDDD
                    }
                 }
                 
                 
                 
                 for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
            	{
            		
            		double &A0 = i->FastGetSolutionStepValue(NODAL_AREA);

            		if (i->FastGetSolutionStepValue(NODAL_AREA)==0)
            		{
            		double pp=0.0;            		
            		
            		i->FastGetSolutionStepValue(DELTA_SIGMA_XX) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_XY) = 0.0;        
                               
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YX) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YY) = 0.0;
                               
			//HISTORICAL CAUCHY STRESS TENSOR 
			i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XX) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XY) = 0.0;        
                               
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YX) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YY) = 0.0;
                               
            		//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");
            		}
            		else{
                	i->FastGetSolutionStepValue(DELTA_SIGMA_XX) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));
                       i->FastGetSolutionStepValue(DELTA_SIGMA_XY) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));        
                               
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YX) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YY) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));
                               
			//HISTORICAL CAUCHY STRESS TENSOR 
			i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XX) += i->FastGetSolutionStepValue(DELTA_SIGMA_XX);
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XY) += i->FastGetSolutionStepValue(DELTA_SIGMA_XY);        
                               
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YX) += i->FastGetSolutionStepValue(DELTA_SIGMA_YX);
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YY) += i->FastGetSolutionStepValue(DELTA_SIGMA_YY);
                               
                       
                       if(i->FastGetSolutionStepValue(IS_INTERFACE) == 0.0)
                       
                       {
                       i->FastGetSolutionStepValue(DELTA_SIGMA_XX) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_XY) = 0.0;        
                               
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YX) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YY) = 0.0;
                               
			//HISTORICAL CAUCHY STRESS TENSOR 
			i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XX) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XY) = 0.0;        
                               
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YX) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YY) = 0.0;
                               
                       }

		}
			
         	  }
         	  
         	  
         	  
	      }//
           }//end 2D
	else if (mdomain_size==3)
           {
           
           double Area;
           array_1d<double, 4> N;
           BoundedMatrix<double, 4, 3> msDN_Dx;
           
        
           Matrix dummy=ZeroMatrix(3,3);
                     
           BoundedMatrix<double,3,3> Sigma=ZeroMatrix(3,3);
	   //BoundedMatrix<double,2,2> HistoricalCauchyStress=ZeroMatrix(2,2);

	   BoundedMatrix<double,3,3> CauchyStress=ZeroMatrix(3,3);
	    
           ///////////////////////
           //first, we initialize CAUCHY_STRESS_TENSOR to zero
	   if (proc_info[TIME]==0.0)
 	      {
              for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ;
                  im != mr_model_part.ElementsEnd() ; ++im)
                  {
                  //GeometryUtils::CalculateGeometryData(im->GetGeometry(), msDN_Dx, N, Area);
                  //if(im->GetGeometry()[0].Is(STRUCTURE) && im->GetGeometry()[1].Is(STRUCTURE) && im->GetGeometry()[2].Is(STRUCTURE) && im->GetGeometry()[3].Is(STRUCTURE))
                      im->SetValue(CAUCHY_STRESS_TENSOR, dummy);
                      //im->SetValue(GREEN_LAGRANGE_STRAIN_TENSOR, dummy);
                      //im->SetValue(PK2_STRESS_TENSOR, dummy);
                      //KRATOS_WATCH(CAUCHY_STRESS_TENSOR)
                  }
                  
                  for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
            	{
                
                i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XX)=0.0;
                i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XY)=0.0;
                i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ)=0.0;

                i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YX)=0.0;
                i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YY)=0.0;
                i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ)=0.0;

                i->FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX)=0.0;
                i->FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY)=0.0;
                i->FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ)=0.0;

		i->FastGetSolutionStepValue(DELTA_SIGMA_XX)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_XY)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_XZ)=0.0;

                i->FastGetSolutionStepValue(DELTA_SIGMA_YX)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_YY)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_YZ)=0.0;

                i->FastGetSolutionStepValue(DELTA_SIGMA_ZX)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_ZY)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_ZZ)=0.0;
                

          	}
          	
	      }
	   //and now we actually compute it
	   else
	     {          
	     
             double dummy_aux=0.0;   
             
            for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
            	{
                i->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
                
                i->FastGetSolutionStepValue(DELTA_SIGMA_XX)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_XY)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_XZ)=0.0;

                i->FastGetSolutionStepValue(DELTA_SIGMA_YX)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_YY)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_YZ)=0.0;

                i->FastGetSolutionStepValue(DELTA_SIGMA_ZX)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_ZY)=0.0;
                i->FastGetSolutionStepValue(DELTA_SIGMA_ZZ)=0.0;
                
          	}
                    
             /*for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ; im != mr_model_part.ElementsEnd() ; ++im)
                 {              
	     		im->Calculate(NODAL_AREA,dummy_aux,proc_info); 
	     	}*/
             for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ; im != mr_model_part.ElementsEnd() ; ++im)
                 {
                      
                         
                 if(im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && im->GetGeometry()[3].FastGetSolutionStepValue(IS_INTERFACE) == 1.0)
                 {
                    GeometryUtils::CalculateGeometryData(im->GetGeometry(), msDN_Dx, N, Area);
                    
                    Geometry< Node < 3 > >& geom = im->GetGeometry();
									//int nn=0;
			im->Calculate(CAUCHY_STRESS_TENSOR,dummy,proc_info);
			//CauchyStress=im->FastGetSolutionStepValue(CAUCHY_STRESS_TENSOR);
			CauchyStress=im->GetValue(CAUCHY_STRESS_TENSOR);
			
                    	for (unsigned int k = 0; k < geom.size(); k++)
			{
                    		//TENSOR PROYECTED IN THE NODES!
                    		geom[k].SetLock();
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_XX) += N[k] * CauchyStress(0,0) * Area;
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_XY) += N[k] * CauchyStress(0,1) * Area;        
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_XZ) += N[k] * CauchyStress(0,2) * Area;
                               
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_YX) += N[k] * CauchyStress(1,0) * Area;
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_YY) += N[k] * CauchyStress(1,1) * Area;
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_YZ) += N[k] * CauchyStress(1,2) * Area;
                               
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_ZX) += N[k] * CauchyStress(2,0) * Area;
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_ZY) += N[k] * CauchyStress(2,1) * Area;
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_ZZ) += N[k] * CauchyStress(2,2) * Area;
                               
                               //DDDDDDDDDDDDDDDDDDD

                		// So it is safe to write in the node in OpenMP
                		if(Area==0){KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");}
                		double &area=geom[k].FastGetSolutionStepValue(NODAL_AREA); 
                		area +=Area * N[k];
				if (area==0) KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","")
                		geom[k].UnSetLock(); // Free the node for other threads
                               //DDDDDDDDDDDDDDDDDDD
                    }
                 }
	      }
	      
	      for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
            	{
            		//KRATOS_WATCH("#########################################")
            		//KRATOS_WATCH("#########################################")
            	
            		//KRATOS_WATCH(i->FastGetSolutionStepValue(NODAL_AREA))
            		
            		double &A0 = i->FastGetSolutionStepValue(NODAL_AREA);

            		if (i->FastGetSolutionStepValue(NODAL_AREA)==0)
            		{
            		double pp=0.0;            		
            		//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");
            		
            		i->FastGetSolutionStepValue(DELTA_SIGMA_XX) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_XY) = 0.0;        
                       i->FastGetSolutionStepValue(DELTA_SIGMA_XZ) = 0.0;
                               
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YX) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YY) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YZ) = 0.0;
                               
                       i->FastGetSolutionStepValue(DELTA_SIGMA_ZX) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_ZY) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_ZZ) = 0.0;
                               
			//HISTORICAL CAUCHY STRESS TENSOR 
			i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XX) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XY) = 0.0;        
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ) = 0.0;
                               
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YX) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YY) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ) = 0.0;
                               
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ) = 0.0;	
            		
            		//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");
            		}
            		else{
                	i->FastGetSolutionStepValue(DELTA_SIGMA_XX) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));
                       i->FastGetSolutionStepValue(DELTA_SIGMA_XY) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));        
                       i->FastGetSolutionStepValue(DELTA_SIGMA_XZ) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));
                               
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YX) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YY) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YZ) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));
                               
                       i->FastGetSolutionStepValue(DELTA_SIGMA_ZX) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));
                       i->FastGetSolutionStepValue(DELTA_SIGMA_ZY) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));
                       i->FastGetSolutionStepValue(DELTA_SIGMA_ZZ) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));
                               
			//HISTORICAL CAUCHY STRESS TENSOR 
			i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XX) += i->FastGetSolutionStepValue(DELTA_SIGMA_XX);
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XY) += i->FastGetSolutionStepValue(DELTA_SIGMA_XY);        
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ) += i->FastGetSolutionStepValue(DELTA_SIGMA_XZ);
                               
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YX) += i->FastGetSolutionStepValue(DELTA_SIGMA_YX);
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YY) += i->FastGetSolutionStepValue(DELTA_SIGMA_YY);
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ) += i->FastGetSolutionStepValue(DELTA_SIGMA_YZ);
                               
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX) += i->FastGetSolutionStepValue(DELTA_SIGMA_ZX);
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY) += i->FastGetSolutionStepValue(DELTA_SIGMA_ZY);
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ) += i->FastGetSolutionStepValue(DELTA_SIGMA_ZZ);		
                       
                       
                       if(i->FastGetSolutionStepValue(IS_INTERFACE) == 0.0)
                       
                       {
                       i->FastGetSolutionStepValue(DELTA_SIGMA_XX) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_XY) = 0.0;        
                       i->FastGetSolutionStepValue(DELTA_SIGMA_XZ) = 0.0;
                               
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YX) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YY) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YZ) = 0.0;
                               
                       i->FastGetSolutionStepValue(DELTA_SIGMA_ZX) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_ZY) = 0.0;
                       i->FastGetSolutionStepValue(DELTA_SIGMA_ZZ) = 0.0;
                               
			//HISTORICAL CAUCHY STRESS TENSOR 
			i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XX) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XY) = 0.0;        
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XZ) = 0.0;
                               
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YX) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YY) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YZ) = 0.0;
                               
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_ZX) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_ZY) = 0.0;
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_ZZ) = 0.0;	
                       
                       
                       }

			/*KRATOS_WATCH(i->FastGetSolutionStepValue(NODAL_AREA))

			KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_XX))
			KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_XY))
			KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_XZ))

			KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_YX))
			KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_YY))
			KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_YZ))

			KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_ZX))
			KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_ZY))
			KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_ZZ))*/
			}
			
         	  }
	}
        //KRATOS_WATCH("Executed of Cauchy stress tensor computation of the hypoelastic element");
        }
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
        return "HypoelasticStressCalculateProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "HypoelasticStressCalculateProcess";
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
    ModelPart& mr_model_part;
    unsigned int mdomain_size;


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
//		HypoelasticStressCalculateProcess& operator=(HypoelasticStressCalculateProcess const& rOther);

    /// Copy constructor.
//		HypoelasticStressCalculateProcess(HypoelasticStressCalculateProcess const& rOther);


    ///@}

}; // Class HypoelasticStressCalculateProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  HypoelasticStressCalculateProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const HypoelasticStressCalculateProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_HYPOELASTIC_CALCULATE_PROCESS_INCLUDED  defined 



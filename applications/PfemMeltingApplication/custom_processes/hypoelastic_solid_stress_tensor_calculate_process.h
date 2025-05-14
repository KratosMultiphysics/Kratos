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
//#include "custom_elements/qfluid.h"


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

           BoundedMatrix<double,2,2> Sigma=ZeroMatrix(2,2);
           
	   BoundedMatrix<double,2,2> Sigma_old=ZeroMatrix(2,2);
  	   BoundedMatrix<double,2,2> Sigma_star;
           

	   BoundedMatrix<double,2,2> CauchyStress=ZeroMatrix(2,2);
	   
   	   BoundedMatrix<double,2,2> CauchyStress_old=ZeroMatrix(2,2);


/*	   if (proc_info[TIME]==0.0)
 	      {
              for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ;
                  im != mr_model_part.ElementsEnd() ; ++im)
                  {
    	          //IN A MONOLITHIC FLUID-SOLID MODEL WE WANT TO EXCEUTE THIS FUNCTION ONLY FOR THE SOLID ELEMENTS
                  //if(im->GetGeometry()[0].Is(STRUCTURE) && im->GetGeometry()[1].Is(STRUCTURE) && im->GetGeometry()[2].Is(STRUCTURE))
                      im->SetValue(CAUCHY_STRESS_TENSOR, dummy);
                      im->Calculate(CAUCHY_STRESS_TENSOR,dummy,proc_info);
                      CauchyStress=im->GetValue(CAUCHY_STRESS_TENSOR);
                      im->SetValue(CAUCHY_STRESS_TENSOR, dummy);
                  }
             }*/
/*

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

KRATOS_THROW_ERROR(std::logic_error,"pressure calculation jjjjjjjjj not implemented","");
	      }
           else{ */

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
                 dummy=ZeroMatrix(2,2);
                 
                 //KRATOS_WATCH("PRIMERO")
                 //KRATOS_WATCH(im->GetValue(CAUCHY_STRESS_TENSOR,(0,0)));
                 //double pepe1=im->GetValue(CAUCHY_STRESS_TENSOR)[0,0];
                 //im->GetValue(GREEN_LAGRANGE_STRAIN_TENSOR)=im->GetValue(CAUCHY_STRESS_TENSOR);
                 
                 
                 //Sigma=im->GetValue(CAUCHY_STRESS_TENSOR,0)
                 //double pepe1=im->GetValue(CAUCHY_STRESS_TENSOR)(0,0);
                 //KRATOS_WATCH(pepe1)
                 //Sigma(0,0)=im->GetValue(CAUCHY_STRESS_TENSOR)(0,0);
                 //Sigma(0,1)=im->GetValue(CAUCHY_STRESS_TENSOR)(0,1);
                 //Sigma(1,0)=im->GetValue(CAUCHY_STRESS_TENSOR)(1,0);
                 //Sigma(1,1)=im->GetValue(CAUCHY_STRESS_TENSOR)(1,1);
                 //KRATOS_WATCH("jhkhkjhkjkjhjkhk")
                 //KRATOS_WATCH(Sigma)
                 
                 //im->Calculate(CAUCHY_STRESS_TENSOR,dummy,proc_info);
                 //KRATOS_WATCH("SEGUNDO")
                 //KRATOS_WATCH(dummy)
                 
                 //KRATOS_WATCH(im->GetValue(CAUCHY_STRESS_TENSOR))
                 //KRATOS_WATCH(im->GetValue(GREEN_LAGRANGE_STRAIN_TENSOR))
                 //double pepe=dummy(0,0);
                 
                 
                 /*Sigma(0,0)=dummy(0,0);
                 Sigma(0,1)=dummy(0,1);
                 Sigma(1,0)=dummy(1,0);
                 Sigma(1,1)=dummy(1,1);*/
                 //KRATOS_WATCH("gggggggggggggggggggggggggggggggggggggggggggg")
                 //im->GetValue(PK2_STRESS_TENSOR) = im->GetValue(CAUCHY_STRESS_TENSOR) + im->GetValue(GREEN_LAGRANGE_STRAIN_TENSOR);
                 //KRATOS_WATCH("gggggggggggggggggggggggggggggggggggggggggggg")
                 //im->GetValue(CAUCHY_STRESS_TENSOR) = im->GetValue(PK2_STRESS_TENSOR);
                       //im->SetValue(CAUCHY_STRESS_TENSOR, dummy);
                  	
                 CauchyStress_old=im->GetValue(CAUCHY_STRESS_TENSOR);
                 
                 
                 //KRATOS_WATCH("<------------------------>")
                 //KRATOS_WATCH("<------------------------>")
                 //KRATOS_WATCH("<------------------------>")
                 //KRATOS_WATCH("INICIO")
                 
                 Sigma_old(0,0) = CauchyStress_old(0,0);
                 Sigma_old(0,1) = CauchyStress_old(0,1);
                 Sigma_old(1,0) = CauchyStress_old(1,0);
                 Sigma_old(1,1) = CauchyStress_old(1,1);
                 
                 
                 //KRATOS_WATCH("<------------------------>")
                 //KRATOS_WATCH("<------------------------>")
                 //KRATOS_WATCH("<------------------------>")
                 //KRATOS_WATCH("TENSION DEL PASO ANTERIOR")
                 
                 //KRATOS_WATCH(Sigma_old(0,0))	
		  //KRATOS_WATCH(Sigma_old(0,1))	
		  //KRATOS_WATCH(Sigma_old(1,0))	
		  //KRATOS_WATCH(Sigma_old(1,1))	
		     
                 
                 im->Calculate(CAUCHY_STRESS_TENSOR,dummy,proc_info);
                 //CauchyStress_old=im->GetValue(CAUCHY_STRESS_TENSOR);
                 //KRATOS_WATCH(CauchyStress_old)
                 //KRATOS_THROW_ERROR(std::logic_error,"pressure calculation jjjjjjjjj not implemented","");
 
                 CauchyStress=dummy;
                     
                 //CauchyStress(0,0)=dummy(0,0);
                 //CauchyStress(0,0)=dummy(0,0);
                 //CauchyStress(0,0)=dummy(0,0);
                 //CauchyStress(0,0)=dummy(0,0);
                 
                 
                 //KRATOS_WATCH("TENSOR INCREMENTO")
                 Sigma(0,0)=im->GetValue(CAUCHY_STRESS_TENSOR)(0,0);  
                 Sigma(0,1)=im->GetValue(CAUCHY_STRESS_TENSOR)(0,1); 
                 Sigma(1,0)=im->GetValue(CAUCHY_STRESS_TENSOR)(1,0); 
                 Sigma(1,1)=im->GetValue(CAUCHY_STRESS_TENSOR)(1,1); 
                    
                     
                //KRATOS_WATCH("TERCERO")
                //KRATOS_WATCH(CauchyStress)

		 //KRATOS_WATCH("CUARTO")
		    
		 /*Sigma(0,0) += CauchyStress(0,0);
                 Sigma(0,1) += CauchyStress(0,1);
                 Sigma(1,0) += CauchyStress(1,0);
                 Sigma(1,1) += CauchyStress(1,1);*/
                 
		 //KRATOS_WATCH(Sigma(0,0))	
		 //KRATOS_WATCH(Sigma(0,1))	
		 //KRATOS_WATCH(Sigma(1,0))	
		 //KRATOS_WATCH(Sigma(1,1))	
		                         
                //KRATOS_WATCH("INCREMENTO DEL TENSORRRRRRRRRRRRRRRR")
                //KRATOS_WATCH(dummy)
                //im->SetValue(CAUCHY_STRESS_TENSOR, dummy);
                //KRATOS_WATCH(im->GetValue(CAUCHY_STRESS_TENSOR))
                    
                    
                //KRATOS_WATCH(dummy)
                //KRATOS_WATCH(CauchyStress)
                    Sigma(0,0)+=Sigma_old(0,0);  
                    Sigma(0,1)+=Sigma_old(0,1); 
                    Sigma(1,0)+=Sigma_old(1,0); 
                    Sigma(1,1)+=Sigma_old(1,1); 


		     //KRATOS_WATCH("TENSOR FINALLLLLLLLLLLLLLLLLLLLL")
		     //KRATOS_WATCH(Sigma(0,0))	
		     //KRATOS_WATCH(Sigma(0,1))	
		     //KRATOS_WATCH(Sigma(1,0))	
		     //KRATOS_WATCH(Sigma(1,1))	

                    
                    im->GetValue(CAUCHY_STRESS_TENSOR)(0,0)=Sigma(0,0);  
                    im->GetValue(CAUCHY_STRESS_TENSOR)(0,1)=Sigma(0,1);  
                    im->GetValue(CAUCHY_STRESS_TENSOR)(1,0)=Sigma(1,0);  
                    im->GetValue(CAUCHY_STRESS_TENSOR)(1,1)=Sigma(1,1);                      
                    
                    Sigma(0,0)=im->GetValue(CAUCHY_STRESS_TENSOR)(0,0);  
                    Sigma(0,1)=im->GetValue(CAUCHY_STRESS_TENSOR)(0,1); 
                    Sigma(1,0)=im->GetValue(CAUCHY_STRESS_TENSOR)(1,0); 
                    Sigma(1,1)=im->GetValue(CAUCHY_STRESS_TENSOR)(1,1); 


		     //KRATOS_WATCH("TENSOR FINALLLLLLLLLLLLLLLL")
		     //KRATOS_WATCH(Sigma(0,0))	
		     //KRATOS_WATCH(Sigma(0,1))	
		     //KRATOS_WATCH(Sigma(1,0))	
		     //KRATOS_WATCH(Sigma(1,1))	

			//Sigma_star=im->GetValue(GREEN_LAGRANGE_STRAIN_TENSOR);
		     	
		     	
		     	//KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TOTAL_CAUCHY_STRESS )
		     	
	
		     	
		     //im->GetValue(GREEN_LAGRANGE_STRAIN_TENSOR)(0,0)=100.0;  
		     //im->SetValuesOnIntegrationPoints( CAUCHY_STRESS_TENSOR, Sigma, proc_info);
                    /*Sigma_star(0,0)=im->GetValue(GREEN_LAGRANGE_STRAIN_TENSOR)(0,0); 
                    Sigma_star(0,1)=im->GetValue(GREEN_LAGRANGE_STRAIN_TENSOR)(0,1); 
                    Sigma_star(1,0)=im->GetValue(GREEN_LAGRANGE_STRAIN_TENSOR)(1,0); 
                    Sigma_star(1,1)=im->GetValue(GREEN_LAGRANGE_STRAIN_TENSOR)(1,1); */                                       
                    
                    
                    Sigma_star(0,0)=100; 
                    //im->GetValue(PK2_STRESS_TENSOR)(0,0)=100.0; 
                    //Sigma_star(0,0)=im->GetValue(GREEN_LAGRANGE_STRAIN_TENSOR)(0,0); 
                    //im->SetValue(PK2_STRESS_TENSOR, dummy);
                    
                    //KRATOS_WATCH(Sigma_star)	
		     /*KRATOS_WATCH(Sigma_star(0,1))	
		     KRATOS_WATCH(Sigma_star(1,0))	
		     KRATOS_WATCH(Sigma_star(1,1))*/	
                    
                    
                    //KRATOS_WATCH("<------------------------>")
                    //KRATOS_WATCH("<------------------------>")
                    //KRATOS_WATCH("<------------------------>")
                  
                    //KRATOS_WATCH("FINAL")
                    
                    //KRATOS_THROW_ERROR(std::logic_error,"pressure calculation jjjjjjjjj not implemented","");
                    
   

                    GeometryUtils::CalculateGeometryData(im->GetGeometry(), msDN_Dx, N, Area);
                    Geometry< Node >& geom = im->GetGeometry();
                    
                    for (unsigned int k = 0; k < geom.size(); k++)
			{
                    		geom[k].SetLock();
                               /*geom[k].FastGetSolutionStepValue(DELTA_SIGMA_XX) += N[k] * dummy(0,0) * Area;
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_XY) += N[k] * dummy(0,1) * Area;

                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_YX) += N[k] * dummy(1,0) * Area;
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_YY) += N[k] * dummy(1,1) * Area;*/
                               
                               
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_XX) += dummy(0,0) * Area;
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_XY) += dummy(0,1) * Area;

                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_YX) += dummy(1,0) * Area;
                               geom[k].FastGetSolutionStepValue(DELTA_SIGMA_YY) += dummy(1,1) * Area;
                               
                               
				//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation jjjjjjjjj not implemented","");
                		// So it is safe to write in the node in OpenMP
                		if(Area==0){KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");}
                		double &area=geom[k].FastGetSolutionStepValue(NODAL_AREA);
                		
                		/*area +=Area * N[k];*/
               		area +=Area;
                		//KRATOS_WATCH(area)
                		//else area +=0.0;
				//if (area==0) KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","")
                		geom[k].UnSetLock(); // Free the node for other threads
                    }
                    
                    //}
                 }

KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");

                 for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
            	{

            		if (i->FastGetSolutionStepValue(NODAL_AREA)>0){

            		//KRATOS_WATCH(i->FastGetSolutionStepValue(NODAL_AREA))
                	i->FastGetSolutionStepValue(DELTA_SIGMA_XX) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));
                       i->FastGetSolutionStepValue(DELTA_SIGMA_XY) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));

                       i->FastGetSolutionStepValue(DELTA_SIGMA_YX) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));
                       i->FastGetSolutionStepValue(DELTA_SIGMA_YY) *= (1.0/i->FastGetSolutionStepValue(NODAL_AREA));

			//HISTORICAL CAUCHY STRESS TENSOR
			i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XX) += i->FastGetSolutionStepValue(DELTA_SIGMA_XX);
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XY) += i->FastGetSolutionStepValue(DELTA_SIGMA_XY);

                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YX) += i->FastGetSolutionStepValue(DELTA_SIGMA_YX);
                       i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YY) += i->FastGetSolutionStepValue(DELTA_SIGMA_YY);



			}

         	  }
         	  
         	  //KRATOS_WATCH("<------------------------>")
         	  //KRATOS_WATCH("<------------------------>")
         	  //KRATOS_WATCH("<------------------------>")
         	           	           	  
         	  
         	for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
            	{

		//KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_XX))
		//KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_XY))
		//KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_YX))
		//KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_YY))
		

		//KRATOS_WATCH(i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XX))
		//KRATOS_WATCH(i->FastGetSolutionStepValue(HISTORICAL_SIGMA_XY))
		//KRATOS_WATCH(i->FastGetSolutionStepValue(HISTORICAL_SIGMA_YX))
		//KRATOS_WATCH(i->FastGetSolutionStepValue(DELTA_SIGMA_YY))

		

}  

		//}

           }//end 2D


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

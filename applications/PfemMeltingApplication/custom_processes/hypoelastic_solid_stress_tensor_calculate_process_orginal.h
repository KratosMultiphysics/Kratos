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

             bool hypo=false; 
             for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ; im != mr_model_part.ElementsEnd() ; ++im)
                 {
	         //IN A MONOLITHIC FLUID-SOLID MODEL WE WANT TO EXCEUTE THIS FUNCTION ONLY FOR THE SOLID ELEMENTS
                 //if(im->GetGeometry()[0].Is(STRUCTURE) && im->GetGeometry()[1].Is(STRUCTURE) && im->GetGeometry()[2].Is(STRUCTURE))
                 //if(im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 )
                 hypo=false;
                 dummy=ZeroMatrix(2,2);
                 //im->SetValue(CAUCHY_STRESS_TENSOR, dummy);
                 if(im->GetGeometry()[0].FastGetSolutionStepValue(IS_SOLID) == 1.0 && im->GetGeometry()[1].FastGetSolutionStepValue(IS_SOLID) == 1.0 && im->GetGeometry()[2].FastGetSolutionStepValue(IS_SOLID) == 1.0 )
                  {  
                  	im->Calculate(CAUCHY_STRESS_TENSOR,dummy,proc_info);
                       //im->SetValue(CAUCHY_STRESS_TENSOR, dummy);
                  	
			hypo=true;
			
		   }
//	            if(hypo==false) 
                    //CauchyStress=im->GetValue(CAUCHY_STRESS_TENSOR);
                    CauchyStress=dummy;  

		     //KRATOS_WATCH(CauchyStress)
                    GeometryUtils::CalculateGeometryData(im->GetGeometry(), msDN_Dx, N, Area);
                    Geometry< Node < 3 > >& geom = im->GetGeometry();
                    
                    if(hypo==true)
                    {
                    for (unsigned int k = 0; k < geom.size(); k++)
			{
                    		//TENSOR PROYECTED IN THE NODES!
                    		//KRATOS_WATCH(CauchyStress)
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
                		//else area +=0.0;
				//if (area==0) KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","")
                		geom[k].UnSetLock(); // Free the node for other threads
                               //DDDDDDDDDDDDDDDDDDD
                    }
                    
                    }
                 }



                 for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
            	{

            		/*if (i->FastGetSolutionStepValue(NODAL_AREA)==0)
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
			//KRATOS_WATCH(i->FastGetSolutionStepValue(NODAL_AREA))
            		//KRATOS_THROW_ERROR(std::logic_error,"pressure calculation 3D not implemented","");
            		}
            		else{*/
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


                     /*  if(i->FastGetSolutionStepValue(IS_SOLID) == 0.0)

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

                       }*/

		}

         	  }



	      }//
           }//end 2D
	else if (mdomain_size==3)
           {

           double Area;
           array_1d<double, 4> N;
           BoundedMatrix<double, 4, 3> msDN_Dx;

			}

         	 // }
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



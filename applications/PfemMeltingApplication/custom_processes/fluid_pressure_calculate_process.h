//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//


#if !defined(KRATOS_FLUID_PRESSURE_CALCULATE_PROCESS_INCLUDED )
#define  KRATOS_FLUID_PRESSURE_CALCULATE_PROCESS_INCLUDED



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

class FluidPressureCalculateProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of HypoelasticStressCalculateProcess
    KRATOS_CLASS_POINTER_DEFINITION(FluidPressureCalculateProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FluidPressureCalculateProcess(ModelPart& model_part, unsigned int domain_size)
        : mr_model_part(model_part),mdomain_size(domain_size)
    {
    }

    /// Destructor.
    ~FluidPressureCalculateProcess() override
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
           //Matrix dummy=ZeroMatrix(2,2);
           double Area;
           double dummy;
           array_1d<double, 3> N;

           BoundedMatrix<double, 3, 2> DN_DX;
           

           //THIS SHOULD BE EXECUTED ONLY FOR THOSE ELEMENTS OF THE MONOLITHIC MODEL THAT ARE IDENTIFIED TO BELONG TO THE SOLID domain
	   //before the first step we initialize Cauchy stress to zero
//	   if (proc_info[TIME]==0.0)
// 	      {
//                for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
//            	{
//                i->FastGetSolutionStepValue(PRESSUREAUX)=0.0;
//                i->FastGetSolutionStepValue(NODAL_VOLUME) = 0.0;
//          	}
//	      }
	   //and now we actually compute it
//	   else
//	     {
	        for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
            	{
                i->FastGetSolutionStepValue(PRESSUREAUX)=0.0;
                i->FastGetSolutionStepValue(NODAL_VOLUME) = 0.0;
          	}
	     //}



	     int solid_nodes = 0;
             for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ; im != mr_model_part.ElementsEnd() ; ++im)
                 {
	         //IN A MONOLITHIC FLUID-SOLID MODEL WE WANT TO EXCEUTE THIS FUNCTION ONLY FOR THE SOLID ELEMENTS
                 //if(im->GetGeometry()[0].Is(STRUCTURE) && im->GetGeometry()[1].Is(STRUCTURE) && im->GetGeometry()[2].Is(STRUCTURE))
                 //if(im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 && im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE) == 1.0 )

                    solid_nodes = 0;
                 	
                 		Geometry< Node >& geom = im->GetGeometry();

				for(unsigned int i2 = 0; i2 < geom.size(); i2++)
				{
					if(geom[i2].FastGetSolutionStepValue(IS_SOLID) == 1) solid_nodes += 1;
					//KRATOS_WATCH("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

					//KRATOS_WATCH("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")										
					//KRATOS_WATCH(solid_nodes)
				}	    
                    
                	if (solid_nodes<3){
                               //KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl; 
                    		im->Calculate(NODAL_VOLUME,dummy,proc_info);
				    for (unsigned int i = 0; i < 3; ++i)
				    {

					geom[i].SetLock(); // So it is safe to write in the node in OpenMP
					geom[i].FastGetSolutionStepValue(NODAL_VOLUME) += Area * N[i];
					//KRATOS_WATCH(Area)
					//KRATOS_WATCH(N)
					geom[i].UnSetLock(); // Free the node for other threads
				    }
            
		     		GeometryUtils::CalculateGeometryData(im->GetGeometry(), DN_DX, N, Area);
                    		Geometry< Node >& geom = im->GetGeometry();
                    
                    		const double dt = mr_model_part.GetProcessInfo()[DELTA_TIME];
                            

  		        	double K=-0.333333333333333333*(geom[0].FastGetSolutionStepValue(BULK_MODULUS)+geom[1].FastGetSolutionStepValue(BULK_MODULUS)+geom[2].FastGetSolutionStepValue(BULK_MODULUS));
                       	
                       	//KRATOS_WATCH(K)
                       	
                       	//KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl;                	
                       	//K=-100000.0;
                       	K=-25000.0;
                       	//K=-500.0;

				//calculate the divergence
				const array_1d<double,3>& v0 = geom[0].FastGetSolutionStepValue(VELOCITY);
				const array_1d<double,3>& v1 = geom[1].FastGetSolutionStepValue(VELOCITY);
				const array_1d<double,3>& v2 = geom[2].FastGetSolutionStepValue(VELOCITY);

				double div_v = DN_DX(0,0)*v0[0] + DN_DX(0,1)*v0[1] + DN_DX(1,0)*v1[0] + DN_DX(1,1)*v1[1] + DN_DX(2,0)*v2[0] + DN_DX(2,1)*v2[1];
				double dp_el = K * dt * div_v * Area;
		               //dp_el = K * div_v * Area;    
				array_1d<double,3> pn = ZeroVector(3); //dimension = number of nodes
				pn[0] = geom[0].FastGetSolutionStepValue(PRESSURE,1);
				pn[1] = geom[1].FastGetSolutionStepValue(PRESSURE,1);
				pn[2] = geom[2].FastGetSolutionStepValue(PRESSURE,1);
				double diag_term = Area/6.0;
				double out_term = Area/12.0;


				geom[0].FastGetSolutionStepValue(PRESSUREAUX) += dp_el *0.333333333333333333 + pn[0] * 0.333333333333333333 * Area;
				geom[1].FastGetSolutionStepValue(PRESSUREAUX) += dp_el *0.333333333333333333 + pn[1] * 0.333333333333333333 * Area;
				geom[2].FastGetSolutionStepValue(PRESSUREAUX) += dp_el *0.333333333333333333 + pn[2] * 0.333333333333333333 * Area;
				}
		}
                    
                    for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); i!=mr_model_part.NodesEnd(); i++)
            	{
            	
            	 const double& ar= i->FastGetSolutionStepValue(NODAL_VOLUME);
                //KRATOS_WATCH(i->FastGetSolutionStepValue(NODAL_VOLUME))


            	 double aux=ar;
            	 if(aux<0.00000000000000000000000000001) aux=1.0;
            	 //KRATOS_WATCH(aux)
            	 //aux=1.0;
            	 i->FastGetSolutionStepValue(PRESSUREAUX)/=aux;
                //i->FastGetSolutionStepValue(PRESSURE,1)=i->FastGetSolutionStepValue(PRESSUREAUX);
                i->FastGetSolutionStepValue(PRESSURE)=i->FastGetSolutionStepValue(PRESSUREAUX);
                //if(i->FastGetSolutionStepValue(IS_FREE_SURFACE)==1) i->FastGetSolutionStepValue(PRESSURE)=0.0;    
                 }
                 
                 
            //KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl; 




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
        return "FluidPressureCalculateProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FluidPressureCalculateProcess";
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
//		FluidPressureCalculateProcess& operator=(HypoelasticStressCalculateProcess const& rOther);

    /// Copy constructor.
//		FluidPressureCalculateProcess(HypoelasticStressCalculateProcess const& rOther);


    ///@}

}; // Class FluidPressureCalculateProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FluidPressureCalculateProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FluidPressureCalculateProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FLUID_CALCULATE_PROCESS_INCLUDED  defined



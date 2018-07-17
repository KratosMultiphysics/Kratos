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
#include "ULF_application.h"


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
        Matrix dummy=ZeroMatrix(2,2);
	
	//THIS SHOULD BE EXECUTED ONLY FOR THOSE ELEMENTS OF THE MONOLITHIC MODEL THAT ARE IDENTIFIED TO BELONG TO THE SOLID domain
	//THIS CAN BE DONE BY USING SOME FLAG... TO DO... now it is applied to all elements
        
	//before the first step we initialize Cauchy stress to zero
	//KRATOS_WATCH(proc_info[TIME])
	
	if (proc_info[TIME]==0.0)
	{
          for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ;
                im != mr_model_part.ElementsEnd() ; ++im)
          {
	  //IN A MONOLITHIC FLUID-SOLID MODEL WE WANT TO EXCEUTE THIS FUNCTION ONLY FOR THE SOLID ELEMENTS
          if(im->GetGeometry()[0].Is(STRUCTURE) && im->GetGeometry()[0].Is(STRUCTURE) && im->GetGeometry()[0].Is(STRUCTURE))
            im->SetValue(CAUCHY_STRESS_TENSOR, dummy);
          }
	}
	//and now we actually compute it
	else
	{          
          for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ;
                im != mr_model_part.ElementsEnd() ; ++im)
          {
	  //IN A MONOLITHIC FLUID-SOLID MODEL WE WANT TO EXCEUTE THIS FUNCTION ONLY FOR THE SOLID ELEMENTS
          if(im->GetGeometry()[0].Is(STRUCTURE) && im->GetGeometry()[0].Is(STRUCTURE) && im->GetGeometry()[0].Is(STRUCTURE))
            im->Calculate(CAUCHY_STRESS_TENSOR,dummy,proc_info);
          }
	}
        KRATOS_WATCH("Executed of Cauchy stress tensor computation of the hypoelastic element");
        

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
    double m_min_h;
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



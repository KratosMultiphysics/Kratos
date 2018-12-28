//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti.
//


//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-10-23 12:50:01 $
//   Revision:            $Revision: 1.4 $
//
//  this process save structural elements in a separate list

#if !defined(KRATOS_ULF_APPLY_BC_PROCESS_INCLUDED )
#define  KRATOS_ULF_APPLY_BC_PROCESS_INCLUDED



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
//#include "custom_utilities/geometry_utilities2D.h"
//#include "custom_elements/updated_lagrangian_fluid.h"


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

  class Pfem2ApplyBCProcess
    : public Process
  {
  public:
    ///@name Type Definitions
  ///@{

  /// Pointer definition of PushStructureProcess
  KRATOS_CLASS_POINTER_DEFINITION(Pfem2ApplyBCProcess);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
 Pfem2ApplyBCProcess(ModelPart& model_part)
      : mr_model_part(model_part)
    {
    }

    /// Destructor.
    virtual ~Pfem2ApplyBCProcess()
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

    virtual void Execute() override
    {
      KRATOS_TRY
        for(ModelPart::NodesContainerType::const_iterator in = mr_model_part.NodesBegin(); in!=mr_model_part.NodesEnd(); in++)
	  {
            //marking wet nodes
            /*					if(in->FastGetSolutionStepValue(IS_STRUCTURE) )
						if( (in->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0)
						in->FastGetSolutionStepValue(IS_FLUID) = 1.0;
						else //it is not anymore of fluid
						in->FastGetSolutionStepValue(IS_FLUID) = 0.0;
						//marking as free surface the lonely nodes
								else
            */						if( (in->GetValue(NEIGHBOUR_ELEMENTS)).size() == 0)
	    in->FastGetSolutionStepValue(IS_BOUNDARY) = 1.0;
        }

        //identify the free surface
        for(ModelPart::NodeIterator i = mr_model_part.NodesBegin() ;
	    i != mr_model_part.NodesEnd() ; ++i)
	  {
            //reset the free surface
            i->FastGetSolutionStepValue(IS_FREE_SURFACE) = 0;
            //identify the free surface and fix the pressure accordingly
            if( i->FastGetSolutionStepValue(IS_BOUNDARY) != 0
		&&
		i->FastGetSolutionStepValue(IS_STRUCTURE) == 0)
	      {
                i->FastGetSolutionStepValue(IS_FREE_SURFACE) = 1;
	      }
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
    virtual std::string Info() const override
    {
      return "Pfem2ApplyBCProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "Pfem2ApplyBCProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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
    //		UlfApplyBCProcess& operator=(UlfApplyBCProcess const& rOther);

    /// Copy constructor.
    //		UlfApplyBCProcess(UlfApplyBCProcess const& rOther);


    ///@}

}; // Class UlfApplyBCProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
 inline std::istream& operator >> (std::istream& rIStream,
				   Pfem2ApplyBCProcess& rThis);

 /// output stream function
 inline std::ostream& operator << (std::ostream& rOStream,
				   const Pfem2ApplyBCProcess& rThis)
 {
   rThis.PrintInfo(rOStream);
   rOStream << std::endl;
   rThis.PrintData(rOStream);

   return rOStream;
 }
 ///@}


}  // namespace Kratos.

#endif // KRATOS_ULF_APPLY_BC_PROCESS_INCLUDED  defined

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Alex Jarauta
//  Co-author  :     Elaf Mahrous


#if !defined(ASSIGN_SURFACE_TENSION2D_CONDITION_INCLUDED )
#define  ASSIGN_SURFACE_TENSION2D_CONDITION_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "processes/process.h"
#include "includes/condition.h"
#include "includes/element.h"
#include "ULF_application.h"
#include "custom_conditions/Surface_Tension2D.h"



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
	assign surface tension condition for water droplets


*/


  class AssignSurfaceTensionConditions
  : public Process
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignSurfaceTensionConditions);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AssignSurfaceTensionConditions()
    {	
    }

    /// Destructor.
    ~AssignSurfaceTensionConditions() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    //	void operator()()
    //	{
    //		MergeParts();
    //	}


    ///@}
    ///@name Operations
    ///@{
    
    void AssignSurfaceTensionConditions2D(ModelPart& ThisModelPart)
    {
        KRATOS_TRY
        
        // This is the code written to apply the surface tension condition looping over nodes
        
        Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);
	int id = ThisModelPart.Conditions().size();
	for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; ++im)
	{
	  if (im->FastGetSolutionStepValue(IS_FREE_SURFACE) != 0.0 || im->FastGetSolutionStepValue(TRIPLE_POINT) != 0.0)
	  {
	    Condition::NodesArrayType temp;
	    temp.reserve(1);
	    temp.push_back(*(im.base()));
	    Condition::Pointer p_cond = (KratosComponents<Condition>::Get("SurfaceTension2D")).Create(id, temp, properties);
	    ThisModelPart.Conditions().push_back(p_cond);
	    id++;
	  }
	}	
	ThisModelPart.Conditions().Sort();	
        KRATOS_WATCH(ThisModelPart.Conditions().size())

	
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
        return "AssignSurfaceTensionConditions";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignSurfaceTensionConditions";
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
//		AssignSurfaceTensionConditions& operator=(AssignSurfaceTensionConditions const& rOther);

    /// Copy constructor.
//		AssignSurfaceTensionConditions(AssignSurfaceTensionConditions const& rOther);


    ///@}

}; // Class AssignSurfaceTensionConditions

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignSurfaceTensionConditions& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignSurfaceTensionConditions& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}



}  // namespace Kratos.

#endif // KRATOS_CREATE_INTERFACE_CONDITIONS_PROCESS_INCLUDED  defined 



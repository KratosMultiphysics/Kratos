//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-10-02 10:47:21 $
//   Revision:            $Revision: 1.8 $
//
//  this process save structural elements in a separate list

#if !defined(KRATOS_SAVE_EXTERNAL_SURFACE_PROCESS_INCLUDED )
#define  KRATOS_SAVE_EXTERNAL_SURFACE_PROCESS_INCLUDED



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

class SaveExternalSurfaceProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(SaveExternalSurfaceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SaveExternalSurfaceProcess()
    {
    }

    /// Destructor.
    virtual ~SaveExternalSurfaceProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

//		void operator()()
//		{
//			SaveStructure();
//		}


    ///@}
    ///@name Operations
    ///@{

    void  SaveExternal(ModelPart& interface_model_part, ModelPart& external_interface_model_part)
    {
        KRATOS_TRY

        //number of structure nodes
        KRATOS_WATCH("SAVING EXTERNAL INTERFACE MODEL PART FOR EMBEDDED COMPUTATIONS")
        for(ModelPart::ElementsContainerType::iterator im = interface_model_part.ElementsBegin() ;
                im != interface_model_part.ElementsEnd() ; ++im)
        {
            //PointerVector<Element> struct_elements_list;
            //check number of structure nodes
            unsigned int n_internal=0; 

            for (unsigned int i=0; i<im->GetGeometry().size(); i++)
            {
                n_internal += int( im->GetGeometry()[i].FastGetSolutionStepValue(FLAG_VARIABLE) );
            }

            if(n_internal==0)
            {
                external_interface_model_part.Elements().push_back(*(im.base()));
            }

        }
	for(ModelPart::NodesContainerType::iterator in = interface_model_part.NodesBegin() ;
                in != interface_model_part.NodesEnd() ; ++in)
        {
		if (int( in->FastGetSolutionStepValue(FLAG_VARIABLE)==0))
			 external_interface_model_part.Nodes().push_back(*(in.base()));
  	
	}
        //WE HAVE TO COPY THE ProcessInfo pointer to the new part, otherwise it is empty
        external_interface_model_part.SetProcessInfo(interface_model_part.pGetProcessInfo());
	external_interface_model_part.SetProperties(interface_model_part.pProperties());

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
    virtual std::string Info() const
    {
        return "SaveExternalSurfaceProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SaveExternalSurfaceProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
    //ModelPart& mr_fluid_model_part;
    //ModelPart& mr_structure_model_part;

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
//		SaveExternalSurfaceProcess& operator=(SaveExternalSurfaceProcess const& rOther);

    /// Copy constructor.
//		SaveExternalSurfaceProcess(SaveExternalSurfaceProcess const& rOther);


    ///@}

}; // Class SaveExternalSurfaceProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SaveExternalSurfaceProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SaveExternalSurfaceProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SAVE_EXTERNAL_SURFACE_PROCESS_INCLUDED  defined 



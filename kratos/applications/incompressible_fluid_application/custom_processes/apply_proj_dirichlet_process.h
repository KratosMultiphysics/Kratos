//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-05-27 12:29:23 $
//   Revision:            $Revision: 1.1 $
//
//  this process save structural elements in a separate list

#if !defined(KRATOS_PROJ_DIRICHLET_PROCESS_INCLUDED )
#define  KRATOS_PROJ_DIRICHLET_PROCESS_INCLUDED


//This process loops over the destination model part elements, and fixes the velocity (applies Dirichlet) to the boundary nodes(on the fictitious side) of the elements, that belong to the interface



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"
#include "incompressible_fluid_application.h"
#include "includes/deprecated_variables.h"

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

class ApplyProjDirichletProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyProjDirichletProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ApplyProjDirichletProcess()
    {
    }

    /// Destructor.
    virtual ~ApplyProjDirichletProcess()
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

void ApplyProjDirichlet(ModelPart& full_model_part)
    {
        KRATOS_TRY
        unsigned int n_int;
        //unsigned int n_fixed;
        
        //first we remove the Dirichlet conditions from the nodes that were defining the interface in the previous step:
        for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ;
                im != full_model_part.ElementsEnd() ; ++im)
        {
            n_int=0.0;
	    //n_fixed=0.0;

	    for (unsigned int i=0; i<im->GetGeometry().size();i++)
		n_int+=im->GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE);

	if (n_int<im->GetGeometry().size() && n_int>0)
	{
            //iterate over the of nodes in the element (interface element)
            for (unsigned int i=0; i<im->GetGeometry().size(); i++)
            {
                n_int=im->GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE);
		//n_fixed=im->GetGeometry()[i].IsFixed(AUX_VEL_X);
		
                // if the node is lying on fictitious side (or its a BAD node with a fixed AUX_VEL) - apply the project Dirichlet condition
                if (n_int==1.0)// || n_fixed==1.0)		
                    //if (((ic->GetGeometry()[i]).GetDof(AUX_VEL_X)).IsFixed())
                {
                    im->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)=im->GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL);
		    //im->GetGeometry()[i].FastGetSolutionStepValue(FRACT_VEL)=im->GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL);		

                    im->GetGeometry()[i].Fix(VELOCITY_X);
                    im->GetGeometry()[i].Fix(VELOCITY_Y);
                    im->GetGeometry()[i].Fix(VELOCITY_Z);
		    /*
                    im->GetGeometry()[i].Fix(FRACT_VEL_X);
                    im->GetGeometry()[i].Fix(FRACT_VEL_Y);
                    im->GetGeometry()[i].Fix(FRACT_VEL_Z);
		    */
                }
            }
	}//end if  (n_int<3 && n_int>0)

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
    virtual std::string Info() const
    {
        return "ApplyProjDirichletProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyProjDirichletProcess";
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
//		ApplyProjDirichletProcess& operator=(ApplyProjDirichletProcess const& rOther);

    /// Copy constructor.
//		ApplyProjDirichletProcess(ApplyProjDirichletProcess const& rOther);


    ///@}

}; // Class ApplyProjDirichletProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyProjDirichletProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyProjDirichletProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_PROJ_DIRICHLET_PROCESS_INCLUDED  defined 




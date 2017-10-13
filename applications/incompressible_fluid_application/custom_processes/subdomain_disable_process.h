//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-05-27 12:29:23 $
//   Revision:            $Revision: 1.1 $
//
//  this process save structural elements in a separate list

#if !defined(KRATOS_SUBDOMAIN_DISABLE_PROCESS_INCLUDED )
#define  KRATOS_SUBDOMAIN_DISABLE_PROCESS_INCLUDED


//This process permits one to "REMOVE" the elements that have all nodes "DISABLED" (flag) from the model part
// this reduced model part will be saved

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
#include "processes/find_nodal_neighbours_process.h"
#include "incompressible_fluid_application.h"

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

class SubdomainDisableProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(SubdomainDisableProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SubdomainDisableProcess()
    {
    }

    /// Destructor.
    virtual ~SubdomainDisableProcess()
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

    void SaveReducedPart(ModelPart& full_model_part, ModelPart& reduced_model_part)
    {
        KRATOS_TRY        
       
	/////////////////////////////////////////////////////////////////////////
	//clear reduced_model_part
        reduced_model_part.Conditions().clear();
        reduced_model_part.Elements().clear();
        reduced_model_part.Nodes().clear();
      
	//change the name of the var to another one  - otherwise confusing
	for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; in != full_model_part.NodesEnd() ; ++in)
	{	  	
	in->FastGetSolutionStepValue(DISABLE)=false;
	}

	int n_nodes=0;

        for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ;
                im != full_model_part.ElementsEnd() ; ++im)
        {
	    int n_int=0;
	    int n_nodes=im->GetGeometry().size();

	    for (int i=0; i<n_nodes; i++)
		n_int+=im->GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE);   

	    //KRATOS_WATCH(n_int);
	    //if the element has at least one non-interface node, save the element and mark all its node with a disable flag (change flag name.. coz in reality its "enabled")
            if (n_int<n_nodes)
            {
                reduced_model_part.AddElement(*(im.base())); 
		for (int i=0;i<n_nodes;i++)
			{
			im->GetGeometry()[i].FastGetSolutionStepValue(DISABLE)=true;
			}               
            }
	    //KRATOS_WATCH(n_int)
	    //KRATOS_WATCH(n_nodes)
            if (n_int>n_nodes)
                KRATOS_THROW_ERROR(std::logic_error,  "Number of DISABLE flags cant exceed number of the element nodes.... " , "");

        }

	for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; in != full_model_part.NodesEnd() ; ++in)
        {

            if (in->FastGetSolutionStepValue(DISABLE)==true)
            {
                reduced_model_part.AddNode(*(in.base()));       
	    }
	}
       
        for(ModelPart::PropertiesContainerType::iterator i_properties = full_model_part.PropertiesBegin() ;
                i_properties != full_model_part.PropertiesEnd() ; ++i_properties)
        {
            reduced_model_part.AddProperties(*(i_properties.base()));

        }
	
	//find neighbors within reduced model part
	if (n_nodes==3)
		{
		FindNodalNeighboursProcess N_FINDER=FindNodalNeighboursProcess(reduced_model_part, 9, 20);
		N_FINDER.Execute();
		}
	else if (n_nodes==4)
		{
		FindNodalNeighboursProcess N_FINDER=FindNodalNeighboursProcess(reduced_model_part, 20, 30);
		N_FINDER.Execute();
		}



	//THIS PART EXCLUDES FROM THE SYSTEM SOLUTION (FOR AUX_VEL) the NODES WHOSE ALL NEIGHBORS ARE FIXED	
	for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ;
                in != full_model_part.NodesEnd() ; ++in)
        {

            int n_int=in->FastGetSolutionStepValue(IS_INTERFACE);

            if (n_int==1.0)
            {
		WeakPointerVector< Node<3> >& neighb_nodes = in->GetValue(NEIGHBOUR_NODES);
                int n_neighbor_nodes=neighb_nodes.size();
		int n_fixed=0;
		array_1d<double,3> av_vel=ZeroVector(3);
		for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin(); i != neighb_nodes.end(); i++)
                {
                    n_fixed+=	i->IsFixed(AUX_VEL_X);
		    av_vel+=i->FastGetSolutionStepValue(AUX_VEL);
                }
		if (n_neighbor_nodes==n_fixed)
			{
			av_vel/=n_neighbor_nodes;
			in->FastGetSolutionStepValue(AUX_VEL)=av_vel;
			in->Fix(AUX_VEL_X);
			in->Fix(AUX_VEL_Y);
			in->Fix(AUX_VEL_Z);
			}
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
    virtual std::string Info() const
    {
        return "SubdomainDisableProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SubdomainDisableProcess";
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
//		SubdomainDisableProcess& operator=(SubdomainDisableProcessconst& rOther);

    /// Copy constructor.
//		SubdomainDisableProcess(SubdomainDisableProcessconst& rOther);


    ///@}

}; // Class SubdomainDisableProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SubdomainDisableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SubdomainDisableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SUBDOMAIN_DISABLE_PROCESS_INCLUDED  defined 



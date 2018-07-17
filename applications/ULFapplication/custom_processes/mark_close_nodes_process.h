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

#if !defined(KRATOS_MARK_CLOSE_NODES_PROCESS_INCLUDED )
#define  KRATOS_MARK_CLOSE_NODES_PROCESS_INCLUDED



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
#include "custom_elements/updated_lagrangian_fluid.h"
#include "custom_elements/updated_lagrangian_fluid3D.h"
#include "custom_elements/updated_lagrangian_fluid_inc.h"
#include "custom_elements/updated_lagrangian_fluid3D_inc.h"


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

class MarkCloseNodesProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(MarkCloseNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MarkCloseNodesProcess(ModelPart& model_part)
        : mr_model_part(model_part)
    {
    }

    /// Destructor.
    ~MarkCloseNodesProcess() override
    {
    }


    ///@}
    ///@name Operators
    ///@{




    ///@}
    ///@name Operations
    ///@{

    void MarkCloseNodes(const double admissible_distance_factor)
    {
        KRATOS_TRY

        double fact2 = admissible_distance_factor*admissible_distance_factor;
        for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin(); in!=mr_model_part.NodesEnd(); in++)
        {
            if(in->FastGetSolutionStepValue(IS_STRUCTURE) == 0) //if it is not a wall node i can erase
            {
                double hnode2 = in->FastGetSolutionStepValue(NODAL_H);
                hnode2 *= hnode2; //take the square

                //loop on neighbours and erase if they are too close
                for( WeakPointerVector< Node<3> >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
                        i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
                {
                    if(i->Is(TO_ERASE) == false) //we can erase the current node only if the neighb is not to be erased
                    {
                        double dx = i->X() - in->X();
                        double dy = i->Y() - in->Y();
                        double dz = i->Z() - in->Z();

                        double dist2 = dx*dx + dy*dy + dz*dz;

                        if(dist2 < fact2 *  hnode2)
                            in->Set(TO_ERASE, true);
                    }
                }
            }
        }
        /*
        this was my old version. now Riccardos version is implemented
        for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin() ;
        		in != mr_model_part.NodesEnd() ; ++in)
        {

        		const double& X0 = in->X();		const double& Y0 = in->Y();
        		KRATOS_WATCH("ENTERED MARKING CLOSE NODES FUCNTION!");

        		for( WeakPointerVector< Node<3> >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
        							i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
        		{
        		const double& X1 = i->X();	const double& Y1 = i->Y();
        		const double& dist_sq = (X1-X0)*(X1-X0)+(Y1-Y0)*(Y1-Y0);
        		//if (dist_sq<(i->GetValue(NODAL_H))*(i->GetValue(NODAL_H)) && in->GetId()>i->GetId())
        		if (dist_sq<0.005 && in->GetId()>i->GetId())
        			{
        			if (i->FastGetSolutionStepValue(IS_STRUCTURE)==false)
        				{
        				i->Is(TO_ERASE)= true;
        				KRATOS_WATCH("ERASING NODE!!!!!!");
        				KRATOS_WATCH(in->GetId());
        				}

        			}

        		}


        }
        */

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
        return "MarkCloseNodesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MarkCloseNodesProcess";
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
//		MarkCloseNodesProcess& operator=(MarkCloseNodesProcess const& rOther);

    /// Copy constructor.
//		MarkCloseNodesProcess(MarkCloseNodesProcess const& rOther);


    ///@}

}; // Class MarkCloseNodesProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MarkCloseNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MarkCloseNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MARK_CLOSE_NODES_PROCESS_INCLUDED  defined 



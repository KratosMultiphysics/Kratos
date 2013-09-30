//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_FRONT_MEETING_RECOGNITION_PROCESS_H_INCLUDED )
#define  KRATOS_FRONT_MEETING_RECOGNITION_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "processes/find_nodal_neighbours_process.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"


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

/// This Process takes a mesh and extract out the individual bodies by analysing the connectivity.
/** .
*/
class FrontMeetingRecognitionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FrontMeetingRecognitionProcess
    KRATOS_CLASS_POINTER_DEFINITION(FrontMeetingRecognitionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FrontMeetingRecognitionProcess(ModelPart& rModelPart) : mrModelPart(rModelPart)
    {
        FindNodalNeighboursProcess find_nodal_neighbours_process(mrModelPart);
        find_nodal_neighbours_process.Execute();
		mNeighbourNodes.resize(mrModelPart.NumberOfNodes());

		//int size = 0;
		
		 ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
		for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
        {
			mNeighbourNodes[i_node->Id() - 1] = i_node->GetValue(NEIGHBOUR_NODES);
			i_node->GetSolutionStepValue(FRONT_MEETING) = 0;
			//size += mNeighbourNodes[i_node->Id() - 1].size();
        }
		//std::cout << "average number of neighbours : " << size / mrModelPart.NumberOfNodes() << std::endl;
   }

    /// Copy constructor.
    FrontMeetingRecognitionProcess(FrontMeetingRecognitionProcess const& rOther) : mrModelPart(rOther.mrModelPart)
    {
    }

    /// Destructor.
    virtual ~FrontMeetingRecognitionProcess() {}


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

    virtual void Execute()
    {
        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();

		for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
			i_node->Set(NOT_VISITED);

        SelectFrontNodes();

		int front_meeting_count = 0;
        for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
        {
            if(i_node->Is(SELECTED))
					front_meeting_count += FindFrontMeeting(i_node);
        }
		if(front_meeting_count){
						std::cout<< std::endl;
						std::cout<< "#######################                   #########################" << std::endl;
						std::cout<< "############################         ##############################" << std::endl;
						std::cout<< "###############################   #################################" << std::endl;
						std::cout<< "################################ ##################################" << std::endl;
						std::cout<< "################################ ##################################" << std::endl;
						std::cout<< "###############################   #################################" << std::endl;
						std::cout<< "############################         ##############################" << std::endl;
						std::cout<< "#######################                   #########################" << std::endl;
				std::cout << front_meeting_count << " front meeting point founded!!" << std::endl;
		}
		//std::cout << current_body_id << " clusters found" << std::endl;
    }

	void SelectFrontNodes(void)
	{
        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();

		for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
            i_node->Set(NOT_SELECTED);

        for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
        {
            if(i_node->GetSolutionStepValue(DISTANCE) > 0.00)  // if the node is air
			{
				WeakPointerVector<Node<3> >& r_neighbour_nodes = mNeighbourNodes[i_node->Id() - 1];
				for(WeakPointerVector<Node<3> >::iterator i_neighbour_node = r_neighbour_nodes.begin() ; i_neighbour_node != r_neighbour_nodes.end() ; i_neighbour_node++)
				{
					if(i_neighbour_node->GetSolutionStepValue(DISTANCE) < 0.00)
					{
						i_node->Set(SELECTED);
						break;
					}
				}
			}
        }
	}

    int FindFrontMeeting(ModelPart::NodesContainerType::iterator iNode)
    {
        ModelPart::NodesContainerType front_nodes;
        WeakPointerVector<Node<3> >& r_neighbour_nodes = mNeighbourNodes[iNode->Id() - 1];

//		static int counter = 0;
		bool is_candidate = true;
        for(WeakPointerVector<Node<3> >::iterator i_neighbour_node = r_neighbour_nodes.begin() ; i_neighbour_node != r_neighbour_nodes.end() ; i_neighbour_node++)
        {
            if(i_neighbour_node->Is(NOT_SELECTED)) // Is not on the front
            {
				if(i_neighbour_node->GetSolutionStepValue(DISTANCE) > 0.00) // The neighbour is air
				{
					is_candidate = false;
					break;
				}
				else
				{
					front_nodes.push_back((*(i_neighbour_node.base())).lock());
				}
            }
        }

		
		
		if(is_candidate)
		{
//			std::size_t front_nodes_size = front_nodes.size();
			for(ModelPart::NodesContainerType::iterator i_node = front_nodes.begin() ; i_node != front_nodes.end() ; i_node++)
				for(ModelPart::NodesContainerType::iterator j_node = i_node + 1 ; j_node != front_nodes.end() ; j_node++)
				{
					array_1d<double, 3> vi = i_node->GetSolutionStepValue(VELOCITY);
					array_1d<double, 3> vj = j_node->GetSolutionStepValue(VELOCITY);
					double cos_alpha = inner_prod(vi,vj) / (norm_2(vi) * norm_2(vj));
//					if(cos_alpha < -0.7 ) // angle > 135
					if(cos_alpha < -0.1 ) // angle > 90
					{
						iNode->GetSolutionStepValue(FRONT_MEETING) = 1;
						return 1;
					}
				}
				
		}
	
		return 0;
	
	}

 
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "FrontMeetingRecognitionProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

    ModelPart& mrModelPart;
	std::vector<WeakPointerVector<Node<3> > > mNeighbourNodes;

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
    FrontMeetingRecognitionProcess& operator=(FrontMeetingRecognitionProcess const& rOther)
    {
        return *this;
    }


    ///@}

}; // Class FrontMeetingRecognitionProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FrontMeetingRecognitionProcess& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FrontMeetingRecognitionProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FRONT_MEETING_RECOGNITION_PROCESS_H_INCLUDED  defined 



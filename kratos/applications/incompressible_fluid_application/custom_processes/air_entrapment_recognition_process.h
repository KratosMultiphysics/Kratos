//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_AIR_ENTRAPMENT_RECOGNITION_PROCESS_H_INCLUDED )
#define  KRATOS_AIR_ENTRAPMENT_RECOGNITION_PROCESS_H_INCLUDED



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
class AirEntrapmentRecognitionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AirEntrapmentRecognitionProcess
    KRATOS_CLASS_POINTER_DEFINITION(AirEntrapmentRecognitionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AirEntrapmentRecognitionProcess(ModelPart& rModelPart) : mrModelPart(rModelPart)
    {
        FindNodalNeighboursProcess find_nodal_neighbours_process(mrModelPart);
        find_nodal_neighbours_process.Execute();
		mNeighbourNodes.resize(mrModelPart.NumberOfNodes());

		//int size = 0;
		
		 ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
		for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
        {
			mNeighbourNodes[i_node->Id() - 1] = i_node->GetValue(NEIGHBOUR_NODES);
			i_node->GetSolutionStepValue(LAST_AIR) = 0;
			//size += mNeighbourNodes[i_node->Id() - 1].size();
        }
		//std::cout << "average number of neighbours : " << size / mrModelPart.NumberOfNodes() << std::endl;
   }

    /// Copy constructor.
    AirEntrapmentRecognitionProcess(AirEntrapmentRecognitionProcess const& rOther) : mrModelPart(rOther.mrModelPart)
    {
    }

    /// Destructor.
    virtual ~AirEntrapmentRecognitionProcess() {}


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

        //int current_body_id = 1;
        for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
        {
            if(i_node->GetSolutionStepValue(DISTANCE) > 0.00)  // if the node is air
				if(i_node->IsNot(VISITED))  // if is not visited
					FindAirEntrapment(i_node);
        }
		//std::cout << current_body_id << " clusters found" << std::endl;
    }

    void FindAirEntrapment(ModelPart::NodesContainerType::iterator iNode)
    {
        ModelPart::NodesContainerType front_nodes;
        WeakPointerVector<Node<3> >& r_neighbour_nodes = mNeighbourNodes[iNode->Id() - 1];
        ModelPart::NodesContainerType buble_nodes;

        for(WeakPointerVector<Node<3> >::iterator i_neighbour_node = r_neighbour_nodes.begin() ; i_neighbour_node != r_neighbour_nodes.end() ; i_neighbour_node++)
        {
            if(i_neighbour_node->GetSolutionStepValue(DISTANCE) > 0.00) // The neighbour is air
            {
				if(i_neighbour_node->IsNot(VISITED))
				{
					i_neighbour_node->Set(VISITED);
					buble_nodes.push_back((*(i_neighbour_node.base())).lock());
					front_nodes.push_back((*(i_neighbour_node.base())).lock());
				}
            }
        }
		
        while(!front_nodes.empty())
        {
            ModelPart::NodesContainerType new_front_nodes;
            for(ModelPart::NodesContainerType::iterator i_node = front_nodes.begin() ; i_node != front_nodes.end() ; i_node++)
            {
				WeakPointerVector<Node<3> >& r_neighbour_nodes = mNeighbourNodes[i_node->Id() - 1];
				for(WeakPointerVector<Node<3> >::iterator i_neighbour_node = r_neighbour_nodes.begin() ; i_neighbour_node != r_neighbour_nodes.end() ; i_neighbour_node++)
				{
					if(i_neighbour_node->GetSolutionStepValue(DISTANCE) > 0.00) // The neighbour is air
					{
						if(i_neighbour_node->IsNot(VISITED))
						{
							i_neighbour_node->Set(VISITED);
							buble_nodes.push_back((*(i_neighbour_node.base())).lock());
							new_front_nodes.push_back((*(i_neighbour_node.base())).lock());
						}
					}
				}
            }
            front_nodes = new_front_nodes;
        }

		int buble_nodes_size = buble_nodes.size();
		if((buble_nodes_size > 3) && ( buble_nodes_size < 8))
		{
			std::cout<< "#######################################################################################" << std::endl;
			std::cout<< "#######################################################################################" << std::endl;
			std::cout<< "#########################################  ############################################" << std::endl;
			std::cout<< "######################################        #########################################" << std::endl;
			std::cout<< "#####################################   buble  ########################################" << std::endl;
			std::cout<< "####################################  with size #######################################" << std::endl;
			std::cout<< "####################################     " << buble_nodes_size << "      #######################################" << std::endl;
			std::cout<< "#####################################          ########################################" << std::endl;
			std::cout<< "######################################        #########################################" << std::endl;
			std::cout<< "#########################################  ############################################" << std::endl;
			std::cout<< "#######################################################################################" << std::endl;
			std::cout<< "#######################################################################################" << std::endl;
			for(ModelPart::NodesContainerType::iterator i_node = buble_nodes.begin() ; i_node != buble_nodes.end() ; i_node++)
				i_node->GetSolutionStepValue(LAST_AIR) = 1;
		}
	
	}

 
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "AirEntrapmentRecognitionProcess";
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
    AirEntrapmentRecognitionProcess& operator=(AirEntrapmentRecognitionProcess const& rOther)
    {
        return *this;
    }


    ///@}

}; // Class AirEntrapmentRecognitionProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AirEntrapmentRecognitionProcess& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AirEntrapmentRecognitionProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_AIR_ENTRAPMENT_RECOGNITION_PROCESS_H_INCLUDED  defined 



//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_SHRINKAGE_POROSITY_CALCULATION_PROCESS_H_INCLUDED )
#define  KRATOS_SHRINKAGE_POROSITY_CALCULATION_PROCESS_H_INCLUDED



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
class ShrinkagePorosityCalculationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ShrinkagePorosityCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(ShrinkagePorosityCalculationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ShrinkagePorosityCalculationProcess(ModelPart& rModelPart) : mrModelPart(rModelPart)
    {
        FindNodalNeighboursProcess find_nodal_neighbours_process(mrModelPart);
        find_nodal_neighbours_process.Execute();
		mNeighbourNodes.resize(mrModelPart.NumberOfNodes());
		mNodalVolumes.resize(mrModelPart.NumberOfNodes(), 0.00);
		mNodalShrinkage.resize(mrModelPart.NumberOfNodes(), 0.00);
//		mClustersShrinkage.resize(1,0.00);
//		mNodeClusterIndex.resize(mrModelPart.NumberOfNodes(), 0);

		int size = 0;
		
		 ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
		for(ModelPart::NodeIterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
        {
			mNeighbourNodes[i_node->Id() - 1] = i_node->GetValue(NEIGHBOUR_NODES);
			size += mNeighbourNodes[i_node->Id() - 1].size();
        }

		for(ModelPart::ElementIterator i_element = mrModelPart.ElementsBegin() ; i_element != mrModelPart.ElementsEnd() ; i_element++)
		{
			Element::GeometryType& geometry = i_element->GetGeometry();
			const double volume = geometry.Volume();
			const int size = geometry.size();
			const double nodal_volume = volume / size;
			for(int i = 0 ; i < size ; i++)
			{
				mNodalVolumes[geometry[i].Id() - 1] += nodal_volume;
			}
		}

		std::cout << "average number of neighbours : " << size / mrModelPart.NumberOfNodes() << std::endl;
   }

    /// Copy constructor.
    ShrinkagePorosityCalculationProcess(ShrinkagePorosityCalculationProcess const& rOther) : mrModelPart(rOther.mrModelPart)
    {
    }

    /// Destructor.
    virtual ~ShrinkagePorosityCalculationProcess() {}


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

		int current_cluster_id = 0;
        for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
        {
            if(i_node->GetSolutionStepValue(SOLID_FRACTION, 1) < 1.00) // if was fluid in previous step
				if(i_node->IsNot(VISITED))  // if is not visited
					CalculateClusterShrinkage(i_node, current_cluster_id++);
        }
		std::cout << current_cluster_id << " clusters found" << std::endl;
    }

   void CalculateClusterShrinkage(ModelPart::NodesContainerType::iterator iNode, int BodyId)
    {
        ModelPart::NodesContainerType front_nodes;
        ModelPart::NodesContainerType cluster_nodes;
        WeakPointerVector<Node<3> >& r_neighbour_nodes = mNeighbourNodes[iNode->Id() - 1];
		double cluster_shrinkage = 0.00;
		const double shrinkage_factor = 0.07;

        for(WeakPointerVector<Node<3> >::iterator i_neighbour_node = r_neighbour_nodes.begin() ; i_neighbour_node != r_neighbour_nodes.end() ; i_neighbour_node++)
        {
            if(i_neighbour_node->GetSolutionStepValue(SOLID_FRACTION, 1) < 1.00) // The neighbour was fluid
            {
				if(i_neighbour_node->IsNot(VISITED))
				{
					i_neighbour_node->Set(VISITED);
					//i_neighbour_node->GetSolutionStepValue(MACRO_POROSITY) = BodyId;
					front_nodes.push_back((*(i_neighbour_node.base())).lock());
					cluster_nodes.push_back((*(i_neighbour_node.base())).lock());
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
					if(i_neighbour_node->GetSolutionStepValue(SOLID_FRACTION, 1) < 1.00) // The neighbour was fluid
					{
						if(i_neighbour_node->IsNot(VISITED))
						{
							i_neighbour_node->Set(VISITED);
							//i_neighbour_node->GetSolutionStepValue(MACRO_POROSITY) = BodyId;
							new_front_nodes.push_back((*(i_neighbour_node.base())).lock());
							cluster_nodes.push_back((*(i_neighbour_node.base())).lock());
						}
					}
				}
            }
            front_nodes = new_front_nodes;
        }
		
        int cluster_fluid_nodes = 0;
        for(ModelPart::NodesContainerType::iterator i_node = cluster_nodes.begin() ; i_node != cluster_nodes.end() ; i_node++)
        {
			if(i_node->GetSolutionStepValue(SOLID_FRACTION) < 1.00) // the node still is fluid
			{
				cluster_fluid_nodes++;
			}
			else // node is solidified and we add its shrinkage to the cluster
			{
				cluster_shrinkage += mNodalVolumes[i_node->Id() - 1] * shrinkage_factor;
			}
		}

		
        for(ModelPart::NodesContainerType::iterator i_node = cluster_nodes.begin() ; i_node != cluster_nodes.end() ; i_node++)
        {
			mNodalShrinkage[i_node->Id() - 1] += cluster_shrinkage;
		}

		if(cluster_fluid_nodes < 8)
		{
			for(ModelPart::NodesContainerType::iterator i_node = cluster_nodes.begin() ; i_node != cluster_nodes.end() ; i_node++)
			{
				i_node->GetSolutionStepValue(MACRO_POROSITY) = mNodalShrinkage[i_node->Id() - 1];
			}
		}


			KRATOS_WATCH(cluster_nodes.size());
		
    }

 //   virtual void Execute1()
 //   {
	//	CalculateCurrentShrinkage();



 //       ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();

 //		for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
	//		i_node->Set(NOT_VISITED);

	//	int current_cluster_id = 0;
 //       for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
 //       {
 //           if(i_node->GetSolutionStepValue(SOLID_FRACTION) < 1.00)
	//			if(i_node->IsNot(VISITED))  // if is not visited
	//				AssignBody(i_node, current_cluster_id++);
 //       }
	//	mClustersShrinkage.resize(current_cluster_id,0.00);
	//	std::cout << current_cluster_id << " clusters found" << std::endl;
 //   }

	//void CalculateCurrentShrinkage()
	//{
	//	const double shrinkage_factor = 0.07;

 //       ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();

 //		for(ModelPart::NodesContainerType::iterator i_node = r_nodes.begin(); i_node!=r_nodes.end(); i_node++)
	//	{
	//		if(i_node->GetSolutionStepValue(SOLID_FRACTION) > 0.99) // is solidified
	//			if(i_node->GetSolutionStepValue(SOLID_FRACTION,1) <= 0.99) // was not solidified in previous step
	//			{
	//				std::size_t cluster_index = mNodeClusterIndex[i_node->Id() - 1];
	//				mClustersShrinkage[cluster_index] += mNodalVolumes[i_node->Id() - 1] * shrinkage_factor;
	//			}
	//	}

	//}


 //   void AssignBody(ModelPart::NodesContainerType::iterator iNode, int BodyId)
 //   {
 //       ModelPart::NodesContainerType front_nodes;
 //       ModelPart::NodesContainerType cluster_nodes;
 //       WeakPointerVector<Node<3> >& r_neighbour_nodes = mNeighbourNodes[iNode->Id() - 1];
	//	double cluster_volume = 0.00;

 //       for(WeakPointerVector<Node<3> >::iterator i_neighbour_node = r_neighbour_nodes.begin() ; i_neighbour_node != r_neighbour_nodes.end() ; i_neighbour_node++)
 //       {
 //           if(i_neighbour_node->GetSolutionStepValue(SOLID_FRACTION) < 1.00) // The neighbour is fluid
 //           {
	//			if(i_neighbour_node->IsNot(VISITED))
	//			{
	//				i_neighbour_node->Set(VISITED);
	//				//i_neighbour_node->GetSolutionStepValue(MACRO_POROSITY) = BodyId;
	//				front_nodes.push_back((*(i_neighbour_node.base())).lock());
	//				cluster_nodes.push_back((*(i_neighbour_node.base())).lock());
	//			}
 //           }
 //       }
 //       while(!front_nodes.empty())
 //       {
 //           ModelPart::NodesContainerType new_front_nodes;
 //           for(ModelPart::NodesContainerType::iterator i_node = front_nodes.begin() ; i_node != front_nodes.end() ; i_node++)
 //           {
	//			cluster_volume += mNodalVolumes[i_node->Id() - 1];
	//			WeakPointerVector<Node<3> >& r_neighbour_nodes = mNeighbourNodes[i_node->Id() - 1];
	//			for(WeakPointerVector<Node<3> >::iterator i_neighbour_node = r_neighbour_nodes.begin() ; i_neighbour_node != r_neighbour_nodes.end() ; i_neighbour_node++)
	//			{
	//				if(i_neighbour_node->GetSolutionStepValue(SOLID_FRACTION) < 1.00) // The neighbour is fluid
	//				{
	//					if(i_neighbour_node->IsNot(VISITED))
	//					{
	//						i_neighbour_node->Set(VISITED);
	//						//i_neighbour_node->GetSolutionStepValue(MACRO_POROSITY) = BodyId;
	//						new_front_nodes.push_back((*(i_neighbour_node.base())).lock());
	//						cluster_nodes.push_back((*(i_neighbour_node.base())).lock());
	//					}
	//				}
	//			}
 //           }
 //           front_nodes = new_front_nodes;
 //       }
	//	
 //           for(ModelPart::NodesContainerType::iterator i_node = cluster_nodes.begin() ; i_node != cluster_nodes.end() ; i_node++)
 //           {
	//			std::size_t cluster_index = mNodeClusterIndex[i_node->Id() - 1];					
	//			mNodalShrinkage[i_node->Id() - 1] += mClustersShrinkage[cluster_index];
	//		}
 //

	//		//if(cluster_nodes.size() < 8)  // cluster is almost solidified
	//		{
	//			for(ModelPart::NodesContainerType::iterator i_node = cluster_nodes.begin() ; i_node != cluster_nodes.end() ; i_node++)
	//			{
	//				std::size_t cluster_index = mNodeClusterIndex[i_node->Id() - 1];					
	//				i_node->GetSolutionStepValue(MACRO_POROSITY) = mNodalShrinkage[i_node->Id() - 1];
	//			}
	//		}

	//		for(ModelPart::NodesContainerType::iterator i_node = cluster_nodes.begin() ; i_node != cluster_nodes.end() ; i_node++)
 //           {
	//			mNodeClusterIndex[i_node->Id() - 1] = BodyId;					
	//		}

	//		KRATOS_WATCH(cluster_nodes.size());
	//	
 //   }



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
        return "ShrinkagePorosityCalculationProcess";
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
	std::vector<double> mNodalVolumes;
	std::vector<double> mNodalShrinkage;
	//std::vector<double> mClustersShrinkage;
	//std::vector<std::size_t> mNodeClusterIndex;

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
    ShrinkagePorosityCalculationProcess& operator=(ShrinkagePorosityCalculationProcess const& rOther)
    {
        return *this;
    }


    ///@}

}; // Class ShrinkagePorosityCalculationProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ShrinkagePorosityCalculationProcess& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ShrinkagePorosityCalculationProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SHRINKAGE_POROSITY_CALCULATION_PROCESS_H_INCLUDED  defined 



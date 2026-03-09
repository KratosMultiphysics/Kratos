//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2007-11-06 12:34:26 $
//   Revision:            $Revision: 1.4 $
//
//  this process save structural elements in a separate list

#if !defined(KRATOS_ADD_WALL_PROCESS_INCLUDED )
#define  KRATOS_ADD_WALL_PROCESS_INCLUDED



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
//#include "custom_elements/updated_lagrangian_fluid.h"
//#include "custom_elements/updated_lagrangian_fluid3D.h"
//#include "custom_elements/updated_lagrangian_fluid_inc.h"
//#include "custom_elements/updated_lagrangian_fluid3D_inc.h"


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

class AddWallProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(AddWallProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AddWallProcess()
    //ModelPart& fluid_model_part, ModelPart& structure_model_part, ModelPart& combined_model_part)
    //: mr_fluid_model_part(fluid_model_part), mr_structure_model_part(structure_model_part), mr_combined_model_part(combined_model_part)
    {
	//KRATOS_WATCH(" INSIDE ADD WALL NODES CONSTRUCTOR") 
    }

    /// Destructor.
    ~AddWallProcess() override
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

    void AddWall(ModelPart& fluid_model_part, ModelPart& wall_model_part)
    {
        KRATOS_TRY	

	//COMPUTING AVERAGE NODAL_H in FLUID
	double av_mesh_size=0.0;
	for(ModelPart::NodesContainerType::iterator in = fluid_model_part.NodesBegin() ;
                in != fluid_model_part.NodesEnd() ; ++in)
        {
	av_mesh_size+=in->FastGetSolutionStepValue(NODAL_H);
	}

	unsigned int n_fluid_nodes=fluid_model_part.Nodes().size();
	av_mesh_size/=n_fluid_nodes;

	unsigned wall_node_id=n_fluid_nodes+1;

	if (av_mesh_size==0)
		KRATOS_THROW_ERROR(std::logic_error,"your wall nodes will have NODAL_H=0","");
	//////////////////////////////////////////////////////////////////////////////////////////
        for(ModelPart::NodesContainerType::iterator in = wall_model_part.NodesBegin() ;
                in != wall_model_part.NodesEnd() ; ++in)
        {
		//set the NODAL_H in the wall node to the average size and then add this node to the fluid model part
		in->FastGetSolutionStepValue(NODAL_H)=av_mesh_size;
		in->SetId(wall_node_id);
		wall_node_id++;
		//ADDING ALL THE SECOND MOULD WALL NODES TO THE MODEL PART
		fluid_model_part.AddNode(*(in.base()),0);
	    
		
        }
	
        //fluid_model_part.AddNodes(wall_model_part.NodesBegin(), wall_model_part.NodesEnd());	
        //sorting and renumbering the fluid elements
        unsigned int id=1;
        for(ModelPart::NodesContainerType::iterator in = fluid_model_part.NodesBegin() ;
                in != fluid_model_part.NodesEnd() ; ++in)
        {
            in->SetId(id);
            id++;
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
    std::string Info() const override
    {
        return "AddWallProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AddWallProcess";
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
    //ModelPart& mr_fluid_model_part;
    //ModelPart& mr_structure_model_part;
    //ModelPart& mr_combined_model_part;

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
//		AddWallProcess& operator=(AddWallProcess const& rOther);

    /// Copy constructor.
//		AddWallProcess(AddWallProcess const& rOther);


    ///@}

}; // Class AddWallProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AddWallProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AddWallProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ADD_WALL_PROCESS_INCLUDED  defined 



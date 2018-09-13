//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2007-11-06 12:34:26 $
//   Revision:            $Revision: 1.4 $
//
//  this process save structural elements in a separate list

#if !defined(KRATOS_REMOVE_SAVE_WALL_PROCESS_INCLUDED )
#define  KRATOS_REMOVE_SAVE_WALL_PROCESS_INCLUDED



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

class RemoveAndSaveWallNodesProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(RemoveAndSaveWallNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RemoveAndSaveWallNodesProcess()
    //ModelPart& fluid_model_part, ModelPart& structure_model_part, ModelPart& combined_model_part)
    //: mr_fluid_model_part(fluid_model_part), mr_structure_model_part(structure_model_part), mr_combined_model_part(combined_model_part)
    {
	//KRATOS_WATCH(" INSIDE REMOVE AND SAVE WALL CONSTRUCTOR") 
    }

    /// Destructor.
    ~RemoveAndSaveWallNodesProcess() override
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

    void RemoveAndSave(ModelPart& fluid_model_part, ModelPart& wall_model_part)
    {
        KRATOS_TRY
	//ModelPart fluid_only_model_part;
       
        for(ModelPart::NodesContainerType::iterator in = fluid_model_part.NodesBegin() ;
                in != fluid_model_part.NodesEnd() ; ++in)
	//in a two stage process I distinguish the second wall by the FLAG_VARIABLE=5
        {
		//second mould nodes are marked with FLAG=5
	if (in->FastGetSolutionStepValue(FLAG_VARIABLE)==5)
		{
		wall_model_part.AddNode(*(in.base()),0);		
		}
        }

	wall_model_part.SetProcessInfo(fluid_model_part.pGetProcessInfo());

	
        //removing second mould nodes (wall nodes) from fluid_model_part
	for(ModelPart::NodesContainerType::iterator in = wall_model_part.NodesBegin() ;
                in != wall_model_part.NodesEnd() ; ++in)
        {
            unsigned int id=in->GetId();
            fluid_model_part.RemoveNode(id);
        }


	unsigned int id=1;
        for(ModelPart::NodesContainerType::iterator in = fluid_model_part.NodesBegin() ;
                in != fluid_model_part.NodesEnd() ; ++in)
        {
            in->SetId(id);
            id++;
        }
	
	//wall_nodes_id  must still be reset when the wall nodes are added
        for(ModelPart::NodesContainerType::iterator in = wall_model_part.NodesBegin() ;
                in != wall_model_part.NodesEnd() ; ++in)
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
        return "RemoveAndSaveWallNodesProcess";
    }

    /// Print information about this object. 
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RemoveAndSaveWallNodesProcess";
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
//		RemoveAndSaveWallNodesProcess& operator=(RemoveAndSaveWallNodesProcess const& rOther);

    /// Copy constructor.
//		RemoveAndSaveWallNodesProcess(RemoveAndSaveWallNodesProcess const& rOther);


    ///@}

}; // Class RemoveAndSaveWallNodesProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  RemoveAndSaveWallNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RemoveAndSaveWallNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REMOVE_SAVE_WALL_PROCESS_INCLUDED  defined 



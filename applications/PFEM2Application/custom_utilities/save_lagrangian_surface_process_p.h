//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti.
//

//  this process saves the boundary of a Lagrangian part for EMBEDDED technique

#if !defined(KRATOS_SAVE_LAGRANGIAN_SURFACE_PROCESS_P_INCLUDED )
#define  KRATOS_SAVE_LAGRANGIAN_SURFACE_PROCESS_P_INCLUDED



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

class SaveLagrangianSurfaceProcess_p
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(SaveLagrangianSurfaceProcess_p);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SaveLagrangianSurfaceProcess_p()
    {
    }

    /// Destructor.
    virtual ~SaveLagrangianSurfaceProcess_p()
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

	//FLUID_MODEL_PART is the whole Lagrangian fluid domain, surface model_part is its boundary
    void  SaveSurfaceConditions_p(ModelPart& lagrangian_model_part, ModelPart& surface_model_part)
    {
        KRATOS_TRY
        surface_model_part.Elements().clear();
        surface_model_part.Conditions().clear();
        surface_model_part.Nodes().clear();

	KRATOS_WATCH("SaveSurfaceConditions");

        //ModelPart::IndexType default_index = 0;

	for(ModelPart::NodesContainerType::iterator in = lagrangian_model_part.NodesBegin() ;
                in != lagrangian_model_part.NodesEnd() ; ++in)
        {
	    if (in->FastGetSolutionStepValue(IS_BOUNDARY)==1.0 && in->GetValue(NEIGHBOUR_CONDITIONS).size()!=0)
		    surface_model_part.Nodes().push_back(*(in.base()));

        }
	int id = 1;
	Properties::Pointer properties = lagrangian_model_part.GetMesh().pGetProperties(1);
        for(ModelPart::ConditionsContainerType::iterator ic = lagrangian_model_part.ConditionsBegin() ; ic != lagrangian_model_part.ConditionsEnd() ; ++ic)
        {
	    //surface_model_part.Conditions().push_back(*(ic.base()));
	    Geometry< Node<3> >& geom = ic->GetGeometry();

            std::vector<std::size_t> NodeIds(3);

	    for(int i=0; i<3; i++) NodeIds[i]= geom[i].Id();

	    surface_model_part.CreateNewElement("Element3D3N",id, NodeIds,0);//,default_index);

	    id = id + 1;
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
    virtual std::string Info() const override
    {
        return "SaveLagrangianSurfaceProcess_p";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SaveLagrangianSurfaceProcess_p";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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
//		SaveLagrangianSurfaceProcess& operator=(SaveLagrangianSurfaceProcess const& rOther);

    /// Copy constructor.
//		SaveLagrangianSurfaceProcess(SaveLagrangianSurfaceProcess const& rOther);


    ///@}

}; // Class SaveLagrangianSurfaceProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SaveLagrangianSurfaceProcess_p& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SaveLagrangianSurfaceProcess_p& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SAVE_LAGRANGIAN_SURFACE_PROCESS_INCLUDED  defined

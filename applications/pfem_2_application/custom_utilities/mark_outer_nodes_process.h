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

#if !defined(KRATOS_MARK_OUTER_NODES_PROCESS_INCLUDED )
#define  KRATOS_MARK_OUTER_NODES_PROCESS_INCLUDED



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
//#include "custom_elements/updated_lagrangian_fluid.h"


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

class MarkOuterNodesProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(MarkOuterNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MarkOuterNodesProcess(ModelPart& model_part)
        : mr_model_part(model_part)
    {
    }

    /// Destructor.
    virtual ~MarkOuterNodesProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    void MarkOuterNodes(const array_1d<double,3>& corner1, const array_1d<double,3>& corner2)
    {
        KRATOS_TRY
        //add a big number to the id of the nodes to be erased
        int n_erased = 0;
        double xmax, xmin, ymax,ymin, zmax, zmin;

        if(corner1[0] > corner2[0])
        {
            xmax = corner1[0];
            xmin = corner2[0];
        }
        else
        {
            xmax = corner2[0];
            xmin = corner1[0];
        }
        if(corner1[1] > corner2[1])
        {
            ymax = corner1[1];
            ymin = corner2[1];
        }
        else
        {
            ymax = corner2[1];
            ymin = corner1[1];
        }
        if(corner1[2] > corner2[2])
        {
            zmax = corner1[2];
            zmin = corner2[2];
        }
        else
        {
            zmax = corner2[2];
            zmin = corner1[2];
        }
        //for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
        for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin() ;
                in != mr_model_part.NodesEnd() ; ++in)
        {
            bool erase = false;
            double& x = in->X();
            double& y = in->Y();
            double& z = in->Z();

            if(x<xmin || x>xmax)		erase = true;
            else if(y<ymin || y>ymax)	erase = true;
            else if(z<zmin || z>zmax)	erase = true;

            if(erase == true)
            {
                n_erased += 1;
                in->Set(TO_ERASE, true);
                std::cout << "erasing outer node at " << in->X() << " " << in->Y() << " " << in->Z() << std::endl;
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
    virtual std::string Info() const override
    {
        return "MarkOuterNodesProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MarkOuterNodesProcess";
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
//		MarkOuterNodesProcess& operator=(MarkOuterNodesProcess const& rOther);

    /// Copy constructor.
//		MarkOuterNodesProcess(MarkOuterNodesProcess const& rOther);


    ///@}

}; // Class MarkOuterNodesProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MarkOuterNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MarkOuterNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MARK_OUTER_NODES_PROCESS_INCLUDED  defined

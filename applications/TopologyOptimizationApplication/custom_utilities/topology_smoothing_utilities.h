//
//   Project Name:        $Project:     Topology_Optimization_Application $
//   Last modified by:    $Author:      Malfavón Farías, Baumgärtner      $
//   Date:                $Date:        May 2016                          $
//   Revision:            $Revision:    0.0                               $
//


#if !defined(KRATOS_TOPOLOGY_SMOOTHING_UTILITIES_H_INCLUDED)
#define  KRATOS_TOPOLOGY_SMOOTHING_UTILITIES_H_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <algorithm>

// External includes
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/process_info.h"

// Application includes
#include "topology_optimization_application.h"
#include "processes/find_conditions_neighbours_process.h" // To find node neighbours using conditions


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

/// Solution utility that smooths a provided meshZz.
/** Detail class definition.

*/

class TopologySmoothingUtilities
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of TopologySmoothingUtilities
    KRATOS_CLASS_POINTER_DEFINITION(TopologySmoothingUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TopologySmoothingUtilities( ModelPart& model_part )
        : mr_model_part(model_part)
    {
    }

    /// Destructor.
    virtual ~TopologySmoothingUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------- SMOOTH EXTRACTED MESH  ------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

    /// Gets the neighbour nodes and applies a Laplacian algorithm to smooth the previously extracted surface mesh

    void SmoothMesh( double iterations )
    {

        KRATOS_TRY;

            std::cout<<"::[Smoothing mesh]::"<<std::endl;

            if (iterations > 0)
            {
                // Locating neighbouring nodes
                FindConditionsNeighboursProcess nodal_finder = FindConditionsNeighboursProcess(mr_model_part, 10, 10);
                nodal_finder.Execute();

                NodesContainerType& rNodes = mr_model_part.Nodes();

                // Repeat smoothing operation for the selected number of iterations
                for(int i = 0; i < iterations; ++i){
                    std::cout<<"  Smoothing iteration number "<< i+1 <<std::endl;

                    for(NodesContainerType::iterator node_i = rNodes.begin(); node_i!=rNodes.end(); node_i++)
                    {
                        // Prepare the needed neighbouring node vector and the needed variables
                        double x_coord = 0.0;
                        double y_coord = 0.0;
                        double z_coord = 0.0;

                        WeakPointerVector< Node<3> >& neighbours = node_i->GetValue(NEIGHBOUR_NODES);

                        // Obtain and sum the X, Y and Z coordinates of all neighbouring nodes
                        for( WeakPointerVector<Node<3> >::iterator neighbour_node = neighbours.begin(); neighbour_node!=neighbours.end(); neighbour_node++)
                        {
                            x_coord += neighbour_node->X();
                            y_coord += neighbour_node->Y();
                            z_coord += neighbour_node->Z();
                        }

                        // Average the new X, Y and Z coordinates and save them temporary (simultaneous version)
                        node_i->SetValue(ELEM_CENTER_X, x_coord /= neighbours.size());
                        node_i->SetValue(ELEM_CENTER_Y, y_coord /= neighbours.size());
                        node_i->SetValue(ELEM_CENTER_Z, z_coord /= neighbours.size());
                    }

                    // Assign the coordinates calculated value into each nodes
                    for(NodesContainerType::iterator node_i = rNodes.begin(); node_i!=rNodes.end(); node_i++)
                    {
                        node_i->X() = node_i->GetValue(ELEM_CENTER_X);
                        node_i->Y() = node_i->GetValue(ELEM_CENTER_Y);
                        node_i->Z() = node_i->GetValue(ELEM_CENTER_Z);
                    }
                }

                std::cout<<"  Surface smoothed succesfully with "<< iterations << " iteration(s)"<<std::endl;
            }
        else
            std::cout<<"  Surface was not smoothed" <<std::endl;


        KRATOS_CATCH("");
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
        return "TopologySmoothingUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TopologySmoothingUtilities";
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
      //TopologySmoothingUtilities& operator=(TopologySmoothingUtilities const& rOther);

    /// Copy constructor.
      //TopologySmoothingUtilities(TopologySmoothingUtilities const& rOther);


    ///@}

}; // Class TopologySmoothingUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* KRATOS_TOPOLOGY_SMOOTHING_UTILITIES_H_INCLUDED */

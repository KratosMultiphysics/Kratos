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
	TopologySmoothingUtilities(  )
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

	void SmoothMesh( ModelPart& mModelPart, double iterations )
	{

		KRATOS_TRY;

		std::cout<<"::[Smoothing mesh]::"<<std::endl;

		Vector smoothed_coordinates;
		smoothed_coordinates.resize(mModelPart.NumberOfNodes() * 3);

		if (iterations > 0)
		{
			// Locating neighbouring nodes
			FindConditionsNeighboursProcess nodal_finder = FindConditionsNeighboursProcess(mModelPart, 10, 10);
			nodal_finder.Execute();

			// Repeat smoothing operation for the selected number of iterations
			for(int i = 0; i < iterations; ++i){
				std::cout<<"  Smoothing iteration number "<< i+1 <<std::endl;

				smoothed_coordinates.clear();

				int itr = 0;
				for(NodesContainerType::iterator node_i = mModelPart.Nodes().begin(); node_i!=mModelPart.Nodes().end(); node_i++)
				{
					// Prepare the needed neighbouring node vector and the needed variables
					smoothed_coordinates[3*itr+0] = 0.0;
					smoothed_coordinates[3*itr+1] = 0.0;
					smoothed_coordinates[3*itr+2] = 0.0;

					WeakPointerVector< Node<3> >& neighbours = node_i->GetValue(NEIGHBOUR_NODES);

					// Obtain and sum the X, Y and Z coordinates of all neighbouring nodes
					for( WeakPointerVector<Node<3> >::iterator neighbour_node = neighbours.begin(); neighbour_node!=neighbours.end(); neighbour_node++)
					{
						smoothed_coordinates[3*itr+0] += neighbour_node->X();
						smoothed_coordinates[3*itr+1] += neighbour_node->Y();
						smoothed_coordinates[3*itr+2] += neighbour_node->Z();
					}

					// Average the new X, Y and Z coordinates and save them temporary (simultaneous version)
					smoothed_coordinates[3*itr+0] /= neighbours.size();
					smoothed_coordinates[3*itr+1] /= neighbours.size();
					smoothed_coordinates[3*itr+2] /= neighbours.size();

					itr++;
				}

				// Assign the coordinates calculated value into each nodes
				itr = 0;
				for(NodesContainerType::iterator node_i = mModelPart.Nodes().begin(); node_i!=mModelPart.Nodes().end(); node_i++)
				{
					node_i->X() = smoothed_coordinates[3*itr+0];
					node_i->Y() = smoothed_coordinates[3*itr+1];
					node_i->Z() = smoothed_coordinates[3*itr+2];

					itr++;
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

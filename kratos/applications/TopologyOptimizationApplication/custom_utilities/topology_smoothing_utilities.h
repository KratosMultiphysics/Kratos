// ==============================================================================
/*
 KratosTopologyOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================
//
//   Project Name:        KratosTopology                        $
//   Last modified by:	  $Author:   daniel.baumgaertner@tum.de $
// 						  $Co-Author: Octaviano Malfavón Farías $
//   Date:                $Date:                    August 2016 $
//   Revision:            $Revision:                        0.0 $
//
// ==============================================================================

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
#include "processes/find_nodal_neighbours_process.h" // To find node neighbours using conditions


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

	/// Smooth mesh by performing a Laplacian smoothing on a given surface mesh
	// Laplacian smoothing modifies a given node-position by an average of the positions of the neighbour-nodes
	void SmoothMesh( ModelPart& mModelPart, double relaxation_factor, double iterations )
	{

		KRATOS_TRY;

		std::cout<<"::[Smoothing mesh]::"<<std::endl;

		Vector smoothed_coordinates;
		smoothed_coordinates.resize(mModelPart.NumberOfNodes() * 3);

		// Start neighbour search process
		FindNodalNeighboursProcess nodal_finder = FindNodalNeighboursProcess(mModelPart, 10, 10);
		nodal_finder.Execute();

		// Perform smoothing several times for
		if (iterations > 0)
		{
			// Repeat smoothing operation for the selected number of iterations
			for(int i = 0; i < iterations; ++i){
				std::cout<<"  Smoothing iteration number "<< i+1 <<std::endl;

				smoothed_coordinates.clear();

				// Note that the Laplacian smoothing is done each taking into account the point to be smoothed
				int itr = 0;
				for(NodesContainerType::iterator node_i = mModelPart.Nodes().begin(); node_i!=mModelPart.Nodes().end(); node_i++)
				{
					// Start averaging by taking into account the point to be smoothed
					smoothed_coordinates[3*itr+0] = node_i->X();
					smoothed_coordinates[3*itr+1] = node_i->Y();
					smoothed_coordinates[3*itr+2] = node_i->Z();

					// Average node position
					WeakPointerVector< Node<3> >& neighbours = node_i->GetValue(NEIGHBOUR_NODES);
					for( WeakPointerVector<Node<3> >::iterator neighbour_node = neighbours.begin(); neighbour_node!=neighbours.end(); neighbour_node++)
					{
						// Obtain and sum the X, Y and Z coordinates of all adjacent nodes
						smoothed_coordinates[3*itr+0] += neighbour_node->X();
						smoothed_coordinates[3*itr+1] += neighbour_node->Y();
						smoothed_coordinates[3*itr+2] += neighbour_node->Z();
					}

					// Average the new X, Y and Z coordinates and save them temporary (simultaneous Laplacian Smoothing)
					smoothed_coordinates[3*itr+0] /= (neighbours.size()+1);
					smoothed_coordinates[3*itr+1] /= (neighbours.size()+1);
					smoothed_coordinates[3*itr+2] /= (neighbours.size()+1);

					itr++;
				}

				// Assign the coordinates calculated value into each nodes
				itr = 0;
				for(NodesContainerType::iterator node_i = mModelPart.Nodes().begin(); node_i!=mModelPart.Nodes().end(); node_i++)
				{
					// Move nodes and overwrite their initial position
					node_i->X() += (smoothed_coordinates[3*itr+0]-node_i->X())*relaxation_factor;
					node_i->Y() += (smoothed_coordinates[3*itr+1]-node_i->Y())*relaxation_factor;
					node_i->Z() += (smoothed_coordinates[3*itr+2]-node_i->Z())*relaxation_factor;
					node_i->X0() = node_i->X();
					node_i->Y0() = node_i->Y();
					node_i->Z0() = node_i->Z();

					itr++;
				}
			}

			std::cout<<"  Surface smoothed successfully with "<< iterations << " iteration(s)"<<std::endl;
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

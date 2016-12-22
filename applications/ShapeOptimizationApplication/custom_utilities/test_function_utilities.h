// ==============================================================================
/*
 KratosShapeOptimizationApplication
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
//   Project Name:        KratosShape                            $
//   Created by:          $Author:    daniel.baumgaertner@tum.de $
//   Date:                $Date:                   December 2016 $
//   Revision:            $Revision:                         0.0 $
//
// ==============================================================================

#ifndef FUNCTION_TESTER_UTILITIES
#define FUNCTION_TESTER_UTILITIES

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "../../kratos/includes/define.h"
#include "../../kratos/processes/process.h"
#include "../../kratos/includes/node.h"
#include "../../kratos/includes/element.h"
#include "../../kratos/includes/model_part.h"
#include "../../kratos/includes/kratos_flags.h"
#include "../../kratos/spatial_containers/spatial_containers.h"
#include "../../kratos/utilities/timer.h"
#include "../../kratos/processes/node_erase_process.h"
#include "../../kratos/utilities/binbased_fast_point_locator.h"
#include "../../kratos/utilities/normal_calculation_utils.h"
#include "../../kratos/spaces/ublas_space.h"
#include "shape_optimization_application.h"

// ==============================================================================

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

*/

class FunctionTester
{
  public:
    ///@name Type Definitions
    ///@{

    // ==========================================================================
    // Type definitions for better reading later
    // ==========================================================================
    typedef array_1d<double, 3> array_3d;
    typedef Node<3> NodeType;
    typedef std::vector<NodeType::Pointer> PointVector;
    typedef std::vector<NodeType::Pointer>::iterator PointIterator;
    typedef std::vector<double> DistanceVector;
    typedef std::vector<double>::iterator DistanceIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    // ==========================================================================
    // Type definitions for linear algebra including sparse systems
    // ==========================================================================
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef typename SparseSpaceType::MatrixType SparseMatrixType;
    typedef typename SparseSpaceType::VectorType VectorType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    /// Pointer definition of FunctionTester
    KRATOS_CLASS_POINTER_DEFINITION(FunctionTester);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FunctionTester(ModelPart &model_part)
        : mr_model_part(model_part)
    {
    }

    /// Destructor.
    virtual ~FunctionTester()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void test_function()
    {
        KRATOS_TRY;

        // Create list of possible nearest neighbors in a KD-Tree
        PointVector list_of_nodes;
        for (ModelPart::NodesContainerType::iterator node_it = mr_model_part.NodesBegin(); node_it != mr_model_part.NodesEnd(); ++node_it)
        {
            NodeType::Pointer pnode = *(node_it.base());
            list_of_nodes.push_back(pnode);
        }

        // Construct KD-Tree
        typedef Bucket< 3, NodeType, PointVector, NodeType::Pointer, PointIterator, DistanceIterator > BucketType;
        typedef Tree< KDTreePartition<BucketType> > tree;
        int bucket_size = 20;
        tree nodes_tree(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);

        // Loop over all integration points of model-part and find corresponding closest neighbors
        for (ModelPart::ElementsContainerType::iterator elem_i = mr_model_part.ElementsBegin(); elem_i != mr_model_part.ElementsEnd(); ++elem_i)
        {
        	// Get geometry information of current element (integration method, integration points, shape function values of integration points)
        	const Element::GeometryType::IntegrationMethod integration_method = elem_i->GetIntegrationMethod();
        	const Element::GeometryType::IntegrationPointsArrayType& integration_points = elem_i->GetGeometry().IntegrationPoints(integration_method);
        	const Matrix& N_container = elem_i->GetGeometry().ShapeFunctionsValues(integration_method);

        	std::cout << "------------------------------" << std::endl;
        	std::cout << "elem_i->Id() = " << elem_i->Id() << std::endl;

        	for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        	{
        		// Compute global coordinates of current integration point and get corresponding weight
        		NodeType::CoordinatesArrayType ip_coordinates = elem_i->GetGeometry().GlobalCoordinates(ip_coordinates, integration_points[PointNumber].Coordinates());
        		NodeType::Pointer point_of_interest = Node < 3 > ::Pointer(new Node<3>(PointNumber, ip_coordinates ));
        		double i_weight = integration_points[PointNumber].Weight();

        		// Search nearest neighbor of current integration point
        		NodeType resulting_nearest_point;
        		NodeType::Pointer nearest_point = nodes_tree.SearchNearestPoint( *point_of_interest );

        		// Get FEM-shape-function-value for current integration point
        		Vector N_FEM_GPi = row( N_container, PointNumber);

        		// Some output
        		KRATOS_WATCH(integration_points[PointNumber]);
        		KRATOS_WATCH(*point_of_interest);
        		KRATOS_WATCH(*nearest_point);
        		KRATOS_WATCH(N_FEM_GPi);
        	}

        }

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------

    // ==============================================================================

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
        return "FunctionTester";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "FunctionTester";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
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

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    ModelPart &mr_model_part;

    // ==============================================================================
    // General working arrays
    // ==============================================================================

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
    //      FunctionTester& operator=(FunctionTester const& rOther);

    /// Copy constructor.
    //      FunctionTester(FunctionTester const& rOther);

    ///@}

}; // Class FunctionTester

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // FUNCTION_TESTER_UTILITIES

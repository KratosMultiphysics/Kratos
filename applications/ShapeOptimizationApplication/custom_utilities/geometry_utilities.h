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

#ifndef GEOMETRY_UTILITIES_H
#define GEOMETRY_UTILITIES_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>      // for std::setprecision

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
#include "../../kratos/utilities/normal_calculation_utils.h"
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

class GeometryUtilities
{
public:
    ///@name Type Definitions
    ///@{

    // ==========================================================================
    // Type definitions for better reading later
    // ==========================================================================
    typedef array_1d<double,3> array_3d;

    /// Pointer definition of GeometryUtilities
    KRATOS_CLASS_POINTER_DEFINITION(GeometryUtilities);

	// Structs needed for operations related to surface extraction
	struct KeyComparor
	{
		bool operator()(const vector<unsigned int>& lhs, const vector<unsigned int>& rhs) const
		{
			if(lhs.size() != rhs.size())
				return false;

			for(unsigned int i=0; i<lhs.size(); i++)
				if(lhs[i] != rhs[i]) return false;

			return true;
		}
	};
	struct KeyHasher
	{
		std::size_t operator()(const vector<int>& k) const
		{
			return boost::hash_range(k.begin(), k.end());
		}
	};

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GeometryUtilities( ModelPart& model_part )
        : mr_model_part(model_part)
    {
        // Set precision for output
        std::cout.precision(12);
    }

    /// Destructor.
    virtual ~GeometryUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================

    void compute_unit_surface_normals()
    {
        KRATOS_TRY;

        // Compute nodal are normal using given Kratos utilities (sets the variable "NORMAL")
        NormalCalculationUtils normal_util = NormalCalculationUtils();
        const unsigned int domain_size = mr_model_part.GetProcessInfo().GetValue(DOMAIN_SIZE);
        normal_util.CalculateOnSimplex(mr_model_part,domain_size);

        // Take into account boundary conditions, normalize area normal and store in respective variable
        for (ModelPart::NodeIterator node_i = mr_model_part.NodesBegin(); node_i != mr_model_part.NodesEnd(); ++node_i)
        {
            // Normalize normal and assign to solution step value
            array_3d area_normal = node_i->FastGetSolutionStepValue(NORMAL);
            array_3d normalized_normal = area_normal / norm_2(area_normal);
            noalias(node_i->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL)) = normalized_normal;
        }

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void project_grad_on_unit_surface_normal( bool constraint_given )
    {
        KRATOS_TRY;

        // We loop over all nodes and compute the part of the sensitivity which is in direction to the surface normal
        for (ModelPart::NodeIterator node_i = mr_model_part.NodesBegin(); node_i != mr_model_part.NodesEnd(); ++node_i)
        {
            // We compute dFdX_n = (dFdX \cdot n) * n
            array_3d node_sens = node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY);
            array_3d node_normal = node_i->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
            double surface_sens = inner_prod(node_sens,node_normal);
            array_3d normal_node_sens = surface_sens * node_normal;

            // Assign resulting sensitivities back to node
            node_i->GetSolutionStepValue(OBJECTIVE_SURFACE_SENSITIVITY) = surface_sens;
            noalias(node_i->FastGetSolutionStepValue(OBJECTIVE_SENSITIVITY)) = normal_node_sens;

            // Repeat for constraint
            if(constraint_given)
            {
                // We compute dFdX_n = (dFdX \cdot n) * n
                node_sens = node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY);
                node_normal = node_i->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
                surface_sens = inner_prod(node_sens,node_normal);
                normal_node_sens =  surface_sens * node_normal;

                // Assign resulting sensitivities back to node
                node_i->GetSolutionStepValue(CONSTRAINT_SURFACE_SENSITIVITY) = surface_sens;
                noalias(node_i->FastGetSolutionStepValue(CONSTRAINT_SENSITIVITY)) = normal_node_sens;
            }
        }

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void extract_surface_nodes( std::string const& NewSubModelPartName )
    {
    	KRATOS_TRY;

    	if(mr_model_part.HasSubModelPart(NewSubModelPartName))
    	{
    		std::cout << "> Specified name for sub-model part already defined. Skipping extraction of surface nodes!" << std::endl;
    		return;
    	}

    	// Create new sub-model part within the given main model part that shall list all surface nodes
    	mr_model_part.CreateSubModelPart(NewSubModelPartName);

    	// Some type-definitions
    	typedef boost::unordered_map<vector<unsigned int>, unsigned int, KeyHasher, KeyComparor > hashmap;

    	// Create map to ask for number of faces for the given set of node ids representing one face in the model part
    	hashmap n_faces_map;

    	// Fill map that counts number of faces for given set of nodes
    	for (ModelPart::ElementIterator itElem = mr_model_part.ElementsBegin(); itElem != mr_model_part.ElementsEnd(); itElem++)
    	{
    		Element::GeometryType::GeometriesArrayType faces = itElem->GetGeometry().Faces();

    		for(unsigned int face=0; face<faces.size(); face++)
    		{
    			// Create vector that stores all node is of current face
    			vector<unsigned int> ids(faces[face].size());

    			// Store node ids
    			for(unsigned int i=0; i<faces[face].size(); i++)
    				ids[i] = faces[face][i].Id();

    			//*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
    			std::sort(ids.begin(), ids.end());

    			// Fill the map
    			n_faces_map[ids] += 1;
    		}
    	}

    	// Vector to store all nodes on surface. Node ids may be listed several times
    	std::vector<std::size_t> temp_surface_node_ids;

    	// Add surface nodes to sub-model part
    	for(typename hashmap::const_iterator it=n_faces_map.begin(); it!=n_faces_map.end(); it++)
    	{
    		// If given node set represents face that is not overlapping with a face of another element, add it as skin element
    		if(it->second == 1)
    		{
    			for(unsigned int i=0; i<it->first.size(); i++)
    				temp_surface_node_ids.push_back(it->first[i]);
    		}
    	}

    	// Add nodes and remove double entries
    	mr_model_part.GetSubModelPart(NewSubModelPartName).AddNodes(temp_surface_node_ids);

    	KRATOS_CATCH("");
    }

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
        return "GeometryUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GeometryUtilities";
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

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
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
//      GeometryUtilities& operator=(GeometryUtilities const& rOther);

    /// Copy constructor.
//      GeometryUtilities(GeometryUtilities const& rOther);


    ///@}

}; // Class GeometryUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // GEOMETRY_UTILITIES_H

// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
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
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "utilities/normal_calculation_utils.h"
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
    GeometryUtilities( ModelPart& modelPart )
        : mrModelPart( modelPart )
    {
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
    void ComputeUnitSurfaceNormals()
    {
        KRATOS_TRY;

        // Compute nodal are normal using given Kratos utilities (sets the variable "NORMAL")
        NormalCalculationUtils normal_util = NormalCalculationUtils();
        const unsigned int domain_size = mrModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE);
        normal_util.CalculateOnSimplex(mrModelPart,domain_size);

        // Take into account boundary conditions, normalize area normal and store in respective variable
        for (ModelPart::NodeIterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); ++node_i)
        {
            // Normalize normal and assign to solution step value
            array_3d& normalized_normal = node_i->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);
            const array_1d<double,3>& area_normal = node_i->FastGetSolutionStepValue(NORMAL);
            noalias(normalized_normal) = area_normal/norm_2(area_normal);
        }

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void ProjectNodalVariableOnUnitSurfaceNormals( const Variable<array_3d> &rNodalVariable )
    {
        KRATOS_TRY;

        // We loop over all nodes and compute the part of the sensitivity which is in direction to the surface normal
        for (ModelPart::NodeIterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); ++node_i)
        {
            array_3d &nodal_variable = node_i->FastGetSolutionStepValue(rNodalVariable);
            array_3d &node_normal = node_i->FastGetSolutionStepValue(NORMALIZED_SURFACE_NORMAL);

            // We compute dFdX_n = (dFdX \cdot n) * n
            double surface_sens = inner_prod(nodal_variable,node_normal);
            nodal_variable = surface_sens * node_normal;
        }

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void UpdateShapeChangeByInputVariable( const Variable<array_3d> &rNodalVariable )
    {
        for (ModelPart::NodeIterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); ++node_i)
        {
            array_3d& variable_update = node_i->FastGetSolutionStepValue(rNodalVariable);                
           
            node_i->X() += variable_update[0];
            node_i->Y() += variable_update[1];
            node_i->Z() += variable_update[2];
            noalias(node_i->FastGetSolutionStepValue(SHAPE_CHANGE)) += variable_update;
        }
    }

    // --------------------------------------------------------------------------
    void ExtractSurfaceNodes( std::string const& NewSubModelPartName )
    {
    	KRATOS_TRY;

    	if(mrModelPart.HasSubModelPart(NewSubModelPartName))
    	{
    		std::cout << "> Specified name for sub-model part already defined. Skipping extraction of surface nodes!" << std::endl;
    		return;
    	}

    	// Create new sub-model part within the given main model part that shall list all surface nodes
    	mrModelPart.CreateSubModelPart(NewSubModelPartName);

    	// Some type-definitions
    	typedef boost::unordered_map<vector<unsigned int>, unsigned int, KeyHasher, KeyComparor > hashmap;

    	// Create map to ask for number of faces for the given set of node ids representing one face in the model part
    	hashmap n_faces_map;

    	// Fill map that counts number of faces for given set of nodes
    	for (ModelPart::ElementIterator itElem = mrModelPart.ElementsBegin(); itElem != mrModelPart.ElementsEnd(); itElem++)
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
    	for(hashmap::const_iterator it=n_faces_map.begin(); it!=n_faces_map.end(); it++)
    	{
    		// If given node set represents face that is not overlapping with a face of another element, add it as skin element
    		if(it->second == 1)
    		{
    			for(unsigned int i=0; i<it->first.size(); i++)
    				temp_surface_node_ids.push_back(it->first[i]);
    		}
    	}

    	// Add nodes and remove double entries
    	mrModelPart.GetSubModelPart(NewSubModelPartName).AddNodes(temp_surface_node_ids);

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
    ModelPart& mrModelPart;

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

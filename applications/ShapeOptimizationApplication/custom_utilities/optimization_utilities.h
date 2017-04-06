// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef OPTIMIZATION_UTILITIES_H
#define OPTIMIZATION_UTILITIES_H

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
#include "../../kratos/utilities/timer.h"
#include "shape_optimization_application.h"
#include "../../kratos/spaces/ublas_space.h"

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

class OptimizationUtilities
{
public:
    ///@name Type Definitions
    ///@{

    typedef array_1d<double,3> array_3d;

    /// Pointer definition of OptimizationUtilities
    KRATOS_CLASS_POINTER_DEFINITION(OptimizationUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    OptimizationUtilities( ModelPart& designSurface, Parameters& optimizationSettings )
        : mrDesignSurface( designSurface ),
          mStepSize( optimizationSettings["line_search"]["step_size"].GetDouble() ),
          mNormalizeSearchDirection( optimizationSettings["line_search"]["normalize_search_direction"].GetBool() )
    {
        // Set precision for output
        std::cout.precision(12);

        // Initialize constraint value
        previous_c = 0.0;

        // Initialize variables
        mCorrectionScaling = 0.1;
    }

    /// Destructor.
    virtual ~OptimizationUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void setPrecisionForOutput()
    {
        std::cout.precision(12);
    }        

    // ==============================================================================
    // General optimization operations
    // ==============================================================================
    void compute_design_update()
    {
        KRATOS_TRY;

        // Computation of update of design variable. Normalization is applied if specified.
        if(mNormalizeSearchDirection)
        {
            // Compute max norm of search direction
            double max_norm_search_dir = 0.0;
            for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
            {
                array_3d search_dir = node_i->FastGetSolutionStepValue(SEARCH_DIRECTION);
                double squared_length = inner_prod(search_dir,search_dir);
                
                if(squared_length>max_norm_search_dir)
                    max_norm_search_dir = squared_length;
            }
            max_norm_search_dir = sqrt(max_norm_search_dir);

            // Compute update
            for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
            {
                array_3d design_update = mStepSize * ( node_i->FastGetSolutionStepValue(SEARCH_DIRECTION)/max_norm_search_dir );
                noalias(node_i->FastGetSolutionStepValue(DESIGN_UPDATE)) = design_update;

                // Sum design updates to obtain control point position
                noalias(node_i->FastGetSolutionStepValue(DESIGN_CHANGE_ABSOLUTE)) += design_update;
            }
        }
        else
        {
            // Compute update
            for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
            {
                array_3d design_update = mStepSize * ( node_i->FastGetSolutionStepValue(SEARCH_DIRECTION) );
                noalias(node_i->FastGetSolutionStepValue(DESIGN_UPDATE)) = design_update;

                // Sum design updates to obtain control point position
                noalias(node_i->FastGetSolutionStepValue(DESIGN_CHANGE_ABSOLUTE)) += design_update;
            }
        }

        KRATOS_CATCH("");
    }

    // ==============================================================================
    // For running unconstrained descent methods
    // ==============================================================================
    void compute_search_direction_steepest_descent()
    {
        KRATOS_TRY;

        // Some output for information
        std::cout << "\n> No constraints given or active. The negative objective gradient is chosen as search direction..." << std::endl;

        // search direction is negative of filtered gradient
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
            node_i->FastGetSolutionStepValue(SEARCH_DIRECTION) = -1.0 * node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY);
        }

        KRATOS_CATCH("");
    }

    // ==============================================================================
    // For running penalized projection method
    // ==============================================================================
    void compute_projected_search_direction( double c )
    {
        KRATOS_TRY;

        // Some output for information
        std::cout << "\n> Constraint is active. Modified search direction on the constraint hyperplane is computed..." << std::endl;

        // Compute norm of constraint gradient
        double norm_2_dCds_i = 0.0;
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
        	array_3d dCds_i = node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);
            norm_2_dCds_i += inner_prod(dCds_i,dCds_i);
        }
       norm_2_dCds_i = sqrt(norm_2_dCds_i);

        // Compute dot product of objective gradient and normalized constraint gradient
        double dot_dFds_dCds = 0.0;
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
        	array_3d dFds_i = node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY);
        	array_3d dCds_i = node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);
            dot_dFds_dCds += inner_prod(dFds_i,(dCds_i / norm_2_dCds_i));
        }

        // Compute and assign projected search direction
        for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        {
        	array_3d dFds_i = node_i->FastGetSolutionStepValue(MAPPED_OBJECTIVE_SENSITIVITY);
        	array_3d dCds_i = node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);

        	array_3d projection_term = dot_dFds_dCds * (dCds_i / norm_2_dCds_i);

            node_i->FastGetSolutionStepValue(SEARCH_DIRECTION) = -1 * (dFds_i - projection_term);
        }

        KRATOS_CATCH("");
    }

    // ==============================================================================
    void correct_projected_search_direction( double c )
    {
        // Check correction necessary
        if(c==0)
         return;

        // Calucation of vector norms 
    	double norm_correction_term = 0.0;
    	double norm_search_direction = 0.0;
    	for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
    	{
    		array_3d correction_term = c * node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);
    		norm_correction_term += inner_prod(correction_term,correction_term);

    		array_3d ds = node_i->FastGetSolutionStepValue(SEARCH_DIRECTION);
    		norm_search_direction += inner_prod(ds,ds);
    	}
    	norm_correction_term = sqrt(norm_correction_term);
    	norm_search_direction = sqrt(norm_search_direction);

    	// Three cases need to be covered

		// 1) In case we have two subsequently decreasing constraint values --> correction is fine --> leave current correction scaling

    	// 2) In case the correction jumps over the constraint (change of sign) --> correction was too big --> reduce
    	if(c*previous_c<0)
    	{
    		mCorrectionScaling *= 0.5;
    		std::cout << "Correction scaling needs to decrease...." << std::endl;
    	}

    	// 3) In case we have subsequently increasing constraint value --> correction was too low --> increase
    	if(std::abs(c)>std::abs(previous_c) && c*previous_c>0)
    	{
    		std::cout << "Correction scaling needs to increase...." << std::endl;
    		mCorrectionScaling = std::min(mCorrectionScaling*2,1.0);
    	}
    	double correction_factor = mCorrectionScaling * norm_search_direction / norm_correction_term;

        // Perform correction
    	for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
    	{
    		array_3d correction_term = correction_factor * c * node_i->FastGetSolutionStepValue(MAPPED_CONSTRAINT_SENSITIVITY);
            // KRATOS_WATCH(correction_term);
            // KRATOS_WATCH(node_i->FastGetSolutionStepValue(SEARCH_DIRECTION));
    		node_i->FastGetSolutionStepValue(SEARCH_DIRECTION) -= correction_term;
    	}

        // Store constraint value for next correction step
        previous_c = c;
    }

    // --------------------------------------------------------------------------
    void set_correction_scaling( double correctionScaling )
    {
        mCorrectionScaling = correctionScaling;
    }

    // --------------------------------------------------------------------------
    double get_correction_scaling()
    {
        return mCorrectionScaling;
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
        return "OptimizationUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "OptimizationUtilities";
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
    ModelPart& mrDesignSurface;
    double mStepSize;
    bool mNormalizeSearchDirection;
    double previous_c;

    // ==============================================================================
    // For running penalized projection method
    // ==============================================================================
    double mCorrectionScaling;

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
//      OptimizationUtilities& operator=(OptimizationUtilities const& rOther);

    /// Copy constructor.
//      OptimizationUtilities(OptimizationUtilities const& rOther);


    ///@}

}; // Class OptimizationUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // OPTIMIZATION_UTILITIES_H

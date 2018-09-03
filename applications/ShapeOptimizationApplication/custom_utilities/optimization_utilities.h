// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
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

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
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
    OptimizationUtilities( ModelPart& designSurface, Parameters optimizationSettings )
        : mrDesignSurface( designSurface ),
          mOptimizationSettings( optimizationSettings )
    {
        // Initialize constraint value
        mConstraintValue = 0.0;
        mPreviousConstraintValue = 0.0;
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

    // ==============================================================================
    // General optimization operations
    // ==============================================================================
    void ComputeControlPointUpdate()
    {
        KRATOS_TRY;

        double step_size = mOptimizationSettings["optimization_algorithm"]["line_search"]["step_size"].GetDouble();
        bool normalize_search_direction = mOptimizationSettings["optimization_algorithm"]["line_search"]["normalize_search_direction"].GetBool();


        // Computation of update of design variable. Normalization is applied if specified.
        if(normalize_search_direction)
        {
            // Compute max norm of search direction
            double max_norm_search_dir = 0.0;
            for (auto & node_i : mrDesignSurface.Nodes())
            {
                array_3d& search_dir = node_i.FastGetSolutionStepValue(SEARCH_DIRECTION);
                double squared_length = inner_prod(search_dir,search_dir);

                if(squared_length>max_norm_search_dir)
                    max_norm_search_dir = squared_length;
            }
            max_norm_search_dir = std::sqrt(max_norm_search_dir);

            // Normalize by max norm
            if(max_norm_search_dir>1e-10)
            {
                for (auto & node_i : mrDesignSurface.Nodes())
                {
                    array_3d normalized_search_direction = node_i.FastGetSolutionStepValue(SEARCH_DIRECTION)/max_norm_search_dir;
                    noalias(node_i.FastGetSolutionStepValue(SEARCH_DIRECTION)) = normalized_search_direction;
                }
            }
            else
                std::cout << "> WARNING: Normalization of search direction by max norm activated but max norm is < 1e-10. Hence normalization is ommited!" << std::endl;
        }

        // Compute update
        for (auto & node_i : mrDesignSurface.Nodes())
            noalias(node_i.FastGetSolutionStepValue(CONTROL_POINT_UPDATE)) = step_size * node_i.FastGetSolutionStepValue(SEARCH_DIRECTION);

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void AddFirstVariableToSecondVariable( const Variable<array_3d> &rFirstVariable, const Variable<array_3d> &rSecondVariable )
    {
        for (auto & node_i : mrDesignSurface.Nodes())
            noalias(node_i.FastGetSolutionStepValue(rSecondVariable)) += node_i.FastGetSolutionStepValue(rFirstVariable);
    }

    // ==============================================================================
    // For running unconstrained descent methods
    // ==============================================================================
    void ComputeSearchDirectionSteepestDescent()
    {
        KRATOS_TRY;

        // Some output for information
        std::cout << "\n> No constraints given or active. The negative objective gradient is chosen as search direction..." << std::endl;

        // search direction is negative of filtered gradient
        for (auto & node_i : mrDesignSurface.Nodes())
        {
            node_i.FastGetSolutionStepValue(SEARCH_DIRECTION) = -1.0 * node_i.FastGetSolutionStepValue(DF1DX_MAPPED);
        }

        KRATOS_CATCH("");
    }

    // ==============================================================================
    // For running penalized projection method
    // ==============================================================================
    void ComputeProjectedSearchDirection()
    {
        KRATOS_TRY;

        // Some output for information
        std::cout << "\n> Constraint is active. Modified search direction on the constraint hyperplane is computed..." << std::endl;

        // Compute norm of constraint gradient
        double norm_2_dCds_i = 0.0;
        for (auto & node_i : mrDesignSurface.Nodes())
        {
        	array_3d& dCds_i = node_i.FastGetSolutionStepValue(DC1DX_MAPPED);
            norm_2_dCds_i += inner_prod(dCds_i,dCds_i);
        }
        norm_2_dCds_i = std::sqrt(norm_2_dCds_i);

        // Avoid division by zero
        if(std::abs(norm_2_dCds_i)<1e-12)
            norm_2_dCds_i = 1.0;

        // Compute dot product of objective gradient and normalized constraint gradient
        double dot_dFds_dCds = 0.0;
        for (auto & node_i : mrDesignSurface.Nodes())
        {
        	array_3d dFds_i = node_i.FastGetSolutionStepValue(DF1DX_MAPPED);
        	array_3d dCds_i = node_i.FastGetSolutionStepValue(DC1DX_MAPPED);
            dot_dFds_dCds += inner_prod(dFds_i,(dCds_i / norm_2_dCds_i));
        }

        // Compute and assign projected search direction
        for (auto & node_i : mrDesignSurface.Nodes())
        {
        	array_3d& dFds_i = node_i.FastGetSolutionStepValue(DF1DX_MAPPED);
        	array_3d& dCds_i = node_i.FastGetSolutionStepValue(DC1DX_MAPPED);

        	array_3d projection_term = dot_dFds_dCds * (dCds_i / norm_2_dCds_i);

            node_i.FastGetSolutionStepValue(SEARCH_DIRECTION) = -1 * (dFds_i - projection_term);
        }

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void CorrectProjectedSearchDirection( double constraint_value )
    {
        mConstraintValue = constraint_value;

        // Check correction necessary
        if(mConstraintValue==0)
         return;

        // Perform correction
        double correction_factor = ComputeCorrectionFactor();
    	for (auto & node_i : mrDesignSurface.Nodes())
    	{
    		array_3d correction_term = correction_factor * mConstraintValue * node_i.FastGetSolutionStepValue(DC1DX_MAPPED);
    		node_i.FastGetSolutionStepValue(SEARCH_DIRECTION) -= correction_term;
    	}

        // Store constraint value for next correction step
        mPreviousConstraintValue = mConstraintValue;
    }

    // --------------------------------------------------------------------------
    double ComputeCorrectionFactor()
    {
    	double norm_correction_term = 0.0;
    	double norm_search_direction = 0.0;
    	for (auto & node_i : mrDesignSurface.Nodes())
    	{
    		array_3d correction_term = mConstraintValue * node_i.FastGetSolutionStepValue(DC1DX_MAPPED);
    		norm_correction_term += inner_prod(correction_term,correction_term);

    		array_3d ds = node_i.FastGetSolutionStepValue(SEARCH_DIRECTION);
    		norm_search_direction += inner_prod(ds,ds);
    	}
    	norm_correction_term = std::sqrt(norm_correction_term);
    	norm_search_direction = std::sqrt(norm_search_direction);
        double correction_scaling = GetCorrectionScaling();

    	return correction_scaling * norm_search_direction / norm_correction_term;
    }

    // --------------------------------------------------------------------------
    double GetCorrectionScaling()
    {
        double correction_scaling = mOptimizationSettings["optimization_algorithm"]["correction_scaling"].GetDouble();
        if(mOptimizationSettings["optimization_algorithm"]["use_adaptive_correction"].GetBool())
        {
            correction_scaling = AdaptCorrectionScaling( correction_scaling );
            mOptimizationSettings["optimization_algorithm"]["correction_scaling"].SetDouble(correction_scaling);
        }
        return correction_scaling;
    }

    // --------------------------------------------------------------------------
    double AdaptCorrectionScaling( double correction_scaling )
    {
    	// Three cases need to be covered
		// 1) In case we have two subsequently decreasing constraint values --> correction is fine --> leave current correction scaling
    	// 2) In case the correction jumps over the constraint (change of sign) --> correction was too big --> reduce
    	if(mConstraintValue*mPreviousConstraintValue<0)
    	{
    		correction_scaling *= 0.5;
    		std::cout << "Correction scaling needs to decrease...." << std::endl;
    	}
    	// 3) In case we have subsequently increasing constraint value --> correction was too low --> increase
    	if(std::abs(mConstraintValue)>std::abs(mPreviousConstraintValue) && mConstraintValue*mPreviousConstraintValue>0)
    	{
    		std::cout << "Correction scaling needs to increase...." << std::endl;
    		correction_scaling = std::min(correction_scaling*2,1.0);
    	}

        return correction_scaling;
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
    Parameters mOptimizationSettings;
    double mConstraintValue;
    double mPreviousConstraintValue;

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

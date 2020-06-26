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
#include "linear_solvers/linear_solver.h"

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
    typedef UblasSpace<double, Matrix, Vector> DenseSpace;

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
        // Initialize member variables for penalized projection
        std::string algorithm_name = optimizationSettings["optimization_algorithm"]["name"].GetString();
        if(algorithm_name == "penalized_projection")
          mCorrectionScaling = optimizationSettings["optimization_algorithm"]["correction_scaling"].GetDouble();
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
    void ComputeControlPointUpdate(const double StepSize)
    {
        KRATOS_TRY;

        // Normalize if specified
        if(mOptimizationSettings["optimization_algorithm"]["line_search"]["normalize_search_direction"].GetBool())
        {
            const double max_norm_search_dir = ComputeMaxNormOfNodalVariable(SEARCH_DIRECTION);
            if(max_norm_search_dir>1e-10)
                for (auto & node_i : mrDesignSurface.Nodes())
                {
                    array_3d& search_dir = node_i.FastGetSolutionStepValue(SEARCH_DIRECTION);
                    search_dir/=max_norm_search_dir;
                }
            else
                KRATOS_WARNING("ShapeOpt::ComputeControlPointUpdate") << "Normalization of search direction by max norm activated but max norm is < 1e-10. Hence normalization is ommited!" << std::endl;
        }

        // Compute update
        for (auto & node_i : mrDesignSurface.Nodes())
            noalias(node_i.FastGetSolutionStepValue(CONTROL_POINT_UPDATE)) = StepSize * node_i.FastGetSolutionStepValue(SEARCH_DIRECTION);

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void AddFirstVariableToSecondVariable( const Variable<array_3d> &rFirstVariable, const Variable<array_3d> &rSecondVariable )
    {
        for (auto & node_i : mrDesignSurface.Nodes())
            noalias(node_i.FastGetSolutionStepValue(rSecondVariable)) += node_i.FastGetSolutionStepValue(rFirstVariable);
    }

    // --------------------------------------------------------------------------
    double ComputeL2NormOfNodalVariable( const Variable<array_3d> &rVariable)
    {
        double l2_norm = 0.0;
        for (auto & node_i : mrDesignSurface.Nodes())
        {
            array_3d& variable_vector = node_i.FastGetSolutionStepValue(rVariable);
            l2_norm += inner_prod(variable_vector,variable_vector);
        }
        return std::sqrt(l2_norm);
    }

    // --------------------------------------------------------------------------
    double ComputeL2NormOfNodalVariable( const Variable<double> &rVariable)
    {
        double l2_norm = 0.0;
        for (auto & node_i : mrDesignSurface.Nodes())
        {
            double &value = node_i.FastGetSolutionStepValue(rVariable);
            l2_norm += value*value;
        }
        return std::sqrt(l2_norm);
    }

    // --------------------------------------------------------------------------
    double ComputeMaxNormOfNodalVariable( const Variable<array_3d> &rVariable)
    {
        double max_norm = 0.0;
        for (auto & node_i : mrDesignSurface.Nodes())
        {
            array_3d& variable_vector = node_i.FastGetSolutionStepValue(rVariable);
            double squared_value = inner_prod(variable_vector,variable_vector);

            max_norm = std::max(squared_value,max_norm);
        }
        return std::sqrt(max_norm);
    }

    // --------------------------------------------------------------------------
    double ComputeMaxNormOfNodalVariable( const Variable<double> &rVariable)
    {
        double max_norm = 0.0;
        for (auto & node_i : mrDesignSurface.Nodes())
        {
            double &value = node_i.FastGetSolutionStepValue(rVariable);
            double squared_value = value*value;

            max_norm = std::max(squared_value,max_norm);
        }
        return std::sqrt(max_norm);
    }

    // ==============================================================================
    // For running unconstrained descent methods
    // ==============================================================================
    void ComputeSearchDirectionSteepestDescent()
    {
        KRATOS_TRY;

        // Some output for information
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("ShapeOpt") << "No constraints given or active. The negative objective gradient is chosen as search direction..." << std::endl;

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
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("ShapeOpt") << "Constraint is active. Modified search direction on the constraint hyperplane is computed..." << std::endl;

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

        if(mOptimizationSettings["optimization_algorithm"]["use_adaptive_correction"].GetBool())
        {
            // Adapt constraint scaling

            // Three cases need to be covered
            // 1) In case we have two subsequently decreasing constraint values --> correction is fine --> leave current correction scaling
            // 2) In case the correction jumps over the constraint (change of sign) --> correction was too big --> reduce
            if(mConstraintValue*mPreviousConstraintValue<0.0)
            {
                mCorrectionScaling *= 0.5;
                KRATOS_INFO("ShapeOpt") << "Correction scaling needs to decrease...." << std::endl;
            }
            // 3) In case we have subsequently increasing constraint value --> correction was too low --> increase
            if(std::abs(mConstraintValue)>std::abs(mPreviousConstraintValue) && mConstraintValue*mPreviousConstraintValue>0)
            {
                KRATOS_INFO("ShapeOpt") << "Correction scaling needs to increase...." << std::endl;
                mCorrectionScaling = std::min(mCorrectionScaling*2,1.0);
            }
        }

    	return mCorrectionScaling * norm_search_direction / norm_correction_term;
    }

    // --------------------------------------------------------------------------
    double GetCorrectionScaling() const
    {
        return mCorrectionScaling;
    }

    /**
     * Assemble the values of the nodal vector variable into a vector
     */
    void AssembleVector(
        Vector& rVector,
        const Variable<array_3d> &rVariable)
    {
        if (rVector.size() != mrDesignSurface.NumberOfNodes()*3){
            rVector.resize(mrDesignSurface.NumberOfNodes()*3);
        }

        int i=0;
        for (auto & node_i : mrDesignSurface.Nodes())
        {
            array_3d& variable_vector = node_i.FastGetSolutionStepValue(rVariable);
            rVector[i*3+0] = variable_vector[0];
            rVector[i*3+1] = variable_vector[1];
            rVector[i*3+2] = variable_vector[2];
            ++i;
        }
    }

    /**
     * Assigns the values of a vector to the nodal vector variables
     */
    void AssignVectorToVariable(
        const Vector& rVector,
        const Variable<array_3d> &rVariable)
    {
        KRATOS_ERROR_IF(rVector.size() != mrDesignSurface.NumberOfNodes()*3)
            << "AssignVectorToVariable: Vector size does not mach number of Nodes!" << std::endl;

        int i=0;
        for (auto & node_i : mrDesignSurface.Nodes())
        {
            array_3d& variable_vector = node_i.FastGetSolutionStepValue(rVariable);
            variable_vector[0] = rVector[i*3+0];
            variable_vector[1] = rVector[i*3+1];
            variable_vector[2] = rVector[i*3+2];
            ++i;
        }
    }

    /**
     * Assemble the values of the nodal vector variables into a dense matrix.
     * One column per variable is created.
     */
    void AssembleMatrix(
        Matrix& rMatrix,
        const std::vector<Variable<array_3d>*>& rVariables
    ) const
    {
        if ((rMatrix.size1() != mrDesignSurface.NumberOfNodes()*3 || rMatrix.size2() !=  rVariables.size())){
            rMatrix.resize(mrDesignSurface.NumberOfNodes()*3, rVariables.size());
        }

        int i=0;
        for (auto & node_i : mrDesignSurface.Nodes())
        {
            int j=0;
            for (Variable<array_3d>* p_variable_j : rVariables)
            {
                const Variable<array_3d>& r_variable_j = *p_variable_j;
                array_3d& variable_vector = node_i.FastGetSolutionStepValue(r_variable_j);
                rMatrix(i*3+0, j) = variable_vector[0];
                rMatrix(i*3+1, j) = variable_vector[1];
                rMatrix(i*3+2, j) = variable_vector[2];
                ++j;
            }
            ++i;
        }
    }


    /**
     * Calculate the projection of the objective gradient into the subspace tangent to
     * the active constraint gradients.
     * In a second step, calculate the restoration move accounting for the current violation of the constraints.
     * Variable naming and implementation based on https://msulaiman.org/onewebmedia/GradProj_2.pdf
     */
    void CalculateProjectedSearchDirectionAndCorrection(
        Vector& rObjectiveGradient,
        Matrix& rConstraintGradients,
        Vector& rConstraintValues,
        LinearSolver<DenseSpace, DenseSpace>& rSolver,
        Vector& rProjectedSearchDirection,
        Vector& rRestoration
        )
    {
        // local variable naming according to https://msulaiman.org/onewebmedia/GradProj_2.pdf
        Vector& nabla_f = rObjectiveGradient;
        Matrix& N = rConstraintGradients;
        Vector& g_a = rConstraintValues;
        Vector& s = rProjectedSearchDirection;
        Vector& c = rRestoration;

        Matrix NTN = prod(trans(N), N);
        Matrix I = IdentityMatrix(N.size2());
        Matrix NTN_inv(NTN.size1(), NTN.size2());

        rSolver.Solve(NTN, NTN_inv, I); // solve with identity to get the inverse

        s = - (nabla_f - prod(N, Vector(prod(NTN_inv, Vector(prod(trans(N), nabla_f))))));

        c = - prod(N, Vector(prod(NTN_inv, g_a)));
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
    double mConstraintValue = 0.0;
    double mPreviousConstraintValue = 0.0;
    double mCorrectionScaling = 1.0;

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

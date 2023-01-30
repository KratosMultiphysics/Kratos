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

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
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

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) OptimizationUtilities
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

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    // General optimization operations
    // ==============================================================================
    static void ComputeControlPointUpdate(ModelPart& rModelPart, const double StepSize, const bool Normalize);

    // --------------------------------------------------------------------------
    static void AddFirstVariableToSecondVariable( ModelPart& rModelPart, const Variable<array_3d> &rFirstVariable, const Variable<array_3d> &rSecondVariable );

    // --------------------------------------------------------------------------
    static double ComputeL2NormOfNodalVariable( ModelPart& rModelPart, const Variable<array_3d> &rVariable);

    // --------------------------------------------------------------------------
    static double ComputeL2NormOfNodalVariable(ModelPart& rModelPart, const Variable<double> &rVariable);

    // --------------------------------------------------------------------------
    static double ComputeMaxNormOfNodalVariable(ModelPart& rModelPart, const Variable<array_3d> &rVariable);

    // --------------------------------------------------------------------------
    static double ComputeMaxNormOfNodalVariable(ModelPart& rModelPart, const Variable<double> &rVariable);

    // ==============================================================================
    // For running unconstrained descent methods
    // ==============================================================================
    static void ComputeSearchDirectionSteepestDescent(ModelPart& rModelPart);

    // ==============================================================================
    // For running penalized projection method
    // ==============================================================================
    static void ComputeProjectedSearchDirection(ModelPart& rModelPart);

    // --------------------------------------------------------------------------
    static double CorrectProjectedSearchDirection(ModelPart& rModelPart, const double PrevConstraintValue, const double ConstraintValue, const double CorrectionScaling, const bool IsAdaptive );

    // --------------------------------------------------------------------------
    static double ComputeCorrectionFactor(ModelPart& rModelPart, const double PrevConstraintValue, const double ConstraintValue, double& CorrectionScaling, const bool IsAdaptive);

    /**
     * Assemble the values of the nodal vector variable into a vector
     */
    static void AssembleVector( ModelPart& rModelPart,
        Vector& rVector,
        const Variable<array_3d> &rVariable);

    /**
     * Assigns the values of a vector to the nodal vector variables
     */
    static void AssignVectorToVariable(ModelPart& rModelPart,
        const Vector& rVector,
        const Variable<array_3d> &rVariable);

    /**
     * Assemble the values of the nodal vector variables into a dense matrix.
     * One column per variable is created.
     */
    static void AssembleMatrix(ModelPart& rModelPart,
        Matrix& rMatrix,
        const std::vector<Variable<array_3d>*>& rVariables);


    /**
     * Calculate the projection of the objective gradient into the subspace tangent to
     * the active constraint gradients.
     * In a second step, calculate the restoration move accounting for the current violation of the constraints.
     * Variable naming and implementation based on https://msulaiman.org/onewebmedia/GradProj_2.pdf
     */
    static void CalculateProjectedSearchDirectionAndCorrection(
        Vector& rObjectiveGradient,
        Matrix& rConstraintGradients,
        Vector& rConstraintValues,
        LinearSolver<DenseSpace, DenseSpace>& rSolver,
        Vector& rProjectedSearchDirection,
        Vector& rRestoration);
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

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) GradientProjectionSolverUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using DenseSpace = UblasSpace<double, Matrix, Vector>;

    ///@}
    ///@name Static operations
    ///@{

    static void CalculateProjectedSearchDirectionAndCorrection(
        Vector& rSearchDirection,
        Vector& rSearchCorrection,
        LinearSolver<DenseSpace, DenseSpace>& rSolver,
        const Vector& rConstraintValues,
        const Vector& rObjectiveGradient,
        const std::vector<Vector>& rConstraintsGradients);

    static void CalculateControlUpdate(
        Vector& rControlUpdate,
        Vector& rSearchDirection,
        Vector& rSearchCorrection,
        const double StepSize,
        const double MaxCorrectionShare);

    ///@}
};

///@}
}
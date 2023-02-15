// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#pragma once

#include <geo_mechanics_application.h>
#include "includes/kernel.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include <solving_strategies/convergencecriterias/mixed_generic_criteria.h>
#include <solving_strategies/strategies/implicit_solving_strategy.h>

namespace Kratos
{
    class KratosGeoApplication
    {

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
    typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;

    public:
        KratosGeoApplication();
        ~KratosGeoApplication() {};

    private:
        Model current_model;
        Kernel kernel;
        KratosGeoMechanicsApplication::Pointer geoApp;
        int echoLevel = 1;

    protected:
        Model& GetModelPointer();
        void ResetModelParts();
        int GetEchoLevel();
        void SetEchoLevel(int level);
        int execute_kratos_calculation(ModelPart& model_part, const std::vector<std::shared_ptr<Process>>& rProcesses,
                                       ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>::Pointer p_solving_strategy,
                                       double time, double delta_time, double number_iterations);

        virtual ConvergenceCriteriaType::Pointer convergence_criteria() = 0;
        virtual LinearSolverType::Pointer solver_settings() = 0;
        virtual ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>::Pointer strategy_settings(ModelPart& model_part) = 0;
    };
}
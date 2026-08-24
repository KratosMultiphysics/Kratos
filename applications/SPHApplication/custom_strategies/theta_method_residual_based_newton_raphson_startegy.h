//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

#pragma once

// System includes
#include <array>

// Project includes
#include "includes/variables.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes
#include "sph_application_variables.h"

namespace Kratos
{

/**
 * @class ThetaResidualBasedNewtonRaphsonStrategy
 * @brief Newton-Raphson strategy for the first-order transient mixed SPH formulation using theta-method as time integrator.
 * @details Velocity and deformation-gradient components are the primary DOFs. After each primary-variable update, the strategy 
 * reconstructs displacement, coordinates, acceleration and the deformation-gradient rates from the theta relations.
 * @remark This class can  only be used with the TotalLagrangianMixedStrainParticle element.
 */
template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class ThetaResidualBasedNewtonRaphsonStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ThetaResidualBasedNewtonRaphsonStrategy);

    using BaseType = ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TSchemeType = typename BaseType::TSchemeType;
    using TConvergenceCriteriaType = typename BaseType::TConvergenceCriteriaType;
    using TBuilderAndSolverType = typename BaseType::TBuilderAndSolverType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    ThetaResidualBasedNewtonRaphsonStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver,
        const int MaxIterations,
        const bool CalculateReactions,
        const bool ReformDofSetAtEachStep,
        const bool MoveMeshFlag)
        : BaseType( 
            rModelPart, 
            pScheme, 
            pConvergenceCriteria, 
            pBuilderAndSolver, 
            MaxIterations, 
            CalculateReactions, 
            ReformDofSetAtEachStep, 
            MoveMeshFlag)
    {
    }

    void Predict() override
    {
        KRATOS_TRY

        BaseType::Predict();
        UpdateDerivedVariables();

        if (this->MoveMeshFlag()) {
            this->MoveMesh();
        }

        KRATOS_CATCH("")
    }

protected:
    void UpdateDatabase(TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb, const bool MoveMesh) override
    {
        KRATOS_TRY

        BaseType::UpdateDatabase(rA, rDx, rb, false);
        UpdateDerivedVariables();

        if (MoveMesh) {
            this->MoveMesh();
        }

        KRATOS_CATCH("")
    }

    void UpdateDerivedVariables()
    {
        auto& r_model_part = this->GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();

        const double delta_time = r_process_info[DELTA_TIME];
        const double theta = r_process_info[TIME_INTEGRATION_THETA];
        const int domain_size = r_process_info[DOMAIN_SIZE];

        const double one_minus_theta = 1.0 - theta;

        const std::array<const Variable<double>*, 9> deformation_gradient_variables = {{
            &DEFORMATION_GRADIENT_XX,
            &DEFORMATION_GRADIENT_YY,
            &DEFORMATION_GRADIENT_XY,
            &DEFORMATION_GRADIENT_YX,
            &DEFORMATION_GRADIENT_ZZ,
            &DEFORMATION_GRADIENT_XZ,
            &DEFORMATION_GRADIENT_YZ,
            &DEFORMATION_GRADIENT_ZX,
            &DEFORMATION_GRADIENT_ZY
        }};
        const std::array<const Variable<double>*, 9> deformation_gradient_rate_variables = {{
            &DEFORMATION_GRADIENT_DOT_XX,
            &DEFORMATION_GRADIENT_DOT_YY,
            &DEFORMATION_GRADIENT_DOT_XY,
            &DEFORMATION_GRADIENT_DOT_YX,
            &DEFORMATION_GRADIENT_DOT_ZZ,
            &DEFORMATION_GRADIENT_DOT_XZ,
            &DEFORMATION_GRADIENT_DOT_YZ,
            &DEFORMATION_GRADIENT_DOT_ZX,
            &DEFORMATION_GRADIENT_DOT_ZY
        }};
        
        const std::size_t number_of_deformation_gradient_components = domain_size == 2 ? 4 : 9;

        for (auto& r_node : r_model_part.Nodes()) {
            
            const auto& r_previous_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT, 1);
            auto& r_current_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            
            const auto& r_previous_velocity = r_node.FastGetSolutionStepValue(VELOCITY, 1);
            const auto& r_current_velocity = r_node.FastGetSolutionStepValue(VELOCITY);
            
            const auto& r_previous_acceleration = r_node.FastGetSolutionStepValue(ACCELERATION, 1);
            auto& r_current_acceleration = r_node.FastGetSolutionStepValue(ACCELERATION);

            // Updating displacement and acceleration using the theta method relations
            for (int component = 0; component < domain_size; ++component) {
                r_current_displacement[component] = r_previous_displacement[component] + delta_time * (one_minus_theta * r_previous_velocity[component] + theta * r_current_velocity[component]);
                r_current_acceleration[component] = ((r_current_velocity[component] - r_previous_velocity[component]) / delta_time - one_minus_theta * r_previous_acceleration[component]) / theta;
            }

            const auto& r_initial_position = r_node.GetInitialPosition();
            auto& r_current_position = r_node.Coordinates();
            
            // Updating the current position using the updated displacement
            for (std::size_t component = 0; component < 3; ++component) {
                r_current_position[component] = r_initial_position[component] + r_current_displacement[component];
            }

            for (std::size_t component = 0; component < number_of_deformation_gradient_components; ++component) {
                const auto& r_deformation_gradient_variable = *deformation_gradient_variables[component];
                const auto& r_deformation_gradient_rate_variable = *deformation_gradient_rate_variables[component];
                
                const double previous_deformation_gradient = r_node.FastGetSolutionStepValue(r_deformation_gradient_variable, 1);
                const double current_deformation_gradient = r_node.FastGetSolutionStepValue(r_deformation_gradient_variable);
                const double previous_deformation_gradient_rate = r_node.FastGetSolutionStepValue(r_deformation_gradient_rate_variable, 1);

                // Updating deformation gradient rate using the theta method relations
                r_node.FastGetSolutionStepValue(r_deformation_gradient_rate_variable) = 
                    ((current_deformation_gradient - previous_deformation_gradient) / delta_time - one_minus_theta * previous_deformation_gradient_rate) / theta;
            }
        }
    }
};

}
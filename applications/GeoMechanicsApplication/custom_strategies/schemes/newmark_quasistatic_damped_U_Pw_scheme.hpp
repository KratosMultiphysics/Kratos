// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

// Application includes
#include "custom_strategies/schemes/newmark_quasistatic_U_Pw_scheme.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class NewmarkQuasistaticDampedUPwScheme : public GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(NewmarkQuasistaticDampedUPwScheme);

    using BaseType              = Scheme<TSparseSpace, TDenseSpace>;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    NewmarkQuasistaticDampedUPwScheme(double beta, double gamma, double theta)
        : GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>(
              {FirstOrderScalarVariable(WATER_PRESSURE, DT_WATER_PRESSURE, DT_PRESSURE_COEFFICIENT)},
              {SecondOrderVectorVariable(DISPLACEMENT), SecondOrderVectorVariable(ROTATION)},
              beta,
              gamma,
              theta)
    {
        // Allocate auxiliary memory
        const auto num_threads = ParallelUtilities::GetNumThreads();
        mDampingMatrix.resize(num_threads);
        mVelocityVector.resize(num_threads);
    }

    void CalculateSystemContributions(Element&                       rCurrentElement,
                                      LocalSystemMatrixType&         LHS_Contribution,
                                      LocalSystemVectorType&         RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      const ProcessInfo&             CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDampingToLHS(LHS_Contribution, mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDampingToRHS(rCurrentElement, RHS_Contribution, mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(Element&                       rCurrentElement,
                                  LocalSystemVectorType&         RHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDampingToRHS(rCurrentElement, RHS_Contribution, mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(Element&                       rCurrentElement,
                                  LocalSystemMatrixType&         LHS_Contribution,
                                  Element::EquationIdVectorType& EquationId,
                                  const ProcessInfo&             CurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution, CurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread], CurrentProcessInfo);

        this->AddDampingToLHS(LHS_Contribution, mDampingMatrix[thread], CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

protected:
    void AddDampingToLHS(LocalSystemMatrixType& LHS_Contribution, LocalSystemMatrixType& C, const ProcessInfo& CurrentProcessInfo)
    {
        // adding damping contribution
        if (C.size1() != 0)
            noalias(LHS_Contribution) += (this->GetGamma() / (this->GetBeta() * this->GetDeltaTime())) * C;
    }

    void AddDampingToRHS(Element&               rCurrentElement,
                         LocalSystemVectorType& RHS_Contribution,
                         LocalSystemMatrixType& C,
                         const ProcessInfo&     CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        // adding damping contribution
        if (C.size1() != 0) {
            rCurrentElement.GetFirstDerivativesVector(mVelocityVector[thread], 0);
            noalias(RHS_Contribution) -= prod(C, mVelocityVector[thread]);
        }
    }

    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            // No accelerations of the displacement/rotation D.O.F. in this quasistatic damped scheme
            this->UpdateVectorFirstTimeDerivative(rNode);

            for (const auto& r_first_order_scalar_variable : this->GetFirstOrderScalarVariables()) {
                this->UpdateScalarTimeDerivative(rNode, r_first_order_scalar_variable.instance,
                                                 r_first_order_scalar_variable.first_time_derivative);
            }
        });

        KRATOS_CATCH("")
    }

private:
    std::vector<Matrix> mDampingMatrix;
    std::vector<Vector> mVelocityVector;
}; // Class NewmarkQuasistaticDampedUPwScheme

} // namespace Kratos
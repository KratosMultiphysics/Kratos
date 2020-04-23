//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Suneth Warnakulasuriya
//                 Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_ALGEBRAIC_FLUX_CORRECTED_SCALAR_STEADY_SCHEME)
#define KRATOS_ALGEBRAIC_FLUX_CORRECTED_SCALAR_STEADY_SCHEME

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_strategies/relaxed_dof_updater.h"
#include "rans_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace>
class AlgebraicFluxCorrectedScalarSteadyScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AlgebraicFluxCorrectedScalarSteadyScheme);

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;

    using DofsArrayType = typename BaseType::DofsArrayType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    AlgebraicFluxCorrectedScalarSteadyScheme(const double RelaxationFactor)
        : BaseType(),
          mRelaxationFactor(RelaxationFactor),
          mrPeriodicIdVar(Variable<int>::StaticObject())
    {
        KRATOS_INFO("AlgebraicFluxCorrectedScalarSteadyScheme")
            << " Using residual based algebraic flux corrected scheme with "
               "relaxation "
               "factor = "
            << std::scientific << mRelaxationFactor << "\n";
    }

    AlgebraicFluxCorrectedScalarSteadyScheme(const double RelaxationFactor,
                                             const Variable<int>& rPeriodicIdVar)
        : BaseType(), mRelaxationFactor(RelaxationFactor), mrPeriodicIdVar(rPeriodicIdVar)
    {
        KRATOS_INFO("AlgebraicFluxCorrectedScalarSteadyScheme")
            << " Using periodic residual based algebraic flux corrected scheme "
               "with relaxation "
               "factor = "
            << std::scientific << mRelaxationFactor << "\n";
    }

    ~AlgebraicFluxCorrectedScalarSteadyScheme() override = default;

    ///@}
    ///@name Operators
    ///@{

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        if (mrPeriodicIdVar != Variable<int>::StaticObject())
        {
            const int number_of_conditions = rModelPart.NumberOfConditions();
#pragma omp parallel for
            for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
            {
                const ModelPart::ConditionType& r_condition =
                    *(rModelPart.ConditionsBegin() + i_cond);
                if (r_condition.Is(PERIODIC))
                {
                    // this only supports 2 noded periodic conditions
                    KRATOS_ERROR_IF(r_condition.GetGeometry().PointsNumber() != 2)
                        << this->Info() << " only supports two noded periodic conditions. Found "
                        << r_condition.Info() << " with "
                        << r_condition.GetGeometry().PointsNumber() << " nodes.\n";

                    const ModelPart::NodeType& r_node_0 = r_condition.GetGeometry()[0];
                    const std::size_t r_node_0_pair_id =
                        r_node_0.FastGetSolutionStepValue(mrPeriodicIdVar);

                    const ModelPart::NodeType& r_node_1 = r_condition.GetGeometry()[1];
                    const std::size_t r_node_1_pair_id =
                        r_node_1.FastGetSolutionStepValue(mrPeriodicIdVar);

                    KRATOS_ERROR_IF(r_node_0_pair_id != r_node_1.Id())
                        << "Periodic condition pair id mismatch in "
                        << mrPeriodicIdVar.Name() << ". [ " << r_node_0_pair_id
                        << " != " << r_node_1.Id() << " ].\n";

                    KRATOS_ERROR_IF(r_node_1_pair_id != r_node_0.Id())
                        << "Periodic condition pair id mismatch in "
                        << mrPeriodicIdVar.Name() << ". [ " << r_node_1_pair_id
                        << " != " << r_node_0.Id() << " ].\n";
                }
            }
        }

        KRATOS_CATCH("");
    }

    void InitializeNonLinIteration(ModelPart& rModelPart,
                                   TSystemMatrixType& A,
                                   TSystemVectorType& Dx,
                                   TSystemVectorType& b) override
    {
        KRATOS_TRY

        ModelPart::NodesContainerType& r_nodes = rModelPart.Nodes();

        VariableUtils variable_utilities;
        variable_utilities.SetHistoricalVariableToZero(
            AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX, r_nodes);
        variable_utilities.SetHistoricalVariableToZero(
            AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX, r_nodes);
        variable_utilities.SetHistoricalVariableToZero(
            AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT, r_nodes);
        variable_utilities.SetHistoricalVariableToZero(
            AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT, r_nodes);

        ModelPart::ElementsContainerType& r_elements = rModelPart.Elements();
        const int number_of_elements = r_elements.size();

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

#pragma omp parallel
        {
            Matrix left_hand_side, artificial_diffusion;
            Vector right_hand_side, values;
            std::vector<IndexType> equation_ids;
#pragma omp for
            for (int i = 0; i < number_of_elements; ++i)
            {
                ModelPart::ElementType& r_element = *(r_elements.begin() + i);
                this->CalculateSystemMatrix(r_element, left_hand_side, r_current_process_info);
                this->CalculateArtificialDiffusionMatrix(artificial_diffusion, left_hand_side);
                r_element.GetValuesVector(values);
                r_element.EquationIdVector(equation_ids, r_current_process_info);

                const int size = artificial_diffusion.size1();

                Vector p_plus = ZeroVector(size);
                Vector p_minus = ZeroVector(size);
                Vector q_plus = ZeroVector(size);
                Vector q_minus = ZeroVector(size);

                Element::GeometryType& r_geometry = r_element.GetGeometry();
                for (int i = 0; i < size; ++i)
                {
                    for (int j = 0; j < size; j++)
                    {
                        if (i != j)
                        {
                            const double f_ij = artificial_diffusion(i, j) *
                                                (values[j] - values[i]);

                            if (left_hand_side(j, i) <= left_hand_side(i, j))
                            {
                                p_plus[i] += std::max(0.0, f_ij);
                                p_minus[i] -= std::max(0.0, -f_ij);
                            }

                            if (equation_ids[i] < equation_ids[j])
                            {
                                q_plus[i] += std::max(0.0, -f_ij);
                                q_minus[i] -= std::max(0.0, f_ij);
                                q_plus[j] += std::max(0.0, f_ij);
                                q_minus[j] -= std::max(0.0, -f_ij);
                            }
                        }
                    }
                }

                for (int i = 0; i < size; ++i)
                {
                    ModelPart::NodeType& r_node = r_geometry[i];
                    r_node.SetLock();
                    r_node.FastGetSolutionStepValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX) +=
                        p_plus[i];
                    r_node.FastGetSolutionStepValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) +=
                        q_plus[i];
                    r_node.FastGetSolutionStepValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX) +=
                        p_minus[i];
                    r_node.FastGetSolutionStepValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) +=
                        q_minus[i];
                    r_node.UnSetLock();
                }
            }
        }

        if (mrPeriodicIdVar != Variable<int>::StaticObject())
        {
            const int number_of_conditions = rModelPart.NumberOfConditions();
#pragma omp parallel for
            for (int i_cond = 0; i_cond < number_of_conditions; ++i_cond)
            {
                ModelPart::ConditionType& r_condition =
                    *(rModelPart.ConditionsBegin() + i_cond);
                if (r_condition.Is(PERIODIC))
                {
                    ModelPart::NodeType& r_node_0 = r_condition.GetGeometry()[0];
                    ModelPart::NodeType& r_node_1 = r_condition.GetGeometry()[1];

                    double p_plus = r_node_0.FastGetSolutionStepValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX);
                    double q_plus = r_node_0.FastGetSolutionStepValue(
                        AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);
                    double p_minus =
                        r_node_0.FastGetSolutionStepValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX);
                    double q_minus = r_node_0.FastGetSolutionStepValue(
                        AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);

                    p_plus += r_node_1.FastGetSolutionStepValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX);
                    q_plus += r_node_1.FastGetSolutionStepValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);
                    p_minus += r_node_1.FastGetSolutionStepValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX);
                    q_minus += r_node_1.FastGetSolutionStepValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);

                    r_node_0.SetLock();
                    r_node_0.FastGetSolutionStepValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX) = p_plus;
                    r_node_0.FastGetSolutionStepValue(
                        AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) = q_plus;
                    r_node_0.FastGetSolutionStepValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX) = p_minus;
                    r_node_0.FastGetSolutionStepValue(
                        AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) = q_minus;
                    r_node_0.UnSetLock();

                    r_node_1.SetLock();
                    r_node_1.FastGetSolutionStepValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX) = p_plus;
                    r_node_1.FastGetSolutionStepValue(
                        AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) = q_plus;
                    r_node_1.FastGetSolutionStepValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX) = p_minus;
                    r_node_1.FastGetSolutionStepValue(
                        AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) = q_minus;
                    r_node_1.UnSetLock();
                }
            }
        }

        Communicator& r_communicator = rModelPart.GetCommunicator();
        r_communicator.AssembleCurrentData(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX);
        r_communicator.AssembleCurrentData(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);
        r_communicator.AssembleCurrentData(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX);
        r_communicator.AssembleCurrentData(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);

        KRATOS_CATCH("")
    }

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                TSystemMatrixType& rA,
                TSystemVectorType& rDx,
                TSystemVectorType& rb) override
    {
        KRATOS_TRY;

        mpDofUpdater->UpdateDofs(rDofSet, rDx, mRelaxationFactor);

        KRATOS_CATCH("");
    }

    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    void CalculateSystemContributions(Element::Pointer rCurrentElement,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      LocalSystemVectorType& RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        rCurrentElement->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentElement->CalculateLocalSystem(
            LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        Matrix SteadyLHS;
        rCurrentElement->CalculateLocalVelocityContribution(
            SteadyLHS, RHS_Contribution, CurrentProcessInfo);
        rCurrentElement->EquationIdVector(EquationId, CurrentProcessInfo);

        if (SteadyLHS.size1() != 0)
            noalias(LHS_Contribution) += SteadyLHS;

        Matrix artificial_diffusion;
        this->CalculateArtificialDiffusionMatrix(artificial_diffusion, LHS_Contribution);

        AddAntiDiffusiveFluxes(RHS_Contribution, LHS_Contribution,
                               *rCurrentElement, artificial_diffusion);
        noalias(LHS_Contribution) += artificial_diffusion;

        Vector U;
        rCurrentElement->GetValuesVector(U);
        noalias(RHS_Contribution) -= prod(LHS_Contribution, U);

        KRATOS_CATCH("");
    }

    void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
                                                LocalSystemMatrixType& LHS_Contribution,
                                                LocalSystemVectorType& RHS_Contribution,
                                                Condition::EquationIdVectorType& EquationId,
                                                ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        rCurrentCondition->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentCondition->CalculateLocalSystem(
            LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        Matrix SteadyLHS;
        rCurrentCondition->CalculateLocalVelocityContribution(
            SteadyLHS, RHS_Contribution, CurrentProcessInfo);
        rCurrentCondition->EquationIdVector(EquationId, CurrentProcessInfo);

        if (SteadyLHS.size1() != 0)
            noalias(LHS_Contribution) += SteadyLHS;

        KRATOS_CATCH("");
    }

    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
                                    LocalSystemVectorType& rRHS_Contribution,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        Matrix LHS_Contribution;
        CalculateSystemContributions(rCurrentElement, LHS_Contribution, rRHS_Contribution,
                                     rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
                                              LocalSystemVectorType& rRHS_Contribution,
                                              Element::EquationIdVectorType& rEquationId,
                                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        Matrix LHS_Contribution;
        Condition_CalculateSystemContributions(rCurrentCondition, LHS_Contribution,
                                               rRHS_Contribution, rEquationId,
                                               rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    ///@}

protected:
    ///@name Protected Operators
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    using DofUpdaterType = RelaxedDofUpdater<TSparseSpace>;
    using DofUpdaterPointerType = typename DofUpdaterType::UniquePointer;

    DofUpdaterPointerType mpDofUpdater = Kratos::make_unique<DofUpdaterType>();

    double mRelaxationFactor;

    const Variable<int>& mrPeriodicIdVar;

    template <typename TItem>
    void CalculateSystemMatrix(TItem& rItem, Matrix& rLeftHandSide, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        Vector rhs;

        rItem.InitializeNonLinearIteration(rCurrentProcessInfo);
        rItem.CalculateLocalSystem(rLeftHandSide, rhs, rCurrentProcessInfo);

        Matrix SteadyLHS;
        rItem.CalculateLocalVelocityContribution(SteadyLHS, rhs, rCurrentProcessInfo);

        if (SteadyLHS.size1() != 0)
            noalias(rLeftHandSide) += SteadyLHS;

        KRATOS_CATCH("");
    }

    void CalculateArtificialDiffusionMatrix(Matrix& rOutput, const Matrix& rInput)
    {
        const int size = rInput.size1();

        if (rOutput.size1() != size || rOutput.size2() != size)
        {
            rOutput.resize(size, size, false);
        }

        rOutput = ZeroMatrix(size, size);

        for (int i = 0; i < size; ++i)
        {
            for (int j = i + 1; j < size; ++j)
            {
                rOutput(i, j) = -std::max(std::max(rInput(i, j), rInput(j, i)), 0.0);
                rOutput(j, i) = rOutput(i, j);
            }
        }

        for (int i = 0; i < size; ++i)
        {
            double value = 0.0;
            for (int j = 0; j < size; ++j)
            {
                value -= rOutput(i, j);
            }
            rOutput(i, i) = value;
        }
    }

    template <typename TItem>
    void AddAntiDiffusiveFluxesToLHS(Vector& rRHS, Matrix& rLHS, const TItem& rItem, const Matrix& rArtificialDiffusion)
    {
        Vector values;
        rItem.GetValuesVector(values);

        const int size = rRHS.size();
        Matrix coeffs = ZeroMatrix(size, size);
        Matrix f = ZeroMatrix(size, size);

        for (int i = 0; i < size; ++i)
        {
            const ModelPart::NodeType& r_node_i = rItem.GetGeometry()[i];
            double r_plus_i{0.0}, r_minus_i{0.0};
            CalculateAntiDiffusiveFluxR(r_plus_i, r_minus_i, r_node_i);

            for (int j = 0; j < size; ++j)
            {
                if (i != j)
                {
                    if (rLHS(j, i) <= rLHS(i, j))
                    {
                        f(i, j) = rArtificialDiffusion(i, j) * (values[j] - values[i]);

                        if (f(i, j) > 0.0)
                        {
                            coeffs(i, j) = r_plus_i;
                        }
                        else if (f(i, j) < 0.0)
                        {
                            coeffs(i, j) = r_minus_i;
                        }
                        else
                        {
                            coeffs(i, j) = 1.0;
                        }
                        coeffs(j, i) = coeffs(i, j);
                    }
                }
            }
        }

        for (int i = 0; i < size; ++i)
        {
            double row_sum = 0.0;
            for (int j = 0; j < size; ++j)
            {
                const double value = (1.0 - coeffs(i, j)) * rArtificialDiffusion(i, j);
                rLHS(i, j) += value;
                row_sum += value;
            }
            rLHS(i, i) -= row_sum;
        }
    }

    template <typename TItem>
    void AddAntiDiffusiveFluxes(Vector& rRHS, Matrix& rLHS, TItem& rItem, const Matrix& rArtificialDiffusion)
    {
        KRATOS_TRY

        Vector values;
        rItem.GetValuesVector(values);

        const int size = rRHS.size();
        Matrix coeffs = ZeroMatrix(size, size);
        Matrix f = ZeroMatrix(size, size);

        for (int i = 0; i < size; ++i)
        {
            const ModelPart::NodeType& r_node_i = rItem.GetGeometry()[i];
            double r_plus_i{0.0}, r_minus_i{0.0};
            CalculateAntiDiffusiveFluxR(r_plus_i, r_minus_i, r_node_i);

            for (int j = 0; j < size; ++j)
            {
                if (i != j)
                {
                    f(i, j) = rArtificialDiffusion(i, j) * (values[j] - values[i]);

                    if (rLHS(j, i) <= rLHS(i, j))
                    {
                        if (f(i, j) > 0.0)
                        {
                            coeffs(i, j) = r_plus_i;
                        }
                        else if (f(i, j) < 0.0)
                        {
                            coeffs(i, j) = r_minus_i;
                        }
                        else
                        {
                            coeffs(i, j) = 1.0;
                        }
                        coeffs(j, i) = coeffs(i, j);
                    }
                }
            }
        }
        Vector temp = ZeroVector(size);
        for (int i = 0; i < size; ++i)
        {
            temp[i] = rRHS[i];
            for (int j = 0; j < size; ++j)
            {
                temp[i] += coeffs(i, j) * f(i, j);
            }
        }

        bool found_negatives = false;
        // for (int i = 0; i < size; ++i)
        // {
        //     if (temp[i] < 0.0 && !rItem.GetGeometry()[i].Is(INLET))
        //     {
        //         found_negatives = true;
        //     }
        // }

        if (!found_negatives)
        {
            noalias(rRHS) = temp;
        }

        KRATOS_CATCH("");
    }

    // template <typename TItem>
    // void AddLumpedMassMatrix(Matrix& rLHS, TItem& rItem, const ProcessInfo& rCurrentProcessInfo) const
    // {
    //     const auto& r_geometry = rItem.GetGeometry();
    //     const int number_of_nodes = r_geometry.PointsNumber();
    //     double residual = 0.0;
    //     rItem.Calculate(RESIDUAL, residual, rCurrentProcessInfo);qu

    //     const double coefficient = CalculatePositivityPreservingMatrix(rLHS);

    //     IdentityMatrix identity_matrix(number_of_nodes);
    //     noalias(rLHS) += identity_matrix * (coefficient * residual);
    // }

    inline double CalculatePositivityPreservingMatrix(const Matrix& rInputMatrix) const
    {
        double coefficient = 0.0;
        for (unsigned int a = 0; a < rInputMatrix.size1(); ++a)
        {
            double row_sum = 0.0;
            for (unsigned int b = 0; b < rInputMatrix.size2(); ++b)
            {
                row_sum += rInputMatrix(a, b);
            }
            coefficient = std::max(coefficient, -row_sum);
        }
        return coefficient;
    }

    void CalculateAntiDiffusiveFluxR(double& rRPlus,
                                     double& rRMinus,
                                     const ModelPart::NodeType& rNode) const
    {
        if (rNode.Is(SLIP) || rNode.Is(INLET))
        {
            rRMinus = 1.0;
            rRPlus = 1.0;
        }
        else
        {
            const double q_plus =
                rNode.FastGetSolutionStepValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);
            const double p_plus =
                rNode.FastGetSolutionStepValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX);
            const double q_minus =
                rNode.FastGetSolutionStepValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);
            const double p_minus =
                rNode.FastGetSolutionStepValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX);

            rRPlus = 1.0;
            if (p_plus > 0.0)
            {
                rRPlus = std::min(1.0, q_plus / p_plus);
            }

            rRMinus = 1.0;
            if (p_minus < 0.0)
            {
                rRMinus = std::min(1.0, q_minus / p_minus);
            }
        }
    }

    ///@}
}; // namespace Kratos

///@}

} // namespace Kratos

#endif /* KRATOS_ALGEBRAIC_FLUX_CORRECTED_SCALAR_STEADY_SCHEME defined */

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_ALGEBRAIC_FLUX_CORRECTED_SCALAR_STEADY_SCHEME)
#define KRATOS_ALGEBRAIC_FLUX_CORRECTED_SCALAR_STEADY_SCHEME

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/openmp_utils.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_strategies/relaxed_dof_updater.h"
#include "rans_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @brief Algebraic flux corrected scalar steady transport scheme.
 *
 * This scheme is based on following publication.
 *
 * D. Kuzmin, Algebraic flux correction for finite element discretizations of coupled systems,
 * Computational Methods for Coupled Problems in Science and Engineering II, CIMNE,
 * Barcelona, (2007), pp. 653–656.
 *
 * This scheme can only be used to solve steady state problems with with elements derrived
 * from ConvectionDiffusionReactionElement.
 *
 * @tparam TSparseSpace                     Sparse space type
 * @tparam TDenseSpace                      Dense space type
 *
 * @see ConvectionDiffusionReactionElement
 */
template <class TSparseSpace, class TDenseSpace>
class AlgebraicFluxCorrectedSteadyScalarScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AlgebraicFluxCorrectedSteadyScalarScheme);

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;

    using DofsArrayType = typename BaseType::DofsArrayType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    AlgebraicFluxCorrectedSteadyScalarScheme(
        const double RelaxationFactor,
        const Flags BoundaryFlags)
    : BaseType(),
        mRelaxationFactor(RelaxationFactor),
        mBoundaryFlags(BoundaryFlags),
        mrPeriodicIdVar(Variable<int>::StaticObject())
    {
        KRATOS_INFO("AlgebraicFluxCorrectedSteadyScalarScheme")
            << " Using residual based algebraic flux corrected scheme with "
               "relaxation "
               "factor = "
            << std::scientific << mRelaxationFactor << "\n";

        mpDofUpdater = Kratos::make_unique<DofUpdaterType>(mRelaxationFactor);
    }

    AlgebraicFluxCorrectedSteadyScalarScheme(
        const double RelaxationFactor,
        const Flags BoundaryFlags,
        const Variable<int>& rPeriodicIdVar)
    : BaseType(),
        mRelaxationFactor(RelaxationFactor),
        mBoundaryFlags(BoundaryFlags),
        mrPeriodicIdVar(rPeriodicIdVar)
    {
        KRATOS_INFO("AlgebraicFluxCorrectedSteadyScalarScheme")
            << " Using periodic residual based algebraic flux corrected scheme "
               "with relaxation "
               "factor = "
            << std::scientific << mRelaxationFactor << "\n";

        mpDofUpdater = Kratos::make_unique<DofUpdaterType>(mRelaxationFactor);
    }

    ~AlgebraicFluxCorrectedSteadyScalarScheme() override = default;

    ///@}
    ///@name Operators
    ///@{

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        BaseType::Initialize(rModelPart);

        block_for_each(rModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
            rNode.SetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX, 0.0);
            rNode.SetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX, 0.0);
            rNode.SetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT, 0.0);
            rNode.SetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT, 0.0);
        });

        if (mrPeriodicIdVar != Variable<int>::StaticObject()) {
            block_for_each(rModelPart.Conditions(), [&](const ModelPart::ConditionType& rCondition) {
                if (rCondition.Is(PERIODIC)) {
                    // this only supports 2 noded periodic conditions
                    KRATOS_ERROR_IF(rCondition.GetGeometry().PointsNumber() != 2)
                        << this->Info() << " only supports two noded periodic conditions. Found "
                        << rCondition.Info() << " with "
                        << rCondition.GetGeometry().PointsNumber() << " nodes.\n";

                    const auto& r_node_0 = rCondition.GetGeometry()[0];
                    const std::size_t r_node_0_pair_id =
                        r_node_0.FastGetSolutionStepValue(mrPeriodicIdVar);

                    const auto& r_node_1 = rCondition.GetGeometry()[1];
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
            });
        }

        // Allocate auxiliary memory.
        const auto num_threads = OpenMPUtils::GetNumThreads();
        mAntiDiffusiveFlux.resize(num_threads);
        mAntiDiffusiveFluxCoefficients.resize(num_threads);
        mValues.resize(num_threads);
        mAuxMatrix.resize(num_threads);

        KRATOS_CATCH("");
    }

    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        auto& r_nodes = rModelPart.Nodes();

        block_for_each(r_nodes, [&](ModelPart::NodeType& rNode) {
            rNode.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX) = 0.0;
            rNode.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX) = 0.0;
            rNode.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) = 0.0;
            rNode.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) = 0.0;
        });

        auto& r_elements = rModelPart.Elements();
        const int number_of_elements = r_elements.size();

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

#pragma omp parallel
        {
            Matrix left_hand_side, artificial_diffusion, aux_matrix;
            Vector right_hand_side, values;
            std::vector<IndexType> equation_ids;
#pragma omp for
            for (int i = 0; i < number_of_elements; ++i) {
                auto& r_element = *(r_elements.begin() + i);
                this->CalculateSystemMatrix<Element>(r_element, left_hand_side,
                                                     right_hand_side, aux_matrix,
                                                     r_current_process_info);
                this->CalculateArtificialDiffusionMatrix(artificial_diffusion, left_hand_side);
                r_element.EquationIdVector(equation_ids, r_current_process_info);
                r_element.GetValuesVector(values);

                const int size = artificial_diffusion.size1();

                Vector p_plus = ZeroVector(size);
                Vector p_minus = ZeroVector(size);
                Vector q_plus = ZeroVector(size);
                Vector q_minus = ZeroVector(size);

                auto& r_geometry = r_element.GetGeometry();
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < size; j++) {
                        if (i != j) {
                            const double f_ij = artificial_diffusion(i, j) *
                                                (values[j] - values[i]);

                            if (left_hand_side(j, i) <= left_hand_side(i, j)) {
                                p_plus[i] += std::max(0.0, f_ij);
                                p_minus[i] -= std::max(0.0, -f_ij);
                            }

                            if (equation_ids[i] < equation_ids[j]) {
                                q_plus[i] += std::max(0.0, -f_ij);
                                q_minus[i] -= std::max(0.0, f_ij);
                                q_plus[j] += std::max(0.0, f_ij);
                                q_minus[j] -= std::max(0.0, -f_ij);
                            }
                        }
                    }
                }

                for (int i = 0; i < size; ++i) {
                    auto& r_node = r_geometry[i];
                    r_node.SetLock();
                    r_node.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX) += p_plus[i];
                    r_node.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) += q_plus[i];
                    r_node.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX) += p_minus[i];
                    r_node.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) += q_minus[i];
                    r_node.UnSetLock();
                }
            }
        }

        if (mrPeriodicIdVar != Variable<int>::StaticObject()) {
            block_for_each(rModelPart.Conditions(), [&](ModelPart::ConditionType& rCondition) {
                if (rCondition.Is(PERIODIC)) {
                    auto& r_node_0 = rCondition.GetGeometry()[0];
                    auto& r_node_1 = rCondition.GetGeometry()[1];

                    double p_plus = r_node_0.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX);
                    double q_plus = r_node_0.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);
                    double p_minus = r_node_0.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX);
                    double q_minus = r_node_0.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);

                    p_plus += r_node_1.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX);
                    q_plus += r_node_1.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);
                    p_minus += r_node_1.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX);
                    q_minus += r_node_1.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);

                    r_node_0.SetLock();
                    r_node_0.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX) = p_plus;
                    r_node_0.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) = q_plus;
                    r_node_0.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX) = p_minus;
                    r_node_0.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) = q_minus;
                    r_node_0.UnSetLock();

                    r_node_1.SetLock();
                    r_node_1.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX) = p_plus;
                    r_node_1.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) = q_plus;
                    r_node_1.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX) = p_minus;
                    r_node_1.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT) = q_minus;
                    r_node_1.UnSetLock();
                }
            });
        }

        Communicator& r_communicator = rModelPart.GetCommunicator();
        r_communicator.AssembleNonHistoricalData(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX);
        r_communicator.AssembleNonHistoricalData(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);
        r_communicator.AssembleNonHistoricalData(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX);
        r_communicator.AssembleNonHistoricalData(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);

        KRATOS_CATCH("")
    }

    void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        KRATOS_TRY;

        mpDofUpdater->UpdateDofs(rDofSet, rDx);

        KRATOS_CATCH("");
    }

    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    void CalculateSystemContributions(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const auto k = OpenMPUtils::ThisThread();

        this->CalculateSystemMatrix<Element>(rElement, rLHS_Contribution, rRHS_Contribution,
                                             mAuxMatrix[k], rCurrentProcessInfo);
        rElement.EquationIdVector(rEquationIdVector, rCurrentProcessInfo);

        this->CalculateArtificialDiffusionMatrix(mAuxMatrix[k], rLHS_Contribution);

        AddAntiDiffusiveFluxes(rRHS_Contribution, rLHS_Contribution, rElement,
                               mAuxMatrix[k]);
        noalias(rLHS_Contribution) += mAuxMatrix[k];

        rElement.GetValuesVector(mValues[k]);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, mValues[k]);

        KRATOS_CATCH("");
    }

    void CalculateSystemContributions(
        Condition& rCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const auto k = OpenMPUtils::ThisThread();

        this->CalculateSystemMatrix<Condition>(rCondition, rLHS_Contribution, rRHS_Contribution,
                                               mAuxMatrix[k], rCurrentProcessInfo);
        rCondition.EquationIdVector(rEquationIdVector, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateRHSContribution(
        Element& rElement,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const auto k = OpenMPUtils::ThisThread();
        CalculateSystemContributions(rElement, mAuxMatrix[k], rRHS_Contribution,
                                     rEquationIdVector, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateRHSContribution(
        Condition& rCondition,
        LocalSystemVectorType& rRHS_Contribution,
        Condition::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const auto k = OpenMPUtils::ThisThread();
        CalculateSystemContributions(rCondition, mAuxMatrix[k], rRHS_Contribution,
                                     rEquationIdVector, rCurrentProcessInfo);

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

    DofUpdaterPointerType mpDofUpdater;

    double mRelaxationFactor;
    const Flags mBoundaryFlags;
    const Variable<int>& mrPeriodicIdVar;

    std::vector<LocalSystemMatrixType> mAuxMatrix;
    std::vector<LocalSystemMatrixType> mAntiDiffusiveFluxCoefficients;
    std::vector<LocalSystemMatrixType> mAntiDiffusiveFlux;
    std::vector<LocalSystemVectorType> mValues;

    /**
     * @brief Common method to calculate Element and Condition system matrices
     *
     * @tparam TItem                Type of item (can be ElementType or ConditionType)
     * @param rItem                 Item instance
     * @param rLeftHandSide         Lefthandside matrix
     * @param rRightHandSide        Righthandside vector
     * @param rAuxMatrix            Auxiliary matrix
     * @param rCurrentProcessInfo   Current process info
     */
    template <typename TItem>
    void CalculateSystemMatrix(
        TItem& rItem,
        LocalSystemMatrixType& rLeftHandSide,
        LocalSystemVectorType& rRightHandSide,
        LocalSystemMatrixType& rAuxMatrix,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rItem.InitializeNonLinearIteration(rCurrentProcessInfo);
        rItem.CalculateLocalSystem(rLeftHandSide, rRightHandSide, rCurrentProcessInfo);
        rItem.CalculateLocalVelocityContribution(rAuxMatrix, rRightHandSide, rCurrentProcessInfo);

        if (rAuxMatrix.size1() != 0) {
            noalias(rLeftHandSide) += rAuxMatrix;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates artificial diffusion matrix for given discretized matrix
     *
     * @param rOutput   Diffusion matrix
     * @param rInput    Input matrix
     */
    void CalculateArtificialDiffusionMatrix(
        Matrix& rOutput,
        const Matrix& rInput)
    {
        const IndexType size = rInput.size1();

        if (rOutput.size1() != size || rOutput.size2() != size) {
            rOutput.resize(size, size, false);
        }

        rOutput = ZeroMatrix(size, size);

        for (IndexType i = 0; i < size; ++i) {
            for (IndexType j = i + 1; j < size; ++j) {
                rOutput(i, j) = -std::max(std::max(rInput(i, j), rInput(j, i)), 0.0);
                rOutput(j, i) = rOutput(i, j);
            }
        }

        for (IndexType i = 0; i < size; ++i) {
            double value = 0.0;
            for (IndexType j = 0; j < size; ++j) {
                value -= rOutput(i, j);
            }
            rOutput(i, i) = value;
        }
    }

    /**
     * @brief Calculates anti-diffusive terms
     *
     * Diffusion calculated by CalculateArtificialDiffusionMatrix alters original problem. Therefore
     * anti-diffusion terms are calculated to cancel diffusion terms where they are not necessary for
     * stabilization of the Convection-Diffusion-Reaction scalar equation.
     *
     * @tparam TItem                Item type (can be ElementType or ConditionType)
     * @param rRHS                  Righthandside vector
     * @param rLHS                  Lefthandside matrix
     * @param rItem                 Item instance
     * @param rArtificialDiffusion  Calculated artificial diffusion
     */
    template <typename TItem>
    void AddAntiDiffusiveFluxes(
        Vector& rRHS,
        const Matrix& rLHS,
        TItem& rItem,
        const Matrix& rArtificialDiffusion)
    {
        KRATOS_TRY

        const auto k = OpenMPUtils::ThisThread();
        const auto size = rRHS.size();

        auto& r_anti_diffusive_flux_coefficients = mAntiDiffusiveFluxCoefficients[k];
        auto& r_anti_diffusive_flux = mAntiDiffusiveFlux[k];
        auto& r_values = mValues[k];

        rItem.GetValuesVector(r_values);
        if (r_anti_diffusive_flux_coefficients.size1() != size ||
            r_anti_diffusive_flux_coefficients.size2() != size) {
            r_anti_diffusive_flux_coefficients.resize(size, size, false);
        }

        if (r_anti_diffusive_flux.size1() != size || r_anti_diffusive_flux.size2() != size) {
            r_anti_diffusive_flux.resize(size, size, false);
        }

        noalias(r_anti_diffusive_flux_coefficients) = ZeroMatrix(size, size);
        noalias(r_anti_diffusive_flux) = ZeroMatrix(size, size);

        for (IndexType i = 0; i < size; ++i) {
            const auto& r_node_i = rItem.GetGeometry()[i];
            double r_plus_i{0.0}, r_minus_i{0.0};
            CalculateAntiDiffusiveFluxR(r_plus_i, r_minus_i, r_node_i);

            for (IndexType j = 0; j < size; ++j) {
                if (i != j) {
                    r_anti_diffusive_flux(i, j) =
                        rArtificialDiffusion(i, j) * (r_values[j] - r_values[i]);

                    if (rLHS(j, i) <= rLHS(i, j)) {
                        if (r_anti_diffusive_flux(i, j) > 0.0) {
                            r_anti_diffusive_flux_coefficients(i, j) = r_plus_i;
                        } else if (r_anti_diffusive_flux(i, j) < 0.0) {
                            r_anti_diffusive_flux_coefficients(i, j) = r_minus_i;
                        } else {
                            r_anti_diffusive_flux_coefficients(i, j) = 1.0;
                        }
                        r_anti_diffusive_flux_coefficients(j, i) =
                            r_anti_diffusive_flux_coefficients(i, j);
                    }
                }
            }
        }

        for (IndexType i = 0; i < size; ++i) {
            for (IndexType j = 0; j < size; ++j) {
                rRHS[i] += r_anti_diffusive_flux_coefficients(i, j) *
                           r_anti_diffusive_flux(i, j);
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates allowed artifical diffusive fluxes
     *
     * @param rRPlus        Allowed positive fluxes
     * @param rRMinus       Allowed negative fluxes
     * @param rNode         Node
     */
    void CalculateAntiDiffusiveFluxR(
        double& rRPlus,
        double& rRMinus,
        const ModelPart::NodeType& rNode) const
    {
        if (rNode.Is(mBoundaryFlags)) {
            rRMinus = 1.0;
            rRPlus = 1.0;
        } else {
            const double q_plus = rNode.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);
            const double p_plus = rNode.GetValue(AFC_POSITIVE_ANTI_DIFFUSIVE_FLUX);
            const double q_minus = rNode.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX_LIMIT);
            const double p_minus = rNode.GetValue(AFC_NEGATIVE_ANTI_DIFFUSIVE_FLUX);

            rRPlus = 1.0;
            if (p_plus > 0.0) {
                rRPlus = std::min(1.0, q_plus / p_plus);
            }

            rRMinus = 1.0;
            if (p_minus < 0.0) {
                rRMinus = std::min(1.0, q_minus / p_minus);
            }
        }
    }

    ///@}
}; // namespace Kratos

///@}

} // namespace Kratos

#endif /* KRATOS_ALGEBRAIC_FLUX_CORRECTED_SCALAR_STEADY_SCHEME defined */

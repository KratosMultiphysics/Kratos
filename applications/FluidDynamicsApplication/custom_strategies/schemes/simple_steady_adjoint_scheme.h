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

#if !defined(KRATOS_SIMPLE_STEADY_ADJOINT_SCHEME_H_INCLUDED)
#define KRATOS_SIMPLE_STEADY_ADJOINT_SCHEME_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// A scheme for steady adjoint equations.
/**
 * Solves the adjoint equations:
 * \f[
 *  \partial_{\mathbf{w}}\mathbf{r}^T \lambda = -\partial_{\mathbf{w}}J^{T}
 * \f]
 *
 * \f$\lambda\f$ is returned by Element::GetFirstDerivativesVector.
 * \f$\partial_{\mathbf{w}}\mathbf{r}^T\f$ is returned by Element::CalculateFirstDerivativesLHS.
 * \f$\partial_{\mathbf{w}}J^{T}\f$ is returned by ResponseFunction::CalculateFirstDerivativesGradient.
 */
template <IndexType TDim, class TSparseSpace, class TDenseSpace>
class SimpleSteadyAdjointScheme : public ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SimpleSteadyAdjointScheme);

    using BaseType = ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using IndexType = std::size_t;

    using NodeType = ModelPart::NodeType;

    static constexpr IndexType TBlockSize = TDim + 1;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    explicit SimpleSteadyAdjointScheme(
        AdjointResponseFunction::Pointer pResponseFunction)
        : BaseType(pResponseFunction)
    {
        // Allocate auxiliary memory.
        const int number_of_threads = OpenMPUtils::GetNumThreads();
        mAuxMatrices.resize(number_of_threads);
    }

    /// Destructor.
    ~SimpleSteadyAdjointScheme() override = default;

    ///@}
    ///@name Operations
    ///@{

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntitySystemContributions(rCurrentElement, rLHS_Contribution, rRHS_Contribution,
                                           rEquationId, rCurrentProcessInfo);
    }

    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityLHSContribution(rCurrentElement, rLHS_Contribution,
                                       rEquationId, rCurrentProcessInfo);
    }

    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntitySystemContributions(rCurrentCondition, rLHS_Contribution, rRHS_Contribution,
                                           rEquationId, rCurrentProcessInfo);
    }

    void CalculateLHSContribution(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityLHSContribution(rCurrentCondition, rLHS_Contribution,
                                       rEquationId, rCurrentProcessInfo);
    }

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    std::vector<Matrix> mAuxMatrices;

    ///@}
    ///@name Private Operations
    ///@{

    template<class TEntityType>
    void CalculateEntitySystemContributions(
        TEntityType& rEntity,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        typename TEntityType::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const auto thread_id = OpenMPUtils::ThisThread();

        const auto& r_const_entity_ref = rEntity;

        CalculateEntityLHSContribution<TEntityType>(
            rEntity, rLHS_Contribution, rEquationId, rCurrentProcessInfo);

        if (rRHS_Contribution.size() != rLHS_Contribution.size1())
            rRHS_Contribution.resize(rLHS_Contribution.size1(), false);

        this->mpResponseFunction->CalculateFirstDerivativesGradient(
            rEntity, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        noalias(rRHS_Contribution) = -rRHS_Contribution;

        // Calculate system contributions in residual form.
        if (rLHS_Contribution.size1() != 0) {
            auto& adjoint_values = this->mAdjointValues[thread_id];
            r_const_entity_ref.GetValuesVector(adjoint_values);
            noalias(rRHS_Contribution) -= prod(rLHS_Contribution, adjoint_values);
        }

        KRATOS_CATCH("");
    }

    template<class TEntityType>
    void CalculateEntityLHSContribution(
        TEntityType& rEntity,
        LocalSystemMatrixType& rLHS,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const auto thread_id = OpenMPUtils::ThisThread();
        auto& aux_matrix = mAuxMatrices[thread_id];

        const auto& r_const_entity_ref = rEntity;
        rEntity.CalculateFirstDerivativesLHS(aux_matrix, rCurrentProcessInfo);
        r_const_entity_ref.EquationIdVector(rEquationId, rCurrentProcessInfo);

        const auto& r_geometry = rEntity.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();

        if (rLHS.size1() != aux_matrix.size1() || rLHS.size2() != aux_matrix.size2()) {
            rLHS.resize(aux_matrix.size1(), aux_matrix.size2(), false);
        }

        rLHS.clear();

        // add residual derivative contributions
        for (IndexType a = 0; a < number_of_nodes; ++a) {
            const auto& r_node = r_geometry[a];
            if (r_node.Is(SLIP)) {
                AddNodalRotatedResidualDerivativeToMatrix(
                    rLHS, aux_matrix, a * TBlockSize, r_node);
            } else {
                AddNodalResidualDerivativeToMatrix(rLHS, aux_matrix, a * TBlockSize);
            }
        }

        KRATOS_CATCH("");
    }

    void AddNodalRotatedResidualDerivativeToMatrix(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const IndexType NodeStartIndex,
        const NodeType& rNode) const
    {
        KRATOS_TRY

        using coordinate_transformation_utils = CoordinateTransformationUtils<Matrix, Vector, double>;

        BoundedVector<double, TDim> residual_derivative, aux_vector;
        BoundedMatrix<double, TDim, TDim> rotation_matrix;
        coordinate_transformation_utils::LocalRotationOperatorPure(rotation_matrix, rNode);

        // add rotated residual derivative contributions
        for (IndexType c = 0; c < rResidualDerivatives.size1(); ++c) {
            // get the residual derivative relevant for node
            FluidCalculationUtilities::ReadSubVector<TDim>(
                residual_derivative, row(rResidualDerivatives, c), NodeStartIndex);

            // rotate residual derivative
            noalias(aux_vector) = prod(rotation_matrix, residual_derivative);

            // add rotated residual derivative to local matrix
            FluidCalculationUtilities::AddSubVector<TDim>(
                rOutput, aux_vector, c, NodeStartIndex);

            // add continuity equation derivatives
            rOutput(c, NodeStartIndex + TDim) +=
                rResidualDerivatives(c, NodeStartIndex + TDim);
        }

        // Apply slip condition in primal scheme makes first momentum dof
        // fixed, making the velocity in the normal direction as rhs.

        // first clear the residual derivative
        for (IndexType c = 0; c < rOutput.size1(); ++c) {
            rOutput(c, NodeStartIndex) = 0.0;
        }

        auto normal = rNode.FastGetSolutionStepValue(NORMAL);
        normal /= norm_2(normal);

        for (IndexType i = 0; i < TDim; ++i) {
            rOutput(NodeStartIndex + i, NodeStartIndex) -= normal[i];
        }

        KRATOS_CATCH("");
    }

    void AddNodalResidualDerivativeToMatrix(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const IndexType NodeStartIndex) const
    {
        KRATOS_TRY

        // add non-rotated residual derivative contributions
        for (IndexType c = 0; c < rResidualDerivatives.size1(); ++c) {
            for (IndexType i = 0; i < TBlockSize; ++i) {
                rOutput(c, NodeStartIndex + i) +=
                    rResidualDerivatives(c, NodeStartIndex + i);
            }
        }

        KRATOS_CATCH("");
    }

    ///@}

}; /* Class SimpleSteadyAdjointScheme */

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_SIMPLE_STEADY_ADJOINT_SCHEME_H_INCLUDED defined */

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
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"

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
template <class TSparseSpace, class TDenseSpace>
class SimpleSteadyAdjointScheme : public ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SimpleSteadyAdjointScheme);

    using BaseType = ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    using IndexType = std::size_t;

    using NodeType = ModelPart::NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    explicit SimpleSteadyAdjointScheme(
        AdjointResponseFunction::Pointer pResponseFunction)
        : ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>(pResponseFunction)
    {
        int num_threads = OpenMPUtils::GetNumThreads();
        mAuxMatrices.resize(num_threads);
    }

    /// Destructor.
    ~SimpleSteadyAdjointScheme() override
    {
    }

    ///@}
    ///@name Operations
    ///@{

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        const IndexType domain_size = rModelPart.GetProcessInfo()[DOMAIN_SIZE];

        // assign domain size specific methods
        if (domain_size == 2) {
            this->mAddNodalRotatedResidualDerivativeToMatrix =
                &SimpleSteadyAdjointScheme::AddNodalRotatedResidualDerivativeToMatrix<2>;
            this->mAddNodalResidualDerivativeToMatrix =
                &SimpleSteadyAdjointScheme::AddNodalResidualDerivativeToMatrix<2>;
        } else if (domain_size == 3) {
            this->mAddNodalRotatedResidualDerivativeToMatrix =
                &SimpleSteadyAdjointScheme::AddNodalRotatedResidualDerivativeToMatrix<3>;
            this->mAddNodalResidualDerivativeToMatrix =
                &SimpleSteadyAdjointScheme::AddNodalResidualDerivativeToMatrix<3>;
        } else {
            KRATOS_ERROR << "Unsupported domain size [ domain_size = " << domain_size
                         << " ].\n";
        }

        BaseType::Initialize(rModelPart);

        // fixing first dof of velocity on slip nodes because when system of equations
        // are rotated according to nodal normal, first dof becomes velocity in the
        // nodal normal direction which is fixed by matrix construction in the
        // primal scheme
        block_for_each(rModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
            if (rNode.Is(SLIP)) {
                rNode.Fix(ADJOINT_FLUID_VECTOR_1_X);
                rNode.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1_X) = 0.0;
            }
        });

        KRATOS_CATCH("");
    }

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
        KRATOS_TRY;

        const auto thread_id = OpenMPUtils::ThisThread();

        CalculateEntityLHSContribution(rCurrentElement, rLHS_Contribution,
                                       mAuxMatrices[thread_id], rEquationId,
                                       rCurrentProcessInfo);

        KRATOS_CATCH("");
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
        KRATOS_TRY;

        const auto thread_id = OpenMPUtils::ThisThread();

        CalculateEntityLHSContribution(rCurrentCondition, rLHS_Contribution,
                                       mAuxMatrices[thread_id], rEquationId,
                                       rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    std::vector<Matrix> mAuxMatrices;

    void (SimpleSteadyAdjointScheme::*mAddNodalRotatedResidualDerivativeToMatrix)(
        Matrix&,
        const Matrix&,
        const IndexType,
        const NodeType&);

    void (SimpleSteadyAdjointScheme::*mAddNodalResidualDerivativeToMatrix)(
        Matrix&, const Matrix&, const IndexType);


    ///@}
    ///@name Private Operations
    ///@{

    template<class TEntityType>
    void CalculateEntityLHSContribution(
        TEntityType& rEntity,
        LocalSystemMatrixType& rRotatedLHS,
        LocalSystemMatrixType& rLHS,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const auto& r_const_entity_ref = rEntity;
        rEntity.CalculateFirstDerivativesLHS(rLHS, rCurrentProcessInfo);

        if (rRotatedLHS.size1() != rLHS.size1() || rRotatedLHS.size2() != rLHS.size2()) {
            rRotatedLHS.resize(rLHS.size1(), rLHS.size2(), false);
        }

        rRotatedLHS.clear();

        auto& r_geometry = rEntity.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();
        const IndexType local_size = rLHS.size1() / number_of_nodes;

        // add residual derivative contributions
        for (IndexType a = 0; a < number_of_nodes; ++a) {
            const auto& r_node = r_geometry[a];
            if (r_node.Is(SLIP)) {
                (this->*(this->mAddNodalRotatedResidualDerivativeToMatrix))(
                    rRotatedLHS, rLHS, a * local_size, r_node);
            } else {
                (this->*(this->mAddNodalResidualDerivativeToMatrix))(
                    rRotatedLHS, rLHS, a * local_size);
            }
        }

        r_const_entity_ref.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

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

        auto& aux_matrix = mAuxMatrices[thread_id];

        const auto& r_const_entity_ref = rEntity;

        CalculateEntityLHSContribution<TEntityType>(rEntity, rLHS_Contribution, aux_matrix, rEquationId, rCurrentProcessInfo);

        if (rRHS_Contribution.size() != rLHS_Contribution.size1())
            rRHS_Contribution.resize(rLHS_Contribution.size1(), false);

        this->mpResponseFunction->CalculateFirstDerivativesGradient(
            rEntity, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        noalias(rRHS_Contribution) = -rRHS_Contribution;

        // Calculate system contributions in residual form.
        r_const_entity_ref.GetFirstDerivativesVector(this->mAdjointValues[thread_id]);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, this->mAdjointValues[thread_id]);

        KRATOS_CATCH("");
    }

    template <unsigned int TDim>
    void AddNodalRotatedResidualDerivativeToMatrix(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const IndexType NodeStartIndex,
        const NodeType& rNode)
    {
        KRATOS_TRY

        using coordinate_transformation_utils = CoordinateTransformationUtils<Matrix, Vector, double>;

        BoundedVector<double, TDim> residual_derivative, aux_vector;

        // get the rotation matrix relevant for rNode
        BoundedMatrix<double, TDim, TDim> rotation_matrix;
        coordinate_transformation_utils::LocalRotationOperatorPure(rotation_matrix, rNode);

        // add rotated residual derivative contributions
        for (IndexType c = 0; c < rResidualDerivatives.size1(); ++c) {
            // get the residual derivative relevant for node
            FluidCalculationUtilities::GetSubVector<TDim>(
                residual_derivative,
                row(rResidualDerivatives, c), NodeStartIndex);

            // rotate residual derivative
            noalias(aux_vector) = prod(rotation_matrix, residual_derivative);

            // add rotated residual derivative to local matrix
            FluidCalculationUtilities::AddSubVectorToMatrix<TDim>(
                rOutput, aux_vector, c, NodeStartIndex);

            // add continuity equation derivatives
            rOutput(c, NodeStartIndex + TDim) +=
                rResidualDerivatives(c, NodeStartIndex + TDim);
        }

        KRATOS_CATCH("");
    }

    template <unsigned int TDim>
    void AddNodalResidualDerivativeToMatrix(
        Matrix& rOutput,
        const Matrix& rResidualDerivatives,
        const IndexType NodeStartIndex)
    {
        KRATOS_TRY

        // add non-rotated residual derivative contributions
        for (IndexType c = 0; c < rResidualDerivatives.size1(); ++c) {
            for (IndexType i = 0; i < TDim + 1; ++i) {
                rOutput(c, NodeStartIndex + i) +=
                    rResidualDerivatives(c, NodeStartIndex + i);
            }
        }

        KRATOS_CATCH("");
    }

    ///@}

}; /* Class SimpleSteadyAdjointScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_SIMPLE_STEADY_ADJOINT_SCHEME_H_INCLUDED defined */

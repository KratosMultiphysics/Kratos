//  KratosAdjointFluidApplication
//
//  License:         BSD License
//                   license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_ADJOINT_BOSSAK_SCHEME)
#define KRATOS_ADJOINT_BOSSAK_SCHEME

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/communicator.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"
#include "solving_strategies/schemes/scheme.h"
#include "containers/variable.h"
#include "solving_strategies/response_functions/response_function.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// A scheme for unsteady adjoint equations using Bossak time discretization.
/**
 * The forward Bossak equations are:
 * \f[
 * \mathbf{M}\dot{\mathbf{w}}^{n-\alpha} = \mathbf{f}(\mathbf{w}^{n};\mathbf{s})
 * \f]
 * \f[
 * \dot{\mathbf{w}}^{n-\alpha}
 * = (1 - \alpha) \dot{\mathbf{w}}^n + \alpha \dot{\mathbf{w}}^{n-1}
 * \f]
 * \f[
 * \dot{\mathbf{w}}^n
 * = \frac{\mathbf{w}^n - \mathbf{w}^{n-1}}{\gamma \Delta t}
 * + \frac{\gamma - 1}{\gamma}\dot{\mathbf{w}}^{n-1}
 * \f]
 *
 * The adjoint Bossak equations are:
 * \f[
 * \frac{1}{\gamma - 1} (\dot{\lambda}^n - \dot{\lambda}^{n+1})
 * + (\partial_{\mathbf{w}^n}\mathbf{f}^n
 * -\partial_{\mathbf{w}^n}(\mathbf{M}^n\dot{\mathbf{w}}^{n-\alpha}))^T\lambda^n
 * = -\partial_{\mathbf{w}^n}J^{nT}
 * \f]
 * \f[
 * \frac{1}{\gamma - 1} \dot{\lambda}^n
 * = \frac{1}{\gamma} \dot{\lambda}^{n+1}
 * - \frac{1 - \alpha}{\gamma \Delta t}M^{nT} \lambda^n
 * - \frac{\alpha}{\gamma \Delta t}M^{(n+1)T} \lambda^{n+1}
 * + \frac{1}{\gamma \Delta t}\partial_{\dot{\mathbf{w}}^n}J^{nT}
 * + \frac{1}{\gamma \Delta t}\partial_{\dot{\mathbf{w}}^n}J^{(n+1)T}
 * \f]
 *
 * with response function
 *\f$J^n=J(\mathbf{w}^n,\dot{\mathbf{w}}^n,\dot{\mathbf{w}}^{n-1};\mathbf{s})\f$.
 */
template <class TSparseSpace, class TDenseSpace>
class AdjointBossakScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointBossakScheme);

    typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::TSystemMatrixType SystemMatrixType;

    typedef typename BaseType::TSystemVectorType SystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    AdjointBossakScheme(Parameters& rParameters, ResponseFunction::Pointer pResponseFunction)
        : Scheme<TSparseSpace, TDenseSpace>()
    {
        KRATOS_TRY;

        Parameters default_params(R"(
        {
            "scheme_type": "bossak",
            "alpha_bossak": -0.3
        })");

        rParameters.ValidateAndAssignDefaults(default_params);

        mAlphaBossak = rParameters["alpha_bossak"].GetDouble();
        mGammaNewmark = 0.5 - mAlphaBossak;
        mInvGamma = 1.0 / mGammaNewmark;
        mInvGammaMinusOne = 1.0 / (mGammaNewmark - 1.0);

        mpResponseFunction = pResponseFunction;

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~AdjointBossakScheme() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        BaseType::Initialize(rModelPart);

        // Allocate auxiliary memory.
        int num_threads = OpenMPUtils::GetNumThreads();
        mAdjointValuesVector.resize(num_threads);
        mAdjointFirstDerivsVector.resize(num_threads);
        mAdjointSecondDerivsVector.resize(num_threads);
        mAdjointAuxVector.resize(num_threads);
        mResponseGradient.resize(num_threads);
        mFirstDerivsResponseGradient.resize(num_threads);
        mSecondDerivsResponseGradient.resize(num_threads);
        mLeftHandSide.resize(num_threads);
        mFirstDerivsLHS.resize(num_threads);
        mSecondDerivsLHS.resize(num_threads);

        // Initialize the adjoint variables to zero (adjoint initial conditions
        // are always zero).
#pragma omp parallel
        {
            ModelPart::NodeIterator nodes_begin;
            ModelPart::NodeIterator nodes_end;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), nodes_begin, nodes_end);
            for (auto it = nodes_begin; it != nodes_end; ++it)
            {
                noalias(it->FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1)) = ADJOINT_FLUID_VECTOR_1.Zero();
                it->FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1) = ADJOINT_FLUID_SCALAR_1.Zero();
                noalias(it->FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3)) = ADJOINT_FLUID_VECTOR_3.Zero();
                noalias(it->FastGetSolutionStepValue(AUX_ADJOINT_FLUID_VECTOR_1)) = AUX_ADJOINT_FLUID_VECTOR_1.Zero();
            }
        }

        mpResponseFunction->Initialize(rModelPart);

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep(ModelPart& rModelPart,
                                SystemMatrixType& rA,
                                SystemVectorType& rDx,
                                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        // Get current time step.
        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        double delta_time = -r_current_process_info[DELTA_TIME]; // DELTA_TIME < 0
        if (delta_time <= 0.0)
            KRATOS_ERROR << "detected for adjoint solution DELTA_TIME >= 0" << std::endl;
        mInvDt = 1.0 / delta_time;

        // Calculate number of neighbour elements for each node.
#pragma omp parallel
        {
            ModelPart::NodeIterator nodes_begin;
            ModelPart::NodeIterator nodes_end;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), nodes_begin, nodes_end);
            for (auto it = nodes_begin; it != nodes_end; ++it)
                it->GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS) = 0.0;
        }

#pragma omp parallel
        {
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(rModelPart.Elements(), elements_begin, elements_end);
            for (auto it = elements_begin; it != elements_end; ++it)
                for (unsigned int i_node = 0; i_node < it->GetGeometry().PointsNumber(); ++i_node)
                {
                    double& r_num_neighbour = it->GetGeometry()[i_node].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
                    #pragma omp atomic
                    r_num_neighbour += 1.0;
                }
        }

        rModelPart.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_ELEMENTS);

        mpResponseFunction->InitializeSolutionStep(rModelPart);

        KRATOS_CATCH("");
    }

    void FinalizeSolutionStep(ModelPart& rModelPart,
                              SystemMatrixType& rA,
                              SystemVectorType& rDx,
                              SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        // Set aux adjoint acceleration to zero before we assemble it.
#pragma omp parallel
        {
            ModelPart::NodeIterator nodes_begin;
            ModelPart::NodeIterator nodes_end;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), nodes_begin, nodes_end);
            for (auto it = nodes_begin; it != nodes_end; ++it)
                noalias(it->FastGetSolutionStepValue(AUX_ADJOINT_FLUID_VECTOR_1)) =
                    AUX_ADJOINT_FLUID_VECTOR_1.Zero();
        }

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const unsigned int domain_size =
            static_cast<unsigned int>(r_current_process_info[DOMAIN_SIZE]);

        // Calculate and store contributions to the adjoint acceleration for
        // the adjoint solution of the next time step.
        // Loop over elements.
#pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(rModelPart.Elements(), elements_begin, elements_end);

            for (auto it = elements_begin; it != elements_end; ++it)
            {
                // Calculate transposed gradient of element residual w.r.t. acceleration.
                it->CalculateSecondDerivativesLHS(mSecondDerivsLHS[k], r_current_process_info);
                mSecondDerivsLHS[k] = mAlphaBossak * mSecondDerivsLHS[k];

                // Calculate transposed gradient of response function on element w.r.t. acceleration.
                mpResponseFunction->CalculateSecondDerivativesGradient(
                    *it, mSecondDerivsLHS[k], mSecondDerivsResponseGradient[k], r_current_process_info);

                // Get adjoint vector.
                it->GetValuesVector(mAdjointValuesVector[k]);

                mAdjointAuxVector[k] =
                    prod(mSecondDerivsLHS[k], mAdjointValuesVector[k]) + mSecondDerivsResponseGradient[k];

                // Assemble.
                unsigned int local_index = 0;
                for (unsigned int i_node = 0; i_node < it->GetGeometry().PointsNumber(); ++i_node)
                {
                    array_1d<double, 3>& r_adjoint_aux_fluid_vector_1 =
                        it->GetGeometry()[i_node].FastGetSolutionStepValue(AUX_ADJOINT_FLUID_VECTOR_1);
                    it->GetGeometry()[i_node].SetLock();
                    for (unsigned int d = 0; d < domain_size; ++d)
                        r_adjoint_aux_fluid_vector_1[d] += mAdjointAuxVector[k][local_index++];
                    it->GetGeometry()[i_node].UnSetLock();
                    ++local_index; // pressure dof
                }
            }
        }
        // Loop over conditions.
// #pragma omp parallel
//         {
//         }

        rModelPart.GetCommunicator().AssembleCurrentData(AUX_ADJOINT_FLUID_VECTOR_1);

        mpResponseFunction->FinalizeSolutionStep(rModelPart);

        KRATOS_CATCH("");
    }

    /// Update the adjoint solution.
    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                SystemMatrixType& rA,
                SystemVectorType& rDx,
                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const unsigned int domain_size =
            static_cast<unsigned int>(r_current_process_info[DOMAIN_SIZE]);
        Communicator& r_comm = rModelPart.GetCommunicator();

        if (r_comm.TotalProcesses() == 1)
        {
            // Update adjoint values.
            int ndofs = static_cast<int>(rDofSet.size());
            #pragma omp parallel for
            for (int i = 0; i < ndofs; ++i)
            {
                typename DofsArrayType::iterator it = rDofSet.begin() + i;
                if (it->IsFree() == true)
                    it->GetSolutionStepValue() +=
                        TSparseSpace::GetValue(rDx, it->EquationId());
            }

            // Assign contributions to adjoint second derivatives that don't 
            // require assembly.
            #pragma omp parallel
            {
                ModelPart::NodeIterator nodes_begin;
                ModelPart::NodeIterator nodes_end;
                OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), nodes_begin, nodes_end);
                for (auto it = nodes_begin; it != nodes_end; ++it)
                {
                    array_1d<double, 3>& r_adjoint_fluid_vector_3 =
                        it->FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3);
                    const array_1d<double, 3>& r_old_adjoint_fluid_vector_3 =
                        it->FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3, 1);
                    const array_1d<double, 3>& r_adjoint_aux_fluid_vector_1 =
                        it->FastGetSolutionStepValue(AUX_ADJOINT_FLUID_VECTOR_1, 1);
                    for (unsigned int d = 0; d < domain_size; ++d)
                        r_adjoint_fluid_vector_3[d] = (mGammaNewmark - 1.0) * mInvGamma *
                            (r_old_adjoint_fluid_vector_3[d] + r_adjoint_aux_fluid_vector_1[d]);
                }
            }
        }
        else
        {
            // Update adjoint values.
            int ndofs = static_cast<int>(rDofSet.size());
            #pragma omp parallel for
            for (int i = 0; i < ndofs; ++i)
            {
                typename DofsArrayType::iterator it = rDofSet.begin() + i;
                if (it->GetSolutionStepValue(PARTITION_INDEX) == r_comm.MyPID())
                    if (it->IsFree() == true)
                        it->GetSolutionStepValue() +=
                            TSparseSpace::GetValue(rDx, it->EquationId());
            }

            // todo: add a function Communicator::SynchronizeDofVariables() to
            // reduce communication here.
            r_comm.SynchronizeNodalSolutionStepsData();

            // Assign contributions to adjoint second derivatives that don't 
            // require assembly.
            #pragma omp parallel
            {
                ModelPart::NodeIterator nodes_begin;
                ModelPart::NodeIterator nodes_end;
                OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), nodes_begin, nodes_end);
                for (auto it = nodes_begin; it != nodes_end; ++it)
                {
                    array_1d<double, 3>& r_adjoint_fluid_vector_3 =
                        it->FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3);

                    // In the end we need to assemble so we only compute this part
                    // on the process that owns the node.
                    if (it->FastGetSolutionStepValue(PARTITION_INDEX) == r_comm.MyPID())
                    {
                        const array_1d<double, 3>& r_old_adjoint_fluid_vector_3 =
                            it->FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3, 1);
                        const array_1d<double, 3>& r_adjoint_aux_fluid_vector_1 =
                            it->FastGetSolutionStepValue(AUX_ADJOINT_FLUID_VECTOR_1, 1);
                        for (unsigned int d = 0; d < domain_size; ++d)
                            r_adjoint_fluid_vector_3[d] = (mGammaNewmark - 1.0) * mInvGamma *
                                (r_old_adjoint_fluid_vector_3[d] + r_adjoint_aux_fluid_vector_1[d]);
                    }
                    else
                    {
                        for (unsigned int d = 0; d < domain_size; ++d)
                            r_adjoint_fluid_vector_3[d] = 0.0;
                    }
                }
            }
        }

        // Add contributions to adjoint second derivatives that require assembly.
        // Loop over elements.
#pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            ModelPart::ElementIterator elements_begin;
            ModelPart::ElementIterator elements_end;
            OpenMPUtils::PartitionedIterators(rModelPart.Elements(), elements_begin, elements_end);
            for (auto it = elements_begin; it != elements_end; ++it)
            {
                // Calculate transposed gradient of element residual w.r.t. second derivatives.
                it->CalculateSecondDerivativesLHS(mSecondDerivsLHS[k], r_current_process_info);
                mSecondDerivsLHS[k] = (1.0 - mAlphaBossak) * mSecondDerivsLHS[k];

                // Calculate transposed gradient of response function on element w.r.t. acceleration.
                mpResponseFunction->CalculateSecondDerivativesGradient(
                    *it, mSecondDerivsLHS[k], mSecondDerivsResponseGradient[k], r_current_process_info);

                // Get adjoint vector.
                it->GetValuesVector(mAdjointValuesVector[k]);

                mAdjointSecondDerivsVector[k] = (mGammaNewmark - 1.0) * mInvGamma *
                    (prod(mSecondDerivsLHS[k], mAdjointValuesVector[k]) + mSecondDerivsResponseGradient[k]);
                
                // Assemble contributions to adjoint acceleration.
                unsigned int local_index = 0;
                for (unsigned int i_node = 0; i_node < it->GetGeometry().PointsNumber(); ++i_node)
                {
                    array_1d<double, 3>& r_adjoint_fluid_vector_3 =
                        it->GetGeometry()[i_node].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3);
                    it->GetGeometry()[i_node].SetLock();
                    for (unsigned int d = 0; d < domain_size; ++d)
                        r_adjoint_fluid_vector_3[d] += mAdjointSecondDerivsVector[k][local_index++];
                    it->GetGeometry()[i_node].UnSetLock();
                    ++local_index; // pressure dof
                }
            }
        }
        // Loop over conditions.
// #pragma omp parallel
//         {
//         }

        rModelPart.GetCommunicator().AssembleCurrentData(ADJOINT_FLUID_VECTOR_3);

        mpResponseFunction->UpdateSensitivities(rModelPart);

        KRATOS_CATCH("");
    }

    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        // Check domain dimension and element.
        const unsigned int working_space_dimension =
            rModelPart.Elements().begin()->WorkingSpaceDimension();

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const unsigned int domain_size =
            static_cast<unsigned int>(r_current_process_info[DOMAIN_SIZE]);
        if (domain_size != 2 && domain_size != 3)
            KRATOS_ERROR << "invalid DOMAIN_SIZE: " << domain_size << std::endl;
        if (domain_size != working_space_dimension)
            KRATOS_ERROR << "DOMAIN_SIZE != WorkingSpaceDimension()" << std::endl;

        if (rModelPart.NodesBegin()->SolutionStepsDataHas(ADJOINT_FLUID_VECTOR_1) == false)
            KRATOS_ERROR << "Nodal solution steps data missing variable: " << ADJOINT_FLUID_VECTOR_1 << std::endl;
        
        if (rModelPart.NodesBegin()->SolutionStepsDataHas(ADJOINT_FLUID_SCALAR_1) == false)
            KRATOS_ERROR << "Nodal solution steps data missing variable: " << ADJOINT_FLUID_SCALAR_1 << std::endl;
        
        if (rModelPart.NodesBegin()->SolutionStepsDataHas(ADJOINT_FLUID_VECTOR_3) == false)
            KRATOS_ERROR << "Nodal solution steps data missing variable: " << ADJOINT_FLUID_VECTOR_3 << std::endl;

        if (rModelPart.NodesBegin()->SolutionStepsDataHas(AUX_ADJOINT_FLUID_VECTOR_1) == false)
            KRATOS_ERROR << "Nodal solution steps data missing variable: " << AUX_ADJOINT_FLUID_VECTOR_1 << std::endl;

        return BaseType::Check(rModelPart); // Check elements and conditions.
        KRATOS_CATCH("");
    }

    /// Calculate residual based element contributions to transient adjoint.
    void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int thread_id = OpenMPUtils::ThisThread();

        // Calculate contribution from old adjoint acceleration.
        pCurrentElement->GetSecondDerivativesVector(rRHS_Contribution, 1);
        const unsigned int domain_size = static_cast<unsigned int>(rCurrentProcessInfo[DOMAIN_SIZE]);
        unsigned int local_index = 0;
        for (unsigned int i_node = 0; i_node < pCurrentElement->GetGeometry().PointsNumber(); ++i_node)
        {
            const array_1d<double, 3>& r_adjoint_aux_fluid_vector_1 =
                        pCurrentElement->GetGeometry()[i_node].FastGetSolutionStepValue(AUX_ADJOINT_FLUID_VECTOR_1, 1);
            double weight = 1.0 / pCurrentElement->GetGeometry()[i_node].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
            for (unsigned int d = 0; d < domain_size; ++d)
            {
                rRHS_Contribution[local_index] = mInvGamma * mInvDt * weight *
                    (mInvGammaMinusOne * rRHS_Contribution[local_index] - r_adjoint_aux_fluid_vector_1[d]);
                ++local_index;
            }
            ++local_index; // pressure dof
        }

        // Calculate transposed gradient of element residual w.r.t. acceleration.
        pCurrentElement->CalculateSecondDerivativesLHS(mSecondDerivsLHS[thread_id], rCurrentProcessInfo);
        mSecondDerivsLHS[thread_id] = (1.0 - mAlphaBossak) * mSecondDerivsLHS[thread_id];

        // Calculate transposed gradient of response function on element w.r.t. acceleration.
        mpResponseFunction->CalculateSecondDerivativesGradient(
            *pCurrentElement, mSecondDerivsLHS[thread_id], mSecondDerivsResponseGradient[thread_id], rCurrentProcessInfo);
        noalias(rRHS_Contribution) -= mInvGamma * mInvDt * mSecondDerivsResponseGradient[thread_id];

        // Calculate transposed gradient of element residual w.r.t. velocity.
        pCurrentElement->CalculateFirstDerivativesLHS(mFirstDerivsLHS[thread_id], rCurrentProcessInfo);

        // Calculate transposed gradient of response function on element w.r.t. velocity.
        mpResponseFunction->CalculateFirstDerivativesGradient(
            *pCurrentElement, mFirstDerivsLHS[thread_id], mFirstDerivsResponseGradient[thread_id], rCurrentProcessInfo);
        noalias(rRHS_Contribution) -= mFirstDerivsResponseGradient[thread_id];

        if (rLHS_Contribution.size1() != mFirstDerivsLHS[thread_id].size1() || rLHS_Contribution.size2() != mFirstDerivsLHS[thread_id].size2())
            rLHS_Contribution.resize(mFirstDerivsLHS[thread_id].size1(), mFirstDerivsLHS[thread_id].size2());
        noalias(rLHS_Contribution) = mFirstDerivsLHS[thread_id] + mInvGamma * mInvDt * mSecondDerivsLHS[thread_id];

        // Calculate system contributions in residual form.
        pCurrentElement->GetValuesVector(mAdjointValuesVector[thread_id]);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, mAdjointValuesVector[thread_id]);

        pCurrentElement->EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Calculate_LHS_Contribution(Element::Pointer pCurrentElement,
                                    LocalSystemMatrixType& rLHS_Contribution,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        LocalSystemVectorType RHS_Contribution;

        RHS_Contribution.resize(rLHS_Contribution.size1(), false);

        CalculateSystemContributions(
            pCurrentElement, rLHS_Contribution, RHS_Contribution, rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    // /// Calculate residual based condition contributions to transient adjoint.
    // void Condition_CalculateSystemContributions(
    //     Condition::Pointer pCurrentCondition,
    //     LocalSystemMatrixType& rLHS_Contribution,
    //     LocalSystemVectorType& rRHS_Contribution,
    //     Condition::EquationIdVectorType& rEquationId,
    //     ProcessInfo& rCurrentProcessInfo) override
    // {
    //     KRATOS_TRY;



    //     KRATOS_CATCH("");
    // }

    // void Condition_Calculate_LHS_Contribution(Condition::Pointer pCurrentCondition,
    //                                           LocalSystemMatrixType& rLHS_Contribution,
    //                                           Condition::EquationIdVectorType& rEquationId,
    //                                           ProcessInfo& rCurrentProcessInfo) override
    // {
    //     KRATOS_TRY;

    //     LocalSystemVectorType RHS_Contribution;

    //     RHS_Contribution.resize(rLHS_Contribution.size1(), false);

    //     this->Condition_CalculateSystemContributions(
    //         pCurrentCondition, rLHS_Contribution, RHS_Contribution, rEquationId, rCurrentProcessInfo);

    //     KRATOS_CATCH("");
    // }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
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

    double mAlphaBossak;
    double mGammaNewmark;
    double mInvDt;
    double mInvGamma;
    double mInvGammaMinusOne;
    ResponseFunction::Pointer mpResponseFunction;
    std::vector<LocalSystemVectorType> mAdjointValuesVector;
    std::vector<LocalSystemVectorType> mAdjointFirstDerivsVector;
    std::vector<LocalSystemVectorType> mAdjointSecondDerivsVector;
    std::vector<LocalSystemVectorType> mAdjointAuxVector;
    std::vector<LocalSystemVectorType> mResponseGradient;
    std::vector<LocalSystemVectorType> mFirstDerivsResponseGradient;
    std::vector<LocalSystemVectorType> mSecondDerivsResponseGradient;
    std::vector<LocalSystemMatrixType> mLeftHandSide;
    std::vector<LocalSystemMatrixType> mFirstDerivsLHS;
    std::vector<LocalSystemMatrixType> mSecondDerivsLHS;

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

    ///@}

}; /* Class Scheme */

///@}

///@name Type Definitions
///@{

///@}

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_ADJOINT_BOSSAK_SCHEME defined */

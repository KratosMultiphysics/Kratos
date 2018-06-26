//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//                   Jordi Cotela, https://github.com/jcotela
//

#if !defined(KRATOS_RESIDUAL_BASED_ADJOINT_BOSSAK_SCHEME_H_INCLUDED)
#define KRATOS_RESIDUAL_BASED_ADJOINT_BOSSAK_SCHEME_H_INCLUDED

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// A scheme for dynamic adjoint equations, using Bossak time integration.
/**
 */
template <class TSparseSpace, class TDenseSpace>
class ResidualBasedAdjointBossakScheme : public ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedAdjointBossakScheme);

    typedef ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::SystemVectorType SystemVectorType;
    typedef typename BaseType::SystemMatrixType SystemMatrixType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ResidualBasedAdjointBossakScheme(Parameters& rParameters, ResponseFunction::Pointer pResponseFunction):
        ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>(pResponseFunction),
        mVelocityUpdateAdjointVariable(VELOCITY),
        mAccelerationUpdateAdjointVariable(VELOCITY),
        mAuxiliaryVariable(VELOCITY)
    {

        Parameters default_parameters(R"({
            "scheme_type": "bossak",
            "alpha_bossak": -0.3,
            "velocity_update_adjoint_variable": "ADJOINT_FLUID_VECTOR_2",
            "acceleration_update_adjoint_variable": "ADJOINT_FLUID_VECTOR_3",
            "auxiliary_variable": "AUX_ADJOINT_FLUID_VECTOR_1"
        })");

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mAlphaBossak = rParameters["alpha_bossak"].GetDouble();
        mBetaNewmark = 0.25 * (1.0 - mAlphaBossak) * (1.0 - mAlphaBossak);
        mGammaNewmark = 0.5 - mAlphaBossak;

        mVelocityUpdateAdjointVariable = GetVariableFromParameters(rParameters, "velocity_update_adjoint_variable");
        mAccelerationUpdateAdjointVariable = GetVariableFromParameters(rParameters, "acceleration_update_adjoint_variable");
        mAuxiliaryVariable = GetVariableFromParameters(rParameters, "auxiliary_variable");
    }

    /// Destructor.
    ~ResidualBasedAdjointBossakScheme() override
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
        mLeftHandSide.resize(num_threads);
        mResponseGradient.resize(num_threads);
        mFirstDerivsLHS.resize(num_threads);
        mFirstDerivsResponseGradient.resize(num_threads);
        mSecondDerivsLHS.resize(num_threads);
        mSecondDerivsResponseGradient.resize(num_threads);
        mAdjointValuesVector.resize(num_threads);
        mAdjointFirstDerivsVector.resize(num_threads);
        mAdjointSecondDerivsVector.resize(num_threads);
        mAdjointAuxiliaryVector.resize(num_threads);

        this->mpResponseFunction->Initialize(rModelPart);

        InitializeNodeNeighbourCount(rModelPart.Nodes());

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep(
                    ModelPart& rModelPart,
                    SystemMatrixType& rA,
                    SystemVectorType& rDx,
                    SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        // Get current time step.
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const ProcessInfo& r_last_process_info = r_current_process_info.GetPreviousSolutionStepInfo(1);

        // Note: solution is backwards in time, but we still want a positive time step
        // (it is the time step in the "forward" Bossak scheme).
        mTimeStep = r_last_process_info.GetValue(TIME) - r_current_process_info.GetValue(TIME);
        KRATOS_ERROR_IF(mTimeStep <= 0.0) << "Backwards in time solution is not decreasing time from last step." << std::endl;

        mInverseDt = 1.0 / mTimeStep;

        CalculateNodeNeighbourCount(rModelPart);

        this->mpResponseFunction->InitializeSolutionStep(rModelPart);

        KRATOS_CATCH("");
    }

    void FinalizeSolutionStep(
                        ModelPart& rModelPart,
                        SystemMatrixType& rA,
                        SystemVectorType& rDx,
                        SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        this->UpdateAuxiliaryVariable(rModelPart);

        this->mpResponseFunction->FinalizeSolutionStep(rModelPart);

        KRATOS_CATCH("");
    }

    void Update(ModelPart& rModelPart,
                        DofsArrayType& rDofSet,
                        SystemMatrixType& rA,
                        SystemVectorType& rDx,
                        SystemVectorType& rb) override
    {
        KRATOS_TRY;

        // Update degrees of freedom: adjoint variables associated to the residual of the physical problem.
        this->mpDofUpdater->UpdateDofs(rDofSet,rDx);

        // Update adjoint variables associated to time integration.
        this->UpdateTimeSchemeAdjoints(rModelPart);

        this->mpResponseFunction->UpdateSensitivities(rModelPart);

        KRATOS_CATCH("");
    }


    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(rModelPart.GetCommunicator().TotalProcesses() != 1) << "MPI version not implemented yet." << std::endl;

        KRATOS_ERROR_IF(rModelPart.NumberOfElements() == 0) << "No elements found in the ModelPart." << std::endl;

        // Check that the element and space dimensions match
        const unsigned int working_space_dimension = rModelPart.ElementsBegin()->WorkingSpaceDimension();
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        const unsigned int domain_size = static_cast<unsigned int>(r_process_info.GetValue(DOMAIN_SIZE));

        KRATOS_ERROR_IF(domain_size != 2 && domain_size != 3) <<
            "invalid DOMAIN_SIZE: " << domain_size << "." << std::endl;
        KRATOS_ERROR_IF(domain_size != working_space_dimension) <<
            "DOMAIN_SIZE " << domain_size << " not equal to the element's Working Space Dimension " <<
            working_space_dimension << "." << std::endl;

        // Check that the required variables are defined
        const Node<3>& r_node = *(rModelPart.NodesBegin());
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mVelocityUpdateAdjointVariable,r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mAccelerationUpdateAdjointVariable,r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mAuxiliaryVariable,r_node);

        return BaseType::Check(rModelPart); // Check elements and conditions.

        KRATOS_CATCH("");
    }

    void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        // Check and resize rLHS and rRHS
        this->CheckAndResizeLocalSystem(pCurrentElement, rLHS_Contribution, rRHS_Contribution);

        // Contribution from variable gradients
        this->CalculateGradientContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        // Contribution from first derivative gradients
        this->CalculateFirstDerivativeContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        // Contribution from second derivative gradients
        this->CalculateSecondDerivativeContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        // Contributions from the previos time step
        this->CalculatePreviousTimeStepContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        // Make the local contribution residual
        this->CalculateResidualLocalContributions(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

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

        CalculateSystemContributions(
            pCurrentElement, rLHS_Contribution, RHS_Contribution, rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Condition_CalculateSystemContributions(Condition::Pointer pCurrentCondition,
                                                        LocalSystemMatrixType& rLHS_Contribution,
                                                        LocalSystemVectorType& rRHS_Contribution,
                                                        Condition::EquationIdVectorType& rEquationId,
                                                        ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    void Condition_Calculate_LHS_Contribution(Condition::Pointer pCurrentCondition,
                                                      LocalSystemMatrixType& rLHS_Contribution,
                                                      Condition::EquationIdVectorType& rEquationId,
                                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

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
    void CheckAndResizeLocalSystem(Element::Pointer pCurrentElement,
                                   LocalSystemMatrixType& rLHS_Contribution,
                                   LocalSystemVectorType& rRHS_Contribution)
    {
        int k = OpenMPUtils::ThisThread();
        auto& r_residual_adjoint = mAdjointValuesVector[k];
        pCurrentElement->GetValuesVector(r_residual_adjoint);

        if (rRHS_Contribution.size() != r_residual_adjoint.size()) {
            rRHS_Contribution.resize(r_residual_adjoint.size(),false);
        }

        if (rLHS_Contribution.size1() != r_residual_adjoint.size()) {
            rLHS_Contribution.resize(r_residual_adjoint.size(),r_residual_adjoint.size(),false);
        }

        this->CheckAndResizeThreadStorage(r_residual_adjoint.size());
    }

    virtual void CalculateGradientContributions(Element::Pointer pCurrentElement,
                                                LocalSystemMatrixType& rLHS_Contribution,
                                                LocalSystemVectorType& rRHS_Contribution,
                                                ProcessInfo& rCurrentProcessInfo)
    {
        int k = OpenMPUtils::ThisThread();
        auto& r_response_function = *(this->mpResponseFunction);

        auto& r_lhs = mLeftHandSide[k];
        auto& r_response_gradient = mResponseGradient[k];
        pCurrentElement->CalculateLeftHandSide(r_lhs,rCurrentProcessInfo);
        r_response_function.CalculateGradient(*pCurrentElement, r_lhs, r_response_gradient, rCurrentProcessInfo);

        noalias(rLHS_Contribution) = r_lhs;
        noalias(rRHS_Contribution) = -1. * r_response_gradient;
    }

    virtual void CalculateFirstDerivativeContributions(Element::Pointer pCurrentElement,
                                                       LocalSystemMatrixType& rLHS_Contribution,
                                                       LocalSystemVectorType& rRHS_Contribution,
                                                       ProcessInfo& rCurrentProcessInfo)
    {
        int k = OpenMPUtils::ThisThread();
        auto& r_response_function = *(this->mpResponseFunction);

        const double first_deriv_coeff = mGammaNewmark / (mBetaNewmark*mTimeStep);
        auto& r_first_lhs = mFirstDerivsLHS[k];
        auto& r_first_response_gradient = mFirstDerivsResponseGradient[k];
        pCurrentElement->CalculateFirstDerivativesLHS(r_first_lhs,rCurrentProcessInfo);
        r_response_function.CalculateFirstDerivativesGradient(*pCurrentElement,r_first_lhs, r_first_response_gradient, rCurrentProcessInfo);

        noalias(rLHS_Contribution) += first_deriv_coeff * r_first_lhs;
        noalias(rRHS_Contribution) -= first_deriv_coeff * r_first_response_gradient;
    }

    virtual void CalculateSecondDerivativeContributions(Element::Pointer pCurrentElement,
                                                        LocalSystemMatrixType& rLHS_Contribution,
                                                        LocalSystemVectorType& rRHS_Contribution,
                                                        ProcessInfo& rCurrentProcessInfo)
    {
        int k = OpenMPUtils::ThisThread();
        auto& r_response_function = *(this->mpResponseFunction);

        const double second_deriv_coeff = mInverseDt * mInverseDt / mBetaNewmark;
        auto& r_second_lhs = mSecondDerivsLHS[k];
        auto& r_second_response_gradient = mSecondDerivsResponseGradient[k];
        pCurrentElement->CalculateSecondDerivativesLHS(r_second_lhs,rCurrentProcessInfo);
        r_second_lhs *= (1.0 - mAlphaBossak);
        r_response_function.CalculateSecondDerivativesGradient(*pCurrentElement,r_second_lhs, r_second_response_gradient, rCurrentProcessInfo);

        noalias(rLHS_Contribution) += second_deriv_coeff * r_second_lhs;
        noalias(rRHS_Contribution) -= second_deriv_coeff * r_second_response_gradient;
    }

    virtual void CalculatePreviousTimeStepContributions(Element::Pointer pCurrentElement,
                                                        LocalSystemMatrixType& rLHS_Contribution,
                                                        LocalSystemVectorType& rRHS_Contribution,
                                                        ProcessInfo& rCurrentProcessInfo)
    {
        const double old_adjoint_velocity_coeff = mInverseDt * (mBetaNewmark - mGammaNewmark * (mGammaNewmark + 0.5)) / (mBetaNewmark*mBetaNewmark);
        const double old_adjoint_acceleration_coeff = -1.0 * (mInverseDt*mInverseDt) * (mGammaNewmark + 0.5) / (mBetaNewmark*mBetaNewmark);
        const double second_deriv_coeff = mInverseDt * mInverseDt / mBetaNewmark;

        const Geometry< Node<3> >& r_geometry = pCurrentElement->GetGeometry();
        const unsigned int domain_size = r_geometry.WorkingSpaceDimension();
        const unsigned int num_nodes = r_geometry.PointsNumber();

        unsigned int local_index = 0;
        for (unsigned int i_node = 0; i_node < num_nodes; ++i_node) {
            auto& r_node = r_geometry[i_node];
            const array_1d<double, 3>& r_aux_adjoint_vector = r_node.FastGetSolutionStepValue(mAuxiliaryVariable,1);
            const array_1d<double, 3>& r_velocity_adjoint = r_node.FastGetSolutionStepValue(mVelocityUpdateAdjointVariable);
            const array_1d<double, 3>& r_acceleration_adjoint = r_node.FastGetSolutionStepValue(mAccelerationUpdateAdjointVariable);
            const double weight = 1.0 / r_node.GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);

            for (unsigned int d = 0; d < domain_size; d++) {
                rRHS_Contribution[local_index] += weight * second_deriv_coeff * r_aux_adjoint_vector[d];
                rRHS_Contribution[local_index] += weight * old_adjoint_velocity_coeff * r_velocity_adjoint[d];
                rRHS_Contribution[local_index] += weight * old_adjoint_acceleration_coeff * r_acceleration_adjoint[d];
                ++local_index;
            }
            ++local_index; // skip continuity equation adjoint rows.
        }
    }

    virtual void CalculateResidualLocalContributions(Element::Pointer pCurrentElement,
                                                     LocalSystemMatrixType& rLHS_Contribution,
                                                     LocalSystemVectorType& rRHS_Contribution,
                                                     ProcessInfo& rCurrentProcessInfo)
    {
        int k = OpenMPUtils::ThisThread();
        auto& r_residual_adjoint = mAdjointValuesVector[k];
        pCurrentElement->GetValuesVector(r_residual_adjoint);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, r_residual_adjoint);
    }

    virtual void InitializeNodeNeighbourCount(ModelPart::NodesContainerType& rNodes)
    {
        // This loop should not be omp parallel
        // The operation is not threadsafe if the value is uninitialized
        for (auto& r_node : rNodes)
            r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS,0.0);
    }

    virtual void CalculateNodeNeighbourCount(ModelPart& rModelPart)
    {
        // Calculate number of neighbour elements for each node.
        const int num_nodes = rModelPart.NumberOfNodes();
        #pragma omp parallel for
        for (int i = 0; i < num_nodes; i++) {
            Node<3>& r_node = *(rModelPart.Nodes().begin()+i);
            r_node.SetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS,0.0);
        }

        const int num_elements = rModelPart.NumberOfElements();
        #pragma omp parallel for
        for (int i = 0; i < num_elements; i++) {
            Element& r_element = *(rModelPart.Elements().begin()+i);
            Geometry<Node<3>>& r_geometry = r_element.GetGeometry();
            for (unsigned int j = 0; j < r_geometry.PointsNumber(); j++) {
                double& r_num_neighbour = r_geometry[j].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
                #pragma omp atomic
                r_num_neighbour += 1.0;
            }
        }

        rModelPart.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_ELEMENTS);
    }

    virtual void UpdateTimeSchemeAdjoints(ModelPart& rModelPart)
    {
        Communicator& r_communicator = rModelPart.GetCommunicator();
        ResponseFunction& r_response_function = *(this->mpResponseFunction);

        const double a22 = 1.0 - mGammaNewmark/mBetaNewmark;
        const double a23 = -1.0 / (mBetaNewmark*mTimeStep);
        const double a32 = (1.0 - 0.5*mGammaNewmark/mBetaNewmark)*mTimeStep;
        const double a33 = (1.0 - 0.5/mBetaNewmark);

        // Process the part that does not require assembly first
        const int number_of_local_nodes = r_communicator.LocalMesh().NumberOfNodes();
        #pragma omp parallel for
        for (int i = 0; i < number_of_local_nodes; i++) {
            Node<3>& r_node = *(r_communicator.LocalMesh().NodesBegin() + i);
            array_1d<double,3>& r_lambda2 = r_node.FastGetSolutionStepValue(mVelocityUpdateAdjointVariable);
            array_1d<double,3>& r_lambda3 = r_node.FastGetSolutionStepValue(mAccelerationUpdateAdjointVariable);
            const array_1d<double,3>& r_lambda2_old = r_node.FastGetSolutionStepValue(mVelocityUpdateAdjointVariable,1);
            const array_1d<double,3>& r_lambda3_old = r_node.FastGetSolutionStepValue(mAccelerationUpdateAdjointVariable,1);
            const array_1d<double, 3>& r_old_aux_adjoint_vector = r_node.FastGetSolutionStepValue(mAuxiliaryVariable,1);

            noalias(r_lambda2) = a22 * r_lambda2_old + a23 * r_lambda3_old;
            noalias(r_lambda3) = a32 * r_lambda2_old + a33 * r_lambda3_old + r_old_aux_adjoint_vector;
        }

        const int number_of_ghost_nodes = r_communicator.GhostMesh().NumberOfNodes();
        #pragma omp parallel for
        for (int i = 0; i < number_of_ghost_nodes; i++) {
            Node<3>& r_node = *(r_communicator.GhostMesh().NodesBegin() + i);
            noalias(r_node.FastGetSolutionStepValue(mVelocityUpdateAdjointVariable)) = mVelocityUpdateAdjointVariable.Zero();
            noalias(r_node.FastGetSolutionStepValue(mAccelerationUpdateAdjointVariable)) = mAccelerationUpdateAdjointVariable.Zero();
        }

        // Loop over elements to assemble the remaining terms
        const int number_of_elements = rModelPart.NumberOfElements();
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        #pragma omp parallel for
        for (int i = 0; i < number_of_elements; i++) {
            Element& r_element = *(rModelPart.ElementsBegin()+i);
            const int k = OpenMPUtils::ThisThread();
            auto& r_lhs = mLeftHandSide[k];
            auto& r_response_gradient = mResponseGradient[k];
            auto& r_first_lhs = mFirstDerivsLHS[k];
            auto& r_first_response_gradient = mFirstDerivsResponseGradient[k];
            auto& r_second_lhs = mSecondDerivsLHS[k];
            auto& r_second_response_gradient = mSecondDerivsResponseGradient[k];
            auto& r_residual_adjoint = mAdjointValuesVector[k];
            auto& r_velocity_adjoint = mAdjointFirstDerivsVector[k];
            auto& r_acceleration_adjoint = mAdjointSecondDerivsVector[k];

            r_element.GetValuesVector(r_residual_adjoint);
            this->CheckAndResizeThreadStorage(r_residual_adjoint.size());

            r_element.CalculateLeftHandSide(r_lhs,r_process_info);
            r_response_function.CalculateGradient(r_element,r_lhs,r_response_gradient,r_process_info);

            r_element.CalculateFirstDerivativesLHS(r_first_lhs,r_process_info);
            r_response_function.CalculateFirstDerivativesGradient(r_element,r_first_lhs,r_first_response_gradient,r_process_info);

            r_element.CalculateSecondDerivativesLHS(r_second_lhs,r_process_info);
            r_second_lhs *= (1.0 - mAlphaBossak);
            r_response_function.CalculateSecondDerivativesGradient(r_element,r_second_lhs,r_second_response_gradient,r_process_info);

            noalias(r_velocity_adjoint) = - r_first_response_gradient - prod(r_first_lhs,r_residual_adjoint);
            noalias(r_acceleration_adjoint) = - r_second_response_gradient - prod(r_second_lhs,r_residual_adjoint);

            // Assemble the contributions to the corresponing nodal unkowns.
            unsigned int local_index = 0;
            Geometry< Node<3> >& r_geometry = r_element.GetGeometry();
            for (unsigned int i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {

                Node<3>& r_node = r_geometry[i_node];
                array_1d<double, 3>& r_adjoint_fluid_vector_2 = r_node.FastGetSolutionStepValue(mVelocityUpdateAdjointVariable);
                array_1d<double, 3>& r_adjoint_fluid_vector_3 = r_node.FastGetSolutionStepValue(mAccelerationUpdateAdjointVariable);

                r_node.SetLock();
                for (unsigned int d = 0; d < r_geometry.WorkingSpaceDimension(); ++d) {

                    r_adjoint_fluid_vector_2[d] += r_velocity_adjoint[local_index];
                    r_adjoint_fluid_vector_3[d] += r_acceleration_adjoint[local_index];
                    ++local_index;
                }
                r_node.UnSetLock();
                ++local_index; // pressure dof
            }
        }

        // Finalize calculation
        r_communicator.AssembleCurrentData(mVelocityUpdateAdjointVariable);
        r_communicator.AssembleCurrentData(mAccelerationUpdateAdjointVariable);
    }


    virtual void UpdateAuxiliaryVariable(ModelPart& rModelPart)
    {
        ResponseFunction& r_response_function = *(this->mpResponseFunction);

        // Process the part that does not require assembly first
        const int number_of_nodes = rModelPart.NumberOfNodes();
        #pragma omp parallel for
        for (int i = 0; i < number_of_nodes; i++) {
            Node<3>& r_node = *(rModelPart.NodesBegin() + i);
            noalias(r_node.FastGetSolutionStepValue(mAuxiliaryVariable)) = mAuxiliaryVariable.Zero();
        }

        // Loop over elements to assemble the remaining terms
        const int number_of_elements = rModelPart.NumberOfElements();
        ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        #pragma omp parallel for
        for (int i = 0; i < number_of_elements; i++) {
            Element& r_element = *(rModelPart.ElementsBegin()+i);
            const int k = OpenMPUtils::ThisThread();
            auto& r_second_lhs = mSecondDerivsLHS[k];
            auto& r_second_response_gradient = mSecondDerivsResponseGradient[k];
            auto& r_residual_adjoint = mAdjointValuesVector[k];
            auto& r_adjoint_auxiliary = mAdjointAuxiliaryVector[k];

            r_element.GetValuesVector(r_residual_adjoint);
            this->CheckAndResizeThreadStorage(r_residual_adjoint.size());

            r_element.CalculateSecondDerivativesLHS(r_second_lhs,r_process_info);
            r_second_lhs *= mAlphaBossak;
            r_response_function.CalculateSecondDerivativesGradient(r_element,r_second_lhs,r_second_response_gradient,r_process_info);

            noalias(r_adjoint_auxiliary) = - prod(r_second_lhs,r_residual_adjoint) - r_second_response_gradient;

            // Assemble the contributions to the corresponing nodal unkowns.
            unsigned int local_index = 0;
            Geometry< Node<3> >& r_geometry = r_element.GetGeometry();
            for (unsigned int i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {

                Node<3>& r_node = r_geometry[i_node];
                array_1d<double, 3>& r_aux_adjoint_fluid_vector_1 = r_node.FastGetSolutionStepValue(mAuxiliaryVariable);

                r_node.SetLock();
                for (unsigned int d = 0; d < r_geometry.WorkingSpaceDimension(); ++d) {

                    r_aux_adjoint_fluid_vector_1[d] += r_adjoint_auxiliary[local_index];
                    ++local_index;
                }
                r_node.UnSetLock();
                ++local_index; // pressure dof
            }
        }

        // Finalize calculation
        rModelPart.GetCommunicator().AssembleCurrentData(mAuxiliaryVariable);
    }

    /// Free memory allocated by this class.
    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

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
    double mBetaNewmark;
    double mGammaNewmark;
    double mTimeStep;
    double mInverseDt;

    std::vector< LocalSystemMatrixType > mLeftHandSide;
    std::vector< LocalSystemVectorType > mResponseGradient;
    std::vector< LocalSystemMatrixType > mFirstDerivsLHS;
    std::vector< LocalSystemVectorType > mFirstDerivsResponseGradient;
    std::vector< LocalSystemMatrixType > mSecondDerivsLHS;
    std::vector< LocalSystemVectorType > mSecondDerivsResponseGradient;
    std::vector< LocalSystemVectorType > mAdjointValuesVector;
    std::vector< LocalSystemVectorType > mAdjointFirstDerivsVector;
    std::vector< LocalSystemVectorType > mAdjointSecondDerivsVector;
    std::vector< LocalSystemVectorType > mAdjointAuxiliaryVector;

    Variable< array_1d<double,3> > mVelocityUpdateAdjointVariable;
    Variable< array_1d<double,3> > mAccelerationUpdateAdjointVariable;
    Variable< array_1d<double,3> > mAuxiliaryVariable;

    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CheckAndResizeThreadStorage(unsigned int SystemSize) {
        const int k = OpenMPUtils::ThisThread();

        if (mLeftHandSide[k].size1() != SystemSize || mLeftHandSide[k].size2() != SystemSize) {
            mLeftHandSide[k].resize(SystemSize, SystemSize, false);
        }

        if (mFirstDerivsLHS[k].size1() != SystemSize || mFirstDerivsLHS[k].size2() != SystemSize) {
            mFirstDerivsLHS[k].resize(SystemSize, SystemSize, false);
        }

        if (mSecondDerivsLHS[k].size1() != SystemSize || mSecondDerivsLHS[k].size2() != SystemSize) {
            mSecondDerivsLHS[k].resize(SystemSize, SystemSize, false);
        }

        if (mResponseGradient[k].size() != SystemSize) {
            mResponseGradient[k].resize(SystemSize,false);
        }

        if (mFirstDerivsResponseGradient[k].size() != SystemSize) {
            mFirstDerivsResponseGradient[k].resize(SystemSize,false);
        }

        if (mSecondDerivsResponseGradient[k].size() != SystemSize) {
            mSecondDerivsResponseGradient[k].resize(SystemSize,false);
        }

        if (mAdjointFirstDerivsVector[k].size() != SystemSize) {
            mAdjointFirstDerivsVector[k].resize(SystemSize,false);
        }

        if (mAdjointSecondDerivsVector[k].size() != SystemSize) {
            mAdjointSecondDerivsVector[k].resize(SystemSize,false);
        }

        if (mAdjointAuxiliaryVector[k].size() != SystemSize) {
            mAdjointAuxiliaryVector[k].resize(SystemSize,false);
        }
    }

    Variable< array_1d<double,3> > GetVariableFromParameters(
        Parameters& rParameters,
        const std::string& rLabel) {
        KRATOS_TRY;

        std::string variable_name = rParameters[rLabel].GetString();
        bool have_variable = KratosComponents< Variable< array_1d<double,3> > >::Has(variable_name);

        KRATOS_ERROR_IF_NOT( have_variable ) << "Variable " << variable_name <<
        " passed as Parameters argument \"" << rLabel << "\" is not an array_1d<double,3> Variable defined in Kratos." << std::endl <<
        "If it is an application variable, you may need to imoirt the application that defines it." << std::endl;

        return KratosComponents< Variable< array_1d<double,3> > >::Get(variable_name);

        KRATOS_CATCH("");
    }

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

}; /* Class ResidualBasedAdjointBossakScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ADJOINT_BOSSAK_SCHEME_H_INCLUDED defined */

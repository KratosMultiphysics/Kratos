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
        ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>(pResponseFunction)
    {

        Parameters default_parameters(R"({
            "scheme_type": "bossak",
            "alpha_bossak": -0.3
        })");

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mAlphaBossak = rParameters["alpha_bossak"].GetDouble();
        mBetaNewmark = 0.25 * (1.0 - mAlphaBossak) * (1.0 - mAlphaBossak);
        mGammaNewmark = 0.5 - mAlphaBossak;
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

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                SystemMatrixType& rA,
                SystemVectorType& rDx,
                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        // Update degrees of freedom: adjoint variables associated to the residual of the physical problem.
        this->UpdateDegreesOfFreedom(rDofSet,rDx);

        // Update additional variables: adjoint variables associated to time integration and auxiliary quantities.
        this->UpdateAdditionalVariables(rModelPart);

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
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_VECTOR_1,r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_SCALAR_1,r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_VECTOR_2,r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_VECTOR_3,r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(AUX_ADJOINT_FLUID_VECTOR_1,r_node);

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

        KRATOS_CATCH("");
    }

    void Calculate_LHS_Contribution(Element::Pointer pCurrentElement,
                                    LocalSystemMatrixType& rLHS_Contribution,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

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

    virtual void UpdateDegreesOfFreedom(
        DofsArrayType& rDofSet,
        SystemVectorType& rDx)
    {
        const int num_dofs = rDofSet.size();

        #pragma omp parallel for
        for (int i = 0; i < num_dofs; i++) {
            typename DofsArrayType::iterator i_dof = rDofSet.begin() + i;
            if (i_dof->IsFree()) {
                i_dof->GetSolutionStepValue() += TSparseSpace::GetValue(rDx, i_dof->EquationId());
            }
        }
    }

    virtual void UpdateAdditionalVariables(ModelPart& rModelPart)
    {
        const ProcessInfo& r_const_process_info = rModelPart.GetProcessInfo();
        Communicator& r_communicator = rModelPart.GetCommunicator();
        ResponseFunction& r_response_function = *(this->mpResponseFunction);

        const double delta_time = r_const_process_info.GetValue(DELTA_TIME);
        const double a22 = 1.0 - mGammaNewmark/mBetaNewmark;
        const double a23 = -1.0 / (mBetaNewmark*delta_time);
        const double a32 = (1.0 - 0.5*mGammaNewmark/mBetaNewmark)*delta_time;
        const double a33 = (1.0 - 0.5/mBetaNewmark);

        // Process the part that does not require assembly first
        const int number_of_local_nodes = r_communicator.LocalMesh().NumberOfNodes();
        #pragma omp parallel for
        for (int i = 0; i < number_of_local_nodes; i++) {
            Node<3>& r_node = *(r_communicator.LocalMesh().NodesBegin() + i);
            array_1d<double,3>& r_lambda2 = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2);
            array_1d<double,3>& r_lambda3 = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3);
            const array_1d<double,3>& r_lambda2_old = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2,1);
            const array_1d<double,3>& r_lambda3_old = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3,1);

            noalias(r_lambda2) = a22 * r_lambda2_old + a23 * r_lambda3_old;
            noalias(r_lambda3) = a32 * r_lambda2_old + a33 * r_lambda3_old;
        }

        const int number_of_ghost_nodes = r_communicator.GhostMesh().NumberOfNodes();
        #pragma omp parallel for
        for (int i = 0; i < number_of_ghost_nodes; i++) {
            Node<3>& r_node = *(r_communicator.GhostMesh().NodesBegin() + i);
            noalias(r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2)) = ADJOINT_FLUID_VECTOR_2.Zero();
            noalias(r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3)) = ADJOINT_FLUID_VECTOR_3.Zero();
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
            auto& r_adjoint_auxiliary = mAdjointAuxiliaryVector[k];

            r_element.CalculateLeftHandSide(r_lhs,r_process_info);
            r_response_function.CalculateGradient(r_element,r_lhs,r_response_gradient,r_process_info);

            r_element.CalculateFirstDerivativesLHS(r_first_lhs,r_process_info);
            r_response_function.CalculateFirstDerivativesGradient(r_element,r_first_lhs,r_first_response_gradient,r_process_info);

            r_element.CalculateSecondDerivativesLHS(r_second_lhs,r_process_info);
            r_second_lhs *= (1.0 - mAlphaBossak);
            r_response_function.CalculateSecondDerivativesGradient(r_element,r_second_lhs,r_second_response_gradient,r_process_info);

            r_element.GetValuesVector(r_residual_adjoint);

            noalias(r_velocity_adjoint) = - r_first_response_gradient - prod(r_first_lhs,r_residual_adjoint);
            noalias(r_acceleration_adjoint) = - r_second_response_gradient - prod(r_second_lhs,r_residual_adjoint);
            noalias(r_adjoint_auxiliary) = (mAlphaBossak/(1.0-mAlphaBossak)) * prod(r_second_lhs,r_residual_adjoint) - r_second_response_gradient;

            // Assemble the contributions to the corresponing nodal unkowns.
            unsigned int local_index = 0;
            Geometry< Node<3> >& r_geometry = r_element.GetGeometry();
            for (unsigned int i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {

                Node<3>& r_node = r_geometry[i_node];
                array_1d<double, 3>& r_adjoint_fluid_vector_2 = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_2);
                array_1d<double, 3>& r_adjoint_fluid_vector_3 = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3);
                array_1d<double, 3>& r_aux_adjoint_fluid_vector_1 = r_node.FastGetSolutionStepValue(AUX_ADJOINT_FLUID_VECTOR_1);
                const array_1d<double, 3>& r_old_aux_adjoint_fluid_vector_1 = r_node.FastGetSolutionStepValue(AUX_ADJOINT_FLUID_VECTOR_1,1);

                r_node.SetLock();
                for (unsigned int d = 0; d < r_geometry.WorkingSpaceDimension(); ++d) {

                    r_adjoint_fluid_vector_2[d] += r_velocity_adjoint[local_index];
                    r_adjoint_fluid_vector_3[d] += r_acceleration_adjoint[local_index] + r_old_aux_adjoint_fluid_vector_1[d];
                    r_aux_adjoint_fluid_vector_1[d] += r_adjoint_auxiliary[local_index];
                    ++local_index;
                }
                r_node.UnSetLock();
                ++local_index; // pressure dof
            }
        }

        // Finalize calculation
        rModelPart.GetCommunicator().AssembleCurrentData(ADJOINT_FLUID_VECTOR_2);
        rModelPart.GetCommunicator().AssembleCurrentData(ADJOINT_FLUID_VECTOR_3);
        rModelPart.GetCommunicator().AssembleCurrentData(AUX_ADJOINT_FLUID_VECTOR_1);

        r_response_function.UpdateSensitivities(rModelPart);
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

}; /* Class ResidualBasedAdjointBossakScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ADJOINT_BOSSAK_SCHEME_H_INCLUDED defined */

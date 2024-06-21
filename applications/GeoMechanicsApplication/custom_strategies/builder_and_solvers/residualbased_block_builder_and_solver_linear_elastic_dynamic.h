// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         BSD License
//                   geo_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi, Aron Noordam
//  Collaborators:   Vicente Mataix
//
//

#pragma once

/* System includes */
#include <unordered_set>

/* External includes */
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "includes/key_hash.h"
#include "includes/kratos_flags.h"
#include "includes/lock_object.h"
#include "includes/model_part.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "spaces/ublas_space.h"
#include "utilities/atomic_utilities.h"
#include "utilities/builtin_timer.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "utilities/timer.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class ResidualBasedBlockBuilderAndSolverWithMassAndDamping
 * @ingroup GeoMechanicsApplication
 * @brief Current class provides an implementation for builder and solving operations, while the
 * global mass and damping matrices are stored.
 * @details When the LHS is build, the global mass and damping matrices are build separately. When
 * building the RHS, the mass and damping matrices are multiplied with respectively the second and
 * first derivative vector to calculate the mass and damping contribution. The RHS is constituted by
 * the unbalanced loads (residual) Degrees of freedom are reordered putting the restrained degrees
 * of freedom at the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information. Calculation of the reactions involves a cost very similar to the calculation of
 * the total residual. This class is intended to be used when the mass and damping matrices are
 * constant throughout the iterations, using this class, when rebuilding the LHS every iteration is
 * not efficient.
 * @tparam TSparseSpace The sparse system considered
 * @tparam TDenseSpace The dense system considered
 * @tparam TLinearSolver The linear solver considered
 * @author Aron Noordam
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class ResidualBasedBlockBuilderAndSolverLinearElasticDynamic
    : public ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the flags
    KRATOS_DEFINE_LOCAL_FLAG(SILENT_WARNINGS);

    /// Definition of the pointer
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverLinearElasticDynamic);

    /// Definition of the base class
    using BaseType = ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// Definition of the classes from the base class
    using TSchemeType = typename BaseType::TSchemeType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;
    using NodesArrayType = typename BaseType::NodesArrayType;
    using ElementsArrayType = typename BaseType::ElementsArrayType;
    using ConditionsArrayType = typename BaseType::ConditionsArrayType;

    /// Additional definitions
    using ElementsContainerType = PointerVectorSet<Element, IndexedObject>;
    using EquationIdVectorType = Element::EquationIdVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    ResidualBasedBlockBuilderAndSolverLinearElasticDynamic() = default;

    /**
     * @brief Constructor. (with parameters)
     */
    explicit ResidualBasedBlockBuilderAndSolverLinearElasticDynamic(typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters)
        : BaseType(pNewLinearSystemSolver)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Constructor.
     */
    explicit ResidualBasedBlockBuilderAndSolverLinearElasticDynamic(typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    ~ResidualBasedBlockBuilderAndSolverLinearElasticDynamic() override = default;

    ///@}
    ///@name Operations
    ///@{


    void InitializeSolutionStep(ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb) override
    {

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        //mOutOfBalanceVector = ZeroVector(BaseType::mEquationSystemSize);
        mPreviousOutOfBalanceVector = ZeroVector(BaseType::mEquationSystemSize);
        mCurrentOutOfBalanceVector = ZeroVector(BaseType::mEquationSystemSize);
        //TSystemVectorType rDx_tot = ZeroVector(BaseType::mEquationSystemSize);

		if (mPreviousExternalForceVector.size() == 0)
		{
			mPreviousExternalForceVector = ZeroVector(BaseType::mEquationSystemSize);
		}
	}

    /**
     * @brief Function to perform the build of the RHS. The vector could be sized as the total
     * number of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rb The RHS vector
     */
    void BuildLHS(typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const auto nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const auto nconditions = static_cast<int>(rModelPart.Conditions().size());

        const ProcessInfo& r_current_process_info               = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator   el_begin   = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        // contributions to the system
        LocalSystemMatrixType lhs_contribution(0, 0);
        LocalSystemMatrixType mass_contribution(0, 0);
        LocalSystemMatrixType damping_contribution(0, 0);

        InitializeDynamicMatrix(mMassMatrix, BaseType::mEquationSystemSize, pScheme, rModelPart);
        InitializeDynamicMatrix(mDampingMatrix, BaseType::mEquationSystemSize, pScheme, rModelPart);

        // vector containing the localization in the system of the different
        // terms
        Element::EquationIdVectorType equation_ids;

        // Assemble all elements
        const auto timer = BuiltinTimer();

#pragma omp parallel firstprivate(nelements, nconditions, lhs_contribution, mass_contribution, \
                                      damping_contribution, equation_ids)
        {
#pragma omp for schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++) {
                auto it_elem = el_begin + k;

                if (it_elem->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateLHSContribution(*it_elem, lhs_contribution, equation_ids, r_current_process_info);

                    // assemble mass and damping matrix
                    it_elem->CalculateMassMatrix(mass_contribution, r_current_process_info);
                    it_elem->CalculateDampingMatrix(damping_contribution, r_current_process_info);

                    if (mass_contribution.size1() != 0) {
                        BaseType::AssembleLHS(mMassMatrix, mass_contribution, equation_ids);
                    }
                    if (damping_contribution.size1() != 0) {
                        BaseType::AssembleLHS(mDampingMatrix, damping_contribution, equation_ids);
                    }

                    // Assemble the elemental contribution
                    BaseType::AssembleLHS(rA, lhs_contribution, equation_ids);
                }
            }

#pragma omp for schedule(guided, 512)
            for (int k = 0; k < nconditions; k++) {
                auto it_cond = cond_begin + k;

                if (it_cond->IsActive()) {
                    // Calculate elemental contribution
                    pScheme->CalculateLHSContribution(*it_cond, lhs_contribution, equation_ids, r_current_process_info);

                    // assemble mass and damping matrix
                    it_cond->CalculateMassMatrix(mass_contribution, r_current_process_info);
                    it_cond->CalculateDampingMatrix(damping_contribution, r_current_process_info);

                    if (mass_contribution.size1() != 0) {
                        BaseType::AssembleLHS(mMassMatrix, mass_contribution, equation_ids);
                    }
                    if (damping_contribution.size1() != 0) {
                        BaseType::AssembleLHS(mDampingMatrix, damping_contribution, equation_ids);
                    }

                    // Assemble the elemental contribution
                    BaseType::AssembleLHS(rA, lhs_contribution, equation_ids);
                }
            }
        }

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", this->GetEchoLevel() >= 1)
            << "Build time: " << timer.ElapsedSeconds() << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic",
                       (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0))
            << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building and solving phase at the same time.
     * @details It is ideally the fastest and safer function to use when it is possible to solve
     * just after building
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void BuildAndSolve(typename TSchemeType::Pointer pScheme,
                       ModelPart&                    rModelPart,
                       TSystemMatrixType&            rA,
                       TSystemVectorType&            rDx,
                       TSystemVectorType&            rb) override
    {
      

        Timer::Start("Build");

        BuildLHS(pScheme, rModelPart, rA);

        BuildRHSNoDirichlet(pScheme, rModelPart, rb);

        Timer::Stop("Build");

        TSystemVectorType dummy_b(rA.size1(), 0.0);
        TSystemVectorType dummy_rDx(rA.size1(), 0.0);

        BaseType::ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);
        BaseType::ApplyDirichletConditions(pScheme, rModelPart, mMassMatrix, dummy_rDx, dummy_b);
        BaseType::ApplyDirichletConditions(pScheme, rModelPart, mDampingMatrix, dummy_rDx, dummy_b);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", this->GetEchoLevel() == 3)
            << "Before the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;

        const auto timer = BuiltinTimer();
        Timer::Start("Solve");

        BaseType::SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        std::copy(mCurrentOutOfBalanceVector.begin(), mCurrentOutOfBalanceVector.end(), mPreviousOutOfBalanceVector.begin());


        Timer::Stop("Solve");
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", this->GetEchoLevel() >= 1)
            << "System solve time: " << timer.ElapsedSeconds() << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", this->GetEchoLevel() == 3)
            << "After the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;


    }

    /**
     * @brief Corresponds to the previews, but the System's matrix is considered already built and
     * only the RHS is built again
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param rA The LHS matrix
     * @param rDx The Unknowns vector
     * @param rb The RHS vector
     */
    void BuildRHSAndSolve(typename TSchemeType::Pointer pScheme,
                          ModelPart&                    rModelPart,
                          TSystemMatrixType&            rA,
                          TSystemVectorType&            rDx,
                          TSystemVectorType&            rb) override
    {
        KRATOS_TRY

        BuildRHS(pScheme, rModelPart, rb);

        // todo check if this is needed
        //BaseType::ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", this->GetEchoLevel() == 3)
            << "Before the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;

        const auto timer = BuiltinTimer();
        Timer::Start("Solve");

        BaseType::SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        std::copy(mCurrentOutOfBalanceVector.begin(), mCurrentOutOfBalanceVector.end(), mPreviousOutOfBalanceVector.begin());

        Timer::Stop("Solve");
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", this->GetEchoLevel() >= 1)
            << "System solve time: " << timer.ElapsedSeconds() << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverLinearElasticDynamic", this->GetEchoLevel() == 3)
            << "After the solution of the system"
            << "\nSystem Matrix = " << rA << "\nUnknowns vector = " << rDx
            << "\nRHS vector = " << rb << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the build of the RHS.
     * @details The vector could be sized as the total number of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void BuildRHS(typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemVectorType& rb) override
    {
        KRATOS_TRY

        Timer::Start("BuildRHS");

        BuildRHSNoDirichlet(pScheme, rModelPart, rb);

        // NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        block_for_each(BaseType::mDofSet, [&](Dof<double>& r_dof) {
            if (r_dof.IsFixed()) {
                const std::size_t i = r_dof.EquationId();
                rb[i]               = 0.0;
            }
        });

        Timer::Stop("BuildRHS");

        KRATOS_CATCH("")
    }


    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
    )
    {

		BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        std::copy(mCurrentExternalForceVector.begin(), mCurrentExternalForceVector.end(), mPreviousExternalForceVector.begin()); 
        

    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                 : "block_builder_and_solver_linear_elastic_dynamic",
            "block_builder"                        : true,
            "diagonal_values_for_dirichlet_dofs"   : "use_max_diagonal",
            "silent_warnings"                      : false
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name() { return "block_builder_and_solver_linear_elastic_dynamic"; }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedBlockBuilderAndSolverLinearElasticDynamic";
    }

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

    void GetFirstAndSecondDerivativeVector(TSystemVectorType& rFirstDerivativeVector,
                                           TSystemVectorType& rSecondDerivativeVector,
                                           ModelPart&         rModelPart)
    {
        block_for_each(rModelPart.Nodes(), [&rFirstDerivativeVector, &rSecondDerivativeVector, this](Node& rNode) {
            if (rNode.IsActive()) {
                GetDerivativesForVariable(DISPLACEMENT_X, rNode, rFirstDerivativeVector, rSecondDerivativeVector);
                GetDerivativesForVariable(DISPLACEMENT_Y, rNode, rFirstDerivativeVector, rSecondDerivativeVector);

                const std::vector<const Variable<double>*> optional_variables = {
                    &ROTATION_X, &ROTATION_Y, &ROTATION_Z, &DISPLACEMENT_Z};

                for (const auto p_variable : optional_variables) {
                    GetDerivativesForOptionalVariable(*p_variable, rNode, rFirstDerivativeVector,
                                                      rSecondDerivativeVector);
                }
            }
        });
    }

    void GetDerivativesForOptionalVariable(const Variable<double>& rVariable,
                                           const Node&             rNode,
                                           TSystemVectorType&      rFirstDerivativeVector,
                                           TSystemVectorType&      rSecondDerivativeVector) const
    {
        if (rNode.HasDofFor(rVariable)) {
            GetDerivativesForVariable(rVariable, rNode, rFirstDerivativeVector, rSecondDerivativeVector);
        }
    }

    void GetDerivativesForVariable(const Variable<double>& rVariable,
                                   const Node&             rNode,
                                   TSystemVectorType&      rFirstDerivativeVector,
                                   TSystemVectorType&      rSecondDerivativeVector) const
    {
        const auto& first_derivative  = rVariable.GetTimeDerivative();
        const auto& second_derivative = first_derivative.GetTimeDerivative();

        const auto equation_id               = rNode.GetDof(rVariable).EquationId();
        rFirstDerivativeVector[equation_id]  = rNode.FastGetSolutionStepValue(first_derivative);
        rSecondDerivativeVector[equation_id] = rNode.FastGetSolutionStepValue(second_derivative);
    }

    void BuildRHSNoDirichlet(typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemVectorType& rb)
    {
        
        // Getting the Elements
        ElementsArrayType& r_elements = rModelPart.Elements();

        // getting the array of the conditions
        ConditionsArrayType& r_conditions = rModelPart.Conditions();

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // contributions to the system
        LocalSystemVectorType local_external_force = LocalSystemVectorType(0);

		TSystemVectorType external_force = TSystemVectorType(BaseType::mEquationSystemSize, 0.0);
        //TSystemVectorType delta_force = ZeroVector(BaseType::mEquationSystemSize);
        // vector containing the localization in the system of the different
        // terms
        Element::EquationIdVectorType equation_ids;

        // assemble all elements

        const int nelements = static_cast<int>(r_elements.size());
#pragma omp parallel firstprivate(nelements, local_external_force, equation_ids)
        {

            local_external_force.resize(0, false);

            // assemble all conditions
            const int nconditions = static_cast<int>(r_conditions.size());
#pragma omp for schedule(guided, 512)
            for (int i = 0; i < nconditions; i++) {
                auto it = r_conditions.begin() + i;
                // If the condition is active
                if (it->IsActive()) {
                    it->CalculateRightHandSide(local_external_force, r_current_process_info);
                    it->EquationIdVector(equation_ids, r_current_process_info);

                    // assemble the elemental contribution
                    BaseType::AssembleRHS(external_force, local_external_force, equation_ids);
                }
            }
        }

     
		mCurrentExternalForceVector = external_force;

        noalias(mCurrentOutOfBalanceVector) = mCurrentExternalForceVector - mPreviousExternalForceVector;

        AddMassAndDampingToRhs(rModelPart, mCurrentOutOfBalanceVector);

        noalias(rb) = mCurrentOutOfBalanceVector - mPreviousOutOfBalanceVector;


    }

private:
    TSystemMatrixType mMassMatrix;
    TSystemMatrixType mDampingMatrix;
    TSystemVectorType mPreviousExternalForceVector;
    TSystemVectorType mCurrentExternalForceVector;

	//TSystemVectorType mCurrentDx;

    TSystemVectorType mPreviousOutOfBalanceVector;
	TSystemVectorType mCurrentOutOfBalanceVector;

	//TSystemVectorType mOutOfBalanceVector;


    void InitializeDynamicMatrix(TSystemMatrixType&            rMatrix,
                                 unsigned int                  MatrixSize,
                                 typename TSchemeType::Pointer pScheme,
                                 ModelPart&                    rModelPart)
    {
        rMatrix.resize(MatrixSize, MatrixSize, false);
        BaseType::ConstructMatrixStructure(pScheme, rMatrix, rModelPart);
        TSparseSpace::SetToZero(rMatrix);
    }

    void CalculateAndAddDynamicContributionToRhs(TSystemVectorType& rSolutionVector,
                                                 TSystemMatrixType& rGlobalMatrix,
                                                 TSystemVectorType& rb)
    {
        TSystemVectorType contribution;
        contribution.resize(BaseType::mEquationSystemSize, false);
        TSparseSpace::SetToZero(contribution);
        TSparseSpace::Mult(rGlobalMatrix, rSolutionVector, contribution);

        TSparseSpace::UnaliasedAdd(rb, 1.0, contribution);
    }

    /**
     * @brief Function to add the mass and damping contribution to the rhs.
     * @details Damping contribution is the dot product of the global damping matrix and the first
     * derivative vector, Mass contribution is the dot product of the global mass matrix and the
     * second derivative vector
     * @param rModelPart The model part of the problem to solve
     * @param rb The RHS vector
     */
    void AddMassAndDampingToRhs(ModelPart& rModelPart, TSystemVectorType& rb)
    {
        // Get first and second derivative vector
        TSystemVectorType first_derivative_vector  = ZeroVector(BaseType::mEquationSystemSize);
        TSystemVectorType second_derivative_vector = ZeroVector(BaseType::mEquationSystemSize);
        GetFirstAndSecondDerivativeVector(first_derivative_vector, second_derivative_vector, rModelPart);

        // calculate and add mass and damping contribution to rhs
        CalculateAndAddDynamicContributionToRhs(second_derivative_vector, mMassMatrix, rb);
        CalculateAndAddDynamicContributionToRhs(first_derivative_vector, mDampingMatrix, rb);
    }

}; /* Class ResidualBasedBlockBuilderAndSolverLinearElasticDynamic */

///@}

///@name Type Definitions
///@{

// Here one should use the KRATOS_CREATE_LOCAL_FLAG, but it does not play nice with template parameters
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
const Kratos::Flags ResidualBasedBlockBuilderAndSolverLinearElasticDynamic<TSparseSpace, TDenseSpace, TLinearSolver>::SILENT_WARNINGS(
    Kratos::Flags::Create(0));

///@}

} /* namespace Kratos.*/

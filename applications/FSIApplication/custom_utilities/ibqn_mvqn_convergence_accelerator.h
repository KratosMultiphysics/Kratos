//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Ruben Zorrilla
//

#if !defined(KRATOS_IBQN_MVQN_CONVERGENCE_ACCELERATOR)
#define  KRATOS_IBQN_MVQN_CONVERGENCE_ACCELERATOR

// System includes

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "utilities/math_utils.h"
#include "utilities/qr_utility.h"
#include "utilities/parallel_utilities.h"

// Application includes

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

/// Forward declaration of standard MVQN
/// This is required since the include ward avoids the inclusion of the standard MVQN
template<class TSparseSpace, class TDenseSpace>
class MVQNFullJacobianConvergenceAccelerator;

/** @brief Interface Block Newton convergence accelerator
 * Interface Block Newton equations convergence accelerator
 * @tparam TSparseSpace Linear algebra sparse space
 * @tparam TDenseSpace Linear algebra dense space
 */
template<class TSparseSpace, class TDenseSpace>
class IBQNMVQNConvergenceAccelerator: public ConvergenceAccelerator<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( IBQNMVQNConvergenceAccelerator );

    typedef ConvergenceAccelerator<TSparseSpace, TDenseSpace>              BaseType;
    typedef typename BaseType::Pointer                              BaseTypePointer;

    typedef typename BaseType::DenseVectorType                           VectorType;
    typedef typename BaseType::DenseVectorPointerType             VectorPointerType;

    typedef typename BaseType::DenseMatrixType                           MatrixType;
    typedef typename BaseType::DenseMatrixPointerType             MatrixPointerType;

    typedef MVQNFullJacobianConvergenceAccelerator<TSparseSpace, TDenseSpace>          MVQNType;
    typedef typename MVQNType::Pointer                                          MVQNPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new IBQNMVQNConvergenceAccelerator object
     * MVQN convergence accelerator Json settings constructor
     * @param rConvAcceleratorParameters Json string encapsulating the settings
     */

    explicit IBQNMVQNConvergenceAccelerator(Parameters rParameters)
    {
        Parameters mvqn_default_parameters(R"({
            "solver_type"            : "IBQN_MVQN",
            "w_0"                    : 0.825,
            "abs_cut_off_tol"        : 1e-8
        })");
        rParameters.ValidateAndAssignDefaults(mvqn_default_parameters);

        mInitialOmega = rParameters["w_0"].GetDouble();
        const double abs_cut_off_tol = rParameters["abs_cut_off_tol"].GetDouble();
        mpConvergenceAcceleratorLeft = Kratos::make_unique<MVQNType>(mInitialOmega, abs_cut_off_tol, true);
        mpConvergenceAcceleratorRight = Kratos::make_unique<MVQNType>(mInitialOmega, abs_cut_off_tol, true);
    }

    /**
     * Copy Constructor.
     */
    IBQNMVQNConvergenceAccelerator(const IBQNMVQNConvergenceAccelerator& rOther) = delete;

    /**
     * Destructor.
     */
    virtual ~IBQNMVQNConvergenceAccelerator(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        KRATOS_TRY;

        mpConvergenceAcceleratorLeft->Initialize();
        mpConvergenceAcceleratorRight->Initialize();

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Initialize the internal iteration counter
     * This method initializes the convergence acceleratior iteration counter at the begining of the step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        mConvergenceAcceleratorIteration = 0;

        mpConvergenceAcceleratorLeft->InitializeSolutionStep();
        mpConvergenceAcceleratorRight->InitializeSolutionStep();

        KRATOS_CATCH( "" );
    }

    void InitializeNonLinearIteration() override
    {
        KRATOS_TRY;

        mpConvergenceAcceleratorLeft->InitializeNonLinearIteration();
        mpConvergenceAcceleratorRight->InitializeNonLinearIteration();

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Performs the solution update
     * This method computes the solution update using a Jacobian approximation.
     * Such Jacobian approximation is computed with the MVQN (MultiVector Quasi-Newton method) algorithm.
     * @param rResidualVector Residual vector from the residual evaluation
     * @param rIterationGuess Current iteration guess to be corrected. Should be initialized to zero outside the convergence accelerator.
     */
    void UpdateSolution(
        const VectorType& rResidualVector,
        VectorType& rIterationGuess) override
    {
        KRATOS_TRY;

        KRATOS_ERROR << "The UpdateSolution() method cannot be called in the interface block Newton case. Use either \'UpdateSolutionLeft\' or \'UpdateSolutionRight\'" << std::endl;

        KRATOS_CATCH( "" );
    }

    //TODO: CREATE AN INTERMEDIATE IBQN BASE CLASS
    virtual void UpdateSolutionRight(
        const VectorType& rForceInputVector,
        const VectorType& rDisplacementOutputVector,
        VectorType& rIterationGuess)
    {
        // Save the uncorrected force vector for the displacement update
        auto p_uncor_disp_aux = Kratos::make_shared<VectorType>(rDisplacementOutputVector);
        std::swap(mpUncorrectedDisplacementVector, p_uncor_disp_aux);

        // Update right inverse Jacobian approximation with the current information
        mpConvergenceAcceleratorRight->UpdateInverseJacobianApproximation(rForceInputVector, rDisplacementOutputVector);

        VectorType right_correction(TSparseSpace::Size(rIterationGuess));
        if (!mFirstRightCorrectionPerformed) {
            // Calculate the first correction as a fixed point relaxation
            TSparseSpace::SetToZero(right_correction);
            TSparseSpace::UnaliasedAdd(right_correction, mInitialOmega, rDisplacementOutputVector - rIterationGuess);

            // Update the first iteration correction flag
            mFirstRightCorrectionPerformed = true;

        } else {
            // Get the interface problem size
            std::size_t problem_size = mpConvergenceAcceleratorRight->GetProblemSize();

            // Get both left and inverse Jacobian
            auto p_inv_jac_right = mpConvergenceAcceleratorRight->pGetInverseJacobianApproximation();
            auto p_inv_jac_left = mpConvergenceAcceleratorLeft->pGetInverseJacobianApproximation();

            // Set the residual of the update problem to be solved
            VectorType aux_RHS = rDisplacementOutputVector - rIterationGuess;
            VectorType aux_right_onto_left(problem_size);
            VectorType force_iteration_update = *mpUncorrectedForceVector - rForceInputVector;
            TSparseSpace::Mult(*p_inv_jac_right, force_iteration_update, aux_right_onto_left);
            TSparseSpace::UnaliasedAdd(aux_RHS, 1.0, aux_right_onto_left);

            // Set the LHS of the update problem to be solved
            MatrixType aux_LHS = IdentityMatrix(problem_size, problem_size);
            IndexPartition<std::size_t>(problem_size).for_each(
                VectorType(problem_size),
                [&aux_LHS, &p_inv_jac_right, &p_inv_jac_left, &problem_size](std::size_t Col, VectorType& rAuxColumnTLS)
            {
                TDenseSpace::GetColumn(Col, *p_inv_jac_left, rAuxColumnTLS);
                for (std::size_t row = 0; row < problem_size; ++row) {
                    aux_LHS(row,Col) -= TDenseSpace::RowDot(row, *p_inv_jac_right, rAuxColumnTLS);
                }
            });

            // Calculate the correction
            // Do the QR decomposition of (I - J_{S}J_{F}) and solve for the force update
            QR<double, row_major> qr_util;
            qr_util.compute(problem_size, problem_size, &(aux_LHS)(0,0));
            qr_util.solve(&(aux_RHS)(0), &(right_correction)(0));
        }

        // Update the iteration guess
        TSparseSpace::UnaliasedAdd(rIterationGuess, 1.0, right_correction);
    }

    //TODO: CREATE AN INTERMEDIATE IBQN BASE CLASS
    virtual void UpdateSolutionLeft(
        const VectorType& rDisplacementInputVector,
        const VectorType& rForceOutputVector,
        VectorType& rIterationGuess)
    {
        KRATOS_TRY;

        // Save the uncorrected force vector for the displacement update
        auto p_uncor_force_aux = Kratos::make_shared<VectorType>(rForceOutputVector);
        std::swap(mpUncorrectedForceVector, p_uncor_force_aux);

        // Update left inverse Jacobian approximation with the current information
        mpConvergenceAcceleratorLeft->UpdateInverseJacobianApproximation(rDisplacementInputVector, rForceOutputVector);

        VectorType left_correction(TSparseSpace::Size(rIterationGuess));
        if (!mFirstLeftCorrectionPerformed) {
            // Do nothing in the first traction correction
            // This means to set as correction the difference between the guess and the previous iteration
            left_correction = rForceOutputVector - rIterationGuess;

            // Update the first iteration correction flag
            mFirstLeftCorrectionPerformed = true;

        } else {
            // Get the interface problem size
            std::size_t problem_size = mpConvergenceAcceleratorLeft->GetProblemSize();

            // Get both left and inverse Jacobian
            auto p_inv_jac_left = mpConvergenceAcceleratorLeft->pGetInverseJacobianApproximation();
            auto p_inv_jac_right = mpConvergenceAcceleratorRight->pGetInverseJacobianApproximation();

            // Set the residual of the update problem to be solved
            VectorType aux_RHS = rForceOutputVector - rIterationGuess;
            VectorType aux_left_onto_right(problem_size);
            VectorType displacement_iteration_update = *mpUncorrectedDisplacementVector - rDisplacementInputVector;
            TSparseSpace::Mult(*p_inv_jac_left, displacement_iteration_update, aux_left_onto_right);
            TSparseSpace::UnaliasedAdd(aux_RHS, 1.0, aux_left_onto_right);

            // Set the LHS of the update problem to be solved
            MatrixType aux_LHS = IdentityMatrix(problem_size, problem_size);
            IndexPartition<std::size_t>(problem_size).for_each(
                VectorType(problem_size),
                [&aux_LHS, &p_inv_jac_left, &p_inv_jac_right, &problem_size](std::size_t Col, VectorType& rAuxColumnTLS)
            {
                TDenseSpace::GetColumn(Col, *p_inv_jac_right, rAuxColumnTLS);
                for (std::size_t row = 0; row < problem_size; ++row) {
                    aux_LHS(row,Col) -= TDenseSpace::RowDot(row, *p_inv_jac_left, rAuxColumnTLS);
                }
            });

            // Calculate the correction
            // Do the QR decomposition of (I - J_{F}J_{S}) and solve for the force update
            QR<double, row_major> qr_util;
            qr_util.compute(problem_size, problem_size, &(aux_LHS)(0,0));
            qr_util.solve(&(aux_RHS)(0), &(left_correction)(0));
        }

        // Update the iteration guess
        TSparseSpace::UnaliasedAdd(rIterationGuess, 1.0, left_correction);

        KRATOS_CATCH("");
    }

    /**
     * @brief Do the MVQN variables update
     * Updates the MVQN iteration variables for the next non-linear iteration
     */
    void FinalizeNonLinearIteration() override
    {
        KRATOS_TRY;

        mConvergenceAcceleratorIteration += 1;

        mpConvergenceAcceleratorLeft->FinalizeNonLinearIteration();
        mpConvergenceAcceleratorRight->FinalizeNonLinearIteration();

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Save the current step Jacobian
     * This method saves the current step Jacobian as previous step Jacobian for the next time step iteration
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        mpConvergenceAcceleratorLeft->FinalizeSolutionStep();
        mpConvergenceAcceleratorRight->FinalizeSolutionStep();

        KRATOS_CATCH( "" );
    }

    void Finalize() override
    {
        KRATOS_TRY;

        mpConvergenceAcceleratorLeft->Finalize();
        mpConvergenceAcceleratorRight->Finalize();

        KRATOS_CATCH( "" );
    }

    void Clear() override
    {
        KRATOS_TRY;

        mpConvergenceAcceleratorLeft->Clear();
        mpConvergenceAcceleratorRight->Clear();

        KRATOS_CATCH( "" );
    }

    bool IsBlockNewton() const override
    {
        return true;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    double mInitialOmega;

    bool mFirstLeftCorrectionPerformed = false; // Indicates that the initial fixed point iteration has been already performed
    bool mFirstRightCorrectionPerformed = false; // Indicates that the initial fixed point iteration has been already performed
    unsigned int mConvergenceAcceleratorIteration = 0; // Convergence accelerator iteration counter

    MVQNPointerType mpConvergenceAcceleratorLeft;
    MVQNPointerType mpConvergenceAcceleratorRight;

    typename BaseType::VectorPointerType mpUncorrectedForceVector;
    typename BaseType::VectorPointerType mpUncorrectedDisplacementVector;


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
    ///@name Serialization
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
}; /* Class IBQNMVQNConvergenceAccelerator */


///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_IBQN_MVQN_CONVERGENCE_ACCELERATOR defined */

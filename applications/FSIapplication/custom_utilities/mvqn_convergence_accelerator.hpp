//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Original author:  Ruben Zorrilla
//

#if !defined(KRATOS_MVQN_CONVERGENCE_ACCELERATOR)
#define  KRATOS_MVQN_CONVERGENCE_ACCELERATOR

/* System includes */

/* External includes */
#include "utilities/math_utils.h"

/* Project includes */
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "convergence_accelerator.hpp"

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

/** @brief MVQN (MultiVectorQuasiNewton method) acceleration scheme
 */
template<class TSpace>
class MVQNFullJacobianConvergenceAccelerator: public ConvergenceAccelerator<TSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( MVQNFullJacobianConvergenceAccelerator );

    typedef ConvergenceAccelerator<TSpace>                                 BaseType;
    typedef typename BaseType::Pointer                              BaseTypePointer;

    typedef typename BaseType::VectorType                                VectorType;
    typedef typename BaseType::VectorPointerType                  VectorPointerType;

    typedef typename BaseType::MatrixType                                MatrixType;
    typedef typename BaseType::MatrixPointerType                  MatrixPointerType;

    //~ typedef typename BaseType::SizeType                                    SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     * Aitken convergence accelerator
     */
    MVQNFullJacobianConvergenceAccelerator( Parameters &rConvAcceleratorParameters )
    {
        Parameters mvqn_default_parameters(R"(
        {
            "solver_type" : "MVQN",
            "w_0"         : 0.825
        }
        )");

        rConvAcceleratorParameters.ValidateAndAssignDefaults(mvqn_default_parameters);

        mOmega_0 = rConvAcceleratorParameters["w_0"].GetDouble();
        mConvergenceAcceleratorIteration = 0;
        mConvergenceAcceleratorFirstCorrectionPerformed = false;
    }

    MVQNFullJacobianConvergenceAccelerator( double rOmegaInitial = 0.825 )
    {
        mOmega_0 = rOmegaInitial;
        mConvergenceAcceleratorIteration = 0;
        mConvergenceAcceleratorFirstCorrectionPerformed = false;
    }

    /**
     * Copy Constructor.
     */
    MVQNFullJacobianConvergenceAccelerator( const MVQNFullJacobianConvergenceAccelerator& rOther )
    {
        mOmega_0 = rOther.mOmega_0;
        mConvergenceAcceleratorIteration = 0;
        mConvergenceAcceleratorFirstCorrectionPerformed = false;
    }

    /**
     * Destructor.
     */
    virtual ~MVQNFullJacobianConvergenceAccelerator
    () {}

    ///@}

    ///@name Operators
    ///@{
    ///@}

    ///@name Operations
    ///@{

    /**
     * Initialize the internal iteration counter
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        mConvergenceAcceleratorIteration = 0;

        KRATOS_CATCH( "" );
    }

    /**
     * Performs the solution update
     * The correction is computed using a Jacobian approximation obtained with the MVQN (MultiVector Quasi-Newton method).
     * @param rResidualVector Residual vector from the residual evaluation
     * @param rIterationGuess Current iteration guess to be corrected. Should be initialized to zero outside the convergence accelerator.
     */
    void UpdateSolution(const VectorType& rResidualVector,
                        VectorType& rIterationGuess) override
    {
        KRATOS_TRY;

        mProblemSize = TSpace::Size(rResidualVector);

        VectorPointerType pAuxResidualVector(new VectorType(rResidualVector));
        VectorPointerType pAuxIterationGuess(new VectorType(rIterationGuess));
        std::swap(mpResidualVector_1, pAuxResidualVector);
        std::swap(mpIterationValue_1, pAuxIterationGuess);

        if (mConvergenceAcceleratorIteration == 0)
        {
            if (mConvergenceAcceleratorFirstCorrectionPerformed == false)
            {
                // The very first correction of the problem is done with a fixed point iteration
                TSpace::UnaliasedAdd(rIterationGuess, mOmega_0, *mpResidualVector_1);

                // Resize the Jacobian approximation
                MatrixPointerType pNewJac_n = MatrixPointerType(new MatrixType(mProblemSize,mProblemSize));
                std::swap(pNewJac_n,mpJac_n);

                // Initialize the Jacobian approximation matrix (exclusively done in the very fist iteration)
                for (unsigned int i = 0; i < mProblemSize; i++)
                {
                    for (unsigned int j = 0; j < mProblemSize; j++)
                    {
                        (*mpJac_n)(i,j) = 0.0;
                    }
                    (*mpJac_n)(i,i) = -1.0;
                }

                mConvergenceAcceleratorFirstCorrectionPerformed = true;
            }
            else
            {
                // Fist step correction is done with the previous step Jacobian
                VectorType AuxVec(mProblemSize);
                TSpace::Mult(*mpJac_n, *mpResidualVector_1, AuxVec);
                TSpace::UnaliasedAdd(rIterationGuess, -1.0, AuxVec);
            }
        }
        else
        {
            if (mConvergenceAcceleratorIteration == 1)
            {
                // Resize the observation matrices in accordance to the problem size
                MatrixPointerType pNewObsMatrixV = MatrixPointerType(new MatrixType(mProblemSize, 1));
                MatrixPointerType pNewObsMatrixW = MatrixPointerType(new MatrixType(mProblemSize, 1));

                std::swap(mpObsMatrixV,pNewObsMatrixV);
                std::swap(mpObsMatrixW,pNewObsMatrixW);

                //~ std::cout << "First observation matrices step fill" << std::endl;

                // First observation matrices fill
                for (unsigned int i = 0; i < mProblemSize; i++)
                {
                    (*mpObsMatrixV)(i,0) = (*mpResidualVector_1)(i) - (*mpResidualVector_0)(i);
                    (*mpObsMatrixW)(i,0) = (*mpIterationValue_1)(i) - (*mpIterationValue_0)(i);
                }
            }
            else
            {
                // Reshape the existent observation matrices
                if (mConvergenceAcceleratorIteration <= mProblemSize)
                {
                    //~ std::cout << "Current iteration info storage: observation matrices resized (column addition)" << std::endl;

                    MatrixPointerType pNewObsMatrixV = MatrixPointerType(new MatrixType(mProblemSize, mConvergenceAcceleratorIteration));
                    MatrixPointerType pNewObsMatrixW = MatrixPointerType(new MatrixType(mProblemSize, mConvergenceAcceleratorIteration));

                    MatrixType &NewObsMatrixV = *pNewObsMatrixV;
                    MatrixType &NewObsMatrixW = *pNewObsMatrixW;

                    // Recover the previous iterations information
                    for (unsigned int i = 0; i < mProblemSize; i++)
                    {
                        for (unsigned int j = 0; j < (mConvergenceAcceleratorIteration-1); j++)
                        {
                            NewObsMatrixV(i,j) = (*mpObsMatrixV)(i,j);
                            NewObsMatrixW(i,j) = (*mpObsMatrixW)(i,j);
                        }
                    }

                    // Fill the attached column with the current iteration information
                    for (unsigned int i = 0; i < mProblemSize; i++)
                    {
                        NewObsMatrixV(i, mConvergenceAcceleratorIteration-1) = (*mpResidualVector_1)(i) - (*mpResidualVector_0)(i);
                        NewObsMatrixW(i, mConvergenceAcceleratorIteration-1) = (*mpIterationValue_1)(i) - (*mpIterationValue_0)(i);
                    }

                    std::swap(mpObsMatrixV,pNewObsMatrixV);
                    std::swap(mpObsMatrixW,pNewObsMatrixW);
                }
                else
                {
                    //~ std::cout << "Current iteration info storage: observation matrices are not resized (oldest column is dropped)" << std::endl;

                    // Observation matrices size are close to the interface DOFs number. Old columns are to be dropped.
                    MatrixPointerType pNewObsMatrixV = MatrixPointerType(new MatrixType(mProblemSize, mProblemSize));
                    MatrixPointerType pNewObsMatrixW = MatrixPointerType(new MatrixType(mProblemSize, mProblemSize));

                    // Drop the oldest column and reorder data
                    for (unsigned int i = 0; i < mProblemSize; i++)
                    {
                        for (unsigned int j = 0; j < (mProblemSize-1); j++)
                        {
                            (*pNewObsMatrixV)(i,j) = (*mpObsMatrixV)(i,j+1);
                            (*pNewObsMatrixW)(i,j) = (*mpObsMatrixW)(i,j+1);
                        }
                    }

                    // Fill the last observation matrices column
                    for (unsigned int i = 0; i < mProblemSize; i++)
                    {
                        (*pNewObsMatrixV)(i, mProblemSize-1) = (*mpResidualVector_1)(i) - (*mpResidualVector_0)(i);
                        (*pNewObsMatrixW)(i, mProblemSize-1) = (*mpIterationValue_1)(i) - (*mpIterationValue_0)(i);
                    }

                    std::swap(mpObsMatrixV,pNewObsMatrixV);
                    std::swap(mpObsMatrixW,pNewObsMatrixW);
                }

            }

            // Compute the jacobian approximation
            MatrixType aux1( mProblemSize, mConvergenceAcceleratorIteration );
            MatrixType aux2( mConvergenceAcceleratorIteration, mConvergenceAcceleratorIteration );
            MatrixType aux2inv( mConvergenceAcceleratorIteration, mConvergenceAcceleratorIteration );
            MatrixType aux3( mProblemSize, mConvergenceAcceleratorIteration );

            //~ std::cout << "Jacobian approximation computation starts..." << std::endl;

            noalias(aux1) = *mpObsMatrixW - prod(*mpJac_n,*mpObsMatrixV); // Note: dense matrix product is not present in the space
            noalias(aux2) = prod(trans(*mpObsMatrixV),*mpObsMatrixV);     // Note: dense matrix product is not present in the space
            bool inversion_successful = InvertMatrix<>(aux2, aux2inv);

            if (inversion_successful == false)
            {
                KRATOS_ERROR << "Matrix inversion error within the Jacobian approximation computation!!";
            }

            noalias(aux3) = prod(aux1,aux2inv); // Note: dense matrix product is not present in the space

            MatrixPointerType pJac_k1 = MatrixPointerType(new MatrixType(mProblemSize, mProblemSize));

            std::swap(mpJac_k1,pJac_k1);

            noalias(*mpJac_k1) = *mpJac_n + prod(aux3,trans(*mpObsMatrixV)); // Note: dense matrix product is not present in the space

            //~ std::cout << "Jacobian approximation computed" << std::endl;

            // Perform the correction
            VectorType AuxVec(mProblemSize);
            TSpace::Mult(*mpJac_k1, *mpResidualVector_1, AuxVec);
            TSpace::UnaliasedAdd(rIterationGuess, -1.0, AuxVec);

        }

        KRATOS_CATCH( "" );
    }

    /**
     * Updates the MVQN iteration values for the next non-linear iteration
     */
    void FinalizeNonLinearIteration() override
    {
        KRATOS_TRY;

        // Variables update
        mpIterationValue_0 = mpIterationValue_1;
        mpResidualVector_0 = mpResidualVector_1;
        mConvergenceAcceleratorIteration += 1;

        KRATOS_CATCH( "" );
    }

    /**
     * Reset the convergence accelerator iterations counter
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        // Update previous time step Jacobian as the last iteration Jacobian.
        mpJac_n = mpJac_k1;

        KRATOS_CATCH( "" );
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

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}

    ///@name Protected member Variables
    ///@{

    double mOmega_0;                                            // Relaxation factor for the initial fixed point iteration
    unsigned int mProblemSize;                                  // Residual to minimize size
    unsigned int mConvergenceAcceleratorIteration;              // Convergence accelerator iteration counter
    bool mConvergenceAcceleratorFirstCorrectionPerformed;       // Indicates that the initial fixed point iteration has been already performed

    VectorPointerType mpResidualVector_0;       // Previous iteration residual vector
    VectorPointerType mpResidualVector_1;       // Current iteration residual vector
    VectorPointerType mpIterationValue_0;       // Previous iteration guess
    VectorPointerType mpIterationValue_1;       // Current iteration guess

    MatrixPointerType mpJac_n;                  // Previous step Jacobian approximation
    MatrixPointerType mpJac_k1;                 // Current iteration Jacobian approximation
    MatrixPointerType mpObsMatrixV;             // Residual increment observation matrix
    MatrixPointerType mpObsMatrixW;             // Solution increment observation matrix

    ///@}

    ///@name Protected Operators
    ///@{
    ///@}

    ///@name Protected Operations
    ///@{

    template<class T>
    bool InvertMatrix(const T& input, T& inverse)
    {
        typedef permutation_matrix<std::size_t> pmatrix;

        // create a working copy of the input
        T A(input);

        // create a permutation matrix for the LU-factorization
        pmatrix pm(A.size1());

        // perform LU-factorization
        int res = lu_factorize(A, pm);
        if (res != 0)
            return false;

        // create identity matrix of "inverse"
        inverse.assign(identity_matrix<double> (A.size1()));

        // backsubstitute to get the inverse
        lu_substitute(A, pm, inverse);

        return true;
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
    ///@{

private:

    ///@name Static Member Variables
    ///@{
    ///@}

    ///@name Member Variables
    ///@{
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

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Un accessible methods
    ///@{
    ///@}

}; /* Class MVQNFullJacobianConvergenceAccelerator */

///@}

///@name Type Definitions
///@{
///@}

///@name Input and output
///@{
///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_MVQN_CONVERGENCE_ACCELERATOR defined */

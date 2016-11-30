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
#include "boost/smart_ptr.hpp"
#include "utilities/math_utils.h"

/* Project includes */
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
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
    MVQNFullJacobianConvergenceAccelerator( double rOmegaInitial = 0.825 )
    {
        mOmega_0 = rOmegaInitial;
        mConvergenceAcceleratorStep = 0;
        mConvergenceAcceleratorIteration = 0;
    }

    /** 
     * Copy Constructor.
     */
    MVQNFullJacobianConvergenceAccelerator( const MVQNFullJacobianConvergenceAccelerator& rOther )
    {
        mOmega_0 = rOther.mOmega_0;
        mConvergenceAcceleratorStep = 0;
        mConvergenceAcceleratorIteration = 0;
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

        mConvergenceAcceleratorStep += 1;
        mConvergenceAcceleratorIteration = 0;   
                
        KRATOS_CATCH( "" );
    }

    /**
     * Performs the solution update 
     * The correction is computed using a Jacobian approximation obtained with the MVQN (MultiVector Quasi-Newton method).
     * @param rResidualVector: Residual vector from the residual evaluation
     * @param rIterationGuess: Current iteration guess to be corrected. Should be initialized to zero outside the convergence accelerator.
     */
    void UpdateSolution(
        const VectorType& rResidualVector,
        VectorType& rIterationGuess) override
    {
        KRATOS_TRY;
      
        mProblemSize = rResidualVector.size();
      
        mResidualVector_1 = rResidualVector;
        mIterationValue_1 = rIterationGuess;
        
        if (mConvergenceAcceleratorIteration == 0)
        {
            if (mConvergenceAcceleratorStep == 1)
            {
                // The very first correction of the problem is done with a point iteration 
                rIterationGuess += mOmega_0*mResidualVector_1;
                
                // Resize the Jacobian approximation 
                MatrixPointerType pNewJac_n = MatrixPointerType(new MatrixType(mProblemSize,mProblemSize));
                mJac_n.swap(*pNewJac_n);
                
                // Initialize the Jacobian approximation matrix (exclusively done in the very fist iteration)
                for (unsigned int i = 0; i < mProblemSize; i++)
                {
                    for (unsigned int j = 0; j < mProblemSize; j++)
                    {
                        mJac_n(i,j) = 0.0;
                    }
                    mJac_n(i,i) = -1.0;
                }
            }
            else
            {
                // Fist step correction is done with the previous step Jacobian
                rIterationGuess -= prod(mJac_n,mResidualVector_1);
            }
        }
        else
        {
            if (mConvergenceAcceleratorIteration == 1)
            {
                // Resize the observation matrices in accordance to the problem size                
                MatrixPointerType pNewObsMatrixV = MatrixPointerType(new MatrixType(mProblemSize, 1));
                MatrixPointerType pNewObsMatrixW = MatrixPointerType(new MatrixType(mProblemSize, 1));
                
                mObsMatrixV.swap(*pNewObsMatrixV);
                mObsMatrixW.swap(*pNewObsMatrixW);
                                                                                                 
                //~ std::cout << "First observation matrices step fill" << std::endl;
                
                // First observation matrices fill
                for (unsigned int i = 0; i < mProblemSize; i++)
                {
                    mObsMatrixV(i,0) = mResidualVector_1(i) - mResidualVector_0(i);
                    mObsMatrixW(i,0) = mIterationValue_1(i) - mIterationValue_0(i);
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
                            NewObsMatrixV(i,j) = (mObsMatrixV)(i,j);
                            NewObsMatrixW(i,j) = (mObsMatrixW)(i,j);
                        }
                    }
                                        
                    // Fill the attached column with the current iteration information
                    for (unsigned int i = 0; i < mProblemSize; i++)
                    {
                        NewObsMatrixV(i, mConvergenceAcceleratorIteration-1) = mResidualVector_1(i) - mResidualVector_0(i);
                        NewObsMatrixW(i, mConvergenceAcceleratorIteration-1) = mIterationValue_1(i) - mIterationValue_0(i);
                    }
                    
                    mObsMatrixV = NewObsMatrixV;
                    mObsMatrixW = NewObsMatrixW;
                }
                else
                {
                    //~ std::cout << "Current iteration info storage: observation matrices are not resized (oldest column is dropped)" << std::endl;
                    
                    // Observation matrices size are close to the interface DOFs number. Old columns are to be dropped.
                    MatrixPointerType pNewObsMatrixV = MatrixPointerType(new MatrixType(mProblemSize, mProblemSize));
                    MatrixPointerType pNewObsMatrixW = MatrixPointerType(new MatrixType(mProblemSize, mProblemSize));
                    
                    MatrixType &NewObsMatrixV = *pNewObsMatrixV;
                    MatrixType &NewObsMatrixW = *pNewObsMatrixW;
                                                        
                    // Drop the oldest column and reorder data
                    for (unsigned int i = 0; i < mProblemSize; i++)
                    {
                        for (unsigned int j = 0; j < (mProblemSize-1); j++)
                        {
                            NewObsMatrixV(i,j) = (mObsMatrixV)(i,j+1);
                            NewObsMatrixW(i,j) = (mObsMatrixW)(i,j+1);
                        }
                    }
                        
                    // Fill the last observation matrices column
                    for (unsigned int i = 0; i < mProblemSize; i++)
                    {
                        NewObsMatrixV(i, mProblemSize-1) = mResidualVector_1(i) - mResidualVector_0(i);
                        NewObsMatrixW(i, mProblemSize-1) = mIterationValue_1(i) - mIterationValue_0(i);
                    }
                                        
                    mObsMatrixV = NewObsMatrixV;
                    mObsMatrixW = NewObsMatrixW;
                }
            }
            
            // Compute the jacobian approximation
            MatrixType aux1( mProblemSize, mConvergenceAcceleratorIteration );
            MatrixType aux2( mConvergenceAcceleratorIteration, mConvergenceAcceleratorIteration );
            MatrixType aux2inv( mConvergenceAcceleratorIteration, mConvergenceAcceleratorIteration );
            MatrixType aux3( mProblemSize, mConvergenceAcceleratorIteration );
            
            //~ std::cout << "Jacobian approximation computation starts..." << std::endl;
            
            
            noalias(aux1) = mObsMatrixW - prod(mJac_n,mObsMatrixV); 
            noalias(aux2) = prod(trans(mObsMatrixV),mObsMatrixV);                        
            bool inversion_successful = InvertMatrix<>(aux2, aux2inv);
            
            if (inversion_successful == false) 
            {
                KRATOS_ERROR << "Matrix inversion error within the Jacobian approximation computation!!";
            }
            
            noalias(aux3) = prod(aux1,aux2inv);
            
            MatrixPointerType pJac_k1 = MatrixPointerType(new MatrixType(mProblemSize, mProblemSize));
            mJac_k1.swap(*pJac_k1);
            
            noalias(mJac_k1) = mJac_n + prod(aux3,trans(mObsMatrixV));
            
            //~ std::cout << "Jacobian approximation computed" << std::endl;
            
            // Perform the correction
            noalias(rIterationGuess) -= prod(mJac_k1,mResidualVector_1);
            
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
        mIterationValue_0 = mIterationValue_1;
        mResidualVector_0 = mResidualVector_1;
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
        noalias(mJac_n) = mJac_k1;

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

    unsigned int mProblemSize;
    unsigned int mConvergenceAcceleratorStep;
    unsigned int mConvergenceAcceleratorIteration;

    double mOmega_0;
    
    VectorType mResidualVector_0;           // Previous iteration residual vector
    VectorType mResidualVector_1;           // Current iteration residual vector
    VectorType mIterationValue_0;           // Previous iteration guess
    VectorType mIterationValue_1;           // Current iteration guess
    
    MatrixType mJac_n;                      // Previous step Jacobian approximation
    MatrixType mJac_k1;                     // Current iteration Jacobian approximation
    MatrixType mObsMatrixV;                 // Residual increment observation matrix
    MatrixType mObsMatrixW;                 // Solution increment observation matrix
    
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

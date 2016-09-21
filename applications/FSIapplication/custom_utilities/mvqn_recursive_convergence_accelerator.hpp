//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Original author:  Ruben Zorrilla Martinez 
//

#if !defined(KRATOS_MVQN_RECURSIVE_CONVERGENCE_ACCELERATOR)
#define  KRATOS_MVQN_RECURSIVE_CONVERGENCE_ACCELERATOR

/* System includes */

/* External includes */
#include "utilities/math_utils.h"
#include <cmath>
#include <numeric>

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

/** @brief Jacobian emulator
 */
template<class TSpace>
class JacobianEmulator
{
public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( JacobianEmulator );
    
    typedef typename TSpace::VectorType                 VectorType;
    typedef typename TSpace::MatrixType                 MatrixType;
    
    
    ///@}

    ///@name Public member Variables
    ///@{

    std::vector<VectorType>             mJacobianObsMatrixV; // Residual increment observation matrix
    std::vector<VectorType>             mJacobianObsMatrixW; // Solution increment observation matrix
    
    ///@}

    ///@name Life Cycle
    ///@{

    /**
     * Old Jacobian pointer constructor.
     * The Jacobian emulator will use information from the previous Jacobian
     */
    JacobianEmulator( boost::shared_ptr<JacobianEmulator> OldJacobianEmulatorPointer )
    {
        mpOldJacobianEmulator = OldJacobianEmulatorPointer;
    }

    /**
     * Old Jacobian pointer constructor with recursive previous Jacobian deleting.
     * The Jacobian emulator will use information from the previous Jacobian
     */
    JacobianEmulator( boost::shared_ptr<JacobianEmulator> OldJacobianEmulatorPointer, unsigned int EmulatorBufferSize )
    {        
        mpOldJacobianEmulator = OldJacobianEmulatorPointer;
                
        // Get the last pointer out of buffer
        JacobianEmulator::Pointer& p = mpOldJacobianEmulator->mpOldJacobianEmulator;
        for(unsigned int i = 0; i < (EmulatorBufferSize-1); i++)
        {
            JacobianEmulator::Pointer& p = p->mpOldJacobianEmulator;
        }
        
        // Destroy the one that goes out from max_steps
        p.reset();
    }

    /**
     * Empty constructor.
     * The Jacobian emulator will consider the identity matrix as previous Jacobian
     */
    JacobianEmulator( )
    {
    }

    /** 
     * Copy Constructor.
     */
    //~ JacobianEmulator( const JacobianEmulator& rOther, int max_steps )
    JacobianEmulator( const JacobianEmulator& rOther )
    {
        mpOldJacobianEmulator = rOther.mpOldJacobianEmulator;
    }

    /** 
     * Destructor.
     */
    virtual ~JacobianEmulator
    () {}

    ///@}
    
    ///@name Operators
    ///@{
    ///@}
    
    ///@name Operations
    ///@{

    /**
     * Projects the approximated Jacobian onto a vector
     * @param rWorkVector: Vector in where the Jacobian is projected
     */
    void ApplyPrevStepJacobian(
        const VectorType& rWorkVector,
        VectorType& rProjectedVector)
    {
        mpOldJacobianEmulator->ApplyJacobian(rWorkVector, rProjectedVector);
    }

    /**
     * Projects the approximated Jacobian onto a vector
     * @param rWorkVector: Vector in where the Jacobian is projected
     */
    void ApplyJacobian(
        const VectorType& rWorkVector,
        VectorType& rProjectedVector)
    {
        KRATOS_TRY;
        
        //~ std::cout << "STARTING ApplyJacobian..." << std::endl;
        
        unsigned int aux_size = mJacobianObsMatrixV.size();
        unsigned int ProblemSize = mJacobianObsMatrixV[0].size();        
        
        VectorType t(aux_size);
        VectorType z(aux_size);
        VectorType y(ProblemSize);
        VectorType w(ProblemSize);
        VectorType v(ProblemSize);
            
        MatrixType auxMat(aux_size, aux_size);
        MatrixType auxMatInv(aux_size, aux_size);
                                        
        for (unsigned int i = 0; i < aux_size; i++)
        {
            auxMat(i,i) = TSpace::Dot(mJacobianObsMatrixV[i], mJacobianObsMatrixV[i]);
                        
            for (unsigned int j = i+1; j < aux_size; j++)
            {
                auxMat(i,j) = TSpace::Dot(mJacobianObsMatrixV[i], mJacobianObsMatrixV[j]);
                auxMat(j,i) = auxMat(i,j);
            }
        }
                
        // CALL THE INVERT MATRIX FROM THE OTHER CLASS    
        bool inversion_successful = InvertMatrix<>(auxMat, auxMatInv);
        if (inversion_successful == false) 
        {
            KRATOS_ERROR << "Matrix inversion error within the Jacobian approximation computation!!";
        }
            
        // t = V_k.T*r_k
        for (unsigned int i = 0; i < aux_size; i++)
        {
            t(i) = TSpace::Dot(mJacobianObsMatrixV[i], rWorkVector);
        }
            
        // z = (V_k.T*V_k)^-1*t
        TSpace::Mult(auxMatInv, t, z);                                        
            
        // TODO: PARALLELIZE THIS OPERATION (it cannot be done with TSpace::Dot beacuse the obs matrices are stored by columns)
        // y = V_k*z - r_k
        for (unsigned int i = 0; i < ProblemSize; i++)
        {
            y(i) = 0.0;
            for (unsigned int j = 0; j < aux_size; j++)
            {
                y(i) += mJacobianObsMatrixV[j][i]*z[j];
            }
        }
        y -= rWorkVector;                                            

        //~ std::cout << "Applying previous step Jacobian..." << std::endl;
        if (mpOldJacobianEmulator == nullptr)
        {
            v = y;
        }
        else
        {
            mpOldJacobianEmulator->ApplyJacobian(y, v);
        }
        //~ std::cout << "Previous step Jacobian applied." << std::endl;
        
        // w = W_k*z
        for (unsigned int i = 0; i < ProblemSize; i++)
        {
            w(i) = 0.0;
            for (unsigned int j = 0; j < aux_size; j++)
            {
                w(i) += mJacobianObsMatrixW[j][i]*z(j);
            }
        }                                 
            
        rProjectedVector = (v - w); 

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
    
    boost::shared_ptr<JacobianEmulator>     mpOldJacobianEmulator;   // Pointer to the old Jacobian            

    ///@}
    
    ///@name Protected Operators
    ///@{
    ///@}
    
    ///@name Protected Operations
    ///@{
    ///@}
        
    /**
     * Utility for matrix inversion
     * @param input: Matrix to invert
     * @param inverse: Inverted matrix
     */ 
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
    
}; /* Class JacobianEmulator */


/** @brief MVQN (MultiVectorQuasiNewton method) acceleration scheme
 */
template<class TSpace>
class MVQNRecursiveJacobianConvergenceAccelerator: public ConvergenceAccelerator<TSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( MVQNRecursiveJacobianConvergenceAccelerator );
    
    typedef ConvergenceAccelerator<TSpace>                                             BaseType;
    typedef typename BaseType::Pointer                                          BaseTypePointer;
    
    typedef typename JacobianEmulator<TSpace>::Pointer              JacobianEmulatorPointerType;
    
    typedef typename BaseType::VectorType                                            VectorType;
    typedef typename BaseType::VectorPointerType                              VectorPointerType;

    typedef typename BaseType::MatrixType                                            MatrixType;
    typedef typename BaseType::MatrixPointerType                              MatrixPointerType;
    
    //~ typedef typename BaseType::SizeType                                    SizeType;   

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     * MVQN convergence accelerator
     */
    MVQNRecursiveJacobianConvergenceAccelerator( double rOmegaInitial = 0.825, unsigned int rJacobianBufferSize = 10 )
    {
        mOmega_0 = rOmegaInitial;
        mJacobianBufferSize = rJacobianBufferSize;
        mConvergenceAcceleratorStep = 0;
        mConvergenceAcceleratorIteration = 0;
    }

    /** 
     * Copy Constructor.
     */
    MVQNRecursiveJacobianConvergenceAccelerator( const MVQNRecursiveJacobianConvergenceAccelerator& rOther )
    {
        mOmega_0 = rOther.mOmega_0;
        mJacobianBufferSize = rOther.mJacobianBufferSize;
        mConvergenceAcceleratorStep = 0;
        mConvergenceAcceleratorIteration = 0;
    }

    /** 
     * Destructor.
     */
    virtual ~MVQNRecursiveJacobianConvergenceAccelerator
    () {}

    ///@}
    
    ///@name Operators
    ///@{
    ///@}
    
    ///@name Operations
    ///@{
    
    //~ /**
     //~ * Construct the initial Jacobian emulator --> IS IT NECESSARY? IF IT IS A MEMBER VAR. IT IS ALREADY EMPTY?¿
     //~ */
    void Initialize() override
    {
        KRATOS_TRY;

        JacobianEmulatorPointerType mpCurrentJacobianEmulatorPointer();   
                        
        KRATOS_CATCH( "" );
    }
    

    /**
     * Initialize the internal iteration counter
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        mConvergenceAcceleratorStep += 1;
        mConvergenceAcceleratorIteration = 0;   
        
        if (mConvergenceAcceleratorStep <= mJacobianBufferSize)
        {
            // Construct the Jacobian emulator
            typename JacobianEmulator<TSpace>::Pointer paux = boost::make_shared < JacobianEmulator<TSpace> > (mpCurrentJacobianEmulatorPointer); 
            mpCurrentJacobianEmulatorPointer.swap(paux);
        }
        else
        {
            // Construct the Jacobian emulator considering the recursive elimination
            typename JacobianEmulator<TSpace>::Pointer paux = boost::make_shared < JacobianEmulator<TSpace> > (mpCurrentJacobianEmulatorPointer, mJacobianBufferSize); 
            mpCurrentJacobianEmulatorPointer.swap(paux);
        }   
                
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
            }
            else
            {
                std::cout << "First step correction" << std::endl;
                VectorType InitialCorrection(mProblemSize);
                
                mpCurrentJacobianEmulatorPointer->ApplyPrevStepJacobian(mResidualVector_1, InitialCorrection);
                rIterationGuess -= InitialCorrection;
            }
        }
        else
        {
            std::cout << "Gathering information" << std::endl;
            
            // Gather the new observation matrices column information
            VectorType newColV(mProblemSize);
            VectorType newColW(mProblemSize);
            
            for (unsigned int i = 0; i < mProblemSize; i++)
            {
                newColV[i] = mResidualVector_1(i) - mResidualVector_0(i);
                newColW[i] = mIterationValue_1(i) - mIterationValue_0(i);
            }
                        
            // Observation matrices information filling
            if (mConvergenceAcceleratorIteration <= mProblemSize)
            {
                // Append the new information to the existent observation matrices
                (mpCurrentJacobianEmulatorPointer)->mJacobianObsMatrixV.push_back(newColV);
                (mpCurrentJacobianEmulatorPointer)->mJacobianObsMatrixW.push_back(newColW);
                
                std::cout << "Observation matrices new information appended" << std::endl;
            }
            else
            {                
                // Observation matrices size are close to the interface DOFs number. Old columns are to be dropped.
                for (unsigned int i = 0; i < (mProblemSize-1); i++)
                {
                    (mpCurrentJacobianEmulatorPointer)->mJacobianObsMatrixV[i] = (mpCurrentJacobianEmulatorPointer)->mJacobianObsMatrixV[i+1];
                    (mpCurrentJacobianEmulatorPointer)->mJacobianObsMatrixW[i] = (mpCurrentJacobianEmulatorPointer)->mJacobianObsMatrixW[i+1];
                }
                
                // Substitute the last column by the new information.
                (mpCurrentJacobianEmulatorPointer)->mJacobianObsMatrixV.back() = newColV;
                (mpCurrentJacobianEmulatorPointer)->mJacobianObsMatrixW.back() = newColW;
                                
                std::cout << "Observation matrices size is kept (oldest column is dropped)" << std::endl;
            }
            
            std::cout << "Jacobian approximation computation starts..." << std::endl;
            
            unsigned int aux_size = mpCurrentJacobianEmulatorPointer->mJacobianObsMatrixV.size();
            
            VectorType t(aux_size);
            VectorType z(aux_size);
            VectorType y(mProblemSize);
            VectorType w(mProblemSize);
            VectorType v(mProblemSize);
            
            MatrixType auxMat(aux_size, aux_size);
            MatrixType auxMatInv(aux_size, aux_size);
            
            for (unsigned int i = 0; i < aux_size; i++)
            {
                auxMat(i,i) = TSpace::Dot(mpCurrentJacobianEmulatorPointer->mJacobianObsMatrixV[i],
                                          mpCurrentJacobianEmulatorPointer->mJacobianObsMatrixV[i]);
                for (unsigned int j = i+1; j < aux_size; j++)
                {
                    auxMat(i,j) = TSpace::Dot(mpCurrentJacobianEmulatorPointer->mJacobianObsMatrixV[i],
                                              mpCurrentJacobianEmulatorPointer->mJacobianObsMatrixV[j]);
                    auxMat(j,i) = auxMat(i,j);
                }
            }
            
            bool inversion_successful = InvertMatrix<>(auxMat, auxMatInv);
                        
            if (inversion_successful == false) 
            {
                KRATOS_ERROR << "Matrix inversion error within the Jacobian approximation computation!!";
            }
                        
            // t = V_k.T*r_k
            for (unsigned int i = 0; i < aux_size; i++)
            {
                t(i) = TSpace::Dot(mpCurrentJacobianEmulatorPointer->mJacobianObsMatrixV[i], mResidualVector_1);
            }
                        
            // z = (V_k.T*V_k)^-1*t
            TSpace::Mult(auxMatInv, t, z);
            
            // TODO: PARALLELIZE THIS OPERATION (it cannot be done with TSpace::Dot beacuse the obs matrices are stored by columns)
            // y = V_k*z - r_k
            for (unsigned int i = 0; i < mProblemSize; i++)
            {
                y(i) = 0.0;
                for (unsigned int j = 0; j < aux_size; j++)
                {
                    y(i) += (mpCurrentJacobianEmulatorPointer->mJacobianObsMatrixV[j][i])*z[j];
                }
            }
            y -= mResidualVector_1;
                        
            std::cout << "ApplyJacobian call..." << std::endl;
            if (mpCurrentJacobianEmulatorPointer == nullptr)
            {
                v = y;
            }
            else
            {
                mpCurrentJacobianEmulatorPointer->ApplyJacobian(y, v);
            }
            std::cout << "ApplyJacobian finished." << std::endl;
            
            // w = W_k*z
            for (unsigned int i = 0; i < mProblemSize; i++)
            {
                w(i) = 0.0;
                for (unsigned int j = 0; j < aux_size; j++)
                {
                    w(i) += (mpCurrentJacobianEmulatorPointer->mJacobianObsMatrixW[j][i])*z(j);
                }
            }              
            
            // Perform the correction
            rIterationGuess += (v - w); 
            
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
    unsigned int mJacobianBufferSize;
    unsigned int mCurrentJacobianBufferSize;
    unsigned int mConvergenceAcceleratorStep;
    unsigned int mConvergenceAcceleratorIteration;

    double mOmega_0;
    
    VectorType mResidualVector_0;                                               // Previous iteration residual vector
    VectorType mResidualVector_1;                                               // Current iteration residual vector
    VectorType mIterationValue_0;                                               // Previous iteration guess
    VectorType mIterationValue_1;                                               // Current iteration guess
    
    JacobianEmulatorPointerType mpCurrentJacobianEmulatorPointer;               // Current step Jacobian approximator pointer
    
    
    //std::vector<VectorType> mObsMatrixV;            // Residual increment observation matrix
    //std::vector<VectorType> mObsMatrixW;            // Solution increment observation matrix
    
    ///@}
    
    ///@name Protected Operators
    ///@{
    ///@}
    
    ///@name Protected Operations
    ///@{
        
    /**
     * Utility for matrix inversion
     * @param input: Matrix to invert
     * @param inverse: Inverted matrix
     */ 
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
        lu_substitute(A, pm, inverse);                   // THE ERROR IS IN HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        return true;
    }

    ///@}
    
    ///@name Protected  Access
    ///@{
    ///@}JacobianEmulator
    
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
    
}; /* Class MVQNRecursiveJacobianConvergenceAccelerator */

///@}


///@}

///@name Type Definitions
///@{
///@}

///@name Input and output
///@{
///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_MVQN_RECURSIVE_CONVERGENCE_ACCELERATOR defined */

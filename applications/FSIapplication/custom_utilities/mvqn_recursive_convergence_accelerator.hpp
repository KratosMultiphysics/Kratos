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
#include "boost/smart_ptr.hpp"
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

/** @brief MVQN (MultiVectorQuasiNewton method) acceleration scheme
 */
template<class TSpace>
class MVQNRecursiveJacobianConvergenceAccelerator: public ConvergenceAccelerator<TSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( MVQNRecursiveJacobianConvergenceAccelerator );
    
    typedef ConvergenceAccelerator<TSpace>                                 BaseType;
    //~ typedef typename RecursiveJacobianCalculator<TSpace>     JacobianCalculatorType;
    
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
     * MVQN convergence accelerator
     */
    MVQNRecursiveJacobianConvergenceAccelerator( double rOmegaInitial = 0.825, unsigned int rJacobianBufferSize = 4 )
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

    /**
     * Initialize
     */
    void Initialize() override
    {
        KRATOS_TRY;

        mBufferObsMatrixV.reserve(mJacobianBufferSize);
        mBufferObsMatrixW.reserve(mJacobianBufferSize);
        
        mCurrentJacobianBufferSize = 0;
                
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
                InitialCorrection = ApproximateJacobianOntoVector(mResidualVector_1, mCurrentJacobianBufferSize);
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
                mObsMatrixV.push_back(newColV);
                mObsMatrixW.push_back(newColW);
                
                std::cout << "Observation matrices new information appended" << std::endl;
            }
            else
            {                
                // Observation matrices size are close to the interface DOFs number. Old columns are to be dropped.
                for (unsigned int i = 0; i < (mProblemSize-1); i++)
                {
                    mObsMatrixV[i] = mObsMatrixV[i+1];
                    mObsMatrixW[i] = mObsMatrixW[i+1];
                }
                
                // Substitute the last column by the new information.
                mObsMatrixV.back() = newColV;
                mObsMatrixW.back() = newColW;
                
                std::cout << "Observation matrices size is kept (oldest column is dropped)" << std::endl;
            }
            
            std::cout << "Jacobian approximation computation starts..." << std::endl;
            
            unsigned int aux_size = mObsMatrixV.size();
            
            VectorType t(aux_size);
            VectorType z(aux_size);
            VectorType y(mProblemSize);
            VectorType w(mProblemSize);
            VectorType v(mProblemSize);
            
            MatrixType auxMat(aux_size, aux_size);
            MatrixType auxMatInv(aux_size, aux_size);
                                        
            for (unsigned int i = 0; i < aux_size; i++)
            {
                auxMat(i,i) = TSpace::Dot(mObsMatrixV[i], mObsMatrixV[i]);
                for (unsigned int j = i+1; j < aux_size; j++)
                {
                    auxMat(i,j) = TSpace::Dot(mObsMatrixV[i], mObsMatrixV[j]);
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
                t(i) = TSpace::Dot(mObsMatrixV[i], mResidualVector_1);
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
                    y(i) += mObsMatrixV[j][i]*z[j];
                }
            }
            y -= mResidualVector_1;                                            

            v = ApproximateJacobianOntoVector(y, mCurrentJacobianBufferSize);
            
            // w = W_k*z
            for (unsigned int i = 0; i < mProblemSize; i++)
            {
                w(i) = 0.0;
                for (unsigned int j = 0; j < aux_size; j++)
                {
                    w(i) += mObsMatrixW[j][i]*z(j);
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
    
    /**
     * Reset the convergence accelerator iterations counter
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;
        
        std::cout << "Observation matrices buffer update starts..." << std::endl;

        // Update the observation matrices buffer
        UpdateJacobianBuffer( mObsMatrixV, mObsMatrixW);
        
        std::cout << "Observation matrices buffer filled." << std::endl;

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
    
    VectorType mResidualVector_0;                   // Previous iteration residual vector
    VectorType mResidualVector_1;                   // Current iteration residual vector
    VectorType mIterationValue_0;                   // Previous iteration guess
    VectorType mIterationValue_1;                   // Current iteration guess
    
    //~ std::vector<VectorType> mJac_n;                 // Previous step Jacobian approximation
    //~ std::vector<VectorType> mJac_k1;                // Current iteration Jacobian approximation
    std::vector<VectorType> mObsMatrixV;            // Residual increment observation matrix
    std::vector<VectorType> mObsMatrixW;            // Solution increment observation matrix
        
    std::vector < std::vector < VectorType > > mBufferObsMatrixV;         // Storage vector for mObsMatrixV matrices
    std::vector < std::vector < VectorType > > mBufferObsMatrixW;         // Storage vector for mObsMatrixW matrices
    
    ///@}
    
    ///@name Protected Operators
    ///@{
    ///@}
    
    ///@name Protected Operations
    ///@{
    
    /**
     * Approximates the previous Jacobian projection onto a vector operation 
     * The operation is computed with a recursive previous step Jacobian approximation using the previous full observation matrices
     * @param y: Vector to be projected on
     * @param y_proj: Projected vector
     */
    VectorType ApproximateJacobianOntoVector(
        const VectorType& rY,
        const unsigned int CurrentJacobianBufferSize)
    {
        KRATOS_TRY;
        
        std::cout << "Previous step Jacobian recursive approximation computation starts... " << std::endl;
        if (CurrentJacobianBufferSize == 0)
        {
            std::cout << "Considering identity matrix as previous step Jacobian approximation." << std::endl;
            return (rY);
        }
        else if (CurrentJacobianBufferSize == 1)
        {
            VectorType solA1, solB1;
            
            std::cout << "Computing previous step Jacobian approximation for buffer size 1..." << std::endl;
            
            solA1 = ApplyJacobianCoefsToVectorA(mBufferObsMatrixV[0], mBufferObsMatrixW[0], rY);
            solB1 = ApplyJacobianCoefsToVectorB(mBufferObsMatrixV[0], mBufferObsMatrixW[0], rY);
            
            std::cout << "Computed previous step Jacobian approximation for buffer size 1." << std::endl;
            
            return (solA1 - solB1);
        }
        else if (CurrentJacobianBufferSize == 2)
        {
            VectorType solA1, solA2, solB1, solB2;
            
            std::cout << "Computing previous step Jacobian approximation for buffer size 2..." << std::endl;
            
            solA1 = ApplyJacobianCoefsToVectorA(mBufferObsMatrixV[1], mBufferObsMatrixW[1], rY);
            //~ solB1 = ApplyJacobianCoefsToVectorB(mBufferObsMatrixV[1], mBufferObsMatrixW[1], rY);
            
            solA2 = ApplyJacobianCoefsToVectorA(mBufferObsMatrixV[0], mBufferObsMatrixV[0], solA1);
            solB2 = ApplyJacobianCoefsToVectorB(mBufferObsMatrixW[0], mBufferObsMatrixW[0], solA1);
            
            std::cout << "Computed previous step Jacobian approximation for buffer size 2." << std::endl;
            
            //~ return (solA2 - solB2 - solB1);
            return (solA2 - solB2);
        }
        else if (CurrentJacobianBufferSize == 3)
        {
            VectorType solA1, solA2, solA3, solB1, solB2, solB3;
            
            std::cout << "Computing previous step Jacobian approximation for buffer size 3..." << std::endl;
            
            //~ solA1 = ApplyJacobianCoefsToVectorA(mBufferObsMatrixV[2], mBufferObsMatrixW[2], rY);
            //~ solB1 = ApplyJacobianCoefsToVectorB(mBufferObsMatrixV[2], mBufferObsMatrixW[2], rY);
            
            solA2 = ApplyJacobianCoefsToVectorA(mBufferObsMatrixV[1], mBufferObsMatrixV[1], rY);
            solB2 = ApplyJacobianCoefsToVectorB(mBufferObsMatrixW[1], mBufferObsMatrixW[1], rY);
            
            solA3 = ApplyJacobianCoefsToVectorA(mBufferObsMatrixV[0], mBufferObsMatrixV[0], solA2);
            solB3 = ApplyJacobianCoefsToVectorB(mBufferObsMatrixW[0], mBufferObsMatrixW[0], solA2);
            
            std::cout << "Computed previous step Jacobian approximation for buffer size 3." << std::endl;
            
            //~ return (solA3 - solB3 - solB2 - solB1);
            return (solA3 - solB3 - solB2);
            
        }
        else if (CurrentJacobianBufferSize == 4)
        {
            VectorType solA1, solA2, solA3, solA4, solB1, solB2, solB3, solB4;
            
            std::cout << "Computing previous step Jacobian approximation for buffer size 4..." << std::endl;
            
            solA1 = ApplyJacobianCoefsToVectorA(mBufferObsMatrixV[3], mBufferObsMatrixW[3], rY);
            //~ solB1 = ApplyJacobianCoefsToVectorB(mBufferObsMatrixV[3], mBufferObsMatrixW[3], rY);
            
            solA2 = ApplyJacobianCoefsToVectorA(mBufferObsMatrixV[2], mBufferObsMatrixV[2], solA1);
            solB2 = ApplyJacobianCoefsToVectorB(mBufferObsMatrixW[2], mBufferObsMatrixW[2], solA1);
            
            solA3 = ApplyJacobianCoefsToVectorA(mBufferObsMatrixV[1], mBufferObsMatrixV[1], solA2);
            solB3 = ApplyJacobianCoefsToVectorB(mBufferObsMatrixW[1], mBufferObsMatrixW[1], solA2);
            
            solA4 = ApplyJacobianCoefsToVectorA(mBufferObsMatrixV[0], mBufferObsMatrixV[0], solA3);
            solB4 = ApplyJacobianCoefsToVectorB(mBufferObsMatrixW[0], mBufferObsMatrixW[0], solA3);
            
            std::cout << "Computed previous step Jacobian approximation for buffer size 4." << std::endl;
            
            //~ return (solA4 - solB4 - solB3 - solB2 - solB1);
            return (solA4 - solB4 - solB3 - solB2);
        }
                    
        KRATOS_CATCH( "" );
    }

    /**
     * Performs some auxiliar operations needed to compute the Jacobian approximation
     * @param ObsMatrixV: Residual increment observation matrix
     * @param ObsMatrixW: Solution increment observation matrix
     * @param WorkVector: Vector to be multiplied
     * @param solA: First term solution
     */     
    VectorType ApplyJacobianCoefsToVectorA(
        const std::vector<VectorType> ObsMatrixV,
        const std::vector<VectorType> ObsMatrixW,
        const VectorType WorkVector)
    {
        KRATOS_TRY;
     
        // TODO: THIS FUNCTION CAN BE USED TO COMPUTE THE CURRENT ITERATION COEFFICIENTS IN THE UpdateSolution FUNCTION OF THE CLASE ABOVE. MOVE TO THIS (avoid repeated code).
        
        std::cout << "Computing recursive approximation coefficient A..." << std::endl;
        
        unsigned int aux_size_1 = ObsMatrixV.size();
        unsigned int aux_size_2 = ObsMatrixV[0].size();
        
        VectorType t(aux_size_1);
        VectorType z(aux_size_1);
        
        MatrixType auxMat(aux_size_1, aux_size_1);
        MatrixType auxMatInv(aux_size_1, aux_size_1);
                                    
        // Note that the TSpace operations cannot be used to compute auxMat since it is stored as a std::vector<VectorType>
        for (unsigned int i = 0; i < aux_size_1; i++)
        {
            auxMat(i,i) = TSpace::Dot(ObsMatrixV[i], ObsMatrixV[i]);
            for (unsigned int j = i+1; j < aux_size_1; j++)
            {
                auxMat(i,j) = TSpace::Dot(ObsMatrixV[i], ObsMatrixV[j]);
                auxMat(j,i) = auxMat(i,j);
            }
        }
        
        bool inversion_successful = InvertMatrix<>(auxMat, auxMatInv);
        if (inversion_successful == false) 
        {
            KRATOS_ERROR << "Matrix inversion error within the Jacobian approximation computation!!";
        }
        
        // t = V_k.T*r_k
        for (unsigned int i = 0; i < aux_size_1; i++)
        {
            t(i) = TSpace::Dot(ObsMatrixV[i], WorkVector);
        }    
        
        // z = (V_k.T*V_k)^-1*t
        TSpace::Mult(auxMatInv, t, z);
        
        // Resize the solA vector
        VectorPointerType psolA = VectorPointerType(new VectorType(aux_size_2));
        VectorType &solA = *psolA;
        
        // solA vector computation
        // solA = V_k*z - r_k
        for (unsigned int i = 0; i < aux_size_2; i++)
        {
            solA(i) = 0.0;
            for (unsigned int j = 0; j < aux_size_1; j++)
            {
                solA(i) += ObsMatrixV[j][i]*z[j];
            }
        }
        solA -= WorkVector; 
        
        return solA;
            
        KRATOS_CATCH( "" );
    }
        
    /**
     * Performs some auxiliar operations needed to compute the Jacobian approximation
     * @param ObsMatrixV: Residual increment observation matrix
     * @param ObsMatrixW: Solution increment observation matrix
     * @param WorkVector: Vector to be multiplied
     * @param solB: Second term solution
     */     
    VectorType ApplyJacobianCoefsToVectorB(
        const std::vector<VectorType> ObsMatrixV,
        const std::vector<VectorType> ObsMatrixW,
        const VectorType WorkVector)
    {
        KRATOS_TRY;
     
        // TODO: THIS FUNCTION CAN BE USED TO COMPUTE THE CURRENT ITERATION COEFFICIENTS IN THE UpdateSolution FUNCTION OF THE CLASE ABOVE. MOVE TO THIS (avoid repeated code).
        
        unsigned int aux_size_1 = ObsMatrixV.size();
        unsigned int aux_size_2 = ObsMatrixV[0].size();
        
        VectorType t(aux_size_1);
        VectorType z(aux_size_1);
        
        MatrixType auxMat(aux_size_1, aux_size_1);
        MatrixType auxMatInv(aux_size_1, aux_size_1);
                                    
        // Note that the TSpace operations cannot be used to compute auxMat since it is stored as a std::vector<VectorType>
        for (unsigned int i = 0; i < aux_size_1; i++)
        {
            for (unsigned int j = 0; j < aux_size_1; j++)
            {
                auxMat(i,j) = TSpace::Dot(ObsMatrixV[i], ObsMatrixV[j]);
            }
        }
        
        bool inversion_successful = InvertMatrix<>(auxMat, auxMatInv);
        if (inversion_successful == false) 
        {
            KRATOS_ERROR << "Matrix inversion error within the Jacobian approximation computation!!";
        }
        
        // t = V_k.T*r_k
        for (unsigned int i = 0; i < aux_size_1; i++)
        {
            t(i) = TSpace::Dot(ObsMatrixV[i], WorkVector);
        }          
        
        // z = (V_k.T*V_k)^-1*t
        TSpace::Mult(auxMatInv, t, z);
        
        // Resize the solA vector
        VectorPointerType psolB = VectorPointerType(new VectorType(aux_size_2));
        VectorType &solB = *psolB;
        
        // solB vector computation
        // solB = W_k*z
        for (unsigned int i = 0; i < aux_size_2; i++)
        {
            solB(i) = 0.0;
            for (unsigned int j = 0; j < aux_size_1; j++)
            {
                solB(i) += ObsMatrixW[j][i]*z(j);
            }
        } 
        
        return solB;  
        
        KRATOS_CATCH( "" );
    }
    
    /**
     * Once convergence has been achieved, store the last observation matrices
     * @param ObsMatrixV: Residual increment observation matrix
     * @param ObsMatrixW: Solution increment observation matrix
     */ 
    void UpdateJacobianBuffer(
        const std::vector<VectorType> ObsMatrixV,
        const std::vector<VectorType> ObsMatrixW)
        {
            KRATOS_TRY;
            
            // Update the current buffer value
            if (mCurrentJacobianBufferSize < mJacobianBufferSize)
            {
                mCurrentJacobianBufferSize += 1;
            }
            
            // Update the buffer content
            if (mJacobianBufferSize != 0)
            {
                if (mBufferObsMatrixV.size() < mCurrentJacobianBufferSize)  // Append the observation matrices
                {
                    mBufferObsMatrixV.push_back(ObsMatrixV);
                    mBufferObsMatrixW.push_back(ObsMatrixW);
                }
                else                                                        // Drop the oldest observation matrices
                {
                    for (unsigned int i = 0; i < (mJacobianBufferSize-1); i++)
                    {
                        mBufferObsMatrixV[i] = mBufferObsMatrixV[i+1];
                        mBufferObsMatrixW[i] = mBufferObsMatrixW[i+1];
                    }
                    mBufferObsMatrixV[mJacobianBufferSize-1] = ObsMatrixV;
                    mBufferObsMatrixW[mJacobianBufferSize-1] = ObsMatrixW;
                }
            }

            KRATOS_CATCH( "" );
        }
        
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
    
}; /* Class MVQNRecursiveJacobianConvergenceAccelerator */

///@}

///@name Type Definitions
///@{
///@}

///@name Input and output
///@{
///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_MVQN_RECURSIVE_CONVERGENCE_ACCELERATOR defined */

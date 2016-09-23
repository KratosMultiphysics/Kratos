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
//~ #include <numeric>

/* Project includes */
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "convergence_accelerator.hpp"

#include "custom_utilities/qr_utility.h"            //QR decomposition utility used in matrix inversion.

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
    //KRATOS_CLASS_POINTER_DEFINITION( JacobianEmulator );
    typedef typename std::unique_ptr< JacobianEmulator<TSpace> > Pointer;
    
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
    JacobianEmulator( JacobianEmulator::Pointer&& OldJacobianEmulatorPointer )
    {
        mpOldJacobianEmulator = std::unique_ptr<JacobianEmulator<TSpace> >(std::move(OldJacobianEmulatorPointer));
    }

    /**
     * Old Jacobian pointer constructor with recursive previous Jacobian deleting.
     * The Jacobian emulator will use information from the previous Jacobian
     */
    JacobianEmulator( JacobianEmulator::Pointer&& OldJacobianEmulatorPointer, unsigned int EmulatorBufferSize )
    {        
        mpOldJacobianEmulator = std::unique_ptr<JacobianEmulator<TSpace> >(std::move(OldJacobianEmulatorPointer));
                
        // Get the last pointer out of buffer
        JacobianEmulator* p = (mpOldJacobianEmulator->mpOldJacobianEmulator).get();
        if(EmulatorBufferSize > 1)
        {
            for(unsigned int i = 1; i < (EmulatorBufferSize); i++)
            {
                if(i == EmulatorBufferSize-1)
                {
                    (p->mpOldJacobianEmulator).release();
                }
                else
                {
                    p = (p->mpOldJacobianEmulator).get();
                }
            }
        }
        else // If Jacobian buffer size equals 1 directly destroy the previous one
        {
            (mpOldJacobianEmulator->mpOldJacobianEmulator).release();
        }
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
        //~ std::cout << "ApplyPrevStepJacobian starts..." << std::endl;
        mpOldJacobianEmulator->ApplyJacobian(rWorkVector, rProjectedVector);
        //~ std::cout << "ApplyPrevStepJacobian finished." << std::endl;
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

        VectorType y(mJacobianObsMatrixV[0].size());
        VectorType w(mJacobianObsMatrixV[0].size());
        VectorType v(mJacobianObsMatrixV[0].size());
        Vector WorkVectorCopy(rWorkVector);
        Vector zQR(mJacobianObsMatrixV.size());
        
        Matrix auxMatQR(mJacobianObsMatrixV[0].size(), mJacobianObsMatrixV.size());
          
        // Loop to store a std::vector<VectorType> type as Matrix type
        for (unsigned int i = 0; i < mJacobianObsMatrixV[0].size(); i++)
        {
            for (unsigned int j = 0; j < mJacobianObsMatrixV.size(); j++)
            {
                auxMatQR(i,j) = mJacobianObsMatrixV[j](i);
            }
        }
        
        // QR decomposition to compute ((V_k.T*V_k)^-1)*V_k.T*r_k
        QR<double, row_major> QR_decomposition;
        QR_decomposition.compute(mJacobianObsMatrixV[0].size(), mJacobianObsMatrixV.size(), &auxMatQR(0,0)); 
        QR_decomposition.solve(&WorkVectorCopy[0], &zQR[0]);
            
        // TODO: PARALLELIZE THIS OPERATION (it cannot be done with TSpace::Dot beacuse the obs matrices are stored by columns)
        // y = V_k*zQR - r_k
        for (unsigned int i = 0; i < mJacobianObsMatrixV[0].size(); i++)
        {
            y(i) = 0.0;
            for (unsigned int j = 0; j < mJacobianObsMatrixV.size(); j++)
            {
                y(i) += mJacobianObsMatrixV[j][i]*zQR[j];
            }
        }
        y -= rWorkVector;

        //~ std::cout << "Applying previous step Jacobian..." << std::endl;
        //~ std::cout << "mpOldJacobianEmulator VALUE IN APPLYJACOBIAN EMULATOR: " << &mpOldJacobianEmulator << std::endl;
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
        for (unsigned int i = 0; i < mJacobianObsMatrixV[0].size(); i++)
        {
            w(i) = 0.0;
            for (unsigned int j = 0; j < mJacobianObsMatrixV.size(); j++)
            {
                w(i) += mJacobianObsMatrixW[j][i]*zQR(j);
            }
        }              
            
        rProjectedVector = (v - w); 

        KRATOS_CATCH( "" );
        
        //~ std::cout << "FINISHED ApplyJacobian..." << std::endl;
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
    
    std::unique_ptr<JacobianEmulator>     mpOldJacobianEmulator;   // Pointer to the old Jacobian            

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
            mpCurrentJacobianEmulatorPointer = std::unique_ptr< JacobianEmulator<TSpace> > (new JacobianEmulator<TSpace>(std::move(mpCurrentJacobianEmulatorPointer))); 
        }
        else
        {
            // Construct the Jacobian emulator considering the recursive elimination
            mpCurrentJacobianEmulatorPointer = std::unique_ptr< JacobianEmulator<TSpace> > (new JacobianEmulator<TSpace>(std::move(mpCurrentJacobianEmulatorPointer), mJacobianBufferSize)); 
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
                //~ KRATOS_WATCH(InitialCorrection)
                
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
            
            VectorType IterationCorrection(mProblemSize);
            mpCurrentJacobianEmulatorPointer->ApplyJacobian(rResidualVector, IterationCorrection);
            
            rIterationGuess += IterationCorrection;            
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

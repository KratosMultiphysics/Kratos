//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Original author:  Ruben Zorrilla
//

#if !defined(KRATOS_MVQN_RECURSIVE_CONVERGENCE_ACCELERATOR)
#define  KRATOS_MVQN_RECURSIVE_CONVERGENCE_ACCELERATOR

/* System includes */

/* External includes */
//~ #include "utilities/math_utils.h"
//~ #include <cmath>
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
    typedef typename std::unique_ptr< JacobianEmulator<TSpace> >        Pointer;

    typedef typename TSpace::VectorType                              VectorType;
    typedef typename TSpace::VectorPointerType                VectorPointerType;

    typedef typename TSpace::MatrixType                              MatrixType;
    typedef typename TSpace::MatrixPointerType                MatrixPointerType;


    ///@}

    ///@name Public member Variables
    ///@{
    ///@}

    ///@name Life Cycle
    ///@{

    /**
     * Old Jacobian pointer constructor.
     * The inverse Jacobian emulator will use information from the previous Jacobian
     */
    JacobianEmulator( Pointer&& OldJacobianEmulatorPointer )
    {
        mpOldJacobianEmulator = std::unique_ptr<JacobianEmulator<TSpace> >(std::move(OldJacobianEmulatorPointer));
    }

    /**
     * Old Jacobian pointer constructor with recursive previous Jacobian deleting.
     * The inverse Jacobian emulator will use information from the previous Jacobian
     */
    JacobianEmulator( Pointer&& OldJacobianEmulatorPointer, const unsigned int EmulatorBufferSize )
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
                    (p->mpOldJacobianEmulator).reset();
                    //~ std::cout << "Out of buffer Jacobian emulator reset." << std::endl;
                }
                else
                {
                    p = (p->mpOldJacobianEmulator).get();
                }
            }
        }
        else // If Jacobian buffer size equals 1 directly destroy the previous one
        {
            (mpOldJacobianEmulator->mpOldJacobianEmulator).reset();
            //~ std::cout << "Out of buffer Jacobian emulator reset." << std::endl;
        }
    }

    /**
     * Empty constructor.
     * The Jacobian emulator will consider minus the identity matrix as previous Jacobian
     */
    JacobianEmulator( )
    {
    }

    /**
     * Copy Constructor.
     */
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
     * Projects the previous step approximated inverse Jacobian onto a vector
     * @param rWorkVector: Vector in where the inverse Jacobian is to be projected
     * @param rProjectedVector: Projected vector output
     */
    // void ApplyPrevStepJacobian(const VectorType& rWorkVector,
    //                            VectorType& rProjectedVector)
    void ApplyPrevStepJacobian(const VectorPointerType pWorkVector,
                               VectorPointerType pProjectedVector)
    {
        // Security check for the empty observation matrices case (when no correction has been done in the previous step)
        if (mpOldJacobianEmulator->mJacobianObsMatrixV.size() != 0)
        {
            // mpOldJacobianEmulator->ApplyJacobian(rWorkVector, rProjectedVector);
            mpOldJacobianEmulator->ApplyJacobian(pWorkVector, pProjectedVector);
        }
        else
        {
            // TSpace::Assign(rProjectedVector, -1.0, rWorkVector); // Consider minus the identity matrix as inverse Jacobian
            TSpace::Assign(*pProjectedVector, -1.0, *pWorkVector); // Consider minus the identity matrix as inverse Jacobian
        }
    }

    /**
     * Projects the approximated inverse Jacobian onto a vector
     * @param rWorkVector: Vector in where the inverse Jacobian is to be projected
     * @param rProjectedVector: Projected vector output
     */
    // void ApplyJacobian(const VectorType& rWorkVector,
    //                    VectorType& rProjectedVector)
    void ApplyJacobian(const VectorPointerType pWorkVector,
                       VectorPointerType pProjectedVector)
    {
        KRATOS_TRY;

        // Security check for the empty observation matrices case (when no correction has been done in the previous step)
        if (mJacobianObsMatrixV.size() == 0)
        {
            if (mpOldJacobianEmulator != nullptr) // If it is available, consider the previous step Jacobian
            {
                // mpOldJacobianEmulator->ApplyJacobian(rWorkVector, rProjectedVector);
                mpOldJacobianEmulator->ApplyJacobian(pWorkVector, pProjectedVector);
            }
            else // When the JacobianEmulator has no PreviousJacobianEmulator consider minus the identity matrix as inverse Jacobian
            {
                // TSpace::Assign(rProjectedVector,-1.0,rWorkVector);
                TSpace::Assign(*pProjectedVector, -1.0, *pWorkVector);
            }
        }
        else
        {
            // VectorType y(TSpace::Size(mJacobianObsMatrixV[0]));
            // VectorType w(TSpace::Size(mJacobianObsMatrixV[0]));
            // VectorType zQR(mJacobianObsMatrixV.size());
            // VectorType WorkVectorCopy(rWorkVector);
            // Matrix auxMatQR(TSpace::Size(mJacobianObsMatrixV[0]), mJacobianObsMatrixV.size());

            VectorPointerType pY(new VectorType(mJacobianObsMatrixV[0]));
            VectorPointerType pW(new VectorType(mJacobianObsMatrixV[0]));
            VectorPointerType pzQR(new VectorType(mJacobianObsMatrixV.size()));
            MatrixPointerType pAuxMatQR(new MatrixType(TSpace::Size(mJacobianObsMatrixV[0]), mJacobianObsMatrixV.size()));

            // Loop to store a std::vector<VectorType> type as Matrix type
            for (unsigned int i = 0; i < TSpace::Size(mJacobianObsMatrixV[0]); i++)
            {
                for (unsigned int j = 0; j < mJacobianObsMatrixV.size(); j++)
                {
                    // auxMatQR(i,j) = mJacobianObsMatrixV[j](i);
                    (*pAuxMatQR)(i,j) = mJacobianObsMatrixV[j](i);
                }
            }

            VectorPointerType pWorkVectorCopy(new VectorType(*pWorkVector));

            // QR decomposition to compute ((V_k.T*V_k)^-1)*V_k.T*r_k
            // TODO: Implement an if(!frozen) to avoid the recomputation of the inverse matrix each time the OldJacobianEmulator is called
            // mQR_decomposition.compute(TSpace::Size(mJacobianObsMatrixV[0]), mJacobianObsMatrixV.size(), &auxMatQR(0,0));
            // mQR_decomposition.solve(&WorkVectorCopy[0], &zQR[0]);
            mQR_decomposition.compute(TSpace::Size(mJacobianObsMatrixV[0]), mJacobianObsMatrixV.size(), &(*pAuxMatQR)(0,0));
            mQR_decomposition.solve(&(*pWorkVectorCopy)(0), &(*pzQR)(0));

            // TODO: PARALLELIZE THIS OPERATION (it cannot be done with TSpace::Dot beacuse the obs matrices are stored by columns)
            // y = V_k*zQR - r_k
            // TSpace::SetToZero(y);
            TSpace::SetToZero(*pY);
            for (unsigned int j = 0; j < mJacobianObsMatrixV.size(); j++)
            {
                // TSpace::UnaliasedAdd(y,zQR(j), mJacobianObsMatrixV[j]);
                TSpace::UnaliasedAdd(*pY, (*pzQR)(j), mJacobianObsMatrixV[j]);
            }
            // TSpace::UnaliasedAdd(y, -1.0, rWorkVector);
            TSpace::UnaliasedAdd(*pY, -1.0, *pWorkVector);

            if (mpOldJacobianEmulator == nullptr)
            {
                // rProjectedVector = y; // Consider minus the identity as previous step Jacobian
                TSpace::Copy(*pY, *pProjectedVector); // Consider minus the identity as previous step Jacobian
            }
            else
            {
                // VectorType y_minus(y);
                // TSpace::Assign(y_minus,-1.0,y);
                VectorPointerType pYminus(new VectorType(*pY));
                TSpace::Assign(*pYminus, -1.0, *pY);
                // mpOldJacobianEmulator->ApplyJacobian(y_minus, rProjectedVector); // The minus comes from the fact that we want to apply r_k - V_k*zQR
                mpOldJacobianEmulator->ApplyJacobian(pYminus, pProjectedVector); // The minus comes from the fact that we want to apply r_k - V_k*zQR
            }

            // w = W_k*z
            // TSpace::SetToZero(w);
            TSpace::SetToZero(*pW);
            for (unsigned int j = 0; j < mJacobianObsMatrixV.size(); j++)
            {
                // TSpace::UnaliasedAdd(w, zQR(j), mJacobianObsMatrixW[j]);
                TSpace::UnaliasedAdd(*pW, (*pzQR)(j), mJacobianObsMatrixW[j]);
            }

            // rProjectedVector += w;
            // TSpace::UnaliasedAdd(*pProjectedVector, 1.0, w);
            TSpace::UnaliasedAdd(*pProjectedVector, 1.0, *pW);

        }

        KRATOS_CATCH( "" );

    }

    /**
    * Appends a new column to the observation matrix V
    * @param newColV: new column to be appended
    */
    void AppendColToV(const VectorType& rNewColV)
    {
        KRATOS_TRY;

        mJacobianObsMatrixV.push_back(rNewColV);

        KRATOS_CATCH( "" );
    }

    /**
    * Appends a new column to the observation matrix W
    * @param newColW: new column to be appended
    */
    void AppendColToW(const VectorType& rNewColW)
    {
        KRATOS_TRY;

        mJacobianObsMatrixW.push_back(rNewColW);

        KRATOS_CATCH( "" );
    }

    /**
    * Drops the oldest column and appends a new column to the observation matrix V
    * @param newColV: new column to be appended
    */
    void DropAndAppendColToV(const VectorType& rNewColV)
    {
        KRATOS_TRY;

        // Observation matrices size are close to the interface DOFs number. Old columns are to be dropped.
        for (unsigned int i = 0; i < (TSpace::Size(mJacobianObsMatrixV[0])-1); i++)
        {
            mJacobianObsMatrixV[i] = mJacobianObsMatrixV[i+1];
        }

        // Substitute the last column by the new information.
        mJacobianObsMatrixV.back() = rNewColV;

        KRATOS_CATCH( "" );
    }

    /**
    * Drops the oldest column and appends a new column to the observation matrix W
    * @param newColW: new column to be appended
    */
    void DropAndAppendColToW(const VectorType& rNewColW)
    {
        KRATOS_TRY;

        // Observation matrices size are close to the interface DOFs number. Old columns are to be dropped.
        for (unsigned int i = 0; i < (TSpace::Size(mJacobianObsMatrixV[0])-1); i++)
        {
            mJacobianObsMatrixW[i] = mJacobianObsMatrixW[i+1];
        }

        // Substitute the last column by the new information.
        mJacobianObsMatrixW.back() = rNewColW;

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
    QR<double, row_major>                 mQR_decomposition;        // QR decomposition object

    Pointer                           mpOldJacobianEmulator;        // Pointer to the old Jacobian

    std::vector<VectorType>             mJacobianObsMatrixV;        // Residual increment observation matrix
    std::vector<VectorType>             mJacobianObsMatrixW;        // Solution increment observation matrix
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
     //~ * Construct the initial inverse Jacobian emulator
     //~ */
    void Initialize() override
    {
        KRATOS_TRY;

        mpCurrentJacobianEmulatorPointer = std::unique_ptr< JacobianEmulator <TSpace> > (new JacobianEmulator<TSpace>());

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
            // Construct the inverse Jacobian emulator
            mpCurrentJacobianEmulatorPointer = std::unique_ptr< JacobianEmulator<TSpace> > (new JacobianEmulator<TSpace>(std::move(mpCurrentJacobianEmulatorPointer)));
        }
        else
        {
            // Construct the inverse Jacobian emulator considering the recursive elimination
            mpCurrentJacobianEmulatorPointer = std::unique_ptr< JacobianEmulator<TSpace> > (new JacobianEmulator<TSpace>(std::move(mpCurrentJacobianEmulatorPointer), mJacobianBufferSize));
        }

        KRATOS_CATCH( "" );
    }

    /**
     * Performs the solution update
     * The correction is computed using an inverse Jacobian approximation obtained with a recursive matrix-free version of the MVQN (MultiVector Quasi-Newton method).
     * @param rResidualVector: Residual vector from the residual evaluation
     * @param rIterationGuess: Current iteration guess to be corrected. Should be initialized outside the convergence accelerator.
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
            if (mConvergenceAcceleratorStep == 1)
            {
                // The very first correction of the problem is done with a fixed point iteration
                // rIterationGuess += mOmega_0*mResidualVector_1;
                TSpace::UnaliasedAdd(rIterationGuess, mOmega_0, *mpResidualVector_1);
            }
            else
            {
                //~ std::cout << "First step correction" << std::endl;
                // VectorType InitialCorrection(mProblemSize);
                VectorPointerType pInitialCorrection(new VectorType(rResidualVector)); // Hack in case the Epetra copy constructor needs to be called

                // The first correction of the current step is done with the previous step inverse Jacobian approximation
                // mpCurrentJacobianEmulatorPointer->ApplyPrevStepJacobian(mResidualVector_1, InitialCorrection);
                mpCurrentJacobianEmulatorPointer->ApplyPrevStepJacobian(mpResidualVector_1, pInitialCorrection);

                // rIterationGuess -= InitialCorrection; // Recall the minus sign coming from the Taylor expansion of the residual (Newton-Raphson)
                TSpace::UnaliasedAdd(rIterationGuess, -1.0, *pInitialCorrection); // Recall the minus sign coming from the Taylor expansion of the residual (Newton-Raphson)
            }
        }
        else
        {
            //~ std::cout << "Gathering information" << std::endl;

            // Gather the new observation matrices column information
            // VectorType newColV(mProblemSize);
            // VectorType newColW(mProblemSize);

            // for (unsigned int i = 0; i < mProblemSize; i++)
            // {
                // newColV[i] = mResidualVector_1(i) - mResidualVector_0(i);
                // newColW[i] = mIterationValue_1(i) - mIterationValue_0(i);
            // }

            VectorPointerType pNewColV(new VectorType(rResidualVector)); // Hack in case the Epetra copy constructor needs to be called
            VectorPointerType pNewColW(new VectorType(rResidualVector)); // Hack in case the Epetra copy constructor needs to be called

            TSpace::ScaleAndAdd(1.0, *mpResidualVector_1, -1.0, *mpResidualVector_0, *pNewColV);
            TSpace::ScaleAndAdd(1.0, *mpIterationValue_1, -1.0, *mpIterationValue_0, *pNewColW);

            // Observation matrices information filling
            if (mConvergenceAcceleratorIteration <= mProblemSize)
            {
                // Append the new information to the existent observation matrices
                // (mpCurrentJacobianEmulatorPointer)->AppendColToV(newColV);
                // (mpCurrentJacobianEmulatorPointer)->AppendColToW(newColW);
                (mpCurrentJacobianEmulatorPointer)->AppendColToV(*pNewColV);
                (mpCurrentJacobianEmulatorPointer)->AppendColToW(*pNewColW);

                //~ std::cout << "Observation matrices new information appended" << std::endl;
            }
            else
            {
                // (mpCurrentJacobianEmulatorPointer)->DropAndAppendColToV(newColV);
                // (mpCurrentJacobianEmulatorPointer)->DropAndAppendColToW(newColW);
                (mpCurrentJacobianEmulatorPointer)->DropAndAppendColToV(*pNewColV);
                (mpCurrentJacobianEmulatorPointer)->DropAndAppendColToW(*pNewColW);

                //~ std::cout << "Observation matrices size is kept (oldest column is dropped)" << std::endl;
            }

            //~ std::cout << "Jacobian approximation computation starts..." << std::endl;

            // Apply the current step inverse Jacobian emulator to the residual vector
            // VectorType IterationCorrection(mProblemSize);
            // mpCurrentJacobianEmulatorPointer->ApplyJacobian(rResidualVector, IterationCorrection);
            VectorPointerType pIterationCorrection(new VectorType(rResidualVector)); // Hack in case the Epetra copy constructor needs to be called
            mpCurrentJacobianEmulatorPointer->ApplyJacobian(mpResidualVector_1, pIterationCorrection);

            // rIterationGuess -= IterationCorrection; // Recall the minus sign coming from the Taylor expansion of the residual (Newton-Raphson)
            TSpace::UnaliasedAdd(rIterationGuess, -1.0, *pIterationCorrection); // Recall the minus sign coming from the Taylor expansion of the residual (Newton-Raphson)
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
        // TSpace::Copy(mIterationValue_1, mIterationValue_0);
        // TSpace::Copy(mResidualVector_1, mResidualVector_0);
        mpIterationValue_0 = mpIterationValue_1;
        mpResidualVector_0 = mpResidualVector_1;

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
    //
    // VectorType mResidualVector_0;                                               // Previous iteration residual vector
    // VectorType mResidualVector_1;                                               // Current iteration residual vector
    // VectorType mIterationValue_0;                                               // Previous iteration guess
    // VectorType mIterationValue_1;                                               // Current iteration guess

    VectorPointerType mpResidualVector_0;                                       // Previous iteration residual vector pointer
    VectorPointerType mpResidualVector_1;                                       // Current iteration residual vector pointer
    VectorPointerType mpIterationValue_0;                                       // Previous iteration guess pointer
    VectorPointerType mpIterationValue_1;                                       // Current iteration guess pointer

    JacobianEmulatorPointerType mpCurrentJacobianEmulatorPointer;               // Current step Jacobian approximator pointer

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

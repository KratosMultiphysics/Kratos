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

/* Project includes */
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/qr_utility.h"
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
    JacobianEmulator( Pointer&& OldJacobianEmulatorPointer ) {
        mpOldJacobianEmulator = std::unique_ptr<JacobianEmulator<TSpace>>(std::move(OldJacobianEmulatorPointer));
    }

    /**
     * Old Jacobian pointer constructor with recursive previous Jacobian deleting.
     * The inverse Jacobian emulator will use information from the previous Jacobian
     */
    JacobianEmulator( Pointer&& OldJacobianEmulatorPointer, const unsigned int EmulatorBufferSize ) {
        mpOldJacobianEmulator = std::unique_ptr<JacobianEmulator<TSpace> >(std::move(OldJacobianEmulatorPointer));

        // Get the last pointer out of buffer
        if(EmulatorBufferSize > 1) {
            JacobianEmulator* p = (mpOldJacobianEmulator->mpOldJacobianEmulator).get();

            for(unsigned int i = 1; i < (EmulatorBufferSize); i++) {
                if(i == EmulatorBufferSize-1) {
                    (p->mpOldJacobianEmulator).reset();
                } else {
                    p = (p->mpOldJacobianEmulator).get();
                }
            }
        } else { // If Jacobian buffer size equals 1 directly destroy the previous one
            (mpOldJacobianEmulator->mpOldJacobianEmulator).reset();
        }
    }

    /**
     * Empty constructor.
     * The Jacobian emulator will consider minus the identity matrix as previous Jacobian
     */
    JacobianEmulator( ) {}

    /**
     * Copy Constructor.
     */
    JacobianEmulator( const JacobianEmulator& rOther ) {
        mpOldJacobianEmulator = rOther.mpOldJacobianEmulator;
    }

    /**
     * Destructor.
     */
    virtual ~JacobianEmulator() {}

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
    void ApplyPrevStepJacobian(const VectorPointerType pWorkVector,
                               VectorPointerType pProjectedVector) {
        // Security check for the empty observation matrices case (when no correction has been done in the previous step)
        if (mpOldJacobianEmulator->mJacobianObsMatrixV.size() != 0) {
            mpOldJacobianEmulator->ApplyJacobian(pWorkVector, pProjectedVector);
        } else {
            TSpace::Assign(*pProjectedVector, -1.0, *pWorkVector); // Consider minus the identity matrix as inverse Jacobian
        }
    }

    /**
     * Projects the approximated inverse Jacobian onto a vector
     * @param rWorkVector: Vector in where the inverse Jacobian is to be projected
     * @param rProjectedVector: Projected vector output
     */
    void ApplyJacobian(const VectorPointerType pWorkVector,
                       VectorPointerType pProjectedVector) {
        KRATOS_TRY;

        // Security check for the empty observation matrices case (when no correction has been done in the previous step)
        if (mJacobianObsMatrixV.size() == 0) {
            if (mpOldJacobianEmulator != nullptr) { // If it is available, consider the previous step Jacobian
                mpOldJacobianEmulator->ApplyJacobian(pWorkVector, pProjectedVector);
            } else { // When the JacobianEmulator has no PreviousJacobianEmulator consider minus the identity matrix as inverse Jacobian
                TSpace::Assign(*pProjectedVector, -1.0, *pWorkVector);
            }
        } else {
            const unsigned int previous_iterations = mJacobianObsMatrixV.size();
            const unsigned int residual_size = TSpace::Size(mJacobianObsMatrixV[0]);

            VectorPointerType pY(new VectorType(residual_size));
            VectorPointerType pW(new VectorType(residual_size));
            VectorPointerType pzQR(new VectorType(previous_iterations));
            MatrixPointerType pAuxMatQR(new MatrixType(residual_size, previous_iterations));

            // Loop to store a std::vector<VectorType> type as Matrix type
            for (unsigned int i = 0; i < residual_size; ++i) {
                for (unsigned int j = 0; j < previous_iterations; ++j) {
                    (*pAuxMatQR)(i,j) = mJacobianObsMatrixV[j](i);
                }
            }

            VectorPointerType pWorkVectorCopy(new VectorType(*pWorkVector));

            // QR decomposition to compute ((V_k.T*V_k)^-1)*V_k.T*r_k
            mQR_decomposition.compute(residual_size, previous_iterations, &(*pAuxMatQR)(0,0));
            mQR_decomposition.solve(&(*pWorkVectorCopy)(0), &(*pzQR)(0));

            TSpace::SetToZero(*pY);
            for (unsigned int j = 0; j < previous_iterations; ++j) {
                TSpace::UnaliasedAdd(*pY, (*pzQR)(j), mJacobianObsMatrixV[j]);
            }

            TSpace::UnaliasedAdd(*pY, -1.0, *pWorkVector);

            if (mpOldJacobianEmulator == nullptr) {
                TSpace::Copy(*pY, *pProjectedVector); // Consider minus the identity as previous step Jacobian
            } else {
                VectorPointerType pYminus(new VectorType(*pY));
                TSpace::Assign(*pYminus, -1.0, *pY);
                mpOldJacobianEmulator->ApplyJacobian(pYminus, pProjectedVector); // The minus comes from the fact that we want to apply r_k - V_k*zQR
            }

            // w = W_k*z
            TSpace::SetToZero(*pW);
            for (unsigned int j = 0; j < previous_iterations; ++j) {
                TSpace::UnaliasedAdd(*pW, (*pzQR)(j), mJacobianObsMatrixW[j]);
            }

            TSpace::UnaliasedAdd(*pProjectedVector, 1.0, *pW);

        }

        KRATOS_CATCH( "" );

    }

    /**
    * Appends a new column to the observation matrix V
    * @param newColV: new column to be appended
    */
    void AppendColToV(const VectorType& rNewColV) {
        KRATOS_TRY;

        mJacobianObsMatrixV.push_back(rNewColV);

        KRATOS_CATCH( "" );
    }

    /**
    * Appends a new column to the observation matrix W
    * @param newColW: new column to be appended
    */
    void AppendColToW(const VectorType& rNewColW) {
        KRATOS_TRY;

        mJacobianObsMatrixW.push_back(rNewColW);

        KRATOS_CATCH( "" );
    }

    /**
    * Drops the oldest column and appends a new column to the observation matrix V
    * @param newColV: new column to be appended
    */
    void DropAndAppendColToV(const VectorType& rNewColV) {
        KRATOS_TRY;

        // Observation matrices size are close to the interface DOFs number. Old columns are to be dropped.
        for (unsigned int i = 0; i < (TSpace::Size(mJacobianObsMatrixV[0])-1); i++) {
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
    void DropAndAppendColToW(const VectorType& rNewColW) {
        KRATOS_TRY;

        // Observation matrices size are close to the interface DOFs number. Old columns are to be dropped.
        for (unsigned int i = 0; i < (TSpace::Size(mJacobianObsMatrixV[0])-1); i++) {
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
    QR<double, row_major>                 mQR_decomposition;    // QR decomposition object

    Pointer                           mpOldJacobianEmulator;    // Pointer to the old Jacobian

    std::vector<VectorType>             mJacobianObsMatrixV;    // Residual increment observation matrix
    std::vector<VectorType>             mJacobianObsMatrixW;    // Solution increment observation matrix

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
class MVQNRecursiveJacobianConvergenceAccelerator: public ConvergenceAccelerator<TSpace> {
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
    MVQNRecursiveJacobianConvergenceAccelerator( Parameters &rConvAcceleratorParameters ) {
        Parameters mvqn_recursive_default_parameters(R"(
        {
            "solver_type"     : "MVQN_recursive",
            "w_0"             : 0.825,
            "buffer_size"     : 10,
            "cut_off_rel_tol" : 1e-8 
        }
        )");

        rConvAcceleratorParameters.ValidateAndAssignDefaults(mvqn_recursive_default_parameters);

        mOmega_0 = rConvAcceleratorParameters["w_0"].GetDouble();
        mJacobianBufferSize = rConvAcceleratorParameters["buffer_size"].GetInt();
        mColumnCutOffRelTol = rConvAcceleratorParameters["cut_off_rel_tol"].GetDouble();;
        mConvergenceAcceleratorStep = 0;
        mConvergenceAcceleratorIteration = 0;
        mConvergenceAcceleratorFirstCorrectionPerformed = false;
    }

    MVQNRecursiveJacobianConvergenceAccelerator( double OmegaInitial = 0.825, unsigned int JacobianBufferSize = 10, double CutOffRelTol = 1e-8 ) {
        mOmega_0 = OmegaInitial;
        mJacobianBufferSize = JacobianBufferSize;
        mColumnCutOffRelTol = CutOffRelTol;
        mConvergenceAcceleratorStep = 0;
        mConvergenceAcceleratorIteration = 0;
        mConvergenceAcceleratorFirstCorrectionPerformed = false;
    }

    /**
     * Copy Constructor.
     */
    MVQNRecursiveJacobianConvergenceAccelerator( const MVQNRecursiveJacobianConvergenceAccelerator& rOther ) {
        mColumnCutOffRelTol = rOther.mColumnCutOffRelTol;
        mOmega_0 = rOther.mOmega_0;
        mJacobianBufferSize = rOther.mJacobianBufferSize;
        mConvergenceAcceleratorStep = 0;
        mConvergenceAcceleratorIteration = 0;
        mConvergenceAcceleratorFirstCorrectionPerformed = false;
    }

    /**
     * Destructor.
     */
    virtual ~MVQNRecursiveJacobianConvergenceAccelerator() {}

    ///@}

    ///@name Operators
    ///@{
    ///@}

    ///@name Operations
    ///@{

    //~ /**
     //~ * Construct the initial inverse Jacobian emulator
     //~ */
    void Initialize() override {
        KRATOS_TRY;

        mpCurrentJacobianEmulatorPointer = std::unique_ptr< JacobianEmulator <TSpace> > (new JacobianEmulator<TSpace>());

        KRATOS_CATCH( "" );
    }


    /**
     * Initialize the internal iteration counter
     */
    void InitializeSolutionStep() override {
        KRATOS_TRY;

        mConvergenceAcceleratorStep += 1;
        mConvergenceAcceleratorIteration = 0;

        if (mConvergenceAcceleratorStep <= mJacobianBufferSize) {
            // Construct the inverse Jacobian emulator
            mpCurrentJacobianEmulatorPointer = std::unique_ptr< JacobianEmulator<TSpace> > (new JacobianEmulator<TSpace>(std::move(mpCurrentJacobianEmulatorPointer)));
        } else {
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
                        VectorType& rIterationGuess) override {
        KRATOS_TRY;

        mProblemSize = TSpace::Size(rResidualVector);

        VectorPointerType pAuxResidualVector(new VectorType(rResidualVector));
        VectorPointerType pAuxIterationGuess(new VectorType(rIterationGuess));
        std::swap(mpResidualVector_1, pAuxResidualVector);
        std::swap(mpIterationValue_1, pAuxIterationGuess);

        if (mConvergenceAcceleratorIteration == 0) {
            if (mConvergenceAcceleratorFirstCorrectionPerformed == false) {
                // The very first correction of the problem is done with a fixed point iteration
                TSpace::UnaliasedAdd(rIterationGuess, mOmega_0, *mpResidualVector_1);

                mConvergenceAcceleratorFirstCorrectionPerformed = true;
            } else {
                VectorPointerType pInitialCorrection(new VectorType(rResidualVector));

                // The first correction of the current step is done with the previous step inverse Jacobian approximation
                mpCurrentJacobianEmulatorPointer->ApplyPrevStepJacobian(mpResidualVector_1, pInitialCorrection);

                TSpace::UnaliasedAdd(rIterationGuess, -1.0, *pInitialCorrection); // Recall the minus sign coming from the Taylor expansion of the residual (Newton-Raphson)
            }
        } else {
            // Gather the new observation matrices column information
            VectorPointerType pNewColV(new VectorType(*mpResidualVector_1));
            VectorPointerType pNewColW(new VectorType(*mpIterationValue_1));

            TSpace::UnaliasedAdd(*pNewColV, -1.0, *mpResidualVector_0); // NewColV = ResidualVector_1 - ResidualVector_0
            TSpace::UnaliasedAdd(*pNewColW, -1.0, *mpIterationValue_0); // NewColW = IterationValue_1 - IterationValue_0

            const double new_col_v_norm = TSpace::TwoNorm(*pNewColV);
            const double new_col_w_norm = TSpace::TwoNorm(*pNewColW);

            // Observation matrices information filling
            if (mConvergenceAcceleratorIteration <= mProblemSize) {
                if (mConvergenceAcceleratorIteration == 1) {
                    // For the 1st iteration, always append the new information to the existent observation matrices
                    (mpCurrentJacobianEmulatorPointer)->AppendColToV(*pNewColV);
                    (mpCurrentJacobianEmulatorPointer)->AppendColToW(*pNewColW);

                    // Set the 1st column as maximum norm value
                    mObsMatrixVMaxNorm = new_col_v_norm;
                    mObsMatrixWMaxNorm = new_col_w_norm;
                } else {
                    // Append the new information to the existent observation matrices acording to the cut off criterion
                    if ((new_col_v_norm > mColumnCutOffRelTol) && (new_col_w_norm > mColumnCutOffRelTol)) {
                        (mpCurrentJacobianEmulatorPointer)->AppendColToV(*pNewColV);
                        (mpCurrentJacobianEmulatorPointer)->AppendColToW(*pNewColW);
                    } else {
                        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
                        std::cout << "WARNING: Current iteration info has not been appended to observation matrices!" << std::endl;
                        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
                    }

                    // Check if the new column norms are larger than the existent ones
                    mObsMatrixVMaxNorm = (mObsMatrixVMaxNorm < new_col_v_norm) ? new_col_v_norm : mObsMatrixVMaxNorm;
                    mObsMatrixWMaxNorm = (mObsMatrixWMaxNorm < new_col_w_norm) ? new_col_w_norm : mObsMatrixWMaxNorm;
                }
            } else {
                // Append the new information to the existent observation matrices acording to the cut off criterion
                if ((new_col_v_norm  > mColumnCutOffRelTol) && (new_col_w_norm > mColumnCutOffRelTol)) {
                    (mpCurrentJacobianEmulatorPointer)->DropAndAppendColToV(*pNewColV);
                    (mpCurrentJacobianEmulatorPointer)->DropAndAppendColToW(*pNewColW);
                } else {
                    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
                    std::cout << "WARNING: Current iteration info has not been appended to observation matrices!" << std::endl;
                    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
                }

                // Check if the new column norms are larger than the existent ones
                mObsMatrixVMaxNorm = (mObsMatrixVMaxNorm < new_col_v_norm) ? new_col_v_norm : mObsMatrixVMaxNorm;
                mObsMatrixWMaxNorm = (mObsMatrixWMaxNorm < new_col_w_norm) ? new_col_w_norm : mObsMatrixWMaxNorm;
            }

            // Apply the current step inverse Jacobian emulator to the residual vector
            VectorPointerType pIterationCorrection(new VectorType(rResidualVector));
            mpCurrentJacobianEmulatorPointer->ApplyJacobian(mpResidualVector_1, pIterationCorrection);
            
            TSpace::UnaliasedAdd(rIterationGuess, -1.0, *pIterationCorrection); // Recall the minus sign coming from the Taylor expansion of the residual (Newton-Raphson)
        }

        KRATOS_CATCH( "" );
    }

    /**
     * Updates the MVQN iteration values for the next non-linear iteration
     */
    void FinalizeNonLinearIteration() override {
        KRATOS_TRY;

        // Variables update
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

    double mColumnCutOffRelTol;     // Relative tolerance for the observation cut off
    double mObsMatrixVMaxNorm;      // Observation matrix V maximum column norm
    double mObsMatrixWMaxNorm;      // Observation matrix W maximum column norm

    double mOmega_0;                                                    // Relaxation factor for the initial fixed point iteration
    unsigned int mProblemSize;                                          // Residual to minimize size
    unsigned int mJacobianBufferSize;                                   // User-defined Jacobian buffer-size
    unsigned int mCurrentJacobianBufferSize;                            // Current Jacobian buffer-size (expected to be less or equal to the user-defined one)
    unsigned int mConvergenceAcceleratorStep;                           // Convergence accelerator steps counter
    unsigned int mConvergenceAcceleratorIteration;                      // Convergence accelerator iteration counter
    bool mConvergenceAcceleratorFirstCorrectionPerformed;               // Indicates that the initial fixed point iteration has been already performed

    VectorPointerType mpResidualVector_0;                               // Previous iteration residual vector pointer
    VectorPointerType mpResidualVector_1;                               // Current iteration residual vector pointer
    VectorPointerType mpIterationValue_0;                               // Previous iteration guess pointer
    VectorPointerType mpIterationValue_1;                               // Current iteration guess pointer

    JacobianEmulatorPointerType mpCurrentJacobianEmulatorPointer;       // Current step Jacobian approximator pointer

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

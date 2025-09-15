//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Original author:  Ruben Zorrilla
//

#pragma once

/* System includes */

/* External includes */
#include "Epetra_SerialDenseSolver.h"

/* Project includes */
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"

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
class TrilinosJacobianEmulator
{
public:

    ///@name Type Definitions
    ///@{
    typedef typename std::unique_ptr< TrilinosJacobianEmulator<TSpace> >         Pointer;

    typedef typename TSpace::VectorType                                       VectorType;
    typedef typename TSpace::VectorPointerType                         VectorPointerType;

    typedef typename TSpace::MatrixType                                       MatrixType;
    typedef typename TSpace::MatrixPointerType                         MatrixPointerType;


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
    TrilinosJacobianEmulator( ModelPart& rInterfaceModelPart,
                              const Epetra_MpiComm& rEpetraCommunicator,
                              Pointer&& OldJacobianEmulatorPointer ) :
        mrInterfaceModelPart(rInterfaceModelPart),
        mrEpetraCommunicator(rEpetraCommunicator)
    {
        mpOldJacobianEmulator = std::unique_ptr<TrilinosJacobianEmulator<TSpace> >(std::move(OldJacobianEmulatorPointer));
    }

    /**
     * Old Jacobian pointer constructor with recursive previous Jacobian deleting.
     * The inverse Jacobian emulator will use information from the previous Jacobian
     */
    TrilinosJacobianEmulator( ModelPart &rInterfaceModelPart,
                              const Epetra_MpiComm &rEpetraCommunicator,
                              Pointer &&OldJacobianEmulatorPointer,
                              const unsigned int EmulatorBufferSize) :
        mrInterfaceModelPart(rInterfaceModelPart),
        mrEpetraCommunicator(rEpetraCommunicator)
    {
        mpOldJacobianEmulator = std::unique_ptr<TrilinosJacobianEmulator<TSpace> >(std::move(OldJacobianEmulatorPointer));

        // Get the last pointer out of buffer
        if(EmulatorBufferSize > 1)
        {
            TrilinosJacobianEmulator* p = (mpOldJacobianEmulator->mpOldJacobianEmulator).get();

            for(unsigned int i = 1; i < (EmulatorBufferSize); i++)
            {
                if(i == EmulatorBufferSize-1)
                {
                    (p->mpOldJacobianEmulator).reset();
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
        }
    }

    /**
     * Empty constructor.
     * The Jacobian emulator will consider minus the identity matrix as previous Jacobian
     */
    TrilinosJacobianEmulator( ModelPart &rInterfaceModelPart,
                              const Epetra_MpiComm &rEpetraCommunicator) :
        mrInterfaceModelPart(rInterfaceModelPart),
        mrEpetraCommunicator(rEpetraCommunicator) {}

    /**
     * Copy Constructor.
     */
    TrilinosJacobianEmulator( const TrilinosJacobianEmulator& rOther )
        : mrInterfaceModelPart(rOther.mrInterfaceModelPart),
          mrEpetraCommunicator(rOther.mrEpetraCommunicator),
          mpOldJacobianEmulator(rOther.mpOldJacobianEmulator)
    {
    }

    /**
     * Destructor.
     */
    virtual ~TrilinosJacobianEmulator
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
    void ApplyPrevStepJacobian(const VectorPointerType pWorkVector,
                               VectorPointerType pProjectedVector)
    {
        // Security check for the empty observation matrices case (when no correction has been done in the previous step)
        if (mpOldJacobianEmulator->mJacobianObsMatrixV.size() != 0)
        {
            mpOldJacobianEmulator->ApplyJacobian(pWorkVector, pProjectedVector);
        }
        else
        {
            TSpace::Assign(*pProjectedVector, -1.0, *pWorkVector); // Consider minus the identity matrix as inverse Jacobian
        }
    }

    /**
     * Projects the approximated inverse Jacobian onto a vector
     * @param rWorkVector: Vector in where the inverse Jacobian is to be projected
     * @param rProjectedVector: Projected vector output
     */
    void ApplyJacobian(const VectorPointerType pWorkVector,
                       VectorPointerType pProjectedVector)
    {
        KRATOS_TRY;

        const unsigned int previous_iterations = mJacobianObsMatrixV.size();

        // Security check for the empty observation matrices case (when no correction has been done in the previous step)
        if (previous_iterations == 0)
        {
            if (mpOldJacobianEmulator != nullptr) // If it is available, consider the previous step Jacobian
            {
                mpOldJacobianEmulator->ApplyJacobian(pWorkVector, pProjectedVector);
            }
            else // When the JacobianEmulator has no PreviousJacobianEmulator consider minus the identity matrix as inverse Jacobian
            {
                TSpace::Assign(*pProjectedVector, -1.0, *pWorkVector);
            }
        }
        else
        {
            // Get the domain size
            unsigned int n_dim = mrInterfaceModelPart.GetProcessInfo()[DOMAIN_SIZE];

            // Construct the interface map
            int NumLocalInterfaceDofs = mrInterfaceModelPart.GetCommunicator().LocalMesh().NumberOfNodes() * n_dim;
            int NumGlobalInterfaceDofs = mrInterfaceModelPart.GetCommunicator().GetDataCommunicator().SumAll(NumLocalInterfaceDofs);
            int IndexBase = 0; // 0 for C-style vectors, 1 for Fortran numbering
            Epetra_Map InterfaceMap(NumGlobalInterfaceDofs, NumLocalInterfaceDofs, IndexBase, mrEpetraCommunicator);

            // Create new vector using given map
            auto pY(new Epetra_FEVector(InterfaceMap));
            auto pW(new Epetra_FEVector(InterfaceMap));

            // Loop to store a std::vector<VectorType> type as Matrix type
            Epetra_SerialDenseMatrix Vtrans_V(previous_iterations, previous_iterations);

            for (unsigned int i=0; i<previous_iterations; ++i)
            {
                Vtrans_V(i,i) = TSpace::Dot(mJacobianObsMatrixV[i], mJacobianObsMatrixV[i]);

                for (unsigned int j=i+1; j<previous_iterations; ++j)
                {
                    Vtrans_V(i,j) = TSpace::Dot(mJacobianObsMatrixV[i], mJacobianObsMatrixV[j]);
                    Vtrans_V(j,i) = Vtrans_V(i,j);
                }
            }

            Epetra_SerialDenseVector Vtrans_r(previous_iterations);
            Epetra_SerialDenseVector zSystemSol(previous_iterations);

            for (unsigned int i=0; i<previous_iterations; ++i)
            {
                Vtrans_r(i) = TSpace::Dot(mJacobianObsMatrixV[i], (*pWorkVector));
            }

            Epetra_SerialDenseSolver EpetraSystemSolver;
            EpetraSystemSolver.SetMatrix(Vtrans_V);
            EpetraSystemSolver.SetVectors(zSystemSol,Vtrans_r);
            EpetraSystemSolver.Solve();

            TSpace::SetToZero(*pY);
            for (unsigned int j = 0; j < previous_iterations; ++j)
            {
                TSpace::UnaliasedAdd(*pY, zSystemSol(j), mJacobianObsMatrixV[j]);
            }

            TSpace::UnaliasedAdd(*pY, -1.0, *pWorkVector);

            if (mpOldJacobianEmulator == nullptr)
            {
                TSpace::Copy(*pY, *pProjectedVector); // Consider minus the identity as previous step Jacobian
            }
            else
            {
                VectorPointerType pYminus(new VectorType(*pY));
                TSpace::Assign(*pYminus, -1.0, *pY);
                mpOldJacobianEmulator->ApplyJacobian(pYminus, pProjectedVector); // The minus comes from the fact that we want to apply r_k - V_k*zQR
            }

            // w = W_k*z
            TSpace::SetToZero(*pW);
            for (unsigned int j = 0; j < previous_iterations; ++j)
            {
                TSpace::UnaliasedAdd(*pW, zSystemSol(j), mJacobianObsMatrixW[j]);
            }

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
    ModelPart&                   mrInterfaceModelPart;        // Interface model part
    const Epetra_MpiComm&              mrEpetraCommunicator;        // Epetra communicator

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

}; /* Class TrilinosJacobianEmulator */


/** @brief MVQN (MultiVectorQuasiNewton method) acceleration scheme
 * Recursive MultiVectorQuasiNewton convergence accelerator. This convergence accelerator
 * is an alternative implementation of the standard MVQN that avoids the storage of the
 * @tparam TSparseSpace Linear algebra sparse space
 * @tparam TDenseSpace Linear algebra dense space
 */
template<class TSparseSpace, class TDenseSpace>
class TrilinosMVQNRecursiveJacobianConvergenceAccelerator: public ConvergenceAccelerator<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( TrilinosMVQNRecursiveJacobianConvergenceAccelerator );

    typedef ConvergenceAccelerator<TSparseSpace,TDenseSpace>                           BaseType;
    typedef typename BaseType::Pointer                                          BaseTypePointer;

    typedef typename TrilinosJacobianEmulator<TSparseSpace>::Pointer JacobianEmulatorPointerType;

    typedef typename TSparseSpace::VectorType                                        VectorType;
    typedef typename TSparseSpace::VectorPointerType                          VectorPointerType;

    typedef typename BaseType::MatrixType                                            MatrixType;
    typedef typename BaseType::MatrixPointerType                              MatrixPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     * MVQN convergence accelerator
     */
    TrilinosMVQNRecursiveJacobianConvergenceAccelerator(
        ModelPart& rInterfaceModelPart,
        const Epetra_MpiComm& rEpetraCommunicator,
        Parameters rConvAcceleratorParameters)
        : mrInterfaceModelPart(rInterfaceModelPart)
        , mrEpetraCommunicator(rEpetraCommunicator)
    {
        Parameters mvqn_recursive_default_parameters(R"(
        {
            "solver_type" : "MVQN_recursive",
            "w_0"         : 0.825,
            "buffer_size" : 10,
            "interface_block_newton" : false
        }
        )");

        rConvAcceleratorParameters.ValidateAndAssignDefaults(mvqn_recursive_default_parameters);

        mProblemSize = 0;
        mOmega_0 = rConvAcceleratorParameters["w_0"].GetDouble();
        mJacobianBufferSize = rConvAcceleratorParameters["buffer_size"].GetInt();
        mConvergenceAcceleratorStep = 0;
        mConvergenceAcceleratorIteration = 0;
        mConvergenceAcceleratorFirstCorrectionPerformed = false;
    }

    TrilinosMVQNRecursiveJacobianConvergenceAccelerator( ModelPart& rInterfaceModelPart,
                                                         const Epetra_MpiComm& rEpetraCommunicator,
                                                         const double OmegaInitial = 0.825,
                                                         const unsigned int JacobianBufferSize = 7 ):
        mrInterfaceModelPart(rInterfaceModelPart),
        mrEpetraCommunicator(rEpetraCommunicator),
        mOmega_0(OmegaInitial),
        mJacobianBufferSize(JacobianBufferSize)
    {
        mProblemSize = 0;
        mConvergenceAcceleratorStep = 0;
        mConvergenceAcceleratorIteration = 0;
        mConvergenceAcceleratorFirstCorrectionPerformed = false;
    }

    /**
     * Copy Constructor.
     */
    TrilinosMVQNRecursiveJacobianConvergenceAccelerator( const TrilinosMVQNRecursiveJacobianConvergenceAccelerator& rOther )
        : mrInterfaceModelPart(rOther.mrBInterfaceModelPart),
          mrEpetraCommunicator(rOther.mrEpetraCommunicator),
          mOmega_0(rOther.mOmega_0),
          mJacobianBufferSize(rOther.mJacobianBufferSize)
    {
        mConvergenceAcceleratorStep = 0;
        mConvergenceAcceleratorIteration = 0;
        mConvergenceAcceleratorFirstCorrectionPerformed = false;
    }

    /**
     * Destructor.
     */
    virtual ~TrilinosMVQNRecursiveJacobianConvergenceAccelerator
    () {}

    ///@}

    ///@name Operators
    ///@{
    ///@}

    ///@name Operations
    ///@{

    /**
     * Construct the initial inverse Jacobian emulator
     */
    void Initialize() override
    {
        KRATOS_TRY;

        mpCurrentJacobianEmulatorPointer = Kratos::make_unique<TrilinosJacobianEmulator<TSparseSpace>>(
            mrInterfaceModelPart,
            mrEpetraCommunicator);

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
            mpCurrentJacobianEmulatorPointer = Kratos::make_unique< TrilinosJacobianEmulator<TSparseSpace>>(
                mrInterfaceModelPart,
                mrEpetraCommunicator,
                std::move(mpCurrentJacobianEmulatorPointer));
        }
        else
        {
            // Construct the inverse Jacobian emulator considering the recursive elimination
            mpCurrentJacobianEmulatorPointer = Kratos::make_unique<TrilinosJacobianEmulator<TSparseSpace>>(
                mrInterfaceModelPart,
                mrEpetraCommunicator,
                std::move(mpCurrentJacobianEmulatorPointer),
                mJacobianBufferSize);
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

        mProblemSize = TSparseSpace::Size(rResidualVector);

        VectorPointerType pAuxResidualVector(new VectorType(rResidualVector));
        VectorPointerType pAuxIterationGuess(new VectorType(rIterationGuess));
        std::swap(mpResidualVector_1, pAuxResidualVector);
        std::swap(mpIterationValue_1, pAuxIterationGuess);

        if (mConvergenceAcceleratorIteration == 0)
        {
            if (mConvergenceAcceleratorFirstCorrectionPerformed == false)
            {
                // The very first correction of the problem is done with a fixed point iteration
                TSparseSpace::UnaliasedAdd(rIterationGuess, mOmega_0, *mpResidualVector_1);

                mConvergenceAcceleratorFirstCorrectionPerformed = true;
            }
            else
            {
                VectorPointerType pInitialCorrection(new VectorType(rResidualVector));

                // The first correction of the current step is done with the previous step inverse Jacobian approximation
                mpCurrentJacobianEmulatorPointer->ApplyPrevStepJacobian(mpResidualVector_1, pInitialCorrection);

                TSparseSpace::UnaliasedAdd(rIterationGuess, -1.0, *pInitialCorrection); // Recall the minus sign coming from the Taylor expansion of the residual (Newton-Raphson)
            }
        }
        else
        {
            // Gather the new observation matrices column information
            VectorPointerType pNewColV(new VectorType(*mpResidualVector_1));
            VectorPointerType pNewColW(new VectorType(*mpIterationValue_1));

            TSparseSpace::UnaliasedAdd(*pNewColV, -1.0, *mpResidualVector_0); // NewColV = ResidualVector_1 - ResidualVector_0
            TSparseSpace::UnaliasedAdd(*pNewColW, -1.0, *mpIterationValue_0); // NewColW = IterationValue_1 - IterationValue_0

            // Observation matrices information filling
            if (mConvergenceAcceleratorIteration <= mProblemSize)
            {
                // Append the new information to the existent observation matrices
                (mpCurrentJacobianEmulatorPointer)->AppendColToV(*pNewColV);
                (mpCurrentJacobianEmulatorPointer)->AppendColToW(*pNewColW);
            }
            else
            {
                (mpCurrentJacobianEmulatorPointer)->DropAndAppendColToV(*pNewColV);
                (mpCurrentJacobianEmulatorPointer)->DropAndAppendColToW(*pNewColW);
            }

            // Apply the current step inverse Jacobian emulator to the residual vector
            VectorPointerType pIterationCorrection(new VectorType(rResidualVector));
            mpCurrentJacobianEmulatorPointer->ApplyJacobian(mpResidualVector_1, pIterationCorrection);

            TSparseSpace::UnaliasedAdd(rIterationGuess, -1.0, *pIterationCorrection); // Recall the minus sign coming from the Taylor expansion of the residual (Newton-Raphson)
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
    ModelPart& mrInterfaceModelPart;                                    // Interface model part
    const Epetra_MpiComm& mrEpetraCommunicator;                         // Epetra communicator
    double mOmega_0;                                                    // Relaxation factor for the initial fixed point iteration
    unsigned int mJacobianBufferSize;                                   // User-defined Jacobian buffer-size

    unsigned int mProblemSize;                                          // Residual to minimize size
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

}; /* Class TrilinosMVQNRecursiveJacobianConvergenceAccelerator */

///@}


///@}

///@name Type Definitions
///@{
///@}

///@name Input and output
///@{
///@}

}  /* namespace Kratos.*/

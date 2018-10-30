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
#include "convergence_accelerator.hpp"
#include "includes/ublas_interface.h"
#include "utilities/svd_utils.h"
#include "utilities/qr_utility.h"


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
        mpOldJacobianEmulator = std::unique_ptr<JacobianEmulator<TSpace>>(std::move(OldJacobianEmulatorPointer));
    }

    /**
     * Old Jacobian pointer constructor with recursive previous Jacobian deleting.
     * The inverse Jacobian emulator will use information from the previous Jacobian
     */
    JacobianEmulator( Pointer&& OldJacobianEmulatorPointer, const unsigned int EmulatorBufferSize )
    {
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
    JacobianEmulator( const JacobianEmulator& rOther )
    {
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
     * Projects the approximated inverse Jacobian onto a vector
     * @param rWorkVector: Vector in where the inverse Jacobian is to be projected
     * @param rProjectedVector: Projected vector output
     */
    void ApplyJacobian(
        const VectorPointerType pWorkVector,
        VectorPointerType pProjectedVector)
    {
        const auto data_cols = this->GetNumberOfDataCols();

        // Security check for the empty observation matrices case (when no correction has been done in the previous step)
        if (data_cols == 0) {
            if (mpOldJacobianEmulator != nullptr) { // If it is available, consider the previous step Jacobian
                mpOldJacobianEmulator->ApplyJacobian(pWorkVector, pProjectedVector);
            } else { // When the JacobianEmulator has no PreviousJacobianEmulator consider minus the identity matrix as inverse Jacobian
                TSpace::Assign(*pProjectedVector, -1.0, *pWorkVector);
            }
        } else {
            const unsigned int residual_size = this->GetResidualSize();

            // Set V and trans(V)
            MatrixPointerType pV = Kratos::make_shared<MatrixType>(residual_size, data_cols);
            MatrixPointerType pVtrans = Kratos::make_shared<MatrixType>(data_cols,residual_size);
            for (std::size_t i = 0; i < data_cols; ++i){
                const auto i_row = mJacobianObsMatrixV[i];
                for (std::size_t j = 0; j < residual_size; ++j){
                    (*pV)(j, i) = i_row(j);
                    (*pVtrans)(i, j) = i_row(j);
                }
            }

            // Set trans(V)*V
            MatrixPointerType pVtransV = Kratos::make_shared<MatrixType>(data_cols, data_cols);
            for (std::size_t i_row = 0; i_row < data_cols; ++i_row){
                for (std::size_t i_col = 0; i_col < data_cols; ++i_col){
                    (*pVtransV)(i_row, i_col) = TSpace::Dot(mJacobianObsMatrixV[i_row],mJacobianObsMatrixV[i_col]);
                }
            }

            // Compute trans(V)*pWorkVector
            VectorPointerType pVtransWorkVect = Kratos::make_shared<VectorType>(data_cols);
            TSpace::Mult(*pVtrans,*pWorkVector,*pVtransWorkVect);

            // Do the QR decomp of trans(V)*V and solve ((trans(V)*V)^-1)*trans(V)*res
            QR<double, row_major> QRUtil;
            QRUtil.compute(data_cols, data_cols, &(*pVtransV)(0,0));
            VectorPointerType pLambda = Kratos::make_shared<VectorType>(data_cols);
            QRUtil.solve(&(*pVtransWorkVect)(0), &(*pLambda)(0));

            // Compute (res - V*Lambda). Note that it is save in pY
            VectorPointerType pY = Kratos::make_shared<VectorType>(residual_size);
            TSpace::Mult(*pV,*pLambda,*pY);
            TSpace::UnaliasedAdd(*pY, -1.0, *pWorkVector);

            // Project over the previous step Jacobian
            if (mpOldJacobianEmulator == nullptr) {
                TSpace::Copy(*pY, *pProjectedVector); // Consider minus the identity as previous step Jacobian
            } else {
                VectorPointerType pYminus(new VectorType(*pY));
                TSpace::Assign(*pYminus, -1.0, *pY);
                mpOldJacobianEmulator->ApplyJacobian(pYminus, pProjectedVector); // The minus comes from the fact that we want to apply r_k - V_k*zQR
            }

            // w = W_k*z
            VectorPointerType pW(new VectorType(residual_size));
            TSpace::SetToZero(*pW);
            for (unsigned int j = 0; j < data_cols; ++j) {
                // TSpace::UnaliasedAdd(*pW, (*pzQR)(j), mJacobianObsMatrixW[j]);
                TSpace::UnaliasedAdd(*pW, (*pLambda)(j), mJacobianObsMatrixW[j]);
            }

            TSpace::UnaliasedAdd(*pProjectedVector, 1.0, *pW);
        }
    }

    /**
    * Appends a two new columns to the observation matrices V and W
    * Then, it checks if the new information columns are linear dependent
    * to the existent ones by computing a QR decomposition and checking
    * the diagonal coefficients of matrix R. If any value is less than
    * the stablished threshold criterion, it is assume that the related
    * column is close to be linear deppendent and is dropped.
    * @param rNewColV new column to be appended to V observation matrix
    * @param rNewColW new column to be appended to W observation matrix
    * @param AbsCutOffEps epsilon value for the absolut cut-off criteria
    * @return Returns true if the columns have been added
    */
    bool AppendDataColumns(
        const VectorType& rNewColV,
        const VectorType& rNewColW,
        const double AbsCutOff = 1e-8)
    {   
        // Add the provided iformation to the observation matrices
        mJacobianObsMatrixV.push_back(rNewColV);
        mJacobianObsMatrixW.push_back(rNewColW);

        // Loop to store a std::vector<VectorType> type as Matrix type
        const std::size_t data_cols = this->GetNumberOfDataCols();

        // Set trans(V)*V
        MatrixPointerType ptransV_V = Kratos::make_shared<MatrixType>(data_cols, data_cols);
        for (std::size_t i_row = 0; i_row < data_cols; ++i_row){
            for (std::size_t i_col = 0; i_col < data_cols; ++i_col){
                (*ptransV_V)(i_row, i_col) = TSpace::Dot(mJacobianObsMatrixV[i_row],mJacobianObsMatrixV[i_col]);
            }
        }

        // Perform the Singular Value Decomposition (SVD) of matrix trans(V)*V
        // SVD decomposition of matrix A yields two orthogonal matrices "u_svd" 
        // and "v_svd" as well as a diagonal matrix "w_svd" cotaining matrix
        // A eigenvalues such that A = u_svd * w_svd * v_svd
        MatrixType u_svd; // Orthogonal matrix (m x m)
        MatrixType w_svd; // Rectangular diagonal matrix (m x n)
        MatrixType v_svd; // Orthogonal matrix (n x n)
        std::string svd_type = "Jacobi"; // SVD decomposition type
        double svd_rel_tol = 1.0e-6; // Relative tolerance of the SVD decomposition (it will be multiplied by the input matrix norm)
        const auto svd_its = SVDUtils<double>::SingularValueDecomposition(*ptransV_V, u_svd, w_svd, v_svd, svd_type, svd_rel_tol);

        // Get the eigenvalues vector. Remember that eigenvalues 
        // of trans(A)*A are equal to the eigenvalues of A^2
        double eig_sum = 0.0;
        std::vector<double> eig_vector(data_cols);
        for (std::size_t i_col = 0; i_col < data_cols; ++i_col){
            const double i_col_eig = std::sqrt(w_svd(i_col,i_col));
            eig_sum += i_col_eig;
            eig_vector[i_col] = i_col_eig;
        }

        // Reorder the eigenvalues vector from max to min and compute the eigenvalues norm
        std::vector<double> eig_vector_ordered(data_cols);
        for (std::size_t i_col = 0; i_col < data_cols; ++i_col){
            double max_eig = 0.0;
            std::size_t max_index;
            for (std::size_t j_col = 0; j_col < eig_vector.size(); ++j_col){
                if (max_eig < eig_vector[j_col]){
                    max_index = j_col;
                    max_eig = eig_vector[j_col];
                }
            }
            eig_vector_ordered[i_col] = max_eig;
            eig_vector.erase(eig_vector.begin() + max_index);
        }

        // Check the representativity of each eigenvalue and its value.
        // If its value is close to zero or out of the representativity
        // the correspondent data columns are dropped from both V and W.
        const double max_eig_V = eig_vector_ordered[0];
        const double min_eig_V = eig_vector_ordered[data_cols - 1];
        if (min_eig_V < AbsCutOff * max_eig_V){
            KRATOS_WARNING("MVQNRecursiveJacobianConvergenceAccelerator") 
                << "Dropping info for eigenvalue " << eig_vector_ordered[data_cols - 1] << " (tolerance " << AbsCutOff * max_eig_V << " )" << std::endl;
            mJacobianObsMatrixV.pop_back();
            mJacobianObsMatrixW.pop_back();
            return false;
        }
        // }

        return true;
    }

    /**
    * Calls the AppendDataColumns() method to add the new data columns
    * to the observation matrices (provided that the new data columns 
    * are not linear dependent to the previous data). Then, if the
    * information has been added, the oldest column is dropped to avoid 
    * the number of data columns become larger than the problem size.
    * @param rNewColV new column to be appended to V observation matrix
    * @param rNewColW new column to be appended to W observation matrix
    * @param CutOffEps epsilon value for the cut-off criteria
    */
    bool DropAndAppendDataColumns(
        const VectorType& rNewColV,
        const VectorType& rNewColW,
        const double AbsCutOffEps = 1e-8)
    {  
        // std::cout << "DropAndAppendDataColumns()" << std::endl;
        const bool info_added = this->AppendDataColumns(rNewColV, rNewColW, AbsCutOffEps);

        // If a new column has been added, drop the oldest data
        if (info_added){
            // Move the current data (oldest data is dropped)
            for (unsigned int i = 0; i < (this->GetResidualSize() - 1); ++i) {
                mJacobianObsMatrixV[i] = mJacobianObsMatrixV[i+1];
                mJacobianObsMatrixW[i] = mJacobianObsMatrixW[i+1];
            }
            // Pop the last term (keep the amount of data columns as the problem size)
            mJacobianObsMatrixV.pop_back();
            mJacobianObsMatrixW.pop_back();
        }

        return info_added;
    }

    /**
     * @brief Get the Number Of Data Cols object
     * This function returns the number of data columns stored.
     * Since the data columns is assumed (and must be) the same in 
     * both observation matrices, it is computed using V matrix.
     * @return std::size_t number of data columns
     */
    inline std::size_t GetNumberOfDataCols() const
    {
        return mJacobianObsMatrixV.size();
    }

    /**
     * @brief Get the Residual Size object
     * This function returns the interface residual size.
     * Since the residual size is assumed (and must be) the same in 
     * every column for both observation matrices, it is computed 
     * using the first column of V matrix.
     * @return std::size_t residual size
     */
    inline std::size_t GetResidualSize() const
    {
        return TSpace::Size(mJacobianObsMatrixV[0]);
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

    Pointer                     mpOldJacobianEmulator;  // Pointer to the old Jacobian

    std::vector<VectorType>       mJacobianObsMatrixV;  // Residual increment observation matrix
    std::vector<VectorType>       mJacobianObsMatrixW;  // Solution increment observation matrix

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
    MVQNRecursiveJacobianConvergenceAccelerator( Parameters &rConvAcceleratorParameters )
    {
        Parameters mvqn_recursive_default_parameters(R"(
        {
            "solver_type"     : "MVQN_recursive",
            "w_0"             : 0.825,
            "buffer_size"     : 10,
            "rel_cut_off_tol" : 1e-2,
            "abs_cut_off_tol" : 1e-8 
        }
        )");

        rConvAcceleratorParameters.ValidateAndAssignDefaults(mvqn_recursive_default_parameters);

        mOmega_0 = rConvAcceleratorParameters["w_0"].GetDouble();
        mAbsCutOff = rConvAcceleratorParameters["abs_cut_off_tol"].GetDouble();
        mJacobianBufferSize = rConvAcceleratorParameters["buffer_size"].GetInt();
        mConvergenceAcceleratorStep = 0;
        mConvergenceAcceleratorIteration = 0;
        mConvergenceAcceleratorFirstCorrectionPerformed = false;
    }

    MVQNRecursiveJacobianConvergenceAccelerator(
        double OmegaInitial = 0.825,
        unsigned int JacobianBufferSize = 10,
        double AbsCutOff = 1e-8)
    {
        mOmega_0 = OmegaInitial;
        mAbsCutOff = AbsCutOff;
        mJacobianBufferSize = JacobianBufferSize;
        mConvergenceAcceleratorStep = 0;
        mConvergenceAcceleratorIteration = 0;
        mConvergenceAcceleratorFirstCorrectionPerformed = false;
    }

    /**
     * Copy Constructor.
     */
    MVQNRecursiveJacobianConvergenceAccelerator( const MVQNRecursiveJacobianConvergenceAccelerator& rOther ) = delete;

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

    /**
     * @brief Initialize the Jacobian emulator
     * This method constructs the very first Jacobian emulator of the simulation
     */
    void Initialize() override
    {
        KRATOS_TRY;

        mpCurrentJacobianEmulatorPointer = std::unique_ptr< JacobianEmulator <TSpace> > (new JacobianEmulator<TSpace>());

        KRATOS_CATCH( "" );
    }

    /**
     * @brief 
     * This method initializes the internal counters and constructs the previous step Jacobian emulator.
     * The Jacobian emulator is recursively build at each time step according to the buffer size.
     */
    void InitializeSolutionStep() override
    {
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
     * @brief Performs the solution update
     * The correction is computed using an inverse Jacobian approximation obtained with a recursive matrix-free version of the MVQN (MultiVector Quasi-Newton method).
     * @param rResidualVector: Residual vector from the residual evaluation
     * @param rIterationGuess: Current iteration guess to be corrected. Should be initialized outside the convergence accelerator.
     */
    void UpdateSolution(
        const VectorType& rResidualVector,
        VectorType& rIterationGuess) override
    {
        KRATOS_TRY;

        const auto problem_size = TSpace::Size(rResidualVector);

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
                mpCurrentJacobianEmulatorPointer->ApplyJacobian(mpResidualVector_1, pInitialCorrection);

                TSpace::UnaliasedAdd(rIterationGuess, -1.0, *pInitialCorrection); // Recall the minus sign coming from the Taylor expansion of the residual (Newton-Raphson)
            }
        } else {
            // Gather the new observation matrices column information
            VectorPointerType pNewColV(new VectorType(*mpResidualVector_1));
            VectorPointerType pNewColW(new VectorType(*mpIterationValue_1));

            TSpace::UnaliasedAdd(*pNewColV, -1.0, *mpResidualVector_0); // NewColV = ResidualVector_1 - ResidualVector_0
            TSpace::UnaliasedAdd(*pNewColW, -1.0, *mpIterationValue_0); // NewColW = IterationValue_1 - IterationValue_0

            // Observation matrices information filling
            bool info_added = false;
            const std::size_t n_data_cols = mpCurrentJacobianEmulatorPointer->GetNumberOfDataCols();
            if (n_data_cols < problem_size) {
                info_added = (mpCurrentJacobianEmulatorPointer)->AppendDataColumns(*pNewColV, *pNewColW, mAbsCutOff);
            } else {
                info_added = (mpCurrentJacobianEmulatorPointer)->DropAndAppendDataColumns(*pNewColV, *pNewColW, mAbsCutOff);
            }
            KRATOS_WARNING_IF("MVQNRecursiveJacobianConvergenceAccelerator", !info_added)  << "Information not added to the observation matrices." << std::endl;

            // Apply the current step inverse Jacobian emulator to the residual vector
            VectorPointerType pIterationCorrection(new VectorType(rResidualVector));
            mpCurrentJacobianEmulatorPointer->ApplyJacobian(mpResidualVector_1, pIterationCorrection);
            
            TSpace::UnaliasedAdd(rIterationGuess, -1.0, *pIterationCorrection); // Recall the minus sign coming from the Taylor expansion of the residual (Newton-Raphson)
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Do the MVQN variables update
     * Updates the MVQN iteration variables for the next non-linear iteration
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

    ///@}
private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    double mOmega_0;                                                    // Relaxation factor for the initial fixed point iteration
    double mAbsCutOff;                                                  // Tolerance for the absolute cut-off criterion
    unsigned int mJacobianBufferSize;                                   // User-defined Jacobian buffer-size
    unsigned int mConvergenceAcceleratorStep;                           // Convergence accelerator steps counter
    unsigned int mConvergenceAcceleratorIteration;                      // Convergence accelerator iteration counter
    bool mConvergenceAcceleratorFirstCorrectionPerformed;               // Indicates that the initial fixed point iteration has been already performed

    VectorPointerType mpResidualVector_0;                               // Previous iteration residual vector pointer
    VectorPointerType mpResidualVector_1;                               // Current iteration residual vector pointer
    VectorPointerType mpIterationValue_0;                               // Previous iteration guess pointer
    VectorPointerType mpIterationValue_1;                               // Current iteration guess pointer

    JacobianEmulatorPointerType mpCurrentJacobianEmulatorPointer;       // Current step Jacobian approximator pointer

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

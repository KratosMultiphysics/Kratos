//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Quirin Aumann
//                   Matthias Ebert
//

#if !defined(MOR_SECOND_ORDER_IRKA_STRATEGY)
#define MOR_SECOND_ORDER_IRKA_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/builtin_timer.h"
// #include "utilities/qr_utility.h"
// #include "utilities/openmp_utils.h"
#include "custom_strategies/custom_strategies/mor_offline_second_order_strategy.hpp"

//default builder and solver
// #include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"

// includes for linear solver factory
#include "factories/linear_solver_factory.h"
#include "includes/kratos_parameters.h"

// include the eigensolver
// #include "custom_solvers/gen_eigensystem_solver.h"
#include "custom_utilities/generalized_eigenvalue_utility.h"
#include "custom_utilities/orthogonalization_utility.hpp"
#include "custom_utilities/complex_sort_utility.hpp"
// #include "../EigenSolversApplication/custom_solvers/eigen_direct_solver.h"

#include "omp.h"

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

/**
 * @class MorSecondOrderIRKAStrategy
 * @ingroup KratosMORApplication
 * @brief This is the Iterative Rational Krylov Algorithm
 * @details This strategy builds a reduction basis according to
 *      Wyatt, Sarah Alice. Issues in interpolatory model reduction:
 *      Inexact solves, second-order systems and DAEs. Diss. Virginia Tech, 2012.
 * @author Matthias Ebert, based on code of Quirin Aumann
 */
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver,
          class TReducedSparseSpace,
          class TReducedDenseSpace,
          bool TUseModalDamping
          >
class MorSecondOrderIRKAStrategy
    : public MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver, TReducedSparseSpace, TReducedDenseSpace >
{

    using complex = std::complex<double>;

  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorSecondOrderIRKAStrategy);

    typedef TUblasSparseSpace<complex> ComplexSparseSpaceType;

    typedef MorOfflineSecondOrderStrategy<TSparseSpace, TDenseSpace, TLinearSolver, TReducedSparseSpace, TReducedDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef TDenseSpace DenseSpaceType;

    typedef TLinearSolver TLinearSolverType;

    typedef typename TDenseSpace::VectorType TDenseVectorType;

    typedef typename TDenseSpace::MatrixType TDenseMatrixType;

    typedef typename TDenseSpace::MatrixPointerType TDenseMatrixPointerType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename ComplexSparseSpaceType::MatrixType ComplexSparseMatrixType;

    typedef typename ComplexSparseSpaceType::MatrixPointerType ComplexSparseMatrixPointerType;

    typedef typename ComplexSparseSpaceType::VectorType ComplexSparseVectorType;

    typedef typename ComplexSparseSpaceType::VectorPointerType ComplexSparseVectorPointerType;

    typedef typename TReducedSparseSpace::MatrixType TReducedSparseMatrixType;

    typedef typename TReducedSparseSpace::VectorType TReducedSparseVectorType;

    typedef typename TReducedDenseSpace::MatrixType TReducedDenseMatrixType;

    typedef typename TReducedDenseSpace::VectorType TReducedDenseVectorType;

    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Constructor for symmetric problems
     * @param rModelPart The model part of the problem
     * @param pScheme The integration schemed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    MorSecondOrderIRKAStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename BaseType::TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename TLinearSolverType::Pointer pNewLinearSolver,
        ComplexVector SamplingPoints,
        size_t MaxIter,
        double Tolerance)
        : BaseType(rModelPart, pScheme, pBuilderAndSolver, pNewLinearSolver, true),
            mMaxIter(MaxIter), mTolerance(Tolerance)
    {
        KRATOS_TRY;

        this->InitializeSamplingPoints(SamplingPoints);

        KRATOS_CATCH("");
    }

    /**
     * Constructor for unsymmetric problems
     * @param rModelPart The model part of the problem
     * @param pScheme The integration schemed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    MorSecondOrderIRKAStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename BaseType::TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename TLinearSolverType::Pointer pNewLinearSolver,
        typename TLinearSolverType::Pointer pNewAdjointLinearSolver,
        ComplexVector SamplingPoints,
        size_t MaxIter,
        double Tolerance)
        : BaseType(rModelPart, pScheme, pBuilderAndSolver, pNewLinearSolver, pNewAdjointLinearSolver, true),
            mMaxIter(MaxIter), mTolerance(Tolerance)
    {
        KRATOS_TRY;

        this->InitializeSamplingPoints(SamplingPoints);

        KRATOS_CATCH("");
    }

    /**
     * Constructor for frequency limited approximation of symmetric problems
     * @param rModelPart The model part of the problem
     * @param pScheme The integration schemed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    MorSecondOrderIRKAStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename BaseType::TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename TLinearSolverType::Pointer pNewLinearSolver,
        ComplexVector SamplingPoints,
        complex LimitLow,
        complex LimitHigh,
        size_t MaxIter,
        double Tolerance)
        : BaseType(rModelPart, pScheme, pBuilderAndSolver, pNewLinearSolver, true),
            mLimitLow(LimitLow), mLimitHigh(LimitHigh), mMaxIter(MaxIter), mTolerance(Tolerance)
    {
        KRATOS_TRY;

        this->InitializeSamplingPoints(SamplingPoints);
        this->mUseFrequencyLimits = true;

        KRATOS_CATCH("");
    }

    /**
     * Constructor for frequency limited approximation of unsymmetric problems
     * @param rModelPart The model part of the problem
     * @param pScheme The integration schemed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    MorSecondOrderIRKAStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename BaseType::TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename TLinearSolverType::Pointer pNewLinearSolver,
        typename TLinearSolverType::Pointer pNewAdjointLinearSolver,
        ComplexVector SamplingPoints,
        complex LimitLow,
        complex LimitHigh,
        size_t MaxIter,
        double Tolerance)
        : BaseType(rModelPart, pScheme, pBuilderAndSolver, pNewLinearSolver, pNewAdjointLinearSolver, true),
            mLimitLow(LimitLow), mLimitHigh(LimitHigh), mMaxIter(MaxIter), mTolerance(Tolerance)
    {
        KRATOS_TRY;

        this->InitializeSamplingPoints(SamplingPoints);
        this->mUseFrequencyLimits = true;

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorSecondOrderIRKAStrategy() override
    {

        if (this->pGetM() != nullptr)
            TSparseSpace::Clear(this->pGetM());
        if (this->pGetK() != nullptr)
            TSparseSpace::Clear(this->pGetK());
        if (this->pGetC() != nullptr)
            TSparseSpace::Clear(this->pGetC());
        // if (this->mp != nullptr)
        //     TSparseSpace::Clear(mpb);


        this->Clear();
    }

    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (this->mSolutionStepIsInitialized == false)
        {
            BaseType::InitializeSolutionStep();

            const size_t reduced_system_size = mSamplingPoints.size();
            const size_t system_size = this->GetSystemSize();

            TReducedDenseSpace::Resize(this->GetKr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetDr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetMr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetBasis(), system_size, reduced_system_size);
            if( !this->SystemIsSymmetric() ) {
                TReducedDenseSpace::Resize(this->GetBasisLeft(), system_size, reduced_system_size);
            }

            this->mSolutionStepIsInitialized = true;
        }


        KRATOS_CATCH("");
    }

    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        return 0;

        KRATOS_CATCH("")

    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY;
        std::cout << "hello! this is where the second order IRKA MOR magic happens" << std::endl;
        typename TSchemeType::Pointer p_scheme = this->GetScheme();
        typename BaseType::TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver();
        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
        complex z(0,1);

        // get the full system size matrices
        TSystemMatrixType& r_K = this->GetK();
        TSystemMatrixType& r_M = this->GetM();
        TSystemMatrixType& r_D = this->GetC();
        TSystemVectorType& r_RHS = this->GetSystemVector();
        TSystemVectorType& r_ov = this->GetOutputVector();
        TReducedDenseMatrixType& r_basis = this->GetBasis();
        TReducedDenseMatrixType& r_basis_l = this->GetBasisLeft();

        // get the reduced matrices
        TReducedDenseMatrixType& r_Kr = this->GetKr();
        TReducedDenseMatrixType& r_Mr = this->GetMr();
        TReducedDenseMatrixType& r_Dr = this->GetDr();

        // set the system sizes
        const size_t system_size = this->GetSystemSize();
        const size_t n_sampling_points = mSamplingPoints.size();
        const size_t reduced_system_size = n_sampling_points;

        // copy mass matrix, damping matrix, and rhs to the complex space
        ComplexSparseMatrixType r_M_tmp(r_M.size1(), r_M.size2());
        noalias(r_M_tmp) = r_M;
        ComplexSparseMatrixType r_D_tmp(r_D.size1(), r_D.size2());
        noalias(r_D_tmp) = r_D;
        ComplexSparseVectorType r_RHS_tmp = ComplexSparseVectorType(r_RHS);
        ComplexSparseVectorType r_ov_tmp = ComplexSparseVectorType(r_ov);

        // create a complex stiffness matrix if required
        TReducedSparseMatrixType r_K_cplx(r_M.size1(), r_M.size2());
        if( TUseModalDamping )
        {
            //cast complex to double (real part). This is a hack for real valued strategies to make it compile. It will not be called
            noalias(r_K_cplx) = r_K + TReducedSparseMatrixType(r_D) * reinterpret_cast<typename TReducedSparseSpace::DataType(&)[2]>(z)[0];

            //set reduced damping matrix to zero
            noalias(r_Dr) = ZeroMatrix(r_Dr.size1(), r_Dr.size2());
        }

        // create solution vector
        ComplexSparseVectorPointerType tmp_dx = ComplexSparseSpaceType::CreateEmptyVectorPointer();
        ComplexSparseVectorType& r_tmp_dx = *tmp_dx;
        ComplexSparseSpaceType::Resize(r_tmp_dx, system_size);
        ComplexSparseSpaceType::Set(r_tmp_dx, 0.0);

        // create dynamic stiffness matrix
        ComplexSparseMatrixType r_kdyn(system_size, system_size);

        if( BaseType::GetEchoLevel() == 7 ) {
            this->template MatrixOutput<SparseSpaceType>(r_K, "K");
            this->template MatrixOutput<SparseSpaceType>(r_D, "D");
            this->template MatrixOutput<SparseSpaceType>(r_M, "M");
            this->template VectorOutput<SparseSpaceType>(r_RHS, "RHS");
            this->template VectorOutput<SparseSpaceType>(r_ov, "o");
        }

        // create loop variables
        size_t iter = 0;
        double error = 1;
        std::vector<double> error_vec(reduced_system_size, 1);
        ComplexVector eigenvalues;
        auto samplingPoints_old = mSamplingPoints;
        double projection_time;
        double construction_time;

        BuiltinTimer irka_overall_time;
        while( iter < mMaxIter && error > mTolerance )
        {
            BuiltinTimer irka_iteration_time;
            BuiltinTimer irka_projection_time;
            BuiltinTimer basis_construction_time;

            for( size_t i=0; i<n_sampling_points/2; ++i )
            {
                // if already converged, skip the expansion point
                if( error_vec[2*i] < mTolerance ) {
                    KRATOS_INFO_IF("\tSkipping expansion points", BaseType::GetEchoLevel() > 1 && rank == 0)
                        << 2*i << "/" << 2*i+1 << " with error " << error_vec[2*i] << std::endl;
                    continue;
                }

                // build dynamic stiffness matrix
                if( TUseModalDamping) {
                    if( this->SystemIsSymmetric() ) {
                        noalias(r_kdyn) = r_K_cplx + std::pow( mSamplingPoints(2*i), 2 ) * r_M_tmp;
                    } else {
                        noalias(r_kdyn) = r_K_cplx - std::pow( mSamplingPoints(2*i), 2 ) * r_M_tmp;
                    }
                } else {
                    if( this->SystemIsSymmetric() ) {
                        noalias(r_kdyn) = r_K + mSamplingPoints(2*i) * r_D_tmp + std::pow( mSamplingPoints(2*i), 2 ) * r_M_tmp;
                    } else {
                        noalias(r_kdyn) = r_K + mSamplingPoints(2*i) * r_D_tmp - std::pow( mSamplingPoints(2*i), 2 ) * r_M_tmp;
                    }
                }

                // solve for current expansion point
                this->mpLinearSolver->Solve( r_kdyn, r_tmp_dx, r_RHS_tmp );

                KRATOS_DEBUG_ERROR_IF(isinf(norm_2(r_tmp_dx))) << "No solution could be obtained (norm infinity)!" << std::endl;

                // update basis
                column(r_basis, 2*i)   = real(r_tmp_dx);
                column(r_basis, 2*i+1) = imag(r_tmp_dx);

                // for unsymmetric systems, solve the adjoint problem and update left basis
                if( !this->SystemIsSymmetric() ) {
                    this->mpAdjointLinearSolver->Solve( r_kdyn, r_tmp_dx, r_ov_tmp);

                    column(r_basis_l, 2*i)   = real(r_tmp_dx);
                    column(r_basis_l, 2*i+1) = imag(r_tmp_dx);
                }
            }
            construction_time = basis_construction_time.ElapsedSeconds();

            // orthogonalize the bases
            OrthogonalizationUtility::OrthogonalizeQR<TReducedDenseSpace>(r_basis);
            if( !this->SystemIsSymmetric() ) {
                OrthogonalizationUtility::OrthogonalizeQR<TReducedDenseSpace>(r_basis_l);
            }

            // project onto reduced space
            if( TUseModalDamping )
            {
                if( this->SystemIsSymmetric() ) {
                    this->template ProjectMatrix<TReducedSparseMatrixType>(r_K_cplx, r_basis, r_Kr);
                    this->template ProjectMatrix<TSystemMatrixType>(r_M, r_basis, r_Mr);
                } else {
                    this->template ProjectMatrix<TReducedSparseMatrixType>(r_K_cplx, r_basis, r_basis_l, r_Kr);
                    this->template ProjectMatrix<TSystemMatrixType>(r_M, r_basis, r_basis_l, r_Mr);
                }
            }
            else
            {
                if( this->SystemIsSymmetric() ) {
                    this->template ProjectMatrix<TSystemMatrixType>(r_K, r_basis, r_Kr);
                    this->template ProjectMatrix<TSystemMatrixType>(r_D, r_basis, r_Dr);
                    this->template ProjectMatrix<TSystemMatrixType>(r_M, r_basis, r_Mr);
                } else {
                    this->template ProjectMatrix<TSystemMatrixType>(r_K, r_basis, r_basis_l, r_Kr);
                    this->template ProjectMatrix<TSystemMatrixType>(r_D, r_basis, r_basis_l, r_Dr);
                    this->template ProjectMatrix<TSystemMatrixType>(r_M, r_basis, r_basis_l, r_Mr);
                }
            }
            projection_time = irka_projection_time.ElapsedSeconds();

            // compute eigenvalues in reduced space
            if( TUseModalDamping )
            {
                // compute generalized eigenvalues
                ComplexVector tmp_eigenvalues;
                GeneralizedEigenvalueUtility::Compute<TReducedDenseSpace>(r_Kr, r_Mr, tmp_eigenvalues);

                std::for_each(tmp_eigenvalues.begin(), tmp_eigenvalues.end(),
                    [](complex &a) {a = std::sqrt(a);});

                eigenvalues.resize(2*reduced_system_size, false);
                noalias(subrange(eigenvalues, 0, reduced_system_size)) = tmp_eigenvalues;
                std::for_each(tmp_eigenvalues.begin(), tmp_eigenvalues.end(),
                    [](complex &a) {a = std::conj(a);});
                noalias(subrange(eigenvalues, reduced_system_size, reduced_system_size*2)) = tmp_eigenvalues;

                // sort by absolute value with negative imaginary part in front (pair first, then stable sort)
                ComplexSortUtility::PairComplexConjugates(eigenvalues);
                std::stable_sort(eigenvalues.begin(), eigenvalues.end(),
                    [](complex z1, complex z2) {return std::abs(z1) < std::abs(z2);});

                this->UpdateSamplingPoints(eigenvalues, reduced_system_size);
            }
            else
            {
                // compute generalized eigenvalues
                GeneralizedEigenvalueUtility::ComputePolynomial<TReducedDenseSpace>(r_Kr, r_Dr, r_Mr, eigenvalues);

                // use mirror images of first r eigenvalues as expansion points
                eigenvalues *= -1.;
                ComplexSortUtility::PairComplexConjugates(eigenvalues);
                this->UpdateSamplingPoints(eigenvalues, reduced_system_size);
            }

            // calculate error
            error = norm_2(mSamplingPoints - samplingPoints_old);
            error /= norm_2(samplingPoints_old);
            for( std::size_t i=0; i<reduced_system_size; ++i ) {
                error_vec[i] = std::abs((mSamplingPoints[i] - samplingPoints_old[i]) / samplingPoints_old[i]);
            }

            // update loop variables
            iter++;
            samplingPoints_old = mSamplingPoints;

            KRATOS_INFO_IF("IRKA error after iteration " + iter, BaseType::GetEchoLevel() > 0 && rank == 0)
                << ": e=" << error << std::endl;
            KRATOS_INFO_IF("IRKA Iteration Solve Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << irka_iteration_time.ElapsedSeconds() << std::endl;
            KRATOS_INFO_IF("\tProjection Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << projection_time << std::endl;
            KRATOS_INFO_IF("\tBasis update Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << construction_time << std::endl;
            KRATOS_INFO_IF("\tCurrent expansion points", BaseType::GetEchoLevel() > 1 && rank == 0)
                << mSamplingPoints << std::endl;
            KRATOS_INFO_IF("\tRespective errors", BaseType::GetEchoLevel() > 1 && rank == 0)
                << error_vec << std::endl;
        }

        KRATOS_INFO_IF("IRKA Complete Solve Time", BaseType::GetEchoLevel() > 0 && rank == 0)
            << irka_overall_time.ElapsedSeconds() << std::endl;

        if( iter == mMaxIter ) {
            this->MaxIterationsExceeded();
        }

		return true;

        KRATOS_CATCH("");
    }

    void FinalizeSolutionStep() override
    {
        // matrices have been projected inside the IRKA loop
        auto& r_RHS = this->GetSystemVector();
        auto& r_force_vector_reduced = this->GetRHSr();
        auto& r_output_vector = this->GetOutputVector();
        auto& r_output_vector_r = this->GetOVr();
        auto& r_basis = this->GetBasis();
        auto& r_basis_l = this->GetBasisLeft();

        const std::size_t reduced_system_size = r_basis.size2();

        r_force_vector_reduced.resize( reduced_system_size, false);
        axpy_prod(r_RHS, r_basis_l, r_force_vector_reduced, true);

        r_output_vector_r.resize( reduced_system_size, false);
        r_output_vector_r = prod( r_output_vector, r_basis );
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ComplexVector GetSamplingPoints()
    {
        return mSamplingPoints;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "MorSecondOrderIRKAStrategy";
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

  private:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

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

    ///@}

  protected:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ComplexVector mSamplingPoints; //vector of currently used sampling points

    bool mUseFrequencyLimits = false; //flag if the approximation should be limited to a specific frequency range
    complex mLimitLow; //lower frequency limit
    complex mLimitHigh; //higher frequency limit

    size_t mMaxIter; //maximum iterations
    double mTolerance; // tolerance

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Initialize sampling points with complex conjugates
     */
    void InitializeSamplingPoints(ComplexVector SamplingPoints)
    {
        KRATOS_TRY;

        // create complex conjugates of the sampling points and sort them
        const size_t n_sampling_points = SamplingPoints.size();
        mSamplingPoints = ComplexZeroVector(2*n_sampling_points);
        for( size_t i=0; i<n_sampling_points; ++i )
        {
            mSamplingPoints(i) = SamplingPoints(i);
            mSamplingPoints(n_sampling_points+i) = std::conj(SamplingPoints(i));
        }
        ComplexSortUtility::PairComplexConjugates(mSamplingPoints);

        KRATOS_ERROR_IF( mMaxIter < 1 ) << "Invalid number of maximal iterations provided" << std::endl;
        KRATOS_ERROR_IF( (mTolerance >= 1.) || (mTolerance < 0.) ) << "Invalid tolerance provided" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * @brief Update the used sampling points from the eigensolver's result
     */
    void UpdateSamplingPoints(vector<complex> eigenvalues, std::size_t reduced_system_size)
    {
        if( this->mUseFrequencyLimits ) {
            std::vector<double> diff(eigenvalues.size());
            std::vector<std::size_t> idx(eigenvalues.size());
            std::iota(idx.begin(), idx.end(), 0);

            double this_abs = 0;
            const double low_abs = std::abs(mLimitLow);
            const double high_abs = std::abs(mLimitHigh);

            // determine which eigenvalues lie inside
            for( std::size_t i=0; i<eigenvalues.size(); ++i ) {
                this_abs = std::real(eigenvalues[i]);
                if( low_abs <= this_abs && this_abs <= high_abs ) {
                    diff[i] = 0;
                } else if( this_abs < low_abs ) {
                    diff[i] = low_abs - this_abs;
                } else if( high_abs < this_abs ) {
                    diff[i] = this_abs - high_abs;
                }
            }

            // sort by smallest difference outside desired range
            std::stable_sort(idx.begin(), idx.end(),
                [&diff](std::size_t i1, std::size_t i2) {return diff[i1] < diff[i2];});

            // update with the first r eigenvalues inside range
            ComplexVector tmp_ev(reduced_system_size);
            for( std::size_t i=0; i<reduced_system_size; ++i ) {
                tmp_ev[i] = eigenvalues[idx[i]];
            }

            KRATOS_WARNING_IF( "Frequency limited IRKA", diff[idx[reduced_system_size]] == 0 ) << "Possibly too small subspace for interval" << std::endl;

            // update sampling points
            noalias(mSamplingPoints) = tmp_ev;

        } else {
            noalias(mSamplingPoints) = subrange(eigenvalues, 0, reduced_system_size);
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods

    ///@{

    /**
     * Copy constructor.
     */

    MorSecondOrderIRKAStrategy(const MorSecondOrderIRKAStrategy &Other){};

    ///@}

}; /* Class MorSecondOrderIRKAStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_SECOND_ORDER_IRKA_STRATEGY  defined */
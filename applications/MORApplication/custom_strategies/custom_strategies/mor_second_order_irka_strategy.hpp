//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(MOR_SECOND_ORDER_IRKA_STRATEGY)
#define MOR_SECOND_ORDER_IRKA_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/builtin_timer.h"
#include "utilities/qr_utility.h"
#include "utilities/openmp_utils.h"
#include "custom_strategies/custom_strategies/mor_offline_second_order_strategy.hpp"

//default builder and solver
#include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"

// includes for linear solver factory
#include "factories/linear_solver_factory.h"
#include "includes/kratos_parameters.h"

// include the eigensolver
// #include "custom_solvers/gen_eigensystem_solver.h"
#include "custom_utilities/generalized_eigenvalue_utility.hpp"
#include "custom_utilities/orthogonalization_utility.hpp"
#include "custom_utilities/complex_sort_utility.hpp"
// #include "../EigenSolversApplication/custom_solvers/eigen_direct_solver.h"

#include "omp.h"

namespace Kratos
{
    using complex = std::complex<double>;

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
 * @ingroup KratosCore
 * @brief This is the Iterative Rational Krylov Algorithm
 * @details This strategy builds the reduced K and M matrices and outputs them
 * @author Matthias Ebert, based on code of Quirin Aumann
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver, //= LinearSolver<TSparseSpace,TDenseSpace>
          class TReducedSparseSpace,
          class TReducedDenseSpace
          >
class MorSecondOrderIRKAStrategy
    // : public MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver >
    : public MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver, TReducedSparseSpace, TReducedDenseSpace >
{
  
  Parameters mParam;

  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorSecondOrderIRKAStrategy);

    typedef TUblasSparseSpace<complex> ComplexSparseSpaceType;
    typedef TUblasDenseSpace<complex> ComplexDenseSpaceType;

    typedef MorOfflineSecondOrderStrategy<TSparseSpace, TDenseSpace, TLinearSolver, TReducedSparseSpace, TReducedDenseSpace> BaseType;

    // typedef LinearSolverFactory<ComplexSparseSpaceType,  ComplexDenseSpaceType> LinearSolverFactoryType;

    // typedef SystemMatrixBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;

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

    // typedef ComplexSparseSpaceType TSolutionSpace;

    typedef typename TReducedSparseSpace::MatrixType TReducedSparseMatrixType;

    typedef typename TReducedSparseSpace::VectorType TReducedSparseVectorType;

    typedef typename TReducedDenseSpace::MatrixType TReducedDenseMatrixType;

    typedef typename TReducedDenseSpace::VectorType TReducedDenseVectorType;

    typedef LinearSolver<ComplexSparseSpaceType, ComplexDenseSpaceType> ComplexLinearSolverType;

    typedef SystemMatrixBuilderAndSolver< ComplexSparseSpaceType, ComplexDenseSpaceType, ComplexLinearSolverType > TComplexBuilderAndSolverType;

    // typedef Preconditioner<SparseSpaceType, DenseSpaceType> PreconditionerType;
    // typedef Reorderer<SparseSpaceType, DenseSpaceType> ReordererType;
    // typedef IterativeSolver<SparseSpaceType, DenseSpaceType, PreconditionerType, ReordererType> IterativeSolverType;

    // typedef LinearSolver<SparseSpaceType, DenseSpaceType> LinearSolverType;
    // typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    // typedef GenEigensystemSolver< SparseSpaceType, LocalSpaceType, PreconditionerType, ReordererType> EigenSystemSolverType;



    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration schemed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    MorSecondOrderIRKAStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename BaseType::TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename TLinearSolverType::Pointer pNewLinearSolver,  // required for intermediate solves
        vector< std::complex<double> > SamplingPoints,              // real part of initial expansion points
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pBuilderAndSolver, pNewLinearSolver, true, MoveMeshFlag)//, 
            // mSamplingPoints(SamplingPoints)
    {
        KRATOS_TRY;

        mTolerance = 1e-12;
        mMaxIter = 20;

        // create complex conjugates of the sampling points
        const size_t n_sampling_points = SamplingPoints.size();
        mSamplingPoints = ComplexZeroVector(2*n_sampling_points);
        for( size_t i=0; i<n_sampling_points; ++i )
        {
            mSamplingPoints(i) = SamplingPoints(i);
            mSamplingPoints(n_sampling_points+i) = std::conj(SamplingPoints(i));
        }
        // subrange(mSamplingPoints, 0, n_sampling_points) = SamplingPoints;
        // subrange(mSamplingPoints, n_sampling_points+1, n_sampling_points) = std::conj(SamplingPoints);

        // KRATOS_WATCH(SamplingPoints)
        // KRATOS_WATCH(mSamplingPoints)
        // Saving the scheme
        // this->SetScheme(pScheme);

        // // Setting up the default builder and solver
        // this->SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer(
        //     new TBuilderAndSolverType(pNewLinearSolver)));

        // // Saving the linear solver        
        // this->SetLinearSolver(pNewLinearSolver);

        // // Set flags to start correctly the calculations
        // mSolutionStepIsInitialized = false;
        // mInitializeWasPerformed = false;

        // // Tells to the builder and solver if the reactions have to be calculated or not
        // this->GetBuilderAndSolver()->SetCalculateReactionsFlag(false);

        // // Tells to the Builder And Solver if the system matrix and vectors need to
        // // be reshaped at each step or not
        // this->GetBuilderAndSolver()->SetReshapeMatrixFlag(false);

        // // Set EchoLevel to the default value (only time is displayed)
        // this->SetEchoLevel(1);

        // // By default the matrices are rebuilt at each iteration
        // this->SetRebuildLevel(0);

        // // Set members
        // // store the real and imaginary parts of the initial expansion points into a complex vector
        // mSamplingPoints.resize(samplingPoints_real.size());
        // for(size_t i=0; i<samplingPoints_real.size(); i++)
        // {
        //     mSamplingPoints(i) = complex( samplingPoints_real(i), samplingPoints_imag(i) );
        // }

        // mpM = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
        // mpK = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
        // mpD = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
        // mpb = ComplexSparseSpaceType::CreateEmptyVectorPointer();

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorSecondOrderIRKAStrategy() override
    {

        // if (mpM != nullptr)
        //     TSolutionSpace::Clear(mpM);
        // if (mpK != nullptr)
        //     TSolutionSpace::Clear(mpK);
        // if (mpD != nullptr)
        //     TSolutionSpace::Clear(mpD);
        // if (mpb != nullptr)
        //     TSolutionSpace::Clear(mpb);


        this->Clear();
    }

    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (this->mSolutionStepIsInitialized == false)
        {
            BaseType::InitializeSolutionStep();

            const std::size_t reduced_system_size = mSamplingPoints.size();
            const unsigned int system_size = this->GetBuilderAndSolver()->GetEquationSystemSize();

            TReducedDenseSpace::Resize(this->GetKr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetDr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetMr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetBasis(), system_size, reduced_system_size);


            this->mSolutionStepIsInitialized = true;
        }


        KRATOS_CATCH("");
    }

    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        std::cout << "CHECK IF SAMPING POINTS ARE COMPLEX CONJUGATES\n";

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

        // get the full system size matrices                
        TSystemMatrixType& r_K = this->GetSystemMatrix();
        TSystemMatrixType& r_M = this->GetMassMatrix();
        TSystemMatrixType& r_D = this->GetDampingMatrix();
        TSystemVectorType& r_RHS = this->GetSystemVector();  // originally RHS, check if results are bad
        TReducedDenseMatrixType& r_basis = this->GetBasis();

        // get the reduced matrices
        TReducedDenseMatrixType& r_Kr = this->GetKr();
        TReducedDenseMatrixType& r_Mr = this->GetMr();
        TReducedDenseMatrixType& r_Dr = this->GetDr();
        TReducedDenseVectorType& r_RHSr = this->GetRHSr();

        // TSystemVectorType tmp(r_K.size1(), 0.0);


        // start timer (after building the full system)
        double start_time = OpenMPUtils::GetCurrentTime();


        // DEBUG: print matrix information 
        //this->EchoInfo(0);


        // set the system sizes
        const unsigned int system_size = p_builder_and_solver->GetEquationSystemSize(); // n
        const std::size_t n_sampling_points = mSamplingPoints.size(); // number of sampling points
        const std::size_t reduced_system_size = n_sampling_points; // r

        std::cout<<"Full system size: "<<system_size<<std::endl;
        std::cout<<"Reduced system size: "<<reduced_system_size<<std::endl;

        // check if the reduced system size is odd (conjugate eigenvalues don't make sense then)
        // if ( reduced_system_size % 2 != 0 ){
        //     KRATOS_TRY;
        //     KRATOS_ERROR << "Reduced system size is odd! Use an even number for the reduced size." << std::endl;
        //     KRATOS_CATCH("")
        // }

        // sort the sampling points after their real part in increasing order
        // conjugate pairs stay together, with the minus element being the first
        // this->sort_conjugate_pairs(mSamplingPoints);
        ComplexSortUtility::PairComplexConjugates(mSamplingPoints);
        KRATOS_WATCH(mSamplingPoints)

        // store the sampling points for error calculations in the convergence loop for convergence check


        ComplexSparseMatrixType r_M_tmp = ComplexSparseMatrixType(r_M);
        ComplexSparseVectorType r_RHS_tmp = ComplexSparseVectorType(r_RHS);
        // r_RHS_tmp = ComplexSparseVectorType(r_RHS);
        // KRATOS_WATCH(*(this->mpLinearSolver))
        // typename TLinearSolverType::Pointer p_compl_solver_par = this->mpLinearSolver;
        // TLinearSolverType bla = TLinearSolverType(r_linear_solver);
        // KRATOS_WATCH(r_linear_solver)
        // KRATOS_WATCH(bla)
        // auto& r_linear_solver = *(this->mpLinearSolver);
        // using actual_type = decltype(*(this->mpLinearSolver));

        // auto xx(r_linear_solver);
        // KRATOS_WATCH(xx)

        ComplexSparseVectorPointerType tmp_dx = ComplexSparseSpaceType::CreateEmptyVectorPointer();
        ComplexSparseVectorType& r_tmp_dx   = *tmp_dx;
        ComplexSparseSpaceType::Resize(r_tmp_dx, system_size); // n x 1
        ComplexSparseSpaceType::Set(r_tmp_dx, 0.0); // set vector to zero

        auto  kdyn = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_kdyn   = *kdyn;
        ComplexSparseSpaceType::Resize(r_kdyn, system_size, system_size); // n x n
        ComplexSparseSpaceType::SetToZero(r_kdyn);

        //initial basis
        for( size_t i=0; i<n_sampling_points/2; ++i )
        {
            r_kdyn = r_D;
            r_kdyn *= mSamplingPoints(2*i);
            r_kdyn += r_K;
            r_kdyn += std::pow( mSamplingPoints(2*i), 2.0 ) * r_M_tmp;

            this->mpLinearSolver->Solve( r_kdyn, r_tmp_dx, r_RHS_tmp); // Ax = b, solve for x

            column(r_basis, 2*i)   = real(r_tmp_dx);
            column(r_basis, 2*i+1) = imag(r_tmp_dx);
        }

        OrthogonalizationUtility::OrthogonalizeQR<TReducedDenseSpace>(r_basis);

        vector<complex> eigenvalues;
        // vector<complex> error_vector(reduced_system_size);
        std::cout << "1\n";
        // noalias(error_vector) = ComplexZeroVector(reduced_system_size);
        std::cout << "2\n";
        size_t iter = 0;
        double error = 1;
        auto samplingPoints_old = mSamplingPoints;
        while( iter < mMaxIter && error > mTolerance )
        {
            BuiltinTimer bla;
            std::cout << "3\n";
            this->ProjectMatrix(r_K, r_basis, r_Kr);
            std::cout << "4\n";
            this->ProjectMatrix(r_D, r_basis, r_Dr);
            std::cout << "5\n";
            this->ProjectMatrix(r_M, r_basis, r_Mr);
            KRATOS_WATCH(bla.ElapsedSeconds())
            // this->ProjectMatrix(this->mpA, this->mpBasis, this->mpAr);

            GeneralizedEigenvalueUtility::Compute<TReducedDenseSpace>(r_Kr, r_Dr, r_Mr, eigenvalues);
            eigenvalues *= -1.;
            ComplexSortUtility::PairComplexConjugates(eigenvalues);
            KRATOS_WATCH(eigenvalues)
            // KRATOS_WATCH(-eigenvalues)

            noalias(mSamplingPoints) = subrange(eigenvalues, 0, reduced_system_size);
            KRATOS_WATCH(mSamplingPoints)
            KRATOS_WATCH(samplingPoints_old)
            error = norm_2(mSamplingPoints - samplingPoints_old);
            error /= norm_2(samplingPoints_old);

            //updated basis
            for( size_t i=0; i<n_sampling_points/2; ++i )
            {
                r_kdyn = r_D;
                r_kdyn *= mSamplingPoints(2*i);
                r_kdyn += r_K;
                r_kdyn += std::pow( mSamplingPoints(2*i), 2.0 ) * r_M_tmp;

                this->mpLinearSolver->Solve( r_kdyn, r_tmp_dx, r_RHS_tmp); // Ax = b, solve for x

                column(r_basis, 2*i)   = real(r_tmp_dx);
                column(r_basis, 2*i+1) = imag(r_tmp_dx);
            }

            OrthogonalizationUtility::OrthogonalizeQR<TReducedDenseSpace>(r_basis);

            iter++;
            samplingPoints_old = mSamplingPoints;
            std::cout << "IRKA error after iteration " << iter << ": e=" << error << "\n";
        }

        // KRATOS_WATCH(r_Vr_dense)
        /*
        // time for parallel solve
        double end_solve_par_first = OpenMPUtils::GetCurrentTime();
            std::cout<<"-- par solve first: "<<end_solve_par_first-start_solve_par_first<<std::endl;


        //orthogonalize V
        mQR_decomposition.compute( system_size, reduced_system_size, &(r_Vr_dense)(0,0) );
        mQR_decomposition.compute_q();

        for(size_t i=0; i < system_size; ++i)
        {
            for(size_t j=0; j < reduced_system_size; ++j)
            {
                r_Vr_sparse(i,j) = mQR_decomposition.Q(i,j);
            }
        }



        // needs to be resized, otherwise segfault
        SparseSpaceType::Resize(r_M_reduced, reduced_system_size, reduced_system_size); // r x r
        SparseSpaceType::Resize(r_K_reduced, reduced_system_size, reduced_system_size); // r x r
        SparseSpaceType::Resize(r_D_reduced, reduced_system_size, reduced_system_size); // r x r
        SparseSpaceType::Resize(r_b_reduced, reduced_system_size);                      // r x 1


        // initialize linearization matrices
        auto L1_ptr = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_L1 = *L1_ptr;
        SparseSpaceType::Resize(r_L1, 2*reduced_system_size, 2*reduced_system_size); // 2r x 2r

        auto L2_ptr = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_L2 = *L2_ptr;
        SparseSpaceType::Resize(r_L2, 2*reduced_system_size, 2*reduced_system_size); // 2r x 2r

        // identity matrix of size r x r
        identity_matrix<double> id_r (reduced_system_size); // r x r

        // build the Eigensolver and overwrite important settings
        // the number of eigenvalues is always double of the size of the reduced system
        // also we don't need to compute eigenvectors 
        mpLinearEigenSolver = new EigenSystemSolverType(mParam); // will be cleared automatically at the end
        mParam["number_of_eigenvalues"].SetInt(reduced_system_size*2); // 2r
        mParam["compute_eigenvectors"].SetBool(false);

        // initialize output for Eigensolver
        TDenseVectorType Eigenvalues;
        TDenseMatrixType Eigenvectors;  //not needed, but need to be set to call the Eigensolver
        vector<complex> lam_2r (2*reduced_system_size);

        // initialize helpers for the convergence check
        vector<double> abs_val_old (reduced_system_size);
        vector<double> abs_val_new (reduced_system_size);


        // initialize helpers for the projections
        auto  V_col = SparseSpaceType::CreateEmptyVectorPointer();
        auto& r_V_col   = *V_col;
        SparseSpaceType::Resize(r_V_col, system_size); // n x 1

        auto  T_col = SparseSpaceType::CreateEmptyVectorPointer();
        auto& r_T_col   = *T_col;
        SparseSpaceType::Resize(r_T_col, system_size); // n x 1

        TSystemVectorType tmp_col(reduced_system_size, 0.0); // r x 1


        // default settings for the iteration loop
        int iter = 1;
        int max_iter = 100;
        double err = 100.0;
        double tol = 1e-8;

        // start iteration loop
        while(iter < max_iter && err > tol)
        {

            // projections onto the reduced space
            for(size_t i=0; i<reduced_system_size; i++){
                r_V_col = column(r_Vr_sparse,i);

                axpy_prod(r_M_tmp, r_V_col, r_T_col);      // M*V  (=T)
                axpy_prod(r_T_col, r_Vr_sparse, tmp_col);  // V' * T
                column(r_M_reduced, i) = tmp_col;          // Mr = V' M V

                axpy_prod(r_K_tmp, r_V_col, r_T_col);      // K*V  (=T)
                axpy_prod(r_T_col, r_Vr_sparse, tmp_col);  // V' * T
                column(r_K_reduced, i) = tmp_col;          // Kr = V' K V

                axpy_prod(r_D_tmp, r_V_col, r_T_col);      // D*V  (=T)
                axpy_prod(r_T_col, r_Vr_sparse, tmp_col);  // V' * T
                column(r_D_reduced, i) = tmp_col;          // Dr = V' D V
            }


            // linearization of the quadratic eigenvalue problem to be solved using a GEP solver
            subrange(r_L1, 0, reduced_system_size, reduced_system_size, 2*reduced_system_size) = id_r;
            subrange(r_L1, reduced_system_size, 2*reduced_system_size, 0, reduced_system_size) = -1.0*r_K_reduced;
            subrange(r_L1, reduced_system_size, 2*reduced_system_size, reduced_system_size, 2*reduced_system_size) = -1.0*r_D_reduced;

            subrange(r_L2, 0, reduced_system_size, 0, reduced_system_size) = id_r;
            subrange(r_L2, reduced_system_size, 2*reduced_system_size, reduced_system_size, 2*reduced_system_size) = r_M_reduced;


            // output of this Solver is real (real parts and imaginary parts alternating)
            mpLinearEigenSolver->Solve(
                r_L1,
                r_L2,
                Eigenvalues,
                Eigenvectors
            );


            //store the Eigenvalues in a complex vector
            #pragma omp parallel for schedule(dynamic)
            for(size_t ii = 0; ii<2*reduced_system_size; ii++)
            {
                lam_2r(ii) = complex( Eigenvalues(2*ii), Eigenvalues(2*ii+1) );
            }


            // mirror the Eigenvalues and sort them after their real part in increasing order
            // conjugate pairs stay together, with the minus element being the first
            lam_2r *= -1;
            this->sort_conjugate_pairs(lam_2r);


            // reduce to r Eigenvalues so that the dimension of the reduced system does not increase
            // choose the r Eigenvalues which are closest to the imaginary axis (done by previous sorting)
            #pragma omp parallel for schedule(dynamic)
            for(size_t ii=0; ii < reduced_system_size; ii++)
            {
                mSamplingPoints(ii) = lam_2r(ii);
            }


            // calculate the error between old and new sampling points
            #pragma omp parallel for schedule(dynamic)           
            for(size_t ii=0 ; ii < reduced_system_size; ii++)
            {
                abs_val_old(ii) = abs(samplingPoints_old(ii));
            }
         
            #pragma omp parallel for schedule(dynamic)
            for(size_t ii=0 ; ii < reduced_system_size; ii++)
            {
                abs_val_new(ii) = abs(mSamplingPoints(ii));
            }

            // sort them to compare only the closest distance between iteration pairs
            sort(abs_val_old.begin(), abs_val_old.end());
            sort(abs_val_new.begin(), abs_val_new.end());

            err = norm_2(abs_val_new - abs_val_old) / norm_2(samplingPoints_old);

            // if the tolerance is already met, break here to avoid unnecessary computations
            if(err < tol) {break;}


            double start_solve_par = OpenMPUtils::GetCurrentTime(); // measure the time for the parallel solve block

            // update the old sampling points
            samplingPoints_old = mSamplingPoints;

            // update V
            #pragma omp parallel
            {

                // initialization have to be repeated because of the underlying fork-join model and scope

                // initialize helper variables for the complex solve
                // each thread will create its own version, so no conflict arises
                auto  tmp_Vr_col_ptr_par = ComplexDenseSpaceType::CreateEmptyVectorPointer();
                auto& r_tmp_Vr_col_par   = *tmp_Vr_col_ptr_par;
                ComplexDenseSpaceType::Resize(r_tmp_Vr_col_par, system_size); // n x 1
                ComplexDenseSpaceType::Set(r_tmp_Vr_col_par, 0.0); // set vector to zero

                auto  tmp_Vn_ptr_par = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
                auto& r_tmp_Vn_par   = *tmp_Vn_ptr_par;
                ComplexSparseSpaceType::Resize(r_tmp_Vn_par, system_size, system_size); // n x n
                TUblasSparseSpace<complex>::SetToZero(r_tmp_Vn_par);

                // initialize a separate solver for each thread
                typename ComplexLinearSolverType::Pointer compl_solver_par = LinearSolverFactoryType().Create(compl_solve_params);
                compl_solver_par->Initialize( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b);

                // start parallel solve
                // dynamic behavior such that threads that finish earlier won't be idle for too long
                #pragma omp for schedule(dynamic)
                for(size_t i=0; i < n_sampling_points/2; ++i)
                {
                    // intermediate result
                    r_tmp_Vn_par = std::pow( mSamplingPoints(2*i), 2.0 ) * r_M + mSamplingPoints(2*i) * r_D + r_K;

                    // solve
                    compl_solver_par->Solve( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b); // Ax = b, solve for x

                    // write the result to the correct positions
                    column(r_Vr_dense, 2*i) = real(r_tmp_Vr_col_par);
                    column(r_Vr_dense, 2*i+1) = imag(r_tmp_Vr_col_par);
                }
            } // end of parallel part

            std::cout<<"parallel part passed - iter "<<iter<<std::endl;

            double end_solve_par = OpenMPUtils::GetCurrentTime();
            std::cout<<"-- par solve: "<<end_solve_par-start_solve_par<<std::endl;


            // orthogonalize V
            mQR_decomposition.compute( system_size, reduced_system_size, &(r_Vr_dense)(0,0) );
            mQR_decomposition.compute_q();

            for(size_t i=0; i < system_size; ++i)
            {
                for(size_t j=0; j < reduced_system_size; ++j)
                {
                    r_Vr_sparse(i,j) = mQR_decomposition.Q(i,j);

                }
            }

            iter++;
        } // - end of loop

        // project the input vector onto the reduced space
        axpy_prod(r_b_tmp, r_Vr_sparse, r_b_reduced);  // br = V' b


        double end_time = OpenMPUtils::GetCurrentTime();
        double duration = end_time - start_time;

        KRATOS_INFO("IRKA") << "Completed in " << duration << " seconds" << " after " << iter << " iterations" << std::endl;
  
        KRATOS_INFO_IF("IRKA", iter >= max_iter) << "Maximum number of iterations reached!"  << std::endl;
  
        std::cout << "MOR offline solve finished" << std::endl;




        // write the matrices into files in order to do post-processing
        // TODO: this is a quick dirty output; build another function/file for the post-processing step

        // store reduced matrices
        std::stringstream matrix_market_m_reduced;
        matrix_market_m_reduced << "M_reduced" << ".mm";
        TDenseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_m_reduced.str()).c_str(), r_M_reduced, false);

        std::stringstream matrix_market_k_reduced;
        matrix_market_k_reduced << "K_reduced" << ".mm";
        TDenseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_k_reduced.str()).c_str(), r_K_reduced, false);

        std::stringstream matrix_market_d_reduced;
        matrix_market_d_reduced << "D_reduced" << ".mm";
        TDenseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_d_reduced.str()).c_str(), r_D_reduced, false);

        std::stringstream matrix_market_b_reduced;
        matrix_market_b_reduced << "b_reduced" << ".mm";
        TDenseSpace::WriteMatrixMarketVector((char *)(matrix_market_b_reduced.str()).c_str(), r_b_reduced);
*/

        // store full system size matrices (double)
/*        
        std::stringstream matrix_market_m_full;
        matrix_market_m_full << "M_full" << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_m_full.str()).c_str(), r_M_tmp, false);
        std::stringstream matrix_market_k_full;
        matrix_market_k_full << "K_full" << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_k_full.str()).c_str(), r_K_tmp, false);
        std::stringstream matrix_market_d_full;
        matrix_market_d_full << "D_full" << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_d_full.str()).c_str(), r_D_tmp, false);
        
*/

/*
        // needed for post processing
        std::stringstream matrix_market_b_full;
        matrix_market_b_full << "b_full" << ".mm";
        TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_b_full.str()).c_str(), r_b_tmp);


        // store the transformation matrix
        std::stringstream matrix_market_vr;
        matrix_market_vr << "Vr" << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_vr.str()).c_str(), r_Vr_sparse, false);

   
        // store the time
        std::ofstream time_file;
        time_file.open("times.txt", std::ios::out | std::ios::app);
        time_file << duration << " seconds\n";
        time_file.close();

*/

		return true;

        KRATOS_CATCH("");
    }

    void FinalizeSolutionStep() override
    {
        std::cout  << "IRKA finalize!!\n";

        auto& r_force_vector_reduced = this->GetRHSr();
        auto& r_output_vector = this->GetOutputVector();
        auto& r_output_vector_r = this->GetOVr();
        auto& r_basis = this->GetBasis();

        const size_t system_size = r_basis.size1();
        const size_t reduced_system_size = r_basis.size2();

        BuiltinTimer system_projection_time;

        r_force_vector_reduced.resize( reduced_system_size, false);
        r_force_vector_reduced = prod( this->GetSystemVector(), r_basis );

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

    vector< complex > mSamplingPoints;
    double mTolerance;
    size_t mMaxIter;
    // QR<double, row_major> mQR_decomposition;
    // EigenSystemSolverType* mpLinearEigenSolver;

    // TSolutionMatrixPointerType mpM; // mass matrix
    // TSolutionMatrixPointerType mpK; // stiffness matrix
    // TSolutionMatrixPointerType mpD; // damping matrix
    // TSolutionVectorPointerType mpb; // RHS vector


    /**
     * @brief Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    * @details Default = false
        - true  : Reforme at each time step
        - false : Form just one (more efficient)
     */
    // bool mReformDofSetAtEachStep;

    // bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

    // bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods


    // tentative sort function, immitating cplxpair of Matlab
    // maybe needs to be optimized, but the size of mSamplingPoints, i.e. the reduced system size, is not so large anyway
    void sort_conjugate_pairs(vector<complex>& v)
    {
        vector<complex> tmp_v = v;
        size_t length_of_v = v.size();
        double min_real_part = v(0).real();
        double related_imag_part = v(0).imag();
        int idx_min_real_part = 0;

        // run through half of the vector for comparison, then it is sorted
        for(size_t k=0; k<length_of_v-1; k+=2)
        {

            //set the initial values
            idx_min_real_part = k;
            min_real_part = v(idx_min_real_part).real();
            related_imag_part = v(idx_min_real_part).imag();
            tmp_v = v;

            // find current min real part in v
            for(size_t i=k; i<length_of_v; i++)
            {
                if(v(i).real() < min_real_part)
                {
                    min_real_part = v(i).real();
                    idx_min_real_part = i;
                    related_imag_part = v(i).imag();
                }
            }

            // swap the complex number of min real part with the first free entry
            v(k) = v(idx_min_real_part);
            v(idx_min_real_part) = tmp_v(k);
            tmp_v = v;

            // search for the related imaginary part
            for(size_t i=k+1; i<length_of_v; i++)
            {
                // compare with "+"" since the conjugate has the opposite sign
                if(abs(v(i).imag() + related_imag_part)<1e-4)
                {
                    // swap with the current first free entry
                    v(k+1) = v(i);
                    v(i) = tmp_v(k+1);

                    // change sign so that the number with a negative imaginary part comes first in the pair
                    if(related_imag_part>0)
                    {
                        v(k)   = std::conj(v(k));
                        v(k+1) = std::conj(v(k+1));
                    }
                    break;
                }
                
            }

        }

    }



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
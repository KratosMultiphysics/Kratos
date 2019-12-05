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
// #include "linear_solvers/iterative_solver.h"
#include "custom_solvers/gen_eigensystem_solver.h"

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
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class MorSecondOrderIRKAStrategy
    : public MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver >
{
  
  Parameters mParam;

  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorSecondOrderIRKAStrategy);

    typedef TUblasSparseSpace<complex> ComplexSparseSpaceType;
    typedef TUblasDenseSpace<complex> ComplexDenseSpaceType;

    typedef MorOfflineSecondOrderStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef LinearSolverFactory<ComplexSparseSpaceType,  ComplexDenseSpaceType> LinearSolverFactoryType;

    typedef SystemMatrixBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef TDenseSpace DenseSpaceType;

    typedef typename TDenseSpace::VectorType TDenseVectorType;

    typedef typename TDenseSpace::MatrixType TDenseMatrixType;

    typedef typename TDenseSpace::MatrixPointerType TDenseMatrixPointerType;

    typedef typename BaseType::TSchemeType TSchemeType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;


    typedef typename ComplexSparseSpaceType::MatrixType TSolutionMatrixType;

    typedef typename ComplexSparseSpaceType::MatrixPointerType TSolutionMatrixPointerType;

    typedef typename ComplexSparseSpaceType::VectorType TSolutionVectorType;

    typedef typename ComplexSparseSpaceType::VectorPointerType TSolutionVectorPointerType;

    typedef ComplexSparseSpaceType TSolutionSpace;

    typedef LinearSolver<ComplexSparseSpaceType, ComplexDenseSpaceType> ComplexLinearSolverType;

    typedef SystemMatrixBuilderAndSolver< ComplexSparseSpaceType, ComplexDenseSpaceType, ComplexLinearSolverType > TComplexBuilderAndSolverType;


    typedef Preconditioner<SparseSpaceType, DenseSpaceType> PreconditionerType;
    typedef Reorderer<SparseSpaceType, DenseSpaceType> ReordererType;
    typedef IterativeSolver<SparseSpaceType, DenseSpaceType, PreconditionerType, ReordererType> IterativeSolverType;

    typedef LinearSolver<SparseSpaceType, DenseSpaceType> LinearSolverType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef GenEigensystemSolver< SparseSpaceType, LocalSpaceType, PreconditionerType, ReordererType> EigenSystemSolverType;



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
        typename TLinearSolver::Pointer pNewLinearSolver,
        Parameters param,
        //typename ComplexLinearSolverType::Pointer pNewComplexLinearSolver,
        //typename TLinearSolver::Pointer pNewLinearEigenSolver,
        vector< double > samplingPoints_real,
        vector< double > samplingPoints_imag,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pNewLinearSolver, MoveMeshFlag), mParam(param)
    {
        KRATOS_TRY;

        // Saving the scheme
        this->SetScheme(pScheme);

        // Setting up the default builder and solver
        this->SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer(
            new TBuilderAndSolverType(pNewLinearSolver)));

        // Saving the linear solver        
        this->SetLinearSolver(pNewLinearSolver);

        // Set flags to start correctly the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        // Tells to the builder and solver if the reactions have to be calculated or not
        this->GetBuilderAndSolver()->SetCalculateReactionsFlag(false);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        this->GetBuilderAndSolver()->SetReshapeMatrixFlag(false);

        // Set EchoLevel to the default value (only time is displayed)
        this->SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(0);

        // Set members
        //TODO: let pybind11 do the conversion
        mSamplingPoints.resize(samplingPoints_real.size());
        for(size_t i=0; i<samplingPoints_real.size(); i++)
        {
            mSamplingPoints(i) = complex( samplingPoints_real(i), samplingPoints_imag(i) );
        }

        //mpComplexLinearSolver = pNewComplexLinearSolver;
        //mpLinearEigenSolver = pNewLinearEigenSolver;

        mpM = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
        mpK = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
        mpD = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
        mpb = ComplexSparseSpaceType::CreateEmptyVectorPointer();


        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorSecondOrderIRKAStrategy() override
    {

        if (mpM != nullptr)
            TSolutionSpace::Clear(mpM);
        if (mpK != nullptr)
            TSolutionSpace::Clear(mpK);
        if (mpD != nullptr)
            TSolutionSpace::Clear(mpD);
        if (mpb != nullptr)
            TSolutionSpace::Clear(mpb);


        this->Clear();
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
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver();

        // initialize the full system size matrices                
        TSystemMatrixType& r_K_tmp = this->GetSystemMatrix();
        TSystemMatrixType& r_M_tmp = this->GetMassMatrix();
        TSystemMatrixType& r_D_tmp = this->GetDampingMatrix();
        TSystemVectorType& r_b_tmp = this->GetSystemVector();  // originally RHS, check if results are bad

        TSystemVectorType tmp(r_K_tmp.size1(), 0.0);

        //define references for the complex matrices
        TSolutionMatrixType& r_K = *mpK;
        TSolutionMatrixType& r_M = *mpM;
        TSolutionMatrixType& r_D = *mpD;
        TSolutionVectorType& r_b = *mpb;


        // build the matrices based on the used model
        p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), r_b_tmp);

        p_builder_and_solver->BuildStiffnessMatrix(p_scheme, BaseType::GetModelPart(), r_K_tmp, tmp);
        p_builder_and_solver->ApplyDirichletConditions(p_scheme, BaseType::GetModelPart(), r_K_tmp, tmp, r_b_tmp);  

        p_builder_and_solver->BuildMassMatrix(p_scheme, BaseType::GetModelPart(), r_M_tmp, tmp);
        p_builder_and_solver->ApplyDirichletConditionsForMassMatrix(p_scheme, BaseType::GetModelPart(), r_M_tmp);

        p_builder_and_solver->BuildDampingMatrix(p_scheme, BaseType::GetModelPart(), r_D_tmp, tmp);
        p_builder_and_solver->ApplyDirichletConditionsForDampingMatrix(p_scheme, BaseType::GetModelPart(), r_D_tmp);


        // transform the real matrices to complex types
        r_K = TSolutionMatrixType(r_K_tmp);
        r_M = TSolutionMatrixType(r_M_tmp);
        r_D = TSolutionMatrixType(r_D_tmp);
        r_b = TSolutionVectorType(r_b_tmp);


        // start timer (after building the full system)
        double start_time = OpenMPUtils::GetCurrentTime();


        // DEBUG: print matrix information 
        //this->EchoInfo(0);


        // set the system sizes
        const unsigned int system_size = p_builder_and_solver->GetEquationSystemSize(); // n
        const std::size_t n_sampling_points = mSamplingPoints.size(); // number of sampling points
        const std::size_t reduced_system_size = n_sampling_points; // r

        // sort the sampling points after their real part in increasing order
        // conjugate pairs stay together, with the minus element being the first
        this->sort_conjugate_pairs(mSamplingPoints);


        // store the sampling points for error calculations in the convergence loop for convergence check
        auto samplingPoints_old = mSamplingPoints;



        // initialize V (=W, due to symmetry in FEM applications)
        auto  Vr_dense_ptr = DenseSpaceType::CreateEmptyMatrixPointer();
        auto& r_Vr_dense   = *Vr_dense_ptr;
        DenseSpaceType::Resize(r_Vr_dense, system_size, reduced_system_size); // n x r
        //need to be dense for the QR-solver; if replaced by Gram Schmidt or similar only use sparse version

        auto& r_Vr_sparse = this->GetBasis();
        SparseSpaceType::Resize(r_Vr_sparse, system_size, reduced_system_size); // n x r
        TSparseSpace::SetToZero(r_Vr_sparse);
        // sparse version for the orthogonalized form

        // initialize helper variables for V
        // only in serial version
        // // auto  tmp_Vn_ptr = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
        // // auto& r_tmp_Vn   = *tmp_Vn_ptr;
        // // ComplexSparseSpaceType::Resize(r_tmp_Vn, system_size, system_size); // n x n
        // // TUblasSparseSpace<complex>::SetToZero(r_tmp_Vn);
        
        // only in serial version
        // // auto  tmp_Vr_col_ptr = ComplexDenseSpaceType::CreateEmptyVectorPointer();
        // // auto& r_tmp_Vr_col   = *tmp_Vr_col_ptr;
        // // ComplexDenseSpaceType::Resize(r_tmp_Vr_col, system_size); // n x 1
        // // ComplexDenseSpaceType::Set(r_tmp_Vr_col, 0.0); // set vector to zero

        // initialize the complex solver
        //mpComplexLinearSolver->Initialize( r_tmp_Vn, r_tmp_Vr_col, r_b);    //check if init also increases speed in loop


        Parameters compl_solve_params(R"(
        {
            "solver_type": "skyline_lu_complex"
        }
        )");

        //std::cout<<"max threads: "<<omp_get_max_threads()<<std::endl;
        
        double start_solve_par_first = OpenMPUtils::GetCurrentTime();
        // build V
        #pragma omp parallel
        {

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
                int tid = omp_get_thread_num();
                #pragma omp critical
                std::cout<<" --- Thread "<<tid<<" sampling point "<<mSamplingPoints(2*i)<<std::endl;
                // intermediate result
                r_tmp_Vn_par = std::pow( mSamplingPoints(2*i), 2.0 ) * r_M + mSamplingPoints(2*i) * r_D + r_K;

                // solve 
                double solve_step_first = OpenMPUtils::GetCurrentTime();
                compl_solver_par->Solve( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b); // Ax = b, solve for x
                double end_solve_step_first = OpenMPUtils::GetCurrentTime();
                #pragma omp critical
                std::cout<<" --- Thread " <<tid<<" solve step: "<<end_solve_step_first-solve_step_first<<std::endl;

                // write the result to the correct positions
                double write_step_first = OpenMPUtils::GetCurrentTime();
                column(r_Vr_dense, 2*i)   = real(r_tmp_Vr_col_par);
                column(r_Vr_dense, 2*i+1) = imag(r_tmp_Vr_col_par);
                double end_write_step_first = OpenMPUtils::GetCurrentTime();
                #pragma omp critical
                std::cout<<" --- Thread " <<tid<<" write step: "<<end_write_step_first-write_step_first<<std::endl;
            }
        } // end of parallel part

        //std::cout<<"parallel part passed"<<std::endl;
        double end_solve_par_first = OpenMPUtils::GetCurrentTime();
            std::cout<<"-- par solve first: "<<end_solve_par_first-start_solve_par_first<<std::endl;


        //orthogonalize V
        // TODO: if time left, replace by some orthogonlization method
        mQR_decomposition.compute( system_size, reduced_system_size, &(r_Vr_dense)(0,0) );
        mQR_decomposition.compute_q();

        //#pragma omp parallel for schedule(dynamic)
        for(size_t i=0; i < system_size; ++i)
        {
            //#pragma omp critical (otherame)
            for(size_t j=0; j < reduced_system_size; ++j)
            {
                r_Vr_sparse(i,j) = mQR_decomposition.Q(i,j);
            }
        }


        // initialize the reduced matrices
        TSystemVectorType& r_b_reduced = this->GetRHSr();
        TSystemMatrixType& r_K_reduced = this->GetKr();
        TSystemMatrixType& r_M_reduced = this->GetMr();
        TSystemMatrixType& r_D_reduced = this->GetDr();

        // needs to be resized, otherwise segfault
        SparseSpaceType::Resize(r_M_reduced, reduced_system_size, reduced_system_size); // r x r
        SparseSpaceType::Resize(r_K_reduced, reduced_system_size, reduced_system_size); // r x r
        SparseSpaceType::Resize(r_D_reduced, reduced_system_size, reduced_system_size); // r x r
        SparseSpaceType::Resize(r_b_reduced, reduced_system_size);                      // r x 1



        // helper for auxiliary products
        //TSystemMatrixType T; //TODO: define actual size for first temp prod (rxn)

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

            double start_projections = OpenMPUtils::GetCurrentTime();


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


            double end_projections = OpenMPUtils::GetCurrentTime();
            std::cout<<"-- projections: "<<end_projections-start_projections<<std::endl;



//            double start_ev_setup = OpenMPUtils::GetCurrentTime();

            // linearization of the quadratic eigenvalue problem to be solved using a GEP solver
            subrange(r_L1, 0, reduced_system_size, reduced_system_size, 2*reduced_system_size) = id_r;
            subrange(r_L1, reduced_system_size, 2*reduced_system_size, 0, reduced_system_size) = -1.0*r_K_reduced;
            subrange(r_L1, reduced_system_size, 2*reduced_system_size, reduced_system_size, 2*reduced_system_size) = -1.0*r_D_reduced;

            subrange(r_L2, 0, reduced_system_size, 0, reduced_system_size) = id_r;
            subrange(r_L2, reduced_system_size, 2*reduced_system_size, reduced_system_size, 2*reduced_system_size) = r_M_reduced;


            // double end_ev_setup = OpenMPUtils::GetCurrentTime();
            // std::cout<<"-- ev setup: "<<end_ev_setup-start_ev_setup<<std::endl;


//            double start_ev_solver = OpenMPUtils::GetCurrentTime();

            // output of this Solver is real (real parts and imaginary parts alternating)
            //test_ev_solver_yay->Solve(
            mpLinearEigenSolver->Solve(
                r_L1,
                r_L2,
                Eigenvalues,
                Eigenvectors
            );

            // double end_ev_solver = OpenMPUtils::GetCurrentTime();
            // std::cout<<"-- ev itself: "<<end_ev_solver-start_ev_solver<<std::endl;

            //store the Eigenvalues in a complex vector
            #pragma omp parallel for schedule(dynamic)
            for(size_t ii = 0; ii<2*reduced_system_size; ii++)
            {
                lam_2r(ii) = complex( Eigenvalues(2*ii), Eigenvalues(2*ii+1) );
            }

    //        double start_conjug = OpenMPUtils::GetCurrentTime();

            // mirror the Eigenvalues and sort them after their real part in increasing order
            // conjugate pairs stay together, with the minus element being the first
            lam_2r *= -1;
            this->sort_conjugate_pairs(lam_2r);

            // double end_conjug = OpenMPUtils::GetCurrentTime();
            // std::cout<<"-- conjugagte sort: "<<end_conjug-start_conjug<<std::endl;


   //         double start_error_calc = OpenMPUtils::GetCurrentTime();

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


            // double end_error_calc = OpenMPUtils::GetCurrentTime();
            // std::cout<<"-- error_calc: "<<end_error_calc-start_error_calc<<std::endl;

            // if the tolerance is already met, break here to avoid unnecessary computations
            if(err < tol) {break;}


            double start_solve_par = OpenMPUtils::GetCurrentTime();

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

    //###################   for serial def
                    // // r_tmp_Vn = std::pow( mSamplingPoints(2*i), 2.0 ) * r_M + mSamplingPoints(2*i) * r_D + r_K;
                    // // mpComplexLinearSolver->Solve( r_tmp_Vn, r_tmp_Vr_col, r_b); // Ax = b, solve for x

                    // // column(r_Vr_dense, 2*i) = real(r_tmp_Vr_col);
                    // // column(r_Vr_dense, 2*i+1) = imag(r_tmp_Vr_col);

                }
            } // end of parallel part

                    std::cout<<"parallel part passed - iter "<<iter<<std::endl;

            double end_solve_par = OpenMPUtils::GetCurrentTime();
            std::cout<<"-- par solve: "<<end_solve_par-start_solve_par<<std::endl;


//            double start_qr_orth = OpenMPUtils::GetCurrentTime();

            // orthogonalize V
            mQR_decomposition.compute( system_size, reduced_system_size, &(r_Vr_dense)(0,0) );
            mQR_decomposition.compute_q();

            //#pragma omp parallel for schedule(dynamic)
            for(size_t i=0; i < system_size; ++i)
            {
                for(size_t j=0; j < reduced_system_size; ++j)
                {
                    r_Vr_sparse(i,j) = mQR_decomposition.Q(i,j);

                }
            }

            // double end_qr_orth = OpenMPUtils::GetCurrentTime();
            // std::cout<<"-- qr orthog.: "<<end_qr_orth-start_qr_orth<<std::endl;

            iter++;
        } // - end of loop

        // project the input vector onto the reduced space
        //r_b_reduced = prod( trans(r_Vr_sparse), r_b_tmp);
        axpy_prod(r_b_tmp, r_Vr_sparse, r_b_reduced);  // br = V' b


        double end_time = OpenMPUtils::GetCurrentTime();
        double duration = end_time - start_time;

        KRATOS_INFO("IRKA") << "Completed in " << duration << " seconds" << " after " << iter << " iterations" << std::endl;
  
        KRATOS_INFO_IF("IRKA", iter >= max_iter) << "Maximum number of iterations reached!"  << std::endl;
  
        //KRATOS_WATCH(r_M_reduced)

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
        

        std::stringstream matrix_market_b_full;
        matrix_market_b_full << "b_full" << ".mm";
        TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_b_full.str()).c_str(), r_b_tmp);


        // store the transformation matrix
        std::stringstream matrix_market_vr;
        matrix_market_vr << "Vr" << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_vr.str()).c_str(), r_Vr_sparse, false);
*/
   
        // store the time
        std::ofstream time_file;
        time_file.open("times.txt", std::ios::out | std::ios::app);
        time_file << duration << " seconds\n";
        time_file.close();




		return true;

        KRATOS_CATCH("");
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
    QR<double, row_major> mQR_decomposition;
    //typename TLinearSolver::Pointer mpLinearEigenSolver;
    EigenSystemSolverType* mpLinearEigenSolver;
    //ComplexLinearSolverType::Pointer mpComplexLinearSolver;

    TSolutionMatrixPointerType mpM; // mass matrix
    TSolutionMatrixPointerType mpK; // stiffness matrix
    TSolutionMatrixPointerType mpD; // damping matrix
    TSolutionVectorPointerType mpb; // RHS vector



    /**
     * @brief Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    * @details Default = false
        - true  : Reforme at each time step
        - false : Form just one (more efficient)
     */
    bool mReformDofSetAtEachStep;

    bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

    bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

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
        // needs to be optimized, but the size of mSamplingPoints, i.e. the reduced system size, is not so large anyway
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

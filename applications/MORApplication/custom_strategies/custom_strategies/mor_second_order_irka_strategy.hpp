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
//#include "includes/kratos_parameters.h"
#include "factories/linear_solver_factory.h"

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
        typename ComplexLinearSolverType::Pointer pNewComplexLinearSolver,
        typename TLinearSolver::Pointer pNewLinearEigenSolver,
        vector< double > samplingPoints_real,
        vector< double > samplingPoints_imag,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pNewLinearSolver, MoveMeshFlag)
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

        mpComplexLinearSolver = pNewComplexLinearSolver;
        mpLinearEigenSolver = pNewLinearEigenSolver;

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

        // initialize helper variables for V
        auto  tmp_Vn_ptr = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_tmp_Vn   = *tmp_Vn_ptr;
        ComplexSparseSpaceType::Resize(r_tmp_Vn, system_size, system_size); // n x n
        TUblasSparseSpace<complex>::SetToZero(r_tmp_Vn);
        
        auto  tmp_Vr_col_ptr = ComplexDenseSpaceType::CreateEmptyVectorPointer();
        auto& r_tmp_Vr_col   = *tmp_Vr_col_ptr;
        ComplexDenseSpaceType::Resize(r_tmp_Vr_col, system_size); // n x 1
        ComplexDenseSpaceType::Set(r_tmp_Vr_col, 0.0); // set vector to zero

        // initialize the complex solver
        mpComplexLinearSolver->Initialize( r_tmp_Vn, r_tmp_Vr_col, r_b);


std::cout<<" -- test factory --"<<std::endl;

Parameters comp_solv_params(R"(
{
    "solver_type": "skyline_lu_complex"
}
)");

//LinearSolverFactoryType::Pointer ls_fact;
//typename ComplexLinearSolverType::Pointer test_solver = ls_fact->Create(comp_solv_params); 
typename ComplexLinearSolverType::Pointer test_solver = LinearSolverFactoryType().Create(comp_solv_params); 
std::cout<<"factory ready"<<std::endl;

//typename ComplexLinearSolverType::Pointer test_solver = LinearSolverFactory<ComplexSparseSpaceType,ComplexDenseSpaceType>::Create(comp_solv_params); 
test_solver->Initialize( r_tmp_Vn, r_tmp_Vr_col, r_b);
std::cout<<"init passed"<<std::endl;
r_tmp_Vn = std::pow( mSamplingPoints(2*0), 2.0 ) * r_M + mSamplingPoints(2*0) * r_D + r_K;
test_solver->Solve( r_tmp_Vn, r_tmp_Vr_col, r_b);
std::cout<<"test solve passed"<<std::endl;

        //typename ComplexLinearSolverType::Pointer solver_arr[4];
        ComplexLinearSolverType solver_arr_direct[4];
        typename ComplexLinearSolverType::Pointer solver_arr[4];
        //typename TComplexBuilderAndSolverType::Pointer solver_arr[4];
        for(int i=0; i<4; i++){
            solver_arr[i] = LinearSolverFactoryType().Create(comp_solv_params); 
            //solver_arr[i] = TComplexBuilderAndSolverType::Pointer(new TComplexBuilderAndSolverType(mpComplexLinearSolver));
            //solver_arr_direct[i] = *mpComplexLinearSolver; // check if really different addresses due to copy constructor
            //solver_arr[i] = &solver_arr_direct[i];
            //solver_arr[i] = ComplexLinearSolverType::Pointer( new ComplexLinearSolverType(mpComplexLinearSolver)); // check if really different addresses due to copy constructor
            std::cout<<"addr solver "<<i<<": "<<solver_arr[i]<<std::endl;   // no... 
            std::cout<<"addrp solver "<<i<<": "<<&solver_arr[i]<<std::endl;
            // // // std::cout<<"addr solver "<<i<<": "<<solver_arr_direct[i]<<std::endl;   // no... 
            // // // std::cout<<"addrp solver "<<i<<": "<<&solver_arr_direct[i]<<std::endl;
        }


        // ComplexLinearSolverType::Pointer test_solver_p = &solver_arr_direct[0]; // geht nicht
        //solver_arr_direct[0].Solve( r_tmp_Vn, r_tmp_Vr_col, r_b);  // calling linear solver base class
        //test_solver_p->Solve( r_tmp_Vn, r_tmp_Vr_col, r_b);
        std::cout<<"passed single call"<<std::endl;

        std::cout<<"max threads: "<<omp_get_max_threads()<<std::endl;
        //omp_set_num_threads(4);

        // build V
        #pragma omp parallel for default(shared) schedule(static,1)
        //#pragma omp parallel for default(shared) schedule(dynamic)
        for(size_t i=0; i < n_sampling_points/2; ++i)
        {

            auto  tmp_Vr_col_ptr_par = ComplexDenseSpaceType::CreateEmptyVectorPointer();
            auto& r_tmp_Vr_col_par   = *tmp_Vr_col_ptr_par;
            ComplexDenseSpaceType::Resize(r_tmp_Vr_col_par, system_size); // n x 1
            ComplexDenseSpaceType::Set(r_tmp_Vr_col_par, 0.0); // set vector to zero

            auto  tmp_Vn_ptr_par = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
            auto& r_tmp_Vn_par   = *tmp_Vn_ptr_par;
            ComplexSparseSpaceType::Resize(r_tmp_Vn_par, system_size, system_size); // n x n
            TUblasSparseSpace<complex>::SetToZero(r_tmp_Vn_par);

            auto  b_par = ComplexSparseSpaceType::CreateEmptyVectorPointer();
            auto& r_b_par   = *b_par;
            ComplexSparseSpaceType::Resize(r_b_par, system_size); // n x 1
            //r_b_par = r_b; // or actually copy, look up
            noalias(r_b_par) = r_b;  // leider auch mit copy derselbe Fehler...

            //typename ComplexLinearSolverType::Pointer tmp_CLsolver;

            r_tmp_Vn_par = std::pow( mSamplingPoints(2*i), 2.0 ) * r_M + mSamplingPoints(2*i) * r_D + r_K;

            //#pragma omp critical (somename)
            //typename ComplexLinearSolverType::Pointer tmp_CLsolver = ComplexLinearSolverType::Pointer(new ComplexLinearSolverType(mpComplexLinearSolver));
            //typename ComplexLinearSolverType::Pointer tmp_CLsolver = mpComplexLinearSolver;
            //#pragma omp critical (buildnew)
            //tmp_CLsolver = mpComplexLinearSolver;
            //tmp_CLsolver->Initialize( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b);
            //tmp_CLsolver = new ComplexLinearSolverType(mpComplexLinearSolver);
            //tmp_CLsolver->Solve( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b); // Ax = b, solve for x
            
            //#pragma omp critical (somename)
            //mpComplexLinearSolver->Solve( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b); // Ax = b, solve for x

            //int tid = omp_get_thread_num();
            //solver_arr[tid]->Solve( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b_par); // Ax = b, solve for x

            //solver_arr_direct[tid].Solve( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b_par); // Ax = b, solve for x
            //solver_arr[tid]->GetLinearSystemSolver()->Solve( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b_par); // Ax = b, solve for x

            typename ComplexLinearSolverType::Pointer test_solver_par = LinearSolverFactoryType().Create(comp_solv_params); 
            test_solver_par->Solve( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b_par); // Ax = b, solve for x


            column(r_Vr_dense, 2*i) = real(r_tmp_Vr_col_par);
            column(r_Vr_dense, 2*i+1) = imag(r_tmp_Vr_col_par);
        }

        std::cout<<"parallel part passed"<<std::endl;


        //orthogonalize V
        // TODO: if time left, replace by some orthogonlization method
        mQR_decomposition.compute( system_size, reduced_system_size, &(r_Vr_dense)(0,0) );
        mQR_decomposition.compute_q();

        //auto& r_tmp_Q = mQR_decomposition.Q;

        //#pragma omp parallel for schedule(dynamic)
        for(size_t i=0; i < system_size; ++i)
        {
            //#pragma omp critical (otherame)
            for(size_t j=0; j < reduced_system_size; ++j)
            {
                r_Vr_sparse(i,j) = mQR_decomposition.Q(i,j);
                //r_Vr_sparse(i,j) = r_tmp_Q(i,j);

            }
        }


        // initialize the reduced matrices
        auto& r_b_reduced = this->GetRHSr();
        auto& r_K_reduced = this->GetKr();
        auto& r_M_reduced = this->GetMr();
        auto& r_D_reduced = this->GetDr();

        // helper for auxiliary products
        TSystemMatrixType T; //TODO: define actual size for first temp prod (rxn)

        // initialize linearization matrices
        auto L1_ptr = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_L1 = *L1_ptr;
        SparseSpaceType::Resize(r_L1, 2*reduced_system_size, 2*reduced_system_size); // 2r x 2r

        auto L2_ptr = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_L2 = *L2_ptr;
        SparseSpaceType::Resize(r_L2, 2*reduced_system_size, 2*reduced_system_size); // 2r x 2r

        // identity matrix of size r x r
        identity_matrix<double> id_r (reduced_system_size); 

        // initialize output for Eigensolver
        TDenseVectorType Eigenvalues;
        TDenseMatrixType Eigenvectors;  //not needed, but need to be set to call the Eigensolver
        vector<complex> lam_2r (2*reduced_system_size);

        // initialize helpers for the convergence check
        vector<double> abs_val_old (reduced_system_size);
        vector<double> abs_val_new (reduced_system_size);


        int iter = 1;
        int max_iter = 100;
        double err = 100.0;
        double tol = 1e-8;
        while(iter < max_iter && err > tol)
        {

            double start_projections = OpenMPUtils::GetCurrentTime();

            // projections onto the reduced space
  //          #pragma omp parallel sections
    //        {
      //          #pragma omp section
        //        {
                    // for(size_t i=0; i<system_size; i++){
                    //     axpy_prod(column(r_M_tmp,i), r_Vr_sparse, subrange(T,0,reduced_system_size,i,i), true);
                    // }

                    double start_axpy = OpenMPUtils::GetCurrentTime();



                        auto  tmp_T_par = SparseSpaceType::CreateEmptyMatrixPointer();
                        auto& r_tmp_T_par   = *tmp_T_par;
                        SparseSpaceType::Resize(r_tmp_T_par, reduced_system_size, system_size); // r x n
                        TSparseSpace::SetToZero(r_tmp_T_par);


                    #pragma omp parallel
                    {


        // //               auto Vr_sparse_par = SparseSpaceType::CreateEmptyMatrixPointer();
        // //               auto& r_Vr_sparse_par = *Vr_sparse_par;
        // // SparseSpaceType::Resize(r_Vr_sparse_par, system_size, reduced_system_size); // n x r
        // // TSparseSpace::SetToZero(r_Vr_sparse_par);

        // //                 noalias(r_Vr_sparse_par) = r_Vr_sparse;
                        //r_Vr_sparse_par = r_Vr_sparse;

                        //const double& test_V = r_Vr_sparse_par(0,0);

                        //#pragma omp critical
                        //std::cout<<"     [V spars par id] = "<<test_V<<std::endl;
                        //std::cout<<"    [V sparse id] = "<<&r_Vr_sparse<<" --- [V sparse par id] = "<<&r_Vr_sparse_par<<std::endl;


                      auto Vr_sparse_par = SparseSpaceType::CreateEmptyMatrixPointer();
                      auto& r_Vr_sparse_par = *Vr_sparse_par;
                        //auto& r_Vr_sparse_par = this->GetBasis();
        SparseSpaceType::Resize(r_Vr_sparse_par, system_size, reduced_system_size); // n x r
        TSparseSpace::SetToZero(r_Vr_sparse_par);

                    #pragma omp critical
                     subrange(r_Vr_sparse_par,0,system_size,0,reduced_system_size) = r_Vr_sparse; 



                        // auto  tmp_T_par = SparseSpaceType::CreateEmptyMatrixPointer();
                        // auto& r_tmp_T_par   = *tmp_T_par;
                        // SparseSpaceType::Resize(r_tmp_T_par, reduced_system_size, system_size); // r x n
                        // TSparseSpace::SetToZero(r_tmp_T_par);    

                        auto  m_col = SparseSpaceType::CreateEmptyVectorPointer();
                        auto& r_m_col   = *m_col;
                        SparseSpaceType::Resize(r_m_col, system_size); // n x 1

                        auto  t_par = SparseSpaceType::CreateEmptyVectorPointer();
                        auto& r_t_par   = *t_par;
                        SparseSpaceType::Resize(r_t_par, reduced_system_size); // r x 1


                    //#pragma omp parallel for shared(r_tmp_T_par)  // shared default anyway, but to be safe...
                    //#pragma omp parallel for schedule(dynamic)
                    #pragma omp for
                    for(size_t i=0; i<system_size; i++){

                        // auto  tmp_T_par = SparseSpaceType::CreateEmptyMatrixPointer();
                        // auto& r_tmp_T_par   = *tmp_T_par;
                        // SparseSpaceType::Resize(r_tmp_T_par, reduced_system_size, system_size); // r x n
                        // TSparseSpace::SetToZero(r_tmp_T_par);

                        // // auto  m_col = SparseSpaceType::CreateEmptyVectorPointer();
                        // // auto& r_m_col   = *m_col;
                        // // SparseSpaceType::Resize(r_m_col, system_size); // n x 1

                        // // auto  t_par = SparseSpaceType::CreateEmptyVectorPointer();
                        // // auto& r_t_par   = *t_par;
                        // // SparseSpaceType::Resize(r_t_par, reduced_system_size); // r x 1


// INEFFICIENT!!!!! -JUST FOR TESTING
        // //               auto Vr_sparse_par = SparseSpaceType::CreateEmptyMatrixPointer();
        // //               auto& r_Vr_sparse_par = *Vr_sparse_par;
        // //                 //auto& r_Vr_sparse_par = this->GetBasis();
        // // SparseSpaceType::Resize(r_Vr_sparse_par, system_size, reduced_system_size); // n x r
        // // TSparseSpace::SetToZero(r_Vr_sparse_par);

        // //             #pragma omp critical
        // //              subrange(r_Vr_sparse_par,0,system_size,0,reduced_system_size) = r_Vr_sparse; 
                    //#pragma omp critical
                    //std::cout<<"     --- copy passed"<<std::endl;




                    // test w = trans(A) * u  im ersten Schritt  - funktioniert
                    //#pragma omp parallel for schedule(dynamic)
                    //for(size_t i=0; i<system_size; i++){
                        //std::cout<<"     test_axpy_00"<<std::endl;
                        r_m_col = column(r_M_tmp,i);
                        //std::cout<<"     test_axpy_01"<<std::endl;
                        //#pragma omp critical (axpy)  //damit passt es, aber gleiches Problem wie beim Solver
                       // axpy_prod(r_m_col,r_Vr_sparse,r_t_par,true);
                        axpy_prod(r_m_col,r_Vr_sparse_par,r_t_par,true);
                       // std::cout<<"     test_axpy_02"<<std::endl;
                        column(r_tmp_T_par,i) = r_t_par;
                       // std::cout<<"     test_axpy_03"<<std::endl;

                    //    if(i==system_size-1){
                    //         #pragma omp barrier
                    //         r_M_reduced = prod( r_tmp_T_par, r_Vr_sparse );
                    //    }
                    }


                    }  // omp parallel
                    r_M_reduced = prod( r_tmp_T_par, r_Vr_sparse );
                    //r_M_reduced = prod( r_tmp_T_par, r_Vr_sparse_par );


                   
                        // // auto  t_par_dense = DenseSpaceType::CreateEmptyVectorPointer();
                        // // auto& r_t_par_dense   = *t_par_dense;
                        // // DenseSpaceType::Resize(r_t_par_dense, reduced_system_size); // r x 1

                        // vector<double> t_par_dense; 
                        // t_par_dense.resize(reduced_system_size);


// //                     for(size_t i=0; i<reduced_system_size; i++){
// //                         r_m_col = column(r_Vr_sparse,i);
// //                         axpy_prod(r_tmp_T_par, r_m_col, r_t_par_dense, true);
// //                         column(r_M_reduced,i) = r_t_par_dense;
// // //                        axpy_prod(r_tmp_T_par, r_m_col, t_par_dense, true);
// //                     }

                    //axpy_prod(r_tmp_T_par, r_M_tmp, r_M_reduced);  // - Speicherfehler



                    //std::cout<<"     test_axpy_04"<<std::endl;

                    // test C = A*B im zweiten Schritt  - Speicherfehler

                    // std::cout<<"     test_axpy_00"<<std::endl;
                    // r_tmp_T_par = prod( trans(r_Vr_sparse), r_M_tmp );
                    // std::cout<<"     test_axpy_01"<<std::endl;
                    // axpy_prod(r_tmp_T_par, r_Vr_sparse, r_M_reduced, true);
                    // std::cout<<"     test_axpy_02"<<std::endl;

                    double end_axpy = OpenMPUtils::GetCurrentTime();
            std::cout<<"       -- axpy: "<<end_axpy-start_axpy<<std::endl;

                    KRATOS_WATCH(r_M_reduced);

                    // different temp vars for the different sizes may be better, improve performance
                    // use OMP utils to thest which function is better, which takes the most time (trans, prod...)

// double start_old_T = OpenMPUtils::GetCurrentTime();

                    // // //  T = prod( trans(r_Vr_sparse), r_M_tmp );
                    // // //  r_M_reduced = prod( T, r_Vr_sparse );

// double end_old_T = OpenMPUtils::GetCurrentTime();
//             std::cout<<"     -- old T: "<<end_old_T-start_old_T<<std::endl;

//                     KRATOS_WATCH(r_M_reduced);

      //          }

     //           #pragma omp section
      //          {
                     T = prod( trans(r_Vr_sparse), r_D_tmp );
                     r_D_reduced = prod( T, r_Vr_sparse );
       //          }

      //           #pragma omp section
       //          {
                     T = prod( trans(r_Vr_sparse), r_K_tmp );
                     r_K_reduced = prod( T, r_Vr_sparse );
      //           }

      //       }


                // // // //    #pragma omp parallel for shared(r_tmp_T_par)  // shared default anyway, but to be safe...
                // // // //     for(size_t i=0; i<system_size; i++){

                // // // //         auto  m_col = SparseSpaceType::CreateEmptyVectorPointer();
                // // // //         auto& r_m_col   = *m_col;
                // // // //         SparseSpaceType::Resize(r_m_col, system_size); // n x 1

                // // // //         auto  t_par = SparseSpaceType::CreateEmptyVectorPointer();
                // // // //         auto& r_t_par   = *t_par;
                // // // //         SparseSpaceType::Resize(r_t_par, reduced_system_size); // r x 1

                // // // //         r_m_col = column(r_D_tmp,i);
                // // // //         axpy_prod(r_m_col,r_Vr_sparse,r_t_par,true);
                // // // //         column(r_tmp_T_par,i) = r_t_par;
                // // // //     }
                // // // //     r_D_reduced = prod( r_tmp_T_par, r_Vr_sparse );



                // // // //    #pragma omp parallel for shared(r_tmp_T_par)  // shared default anyway, but to be safe...
                // // // //     for(size_t i=0; i<system_size; i++){

                // // // //         auto  m_col = SparseSpaceType::CreateEmptyVectorPointer();
                // // // //         auto& r_m_col   = *m_col;
                // // // //         SparseSpaceType::Resize(r_m_col, system_size); // n x 1

                // // // //         auto  t_par = SparseSpaceType::CreateEmptyVectorPointer();
                // // // //         auto& r_t_par   = *t_par;
                // // // //         SparseSpaceType::Resize(r_t_par, reduced_system_size); // r x 1

                // // // //         r_m_col = column(r_K_tmp,i);
                // // // //         axpy_prod(r_m_col,r_Vr_sparse,r_t_par,true);
                // // // //         column(r_tmp_T_par,i) = r_t_par;
                // // // //     }
                // // // //     r_K_reduced = prod( r_tmp_T_par, r_Vr_sparse );





            double end_projections = OpenMPUtils::GetCurrentTime();
            std::cout<<"-- projections: "<<end_projections-start_projections<<std::endl;



            double start_ev_setup = OpenMPUtils::GetCurrentTime();

            // linearization of the quadratic eigenvalue problem to be solved using a GEP solver
            subrange(r_L1, 0, reduced_system_size, reduced_system_size, 2*reduced_system_size) = id_r;
            subrange(r_L1, reduced_system_size, 2*reduced_system_size, 0, reduced_system_size) = -1.0*r_K_reduced;
            subrange(r_L1, reduced_system_size, 2*reduced_system_size, reduced_system_size, 2*reduced_system_size) = -1.0*r_D_reduced;

            subrange(r_L2, 0, reduced_system_size, 0, reduced_system_size) = id_r;
            subrange(r_L2, reduced_system_size, 2*reduced_system_size, reduced_system_size, 2*reduced_system_size) = r_M_reduced;


            //TODO: hard code settings for the Eigensolver needed for IRKA
            //i.e.: 
            //number_of_eigenvalues = 2*reduced_system_size
            //compute_eigenvectors = false

            double end_ev_setup = OpenMPUtils::GetCurrentTime();
            std::cout<<"-- ev setup: "<<end_ev_setup-start_ev_setup<<std::endl;


            double start_ev_solver = OpenMPUtils::GetCurrentTime();

            // output of this Solver is real (real parts and imaginary parts alternating)
            mpLinearEigenSolver->Solve(
                r_L1,
                r_L2,
                Eigenvalues,
                Eigenvectors
            );

            double end_ev_solver = OpenMPUtils::GetCurrentTime();
            std::cout<<"-- ev itself: "<<end_ev_solver-start_ev_solver<<std::endl;

            //store the Eigenvalues in a complex vector
            #pragma omp parallel for schedule(dynamic)
            for(size_t ii = 0; ii<2*reduced_system_size; ii++)
            {
                lam_2r(ii) = complex( Eigenvalues(2*ii), Eigenvalues(2*ii+1) );
            }

            double start_conjug = OpenMPUtils::GetCurrentTime();

            // mirror the Eigenvalues and sort them after their real part in increasing order
            // conjugate pairs stay together, with the minus element being the first
            lam_2r *= -1;
            this->sort_conjugate_pairs(lam_2r);

            double end_conjug = OpenMPUtils::GetCurrentTime();
            std::cout<<"-- conjugagte sort: "<<end_conjug-start_conjug<<std::endl;


            double start_error_calc = OpenMPUtils::GetCurrentTime();

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

            sort(abs_val_old.begin(), abs_val_old.end());
            sort(abs_val_new.begin(), abs_val_new.end());

            err = norm_2(abs_val_new - abs_val_old) / norm_2(samplingPoints_old);


            double end_error_calc = OpenMPUtils::GetCurrentTime();
            std::cout<<"-- error_calc: "<<end_error_calc-start_error_calc<<std::endl;

            // if the tolerance is already met, break here to avoid unnecessary computations
            if(err < tol) {break;}


            double start_solve_par = OpenMPUtils::GetCurrentTime();

            // update the old sampling points
            samplingPoints_old = mSamplingPoints;

            // update V
            //#pragma omp parallel for default(shared) schedule(static,1) num_threads(4)
            #pragma omp parallel for default(shared) schedule(static,1)
            for(size_t i=0; i < n_sampling_points/2; ++i)
            {
                int tid = omp_get_thread_num();
                
                #pragma omp critical (pr)
                std::cout<<"    thread "<<tid<<" starts"<<std::endl;
                


                auto  tmp_Vr_col_ptr_par = ComplexDenseSpaceType::CreateEmptyVectorPointer();
                auto& r_tmp_Vr_col_par   = *tmp_Vr_col_ptr_par;
                ComplexDenseSpaceType::Resize(r_tmp_Vr_col_par, system_size); // n x 1
                ComplexDenseSpaceType::Set(r_tmp_Vr_col_par, 0.0); // set vector to zero

                auto  tmp_Vn_ptr_par = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
                auto& r_tmp_Vn_par   = *tmp_Vn_ptr_par;
                ComplexSparseSpaceType::Resize(r_tmp_Vn_par, system_size, system_size); // n x n
                TUblasSparseSpace<complex>::SetToZero(r_tmp_Vn_par);


                r_tmp_Vn_par = std::pow( mSamplingPoints(2*i), 2.0 ) * r_M + mSamplingPoints(2*i) * r_D + r_K;

                //#pragma omp critical (somename)
                //typename ComplexLinearSolverType::Pointer tmp_CLsolver = mpComplexLinearSolver;
                //tmp_CLsolver->Solve( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b); // Ax = b, solve for x
                //mpComplexLinearSolver->Solve( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b); // Ax = b, solve for x

                // // int tid = omp_get_thread_num();
                // // solver_arr[tid]->Solve( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b); // Ax = b, solve for x

                typename ComplexLinearSolverType::Pointer test_solver_par = LinearSolverFactoryType().Create(comp_solv_params); 
                test_solver_par->Solve( r_tmp_Vn_par, r_tmp_Vr_col_par, r_b); // Ax = b, solve for x

                column(r_Vr_dense, 2*i) = real(r_tmp_Vr_col_par);
                column(r_Vr_dense, 2*i+1) = imag(r_tmp_Vr_col_par);

                #pragma omp critical (pr)
                std::cout<<"    thread "<<tid<<" finished"<<std::endl;
            }

                    std::cout<<"parallel part passed - iter "<<iter<<std::endl;

            double end_solve_par = OpenMPUtils::GetCurrentTime();
            std::cout<<"-- par solve: "<<end_solve_par-start_solve_par<<std::endl;


            double start_qr_orth = OpenMPUtils::GetCurrentTime();

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

            double end_qr_orth = OpenMPUtils::GetCurrentTime();
            std::cout<<"-- qr orthog.: "<<end_qr_orth-start_qr_orth<<std::endl;

            iter++;
        } // - end of loop

        // project the input vector onto the reduced space
        r_b_reduced = prod( trans(r_Vr_sparse), r_b_tmp);

        double end_time = OpenMPUtils::GetCurrentTime();
        double duration = end_time - start_time;

        KRATOS_INFO("IRKA") << "Completed in " << duration << " seconds" << " after " << iter << " iterations" << std::endl;
  
        KRATOS_INFO_IF("IRKA", iter >= max_iter) << "Maximum number of iterations reached!"  << std::endl;
  
        //KRATOS_WATCH(r_M_reduced)

        std::cout << "MOR offline solve finished" << std::endl;




        // write the matrices into files in order to do post-processing
        // TODO: this is a quick dirty output; build another function/file for the post-processing step

        // write reduced matrices
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


        // write full system size matrices (double)
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


        // write the transformation matrix
        std::stringstream matrix_market_vr;
        matrix_market_vr << "Vr" << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_vr.str()).c_str(), r_Vr_sparse, false);

   

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
    typename TLinearSolver::Pointer mpLinearEigenSolver;
    ComplexLinearSolverType::Pointer mpComplexLinearSolver;

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

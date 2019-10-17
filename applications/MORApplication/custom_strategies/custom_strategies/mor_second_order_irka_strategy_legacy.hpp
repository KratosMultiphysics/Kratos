// // // // //    |  /           |
// // // // //    ' /   __| _` | __|  _ \   __|
// // // // //    . \  |   (   | |   (   |\__ \.
// // // // //   _|\_\_|  \__,_|\__|\___/ ____/
// // // // //                   Multi-Physics
// // // // //
// // // // //  License:		 BSD License
// // // // //					 Kratos default license: kratos/license.txt
// // // // //
// // // // //  Main authors:    Riccardo Rossi
// // // // //

// // // // #if !defined(MOR_SECOND_ORDER_IRKA_STRATEGY)
// // // // #define MOR_SECOND_ORDER_IRKA_STRATEGY

// // // // // System includesEeigenvalue problem

// // // // // External includes

// // // // // Project includes
// // // // #include "includes/define.h"
// // // // // #include "solving_strategies/strategies/solving_strategy.h"
// // // // // #include "solving_strategies/convergencecriterias/convergence_criteria.h"
// // // // #include "utilities/builtin_timer.h"
// // // // #include "utilities/qr_utility.h"
// // // // #include "custom_strategies/custom_strategies/mor_offline_second_order_strategy.hpp"

// // // // //default builder and solver
// // // // #include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"

// // // // // TODO: make this include work
// // // // //#include "EigenSolversApplication/costum_solvers/eigensystem_solver.h"
// // // // //#include "../ExternalSolversApplication/external_includes/feast_solver.h"

// // // // //test include eigen matrix
// // // // #include <Eigen/Dense>

// // // // namespace Kratos
// // // // {

// // // // ///@name Kratos Globals
// // // // ///@{

// // // // ///@}
// // // // ///@name Type Definitions
// // // // ///@{

// // // // ///@}

// // // // ///@name  Enum's
// // // // ///@{

// // // // ///@}
// // // // ///@name  Functions
// // // // ///@{

// // // // ///@}
// // // // ///@name Kratos Classes
// // // // ///@{

// // // // /**
// // // //  * @class MorSecondOrderIRKAStrategy
// // // //  * @ingroup KratosCore
// // // //  * @brief This is the Iterative Rational Krylov Algorithm
// // // //  * @details This strategy builds the K and M matrices and outputs them
// // // //  * @author Matthias Ebert, based on code of Aditya Ghantasala and Quirin Aumann
// // // //  */
// // // // template <class TSparseSpace,
// // // //           class TDenseSpace,  // = DenseSpace<double>,
// // // //           class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
// // // //           //class TLinearSolver  // feast? TODO: check if this is necessary
// // // //           //class TTrilLinearSolver
// // // //           >
// // // // class MorSecondOrderIRKAStrategy
// // // //     // : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
// // // //     : public MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver >
// // // // {
// // // //   public:
// // // //     ///@name Type Definitions
// // // //     ///@{
// // // //     // Counted pointer of ClassName
// // // //     KRATOS_CLASS_POINTER_DEFINITION(MorSecondOrderIRKAStrategy);

// // // //     // typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
// // // //     typedef MorOfflineSecondOrderStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

// // // //     typedef SystemMatrixBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;

// // // //     typedef typename BaseType::TDataType TDataType;

// // // //     typedef TSparseSpace SparseSpaceType;

// // // //     typedef TDenseSpace DenseSpaceType;

// // // //     typedef typename TDenseSpace::VectorType TDenseVectorType;

// // // //     typedef typename TDenseSpace::MatrixType TDenseMatrixType;

// // // //     typedef typename TDenseSpace::MatrixPointerType TDenseMatrixPointerType;

// // // //     typedef typename BaseType::TSchemeType TSchemeType;

// // // //     //typedef typename BaseType::DofSetType DofSetType;

// // // //     typedef typename BaseType::DofsArrayType DofsArrayType;

// // // //     typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

// // // //     typedef typename BaseType::TSystemVectorType TSystemVectorType;

// // // //     typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

// // // //     typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

// // // //     typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

// // // //     typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

// // // //     ///@}
// // // //     ///@name Life Cycle

// // // //     ///@{

// // // //     /**
// // // //      * Default constructor
// // // //      * @param rModelPart The model part of the problem
// // // //      * @param pScheme The integration schemed
// // // //      * @param MoveMeshFlag The flag that allows to move the mesh
// // // //      */
// // // //     MorSecondOrderIRKAStrategy(
// // // //         ModelPart& rModelPart,
// // // //         typename TSchemeType::Pointer pScheme,
// // // //         typename TLinearSolver::Pointer pNewLinearSolver,
// // // //         typename TLinearSolver::Pointer pNewLinearEigSolver,
// // // //         vector< double > samplingPoints,
// // // //         bool MoveMeshFlag = false)
// // // //         : BaseType(rModelPart, pScheme, pNewLinearSolver, MoveMeshFlag)  //hier evtl. direkt LU solver //TODO:  
// // // //     {
// // // //         KRATOS_TRY;

// // // //         // Saving the scheme
// // // //         this->SetScheme(pScheme);


// // // //         // Setting up the default builder and solver
// // // //         this->SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer(
// // // //             new TBuilderAndSolverType(pNewLinearSolver)));

// // // //         // Saving the linear solver        
// // // //         this->SetLinearSolver(pNewLinearSolver);

// // // //         // Set flags to correctly start the calculations
// // // //         mSolutionStepIsInitialized = false;
// // // //         mInitializeWasPerformed = false;

// // // //         // Tells to the builder and solver if the reactions have to be calculated or not
// // // //         this->GetBuilderAndSolver()->SetCalculateReactionsFlag(false);

// // // //         // Tells to the Builder And Solver if the system matrix and vectors need to
// // // //         // be reshaped at each step or not
// // // //         this->GetBuilderAndSolver()->SetReshapeMatrixFlag(false);

// // // //         // Set EchoLevel to the default value (only time is displayed)
// // // //         this->SetEchoLevel(1);

// // // //         // By default the matrices are rebuilt at each iteration
// // // //         this->SetRebuildLevel(0);

// // // //         mSamplingPoints = samplingPoints;

// // // //         mpNewLinearEigSolver = pNewLinearEigSolver;

// // // //         KRATOS_CATCH("");
// // // //     }

// // // //     /**
// // // //      * @brief Destructor.
// // // //      * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
// // // //      */
// // // //     ~MorSecondOrderIRKAStrategy() override
// // // //     {
// // // //         this->Clear();
// // // //     }

    

// // // //     //*********************************************************************************
// // // //     /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

// // // //     /**
// // // //      * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
// // // //      */
// // // //     bool SolveSolutionStep() override
// // // //     {
// // // //         KRATOS_TRY;
// // // //         std::cout << "hello! this is where the second order IRKA MOR magic happens" << std::endl;
// // // //         typename TSchemeType::Pointer p_scheme = this->GetScheme();
// // // //         typename TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver(); // LU
// // // //         //TODO: write functions to clear
// // // //         typename TBuilderAndSolverType::Pointer p_builder_and_solver_feast = typename TBuilderAndSolverType::Pointer(new TBuilderAndSolverType(mpNewLinearEigSolver)); // FEAST
        
// // // //         TSystemMatrixType& r_K_size_n = this->GetSystemMatrix();  // n x n
// // // //         TSystemMatrixType& r_M_size_n = this->GetMassMatrix();  // n x n
// // // //         TSystemMatrixType& r_D_size_n = this->GetDampingMatrix();  // n x n
// // // //         TSystemVectorType& r_b_size_n = this->GetSystemVector();  // n x 1

// // // //         TSystemVectorType tmp(r_K_size_n.size1(), 0.0);

// // // //         p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), r_b_size_n);

// // // //         p_builder_and_solver->BuildStiffnessMatrix(p_scheme, BaseType::GetModelPart(), r_K_size_n, tmp);
// // // //         p_builder_and_solver->ApplyDirichletConditions(p_scheme, BaseType::GetModelPart(), r_K_size_n, tmp, r_b_size_n);  

// // // //         p_builder_and_solver->BuildMassMatrix(p_scheme, BaseType::GetModelPart(), r_M_size_n, tmp);
// // // //         p_builder_and_solver->ApplyDirichletConditionsForMassMatrix(p_scheme, BaseType::GetModelPart(), r_M_size_n);

// // // //         p_builder_and_solver->BuildDampingMatrix(p_scheme, BaseType::GetModelPart(), r_D_size_n, tmp);
// // // //         p_builder_and_solver->ApplyDirichletConditionsForDampingMatrix(p_scheme, BaseType::GetModelPart(), r_D_size_n);

// // // //         // EchoInfo(0);
// // // //         const unsigned int system_size = p_builder_and_solver->GetEquationSystemSize();
// // // //         //sampling points
// // // //         KRATOS_WATCH(mSamplingPoints)
// // // //         const std::size_t n_sampling_points = mSamplingPoints.size();
// // // //         const std::size_t reduced_system_size = n_sampling_points; //number of sampling points

// // // //         KRATOS_WATCH(system_size) // n
// // // //         KRATOS_WATCH(reduced_system_size) // r

// // // //         // print all nxn matrices for DEBUG
// // // //         //TODO: print directly, too big for output
        
// // // //         //KRATOS_WATCH(r_M_size_n)
// // // //         //KRATOS_WATCH(r_K_size_n)
// // // //         //KRATOS_WATCH(r_D_size_n) // G
// // // //         //KRATOS_WATCH(r_b_size_n)
        

// // // //         // // // std::cout<<"r_D_size_n:\n";
// // // //         // // // for (unsigned int i=0; i<system_size; i++){
// // // //         // // //     for(unsigned int j=0; j<system_size; j++){
// // // //         // // //         if(r_D_size_n(i,j)!=0){
// // // //         // // //         std::cout<< r_D_size_n(i,j) << " (" << i <<","<< j << ")" << std::endl;
// // // //         // // //         }
// // // //         // // //     }
// // // //         // // // }


// // // //         // initialize V (=W, due to symmetry in FEM application)
// // // //         auto Vrp_col = DenseSpaceType::CreateEmptyVectorPointer();
// // // //         auto& r_Vr_col = *Vrp_col;
// // // //         DenseSpaceType::Resize(r_Vr_col, system_size);  // n x 1
// // // //         DenseSpaceType::Set(r_Vr_col,0.0);

// // // //         auto Vrp = DenseSpaceType::CreateEmptyMatrixPointer();
// // // //         auto& r_Vr = *Vrp;
// // // //         DenseSpaceType::Resize(r_Vr,system_size, reduced_system_size);  // n x r
// // // //         //SparseSpaceType::Set(Vr,0.0);  // only works with vectors


// // // //         // temporal helper
// // // //         auto V_h = SparseSpaceType::CreateEmptyMatrixPointer();
// // // //         auto& r_V_h = *V_h;
// // // //         SparseSpaceType::Resize(r_V_h, system_size, system_size); // n x n


// // // //         // step 2 as described in Wyatt 2012, Alg. 5.3.2     
// // // //         // TODO: put this in a function, since this also present in the while loop   
// // // //         for( size_t i = 0; i < n_sampling_points; ++i )
// // // //         {
// // // //             KRATOS_WATCH( mSamplingPoints(i) )
// // // //             r_V_h = std::pow( mSamplingPoints(i), 2.0 ) * r_M_size_n + mSamplingPoints(i) * r_D_size_n + r_K_size_n;
            
// // // //             p_builder_and_solver->GetLinearSystemSolver()->Solve( r_V_h, r_Vr_col, r_b_size_n );
// // // //             // multiply here with tangential direction b1, ..., br? Or postpone to QR later?
// // // //             column(r_Vr, i) = r_Vr_col;

// // // //         }

// // // //         mQR_decomposition.compute( system_size, reduced_system_size, &(r_Vr)(0,0) );
// // // //         mQR_decomposition.compute_q();

        
// // // //         for( size_t i = 0; i < system_size; ++i )
// // // //         {
// // // //             for( size_t j = 0; j < reduced_system_size; ++j )
// // // //             {
// // // //                 r_Vr(i,j) = mQR_decomposition.Q(i,j);
                
// // // //             }
// // // //         }
// // // //         // step 2 finished



// // // //         // step 4 as described in Wyatt 2012, Alg. 5.3.2  
// // // //         // iterative projection onto the Krylov subspace
// // // //         auto& r_b_reduced = this->GetRHSr(); // r x 1
// // // //         auto& r_K_reduced = this->GetKr(); // r x r
// // // //         auto& r_M_reduced = this->GetMr(); // r x r
// // // //         auto& r_D_reduced = this->GetDr(); // r x r
// // // //         TSystemMatrixType T; // temp
// // // //         vector<double> mSamplingPoints_old; // for convergence check; TODO: make complex (also in offline_strategy)


// // // //         int iter = 4;
// // // //         int err  = 1;
// // // //         // TODO: adapt max iter and tol; additional parameters, settings?
// // // //         while(iter<5 && err > 1e-4){

// // // //             mSamplingPoints_old = mSamplingPoints;

// // // //             // step 4 a)
// // // //             // B_r, reduced right hand side
// // // //             r_b_reduced = prod( trans(r_Vr), r_b_size_n );

// // // //             // K_r, reduced stiffness matrix
// // // //             T = prod( trans( r_Vr ), r_K_size_n );
// // // //             r_K_reduced = prod( T, r_Vr );
            
// // // //             // M_r, reduced mass matrix
// // // //             T = prod( trans( r_Vr ), r_M_size_n );
// // // //             r_M_reduced = prod( T, r_Vr );

// // // //             // D_r, reduced damping matrix
// // // //             T = prod( trans( r_Vr ), r_D_size_n );
// // // //             r_D_reduced = prod( T, r_Vr );
// // // //             // step 4 a) finished


// // // //             // step 4 c) polyeig from Matlab
// // // //             // generalized or FEAST solver? no special quadratic solver in KRATOS as it seems

// // // //             //TODO:  building the feast solve step here;  after everything works, move the definitions etc.
// // // //             TDenseVectorType Eigenvalues;
// // // //             TDenseMatrixType Eigenvectors;

// // // //             // DenseSpaceType::Resize(Eigenvalues, 2*reduced_system_size);
// // // //             // DenseSpaceType::Resize(Eigenvectors,2*reduced_system_size, 2*reduced_system_size);

// // // //             // linearize the EV problem, since FEAST can only deal with the generalized problem
// // // //             // auto L = SparseSpaceType::CreateEmptyMatrixPointer();
// // // //             // auto& r_L = *L;
// // // //             // SparseSpaceType::Resize(r_L, 2*reduced_system_size, 2*reduced_system_size); // 2r x 2r

// // // //             auto L1 = SparseSpaceType::CreateEmptyMatrixPointer();
// // // //             auto& r_L1 = *L1;
// // // //             SparseSpaceType::Resize(r_L1, 2*reduced_system_size, 2*reduced_system_size); // 2r x 2r

// // // //             auto L2 = SparseSpaceType::CreateEmptyMatrixPointer();
// // // //             auto& r_L2 = *L2;
// // // //             SparseSpaceType::Resize(r_L2, 2*reduced_system_size, 2*reduced_system_size); // 2r x 2r

// // // //             // "no type named size_type in class shared_ptr"
// // // //             //subrange(L1, (size_t) 0, reduced_system_size-1, (size_t) 0, reduced_system_size-1) = r_M_reduced;
// // // //             subrange(r_L1, 0, reduced_system_size, 0, reduced_system_size) = r_M_reduced;

// // // //             // "bad alloc" or memory corruption when running the python file...
// // // //             // for(size_t i = reduced_system_size; i < 2*reduced_system_size; ++i){
// // // //             //     r_L1(i,i) = 1.0;
// // // //             // }

// // // //             identity_matrix<double> id_m (reduced_system_size);

// // // //             subrange(r_L1, reduced_system_size, 2*reduced_system_size, reduced_system_size, 2*reduced_system_size) = id_m;


// // // //             subrange(r_L2, 0, reduced_system_size, 0, reduced_system_size) = r_D_reduced;
// // // //             subrange(r_L2, 0, reduced_system_size, reduced_system_size, 2*reduced_system_size) = r_K_reduced;            

// // // //             // for(size_t i = 0; i < reduced_system_size; ++i){
// // // //             //     r_L2(reduced_system_size+i,i) = -1;
// // // //             // }    

// // // //             subrange(r_L2, reduced_system_size, 2*reduced_system_size, 0, reduced_system_size) = -1.0*id_m;

// // // //             //r_L1(reduced_system_size,reduced_system_size) = 1;

// // // //             //r_L = 

// // // //             //KRATOS_WATCH(subrange(r_L1, reduced_system_size-2, reduced_system_size+2, reduced_system_size-2, reduced_system_size+2));
// // // //             //KRATOS_WATCH(subrange(r_L2, reduced_system_size-2, reduced_system_size+2, reduced_system_size-2, reduced_system_size+2));


// // // //             // // // // // Eigen::Matrix3d t_eig_m3;
// // // //             // // // // // t_eig_m3(0,0) = 2;
// // // //             // // // // // t_eig_m3(0,1) = 1;
// // // //             // // // // // t_eig_m3(0,2) = 1;

// // // //             // // // // // t_eig_m3(1,0) = 1;
// // // //             // // // // // t_eig_m3(1,1) = 2;
// // // //             // // // // // t_eig_m3(1,2) = 1;

// // // //             // // // // // t_eig_m3(2,0) = 1;
// // // //             // // // // // t_eig_m3(2,1) = 1;
// // // //             // // // // // t_eig_m3(2,2) = 2;

// // // //             // // // // // Eigen::Matrix3d t_eig_id3; 
// // // //             // // // // // t_eig_id3 = Eigen::Matrix3d::Identity();

// // // //             // // // // // KRATOS_WATCH(t_eig_m3);
// // // //             // // // // // KRATOS_WATCH(t_eig_id3);

// // // //             // // // // // Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eig_test;
// // // //             // // // // // eig_test.compute(t_eig_m3,t_eig_id3);
// // // //             // // // // // KRATOS_WATCH(eig_test.eigenvalues());
// // // //             // // // // // KRATOS_WATCH(eig_test.eigenvectors());

// // // //             // -> same results (also correct according to WolframAlpha)
// // // //             // -> ??? happens in eigensystem_solver.h for #eigenvalues>1 ...
// // // //             // -> need to code another wrapper? 
// // // //             // -> maybe also turn to feast then anyway (complex case included)
// // // //             // possible bug (?): 
// // // //             //               replace "int nc = std::min(2 * nroot, nroot + 8);"   //why??
// // // //             //               by/with "int nc = nroot;"
// // // //             //  then it also works via eigensystem_solver.h with #eigenvalues>1










// // // //             // DenseSpaceType::Resize(Eigenvalues,  3);
// // // //             // DenseSpaceType::Resize(Eigenvectors, 3, 3);
            
// // // //             // // // // // auto L1_test = SparseSpaceType::CreateEmptyMatrixPointer();
// // // //             // // // // // auto& r_L1_test = *L1_test;
// // // //             // // // // // SparseSpaceType::Resize(r_L1_test, 3, 3);
            
// // // //             // // // // // auto L2_test = SparseSpaceType::CreateEmptyMatrixPointer();
// // // //             // // // // // auto& r_L2_test = *L2_test;
// // // //             // // // // // SparseSpaceType::Resize(r_L2_test, 3, 3);

// // // //             // // // // // for(size_t i = 0; i < 3; ++i){
// // // //             // // // // //     r_L2_test(i,i) = 1;
// // // //             // // // // // }

// // // //             // // // // // r_L1_test(0,0) = 2;
// // // //             // // // // // r_L1_test(0,1) = 1;
// // // //             // // // // // r_L1_test(0,2) = 1;

// // // //             // // // // // r_L1_test(1,0) = 1;
// // // //             // // // // // r_L1_test(1,1) = 2;
// // // //             // // // // // r_L1_test(1,2) = 1;

// // // //             // // // // // r_L1_test(2,0) = 1;
// // // //             // // // // // r_L1_test(2,1) = 1;
// // // //             // // // // // r_L1_test(2,2) = 2;

// // // //             // // // // // KRATOS_WATCH(r_L1_test);
// // // //             // // // // // KRATOS_WATCH(r_L2_test);

// // // //             //TODO: issue linearized version: not spd

// // // //             //TDenseVectorType Eigenvalues_im;
// // // // //TODO: hier #ev setzen GetlinearSystemsolver->setnumev
// // // //             p_builder_and_solver_feast->GetLinearSystemSolver()->Solve(
// // // //                 r_L2,
// // // //                 r_L1,
// // // //                 //r_L1_test,
// // // //                 //r_L2_test,
// // // //                 //r_K_reduced,
// // // //                 //r_M_reduced,
// // // //                 Eigenvalues,
// // // //                 //Eigenvalues_im,
// // // //                 Eigenvectors);

// // // //                 KRATOS_WATCH(Eigenvalues);
// // // //                 //KRATOS_WATCH(Eigenvalues_im);
// // // //                 KRATOS_WATCH(Eigenvectors);
                



// // // //             Eigen::Matrix<double,8,8> r_L1_mtype;
// // // //             /*r_L1_mtype << 1.02562,-0.479005,0.0595067,-0.111968,0,0,0,0,
// // // //                           -0.479005,0.485658,-0.21239,0.0710417,0,0,0,0,
// // // //                           0.0595067,-0.21239,0.251714,-0.106319,0,0,0,0,
// // // //                           -0.111968,0.0710417,-0.106319,0.141759,0,0,0,0,
// // // //                           0,0,0,0,1,0,0,0,
// // // //                           0,0,0,0,0,1,0,0,
// // // //                           0,0,0,0,0,0,1,0,
// // // //                           0,0,0,0,0,0,0,1;
// // // //                           */

// // // //             // M in L1
// // // //             r_L1_mtype(0,0) = 1.02562; r_L1_mtype(0,1) = -0.479005; r_L1_mtype(0,2) = 0.0595067; r_L1_mtype(0,3) = -0.111968;
// // // //             r_L1_mtype(1,0) = -0.479005; r_L1_mtype(1,1) = 0.485658; r_L1_mtype(1,2) = -0.21239; r_L1_mtype(1,3) = 0.0710417;
// // // //             r_L1_mtype(2,0) = 0.0595067; r_L1_mtype(2,1) = -0.21239; r_L1_mtype(2,2) = 0.251714; r_L1_mtype(2,3) = -0.106319;
// // // //             r_L1_mtype(3,0) = -0.111968; r_L1_mtype(3,1) = 0.0710417; r_L1_mtype(3,2) = -0.106319; r_L1_mtype(3,3) = 0.141759;

// // // //             // I in L1
// // // //             r_L1_mtype(4,4) = 1.0; r_L1_mtype(4,5) = 0.0;  r_L1_mtype(4,6) = 0.0;  r_L1_mtype(4,7) = 0.0;
// // // //             r_L1_mtype(5,4) = 0.0; r_L1_mtype(5,5) = 1.0;  r_L1_mtype(5,6) = 0.0;  r_L1_mtype(5,7) = 0.0;
// // // //             r_L1_mtype(6,4) = 0.0; r_L1_mtype(6,5) = 0.0;  r_L1_mtype(6,6) = 1.0;  r_L1_mtype(6,7) = 0.0;
// // // //             r_L1_mtype(7,4) = 0.0; r_L1_mtype(7,5) = 0.0;  r_L1_mtype(7,6) = 0.0;  r_L1_mtype(7,7) = 1.0;


// // // //             // Rest 0
// // // //             for(int i=4; i<8; i++){
// // // //                 for(int j=0; j<4; j++){
// // // //                     r_L1_mtype(i,j) = 0.0;
// // // //                 }
// // // //             }

// // // //             for(int i=0; i<4; i++){
// // // //                 for(int j=4; j<8; j++){
// // // //                     r_L1_mtype(i,j) = 0.0;
// // // //                 }
// // // //             }



// // // //             //KRATOS_WATCH(r_L1_mtype);



// // // //             Eigen::Matrix<double,8,8> r_L2_mtype;
// // // //             /*r_L2_mtype << 0,0,0,0,-1,0,0,0,
// // // //                           0,0,0,0,0,-1,0,0,
// // // //                           0,0,0,0,0,0,-1,0,
// // // //                           0,0,0,0,0,0,0,-1,
// // // //                           -1,0,0,0,0,0,0,0,
// // // //                           0,-1,0,0,0,0,0,0,
// // // //                           0,0,-1,0,0,0,0,0,
// // // //                           0,0,0,-1,0,0,0,0;*/

// // // //             // K in L2
// // // //             r_L2_mtype(0,4) = 7980.85; r_L2_mtype(0,5) = 5016.95; r_L2_mtype(0,6) = 1188.86; r_L2_mtype(0,7) = 226.714;
// // // //             r_L2_mtype(1,4) = 5016.95; r_L2_mtype(1,5) = 105123; r_L2_mtype(1,6) = 60796.1; r_L2_mtype(1,7) = 16482.2;
// // // //             r_L2_mtype(2,4) = 1188.86; r_L2_mtype(2,5) = 60796.1; r_L2_mtype(2,6) = 456395; r_L2_mtype(2,7) = 251600;
// // // //             r_L2_mtype(3,4) = 226.714; r_L2_mtype(3,5) = 16482.2; r_L2_mtype(3,6) = 251600; r_L2_mtype(3,7) = 1.07785e+06;

// // // //             // -I in L2
// // // //             r_L2_mtype(4,0) = -1.0; r_L2_mtype(4,1) = 0.0; r_L2_mtype(4,2) = 0.0; r_L2_mtype(4,3) = 0.0;
// // // //             r_L2_mtype(5,0) = 0.0; r_L2_mtype(5,1) = -1.0; r_L2_mtype(5,2) = 0.0; r_L2_mtype(5,3) = 0.0;
// // // //             r_L2_mtype(6,0) = 0.0; r_L2_mtype(6,1) = 0.0; r_L2_mtype(6,2) = -1.0; r_L2_mtype(6,3) = 0.0;
// // // //             r_L2_mtype(7,0) = 0.0; r_L2_mtype(7,1) = 0.0; r_L2_mtype(7,2) = 0.0; r_L2_mtype(7,3) = -1.0;



// // // //             // Rest 0  (C)
// // // //             for(int i=0; i<4; i++){
// // // //                 for(int j=0; j<4; j++){
// // // //                     r_L2_mtype(i,j) = 0.0;
// // // //                 }
// // // //             }

// // // //             for(int i=4; i<8; i++){
// // // //                 for(int j=4; j<8; j++){
// // // //                     r_L2_mtype(i,j) = 0.0;
// // // //                 }
// // // //             }
            
// // // //             //(r_L2_mtype);





// // // //                 //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eig_test2;
// // // //                 Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> eig_test2;
// // // //             eig_test2.compute(r_L2_mtype,r_L1_mtype);
// // // //             //KRATOS_WATCH(eig_test2.eigenvalues());
// // // //             //KRATOS_WATCH(eig_test2.alphas());
// // // //             //KRATOS_WATCH(eig_test2.betas());
// // // //             //KRATOS_WATCH(eig_test2.eigenvectors());

// // // //             //this->AssignVariables(Eigenvalues,Eigenvectors);
// // // //             // // KRATOS_WATCH(Eigenvalues);
// // // //             // // KRATOS_WATCH(Eigenvectors);


// // // //             // auto& r_pa_test_tmp = eig_test2.eigenvectors();
// // // //             // auto& r_pa_test = prod(r_L2_mtype,r_pa_test_tmp);
// // // //             // KRATOS_WATCH(r_pa_test);
// // // //             // auto& r_pb_test_tmp = eig_test2.eigenvectors();
// // // //             // r_pb_test_tmp = prod(r_L1_mtype,r_pb_test_tmp);
// // // //             // auto& r_pb_test_tmp_val = eig_test2.eigenvalues();
// // // //             // auto& r_pb_test = prod(r_pb_test_tmp_val, r_pb_test_tmp);
// // // //             // KRATOS_WATCH(r_pb_test);


// // // //         // print all rxr matrices for DEBUG
// // // //         /*
// // // //         KRATOS_WATCH(r_L1);
// // // //         KRATOS_WATCH(r_L2);
// // // //         KRATOS_WATCH(r_M_reduced)
// // // //         KRATOS_WATCH(r_K_reduced)
// // // //         KRATOS_WATCH(r_D_reduced)
// // // //         KRATOS_WATCH(r_b_reduced)
// // // // */
    


// // // //             // step 4 d) shift selection step
// // // //             // first r eigenvalues
// // // //             mSamplingPoints = subrange(Eigenvalues, 0, reduced_system_size);
// // // //             KRATOS_WATCH(mSamplingPoints)

// // // //             // step 4 e) as described in Wyatt 2012, Alg. 5.3.2     
// // // //             // TODO: put this in a function, since this the same as step 2
// // // //             for( size_t i = 0; i < n_sampling_points; ++i )
// // // //             {
// // // //                 //KRATOS_WATCH( mSamplingPoints(i) )
// // // //                 r_V_h = std::pow( mSamplingPoints(i), 2.0 ) * r_M_size_n + mSamplingPoints(i) * r_D_size_n + r_K_size_n;
                
// // // //                 p_builder_and_solver->GetLinearSystemSolver()->Solve( r_V_h, r_Vr_col, r_b_size_n );
// // // //                 // multiply here with tangential direction b1, ..., br? Or postpone to QR later?
// // // //                 column(r_Vr, i) = r_Vr_col;

// // // //             }

// // // //             mQR_decomposition.compute( system_size, reduced_system_size, &(r_Vr)(0,0) );
// // // //             mQR_decomposition.compute_q();

            
// // // //             for( size_t i = 0; i < system_size; ++i )
// // // //             {
// // // //                 for( size_t j = 0; j < reduced_system_size; ++j )
// // // //                 {
// // // //                     r_Vr(i,j) = mQR_decomposition.Q(i,j);
                    
// // // //                 }
// // // //             }
// // // //             // step 4 e) finished


// // // //             // check convergence
// // // //             // TODO: uncomment; leave it out for the moment, until eigensolver works (i.e. new sampling points are calculated)
// // // //             // err = norm_2(mSamplingPoints - mSamplingPoints_old)/norm_2(mSamplingPoints_old);
// // // //             // KRATOS_WATCH(iter)
// // // //             // KRATOS_WATCH(err)

// // // //             iter++;
// // // //         }


        
        
// // // //         // print all rxr matrices for DEBUG
// // // //                 // KRATOS_WATCH(r_M_reduced)
// // // //                 // KRATOS_WATCH(r_K_reduced)
// // // //                 // KRATOS_WATCH(r_D_reduced)
// // // //                 // KRATOS_WATCH(r_b_reduced)
        
        




// // // //         std::cout << "MOR offline solve finished" << std::endl;
        
// // // // 		return true;

// // // //         KRATOS_CATCH("");
// // // //     }

// // // //     ///@}
// // // //     ///@name Operators
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Operations
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Access
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Input and output
// // // //     ///@{

// // // //     /// Turn back information as a string.
// // // //     virtual std::string Info() const override
// // // //     {
// // // //         return "MorSecondOrderIRKAStrategy";
// // // //     }

// // // //     ///@}
// // // //     ///@name Inquiry
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Friends
// // // //     ///@{

// // // //     ///@}

// // // //   private:
// // // //     ///@name Protected static Member Variables
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Protected member Variables
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Protected Operators
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Protected Operations
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Protected  Access
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Protected Inquiry
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Protected LifeCycle
// // // //     ///@{

// // // //     ///@}

// // // //   protected:
// // // //     ///@name Static Member Variables
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Member Variables
// // // //     ///@{

// // // //     vector< double > mSamplingPoints;
// // // //     QR<double, row_major> mQR_decomposition;
// // // //     typename TLinearSolver::Pointer mpNewLinearEigSolver;

// // // //     /**
// // // //      * @brief Flag telling if it is needed to reform the DofSet at each
// // // //     solution step or if it is possible to form it just once
// // // //     * @details Default = false
// // // //         - true  : Reforme at each time step
// // // //         - false : Form just one (more efficient)
// // // //      */
// // // //     bool mReformDofSetAtEachStep;

// // // //     bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

// // // //     bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

// // // //     ///@}
// // // //     ///@name Private Operators
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Private Operations
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Private  Access
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Private Inquiry
// // // //     ///@{

// // // //     ///@}
// // // //     ///@name Un accessible methods
// // // //     ///@{

// // // //     /**
// // // //      * Copy constructor.
// // // //      */

// // // //     MorSecondOrderIRKAStrategy(const MorSecondOrderIRKAStrategy &Other){};

// // // //     ///@}

// // // // }; /* Class MorSecondOrderIRKAStrategy */

// // // // ///@}

// // // // ///@name Type Definitions
// // // // ///@{

// // // // ///@}

// // // // } /* namespace Kratos. */

// // // // #endif /* MOR_SECOND_ORDER_IRKA_STRATEGY  defined */
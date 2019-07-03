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

// System includesEeigenvalue problem

// External includes

// Project includes
#include "includes/define.h"
// #include "solving_strategies/strategies/solving_strategy.h"
// #include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"
#include "utilities/qr_utility.h"
#include "custom_strategies/custom_strategies/mor_offline_second_order_strategy.hpp"

//default builder and solver
#include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"

//#include "EigenSolversApplication/costum_solvers/eigensystem_solver.h"

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
 * @ingroup KratosCore
 * @brief This is the Iterative Rational Krylov Algorithm
 * @details This strategy builds the K and M matrices and outputs them
 * @author Matthias Ebert, based on code of Aditya Ghantasala and Quirin Aumann
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class MorSecondOrderIRKAStrategy
    // : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    : public MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver >
{
  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorSecondOrderIRKAStrategy);

    // typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef MorOfflineSecondOrderStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef SystemMatrixBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef TDenseSpace DenseSpaceType;

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
        vector< double > samplingPoints,
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

        // Set flags to correctly start the calculations
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

        mSamplingPoints = samplingPoints;

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorSecondOrderIRKAStrategy() override
    {
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
        
        TSystemMatrixType& r_K_size_n = this->GetSystemMatrix();  // n x n
        TSystemMatrixType& r_M_size_n = this->GetMassMatrix();  // n x n
        TSystemMatrixType& r_D_size_n = this->GetDampingMatrix();  // n x n
        TSystemVectorType& r_b_size_n = this->GetSystemVector();  // n x 1

        TSystemVectorType tmp(r_K_size_n.size1(), 0.0);

        p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), r_b_size_n);

        p_builder_and_solver->BuildStiffnessMatrix(p_scheme, BaseType::GetModelPart(), r_K_size_n, tmp);
        p_builder_and_solver->ApplyDirichletConditions(p_scheme, BaseType::GetModelPart(), r_K_size_n, tmp, r_b_size_n);  

        p_builder_and_solver->BuildMassMatrix(p_scheme, BaseType::GetModelPart(), r_M_size_n, tmp);
        p_builder_and_solver->ApplyDirichletConditionsForMassMatrix(p_scheme, BaseType::GetModelPart(), r_M_size_n);

        p_builder_and_solver->BuildDampingMatrix(p_scheme, BaseType::GetModelPart(), r_D_size_n, tmp);
        p_builder_and_solver->ApplyDirichletConditionsForDampingMatrix(p_scheme, BaseType::GetModelPart(), r_D_size_n);

        // EchoInfo(0);
        const unsigned int system_size = p_builder_and_solver->GetEquationSystemSize();
        //sampling points
        KRATOS_WATCH(mSamplingPoints)
        const std::size_t n_sampling_points = mSamplingPoints.size();
        const std::size_t reduced_system_size = n_sampling_points; //number of sampling points

        KRATOS_WATCH(system_size) // n
        KRATOS_WATCH(reduced_system_size) // r


        // initialize V (=W, due to symmetry in FEM application)
        auto Vrp_col = DenseSpaceType::CreateEmptyVectorPointer();
        auto& r_Vr_col = *Vrp_col;
        DenseSpaceType::Resize(r_Vr_col, system_size);  // n x 1
        DenseSpaceType::Set(r_Vr_col,0.0);

        auto Vrp = DenseSpaceType::CreateEmptyMatrixPointer();
        auto& r_Vr = *Vrp;
        DenseSpaceType::Resize(r_Vr,system_size, reduced_system_size);  // n x r
        //SparseSpaceType::Set(Vr,0.0);  // only works with vectors


        //auto Vrp_orth = SparseSpaceType::CreateEmptyMatrixPointer();
        //auto& r_Vr_orth = *Vrp_orth;
        //SparseSpaceType::Resize(r_Vr_orth,system_size, reduced_system_size);  // n x r



        // temporal helper
        auto V_h = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_V_h = *V_h;
        SparseSpaceType::Resize(r_V_h, system_size, system_size); // n x n


        // step 2 as described in Wyatt 2012, Alg. 5.3.2     
        // TODO: put this in a function, since this also present in the while loop   
        for( size_t i = 0; i < n_sampling_points; ++i )
        {
            KRATOS_WATCH( mSamplingPoints(i) )
            r_V_h = std::pow( mSamplingPoints(i), 2.0 ) * r_M_size_n + mSamplingPoints(i) * r_D_size_n + r_K_size_n;
            
            p_builder_and_solver->GetLinearSystemSolver()->Solve( r_V_h, r_Vr_col, r_b_size_n );
            // multiply here with tangential direction b1, ..., br? Or postpone to QR later?
            column(r_Vr, i) = r_Vr_col;

        }

        mQR_decomposition.compute( system_size, reduced_system_size, &(r_Vr)(0,0) );
        mQR_decomposition.compute_q();

        
        for( size_t i = 0; i < system_size; ++i )
        {
            for( size_t j = 0; j < reduced_system_size; ++j )
            {
                r_Vr(i,j) = mQR_decomposition.Q(i,j);
                
            }
        }
        // step 2 finished



        // step 4 as described in Wyatt 2012, Alg. 5.3.2  
        // iterative projection onto the Krylov subspace
        auto& r_b_reduced = this->GetRHSr(); // r x 1
        auto& r_K_reduced = this->GetKr(); // r x r
        auto& r_M_reduced = this->GetMr(); // r x r
        auto& r_D_reduced = this->GetDr(); // r x r
        TSystemMatrixType T; // temp


        int iter = 1;
        int err  = 1;
        // TODO: adapt max iter and tol; additional parameters, settings?
        while(iter<10 && err > 1e-4){
            // step 4 a)
            // B_r, reduced right hand side
            r_b_reduced = prod( trans(r_Vr), r_b_size_n );

            // K_r, reduced stiffness matrix
            T = prod( trans( r_Vr ), r_K_size_n );
            r_K_reduced = prod( T, r_Vr );
            
            // M_r, reduced mass matrix
            T = prod( trans( r_Vr ), r_M_size_n );
            r_M_reduced = prod( T, r_Vr );

            // D_r, reduced damping matrix
            T = prod( trans( r_Vr ), r_D_size_n );
            r_D_reduced = prod( T, r_Vr );
            // step 4 a) finished


            // step 4 c) polyeig from Matlab
            // generalized or FEAST solver? no special quadratic solver in KRATOS as it seems

            // step 4 d) shift selection step
            // first r eigenvalues

            // step 4 e) as described in Wyatt 2012, Alg. 5.3.2     
            // TODO: put this in a function, since this the same as step 2
            for( size_t i = 0; i < n_sampling_points; ++i )
            {
                //KRATOS_WATCH( mSamplingPoints(i) )
                r_V_h = std::pow( mSamplingPoints(i), 2.0 ) * r_M_size_n + mSamplingPoints(i) * r_D_size_n + r_K_size_n;
                
                p_builder_and_solver->GetLinearSystemSolver()->Solve( r_V_h, r_Vr_col, r_b_size_n );
                // multiply here with tangential direction b1, ..., br? Or postpone to QR later?
                column(r_Vr, i) = r_Vr_col;

            }

            mQR_decomposition.compute( system_size, reduced_system_size, &(r_Vr)(0,0) );
            mQR_decomposition.compute_q();

            
            for( size_t i = 0; i < system_size; ++i )
            {
                for( size_t j = 0; j < reduced_system_size; ++j )
                {
                    r_Vr(i,j) = mQR_decomposition.Q(i,j);
                    
                }
            }
            // step 4 e) finished

            iter++;
        }







        
        KRATOS_WATCH(r_M_reduced)











        //p_builder_and_solver->GetLinearSystemSolver()->Solve( mSamplingPoints[0]*mSamplingPoints[0]*M_full + mSamplingPoints[0]*D_full + K_full, Vr_col, b_full );
        //KRATOS_WATCH(Vr_col)
/* 
        auto M_reduced_p = DenseSpaceType::CreateEmptyMatrixPointer();
        auto& M_reduced = *M_reduced_p;
        DenseSpaceType::Resize(M_reduced, reduced_system_size, reduced_system_size);

        auto D_reduced_p = DenseSpaceType::CreateEmptyMatrixPointer();
        auto& D_reduced = *D_reduced_p;
        DenseSpaceType::Resize(D_reduced, reduced_system_size, reduced_system_size);

        auto K_reduced_p = DenseSpaceType::CreateEmptyMatrixPointer();
        auto& K_reduced = *K_reduced_p;
        DenseSpaceType::Resize(M_reduced, reduced_system_size, reduced_system_size);

        // auto M_reduced_p = DenseSpaceType::CreateEmptyMatrixPointer();
        // auto& M_reduced = *M_reduced_p;
        // DenseSpaceType::Resize(M_reduced, reduced_system_size, reduced_system_size);


        auto tmp_transfer_p = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& tmp_transfer = *tmp_transfer_p;
        SparseSpaceType::Resize(tmp_transfer, system_size, system_size);

        for( size_t i = 0; i < n_sampling_points; ++i ){
            double current_s = mSamplingPoints[i];
            tmp_transfer = current_s*current_s*M_full + current_s*D_full + K_full;
            // auto aaa = p_builder_and_solver->GetLinearSystemSolver();
            p_builder_and_solver->GetLinearSystemSolver()->Solve(tmp_transfer, Vr_col, b_full );
            column( Vr, (i) ) = Vr_col;
        }


        unsigned int iter = 1;
        double err = 1.0;
        while((iter < 100)&&(err > 10e-6)){ // (err>tol){
            vector<double> mSamplingPoints_old = mSamplingPoints;

            M_reduced = prod(trans(Vr), prod(M_full,Vr));
            D_reduced = prod(trans(Vr), prod(D_full,Vr));
            K_reduced = prod(trans(Vr), prod(K_full,Vr));

            for( size_t i = 0; i < n_sampling_points; ++i ){
                double current_s = mSamplingPoints[i];
                tmp_transfer = current_s*current_s*M_full + current_s*D_full + K_full;
                // p_builder_and_solver->GetLinearSystemSolver()->Solve(tmp_transfer, Vr_col, b_full );
                column( Vr, (i) ) = Vr_col;
            }

            // check convergence
            err = norm_2(mSamplingPoints - mSamplingPoints_old)/norm_2(mSamplingPoints_old);
            //KRATOS_WATCH(err)
            iter++;

        }

        // M_reduced = prod(prod(trans(Vr), M_full), Vr);
        // D_reduced = prod(trans(Vr), prod(D_full,Vr));
        // K_reduced = prod(trans(Vr), prod(K_full,Vr));
        // C_reduced = C * Vr;


        KRATOS_WATCH(M_reduced)
        */






        //KRATOS_WATCH(Vr)
 /*      
        //initialize sb, As, AAs vectors
        auto s = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rs = *s;
        SparseSpaceType::Resize(rs,system_size);
        SparseSpaceType::Set(rs,0.0);
        auto As = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rAs = *As;
        SparseSpaceType::Resize(rAs,system_size);
        SparseSpaceType::Set(rAs,0.0);
        auto AAs = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rAAs = *AAs;
        SparseSpaceType::Resize(rAAs,system_size);
        SparseSpaceType::Set(rAAs,0.0);

        auto kdyn = SparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_kdyn = *kdyn;
        SparseSpaceType::Resize(r_kdyn, system_size, system_size);

        auto tmp_basis = DenseSpaceType::CreateEmptyMatrixPointer();
        auto& r_tmp_basis = *tmp_basis;
        DenseSpaceType::Resize(r_tmp_basis, system_size, reduced_system_size);

        auto& r_basis = this->GetBasis(); 
        SparseSpaceType::Resize(r_basis, system_size, reduced_system_size);

        TSystemVectorType aux;
        
        for( size_t i = 0; i < n_sampling_points; ++i )
        {
            KRATOS_WATCH( mSamplingPoints(i) )
            r_kdyn = K_full - ( std::pow( mSamplingPoints(i), 2.0 ) * M_full );    // Without Damping  
            //r_kdyn = K_full - ( std::pow( mSamplingPoints(i), 2.0 ) * M_full )-D_full;  With Damping
            
            p_builder_and_solver->GetLinearSystemSolver()->Solve( r_kdyn, rs, b_full );
            aux = prod( M_full, rs );
            p_builder_and_solver->GetLinearSystemSolver()->Solve( r_kdyn, rAs, aux );
            aux = prod( M_full, rAs );
            p_builder_and_solver->GetLinearSystemSolver()->Solve( r_kdyn, rAAs, aux );

            column( r_tmp_basis, (i*3) ) = rs;
            column( r_tmp_basis, (i*3)+1 ) = rAs;
            column( r_tmp_basis, (i*3)+2 ) = rAAs;
        }

        //orthogonalize the basis -> basis_r
        // mQR_decomposition.compute( system_size, 3*n_sampling_points, &(r_basis)(0,0) );
        mQR_decomposition.compute( system_size, 3*n_sampling_points, &(r_tmp_basis)(0,0) );
        mQR_decomposition.compute_q();
        
        for( size_t i = 0; i < system_size; ++i )
        {
            for( size_t j = 0; j < (3*n_sampling_points); ++j )
            {
                r_basis(i,j) = mQR_decomposition.Q(i,j);
                
            }
        }

        // project the system matrices onto the Krylov subspace
        auto& r_force_vector_reduced = this->GetRHSr();
        auto& r_stiffness_matrix_reduced = this->GetKr();
        auto& r_mass_matrix_reduced = this->GetMr();
        auto& r_damping_matrix_reduced = this->GetDr();

        r_force_vector_reduced = prod( b_full, r_basis );

        TSystemMatrixType T = prod( trans( r_basis ), K_full );
        r_stiffness_matrix_reduced = prod( T, r_basis );

        T = prod( trans( r_basis ), M_full );
        r_mass_matrix_reduced = prod( T, r_basis );

        T = prod( trans( r_basis ), D_full );
        r_damping_matrix_reduced = prod( T, r_basis );
        
        KRATOS_WATCH(r_mass_matrix_reduced)
*/
        std::cout << "MOR offline solve finished" << std::endl;
        
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

    vector< double > mSamplingPoints;
    QR<double, row_major> mQR_decomposition;

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
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

// TODO: make this include work
//#include "EigenSolversApplication/costum_solvers/eigensystem_solver.h"
//#include "../ExternalSolversApplication/external_includes/feast_solver.h"

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
          class TLinearSolver//, //= LinearSolver<TSparseSpace,TDenseSpace>
          //class TLinearSolver  // feast? TODO: check if this is necessary
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
        typename TLinearSolver::Pointer pNewLinearEigSolver,
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

        mpNewLinearEigSolver = pNewLinearEigSolver;

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
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver(); // LU
        //TODO: write functions to clear
        typename TBuilderAndSolverType::Pointer p_builder_and_solver_feast = typename TBuilderAndSolverType::Pointer(new TBuilderAndSolverType(mpNewLinearEigSolver)); // FEAST
        
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
        vector<double> mSamplingPoints_old; // for convergence check; TODO: make complex (also in offline_strategy)


        int iter = 1;
        int err  = 1;
        // TODO: adapt max iter and tol; additional parameters, settings?
        while(iter<10 && err > 1e-4){

            mSamplingPoints_old = mSamplingPoints;

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

            //TODO:  building the feast solve step here;  after everything works, move the definitions etc.
            TDenseVectorType Eigenvalues;
            TDenseMatrixType Eigenvectors;


            // linearize the EV problem, since FEAST can only deal with the generalized problem
            // auto L = SparseSpaceType::CreateEmptyMatrixPointer();
            // auto& r_L = *L;
            // SparseSpaceType::Resize(r_L, 2*reduced_system_size, 2*reduced_system_size); // 2r x 2r

            auto L1 = SparseSpaceType::CreateEmptyMatrixPointer();
            auto& r_L1 = *L1;
            SparseSpaceType::Resize(r_L1, 2*reduced_system_size, 2*reduced_system_size); // 2r x 2r

            auto L2 = SparseSpaceType::CreateEmptyMatrixPointer();
            auto& r_L2 = *L2;
            SparseSpaceType::Resize(r_L2, 2*reduced_system_size, 2*reduced_system_size); // 2r x 2r

            // "no type named size_type in class shared_ptr"
            //subrange(L1, (size_t) 0, reduced_system_size-1, (size_t) 0, reduced_system_size-1) = r_M_reduced;
            subrange(r_L1, 0, reduced_system_size-1, 0, reduced_system_size-1) = r_M_reduced;

            // "bad alloc" or memory corruption when running the python file...
            // for(size_t i = reduced_system_size; i < 2*reduced_system_size; ++i){
            //     r_L1(i,i) = 1.0;
            // }

            identity_matrix<double> id_m (reduced_system_size);

            subrange(r_L1, reduced_system_size, 2*reduced_system_size, reduced_system_size, 2*reduced_system_size) = id_m;


            subrange(r_L2, 0, reduced_system_size, 0, reduced_system_size) = r_D_reduced;
            subrange(r_L2, 0, reduced_system_size, reduced_system_size, 2*reduced_system_size) = r_K_reduced;            

            // for(size_t i = 0; i < reduced_system_size; ++i){
            //     r_L2(reduced_system_size+i,i) = -1;
            // }

            subrange(r_L2, reduced_system_size, 2*reduced_system_size, 0, reduced_system_size) = -1.0*id_m;

            //r_L1(reduced_system_size,reduced_system_size) = 1;

            //r_L = 

            //KRATOS_WATCH(subrange(r_L1, reduced_system_size-2, reduced_system_size+2, reduced_system_size-2, reduced_system_size+2));
            //KRATOS_WATCH(subrange(r_L2, reduced_system_size-2, reduced_system_size+2, reduced_system_size-2, reduced_system_size+2));

 
            p_builder_and_solver_feast->GetLinearSystemSolver()->Solve(
                r_L2,
                r_L1,
                //r_K_reduced,
                //r_M_reduced,
                Eigenvalues,
                Eigenvectors);

            //this->AssignVariables(Eigenvalues,Eigenvectors);
            KRATOS_WATCH(Eigenvalues);
            KRATOS_WATCH(Eigenvectors);


            // step 4 d) shift selection step
            // first r eigenvalues
            //mSamplingPoints = subrange(Eigenvalues, 0, reduced_system_size-1);

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


            // check convergence
            // TODO: uncomment; leave it out for the moment, until eigensolver works (i.e. new sampling points are calculated)
            //err = norm_2(mSamplingPoints - mSamplingPoints_old)/norm_2(mSamplingPoints_old);
            //KRATOS_WATCH(iter)
            //KRATOS_WATCH(err)

            iter++;
        }


        
        KRATOS_WATCH(r_M_reduced)





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
    typename TLinearSolver::Pointer mpNewLinearEigSolver;

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
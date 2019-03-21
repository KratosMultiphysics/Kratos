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

#if !defined(MOR_SECOND_ORDER_KRYLOV_STRATEGY)
#define MOR_SECOND_ORDER_KRYLOV_STRATEGY

// System includes

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
 * @class MorSecondOrderKrylovStrategy
 * @ingroup KratosCore
 * @brief This is the linear MOR matrix output strategy
 * @details This strategy builds the K and M matrices and outputs them
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class MorSecondOrderKrylovStrategy
    // : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    : public MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver >
{
  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorSecondOrderKrylovStrategy);

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
    MorSecondOrderKrylovStrategy(
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

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
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
    ~MorSecondOrderKrylovStrategy() override
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
        std::cout << "hello! this is where the second order krylov MOR magic happens" << std::endl;
        typename TSchemeType::Pointer p_scheme = this->GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver();
        
        TSystemMatrixType& rA = this->GetSystemMatrix();
        TSystemMatrixType& rM = this->GetMassMatrix();
        TSystemMatrixType& rS = this->GetDampingMatrix();
        TSystemVectorType& rRHS = this->GetSystemVector();
        SparseSpaceType::Set(rRHS,0.0); //why??

        // std::cout << "system matrices before initialization" << std::endl;
        // KRATOS_WATCH(rA)
        // KRATOS_WATCH(rM)
        // KRATOS_WATCH(rRHS)

        TSystemVectorType tmp(rA.size1(), 0.0);
        //std::cout<<"temp = "<<rA.size1()<<std::endl;

        //TSystemVectorType tmp1(rS.size1(), 0.0);
        //std::cout<<"temp1 = "<<rS.size1()<<std::endl;

        p_builder_and_solver->BuildStiffnessMatrix(p_scheme, BaseType::GetModelPart(), rA, tmp);

        // Applying the Dirichlet Boundary conditions
        p_builder_and_solver->ApplyDirichletConditions(p_scheme, BaseType::GetModelPart(), rA, tmp, rRHS);  

        p_builder_and_solver->BuildMassMatrix(p_scheme, BaseType::GetModelPart(), rM, tmp);

        p_builder_and_solver->ApplyDirichletConditionsForMassMatrix(p_scheme, BaseType::GetModelPart(), rM);

        p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), rRHS);

        p_builder_and_solver->BuildDampingMatrix(p_scheme, BaseType::GetModelPart(), rS, tmp);

        p_builder_and_solver->ApplyDirichletConditionsForDampingMatrix(p_scheme, BaseType::GetModelPart(), rS);
        
        //std::cout<<"Damping matrix is "<<rS<<std::endl;
        
        //std::cout<<"\nDamping matrix"<<std::endl<<rS<<std::endl;
        //std::cout<<"\nMass Matrix"<<std::endl<<rM<<std::endl;
        //std::cout<<"\nStiffness matrix"<<std::endl<<rA<<std::endl;

        // EchoInfo(0);
        const unsigned int system_size = p_builder_and_solver->GetEquationSystemSize();
        //sampling points
        KRATOS_WATCH(mSamplingPoints)
        const std::size_t n_sampling_points = mSamplingPoints.size();
        const std::size_t reduced_system_size = 3 * n_sampling_points;
        
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

        vector< double > aux;
        vector <double > Q_vector;
        
        for( size_t i = 0; i < n_sampling_points; ++i )
        {
            KRATOS_WATCH( mSamplingPoints(i) )
            r_kdyn = rA - ( std::pow( mSamplingPoints(i), 2.0 ) * rM );    // Without Damping  
            //r_kdyn = rA - ( std::pow( mSamplingPoints(i), 2.0 ) * rM )-rS;  With Damping
            
            // KRATOS_WATCH(r_kdyn)
            // KRATOS_WATCH(rRHS)
            // std::cout << "rs vorher" << std::endl;
            // KRATOS_WATCH(rs)
            // KRATOS_WATCH(r_force_vector)
            p_builder_and_solver->GetLinearSystemSolver()->Solve( r_kdyn, rs, rRHS );
            // KRATOS_WATCH(rs)
            aux = prod( rM, rs );
            //KRATOS_WATCH(aux)
            p_builder_and_solver->GetLinearSystemSolver()->Solve( r_kdyn, rAs, aux );
            aux = prod( rM, rAs );
            p_builder_and_solver->GetLinearSystemSolver()->Solve( r_kdyn, rAAs, aux );

            // KRATOS_WATCH(rs)
            // KRATOS_WATCH(rAs)
            // KRATOS_WATCH(rAAs)

            column( r_tmp_basis, (i*3) ) = rs;
            column( r_tmp_basis, (i*3)+1 ) = rAs;
            column( r_tmp_basis, (i*3)+2 ) = rAAs;


            // KRATOS_WATCH(r_tmp_basis)
        }

        //orthogonalize the basis -> basis_r
        // mQR_decomposition.compute( system_size, 3*n_sampling_points, &(r_basis)(0,0) );
        mQR_decomposition.compute( system_size, 3*n_sampling_points, &(r_tmp_basis)(0,0) );
        // std::cout << "yo2" << std::endl;
        // KRATOS_WATCH(r_basis)
        mQR_decomposition.compute_q();

        // auto basis_r = SparseSpaceType::CreateEmptyMatrixPointer();
        // auto& r_basis_r = *mpReducedBasis;
        // SparseSpaceType::Resize(r_basis_r, system_size, 3*n_sampling_points);
        auto& r_basis = this->GetBasis(); 
        // auto& r_basis = *mpBasis;
        SparseSpaceType::Resize(r_basis, system_size, reduced_system_size);
        
        for( size_t i = 0; i < system_size; ++i )
        {
            for( size_t j = 0; j < (3*n_sampling_points); ++j )
            {
                r_basis(i,j) = mQR_decomposition.Q(i,j);
                
            }
        }

        std::cout<<Q_vector.size()<<std::endl;
        // auto r_basis_r = r_basis;
        // auto& r_force_vector_reduced = *mpRHSr;
        // mpForceVectorReduced = ZeroMatrix( 3*n_sampling_points );
        // r_force_vector_reduced = (prod( rRHS, r_basis ));
        // auto& r_stiffness_matrix_reduced = *mpAr;
        // auto& r_mass_matrix_reduced = *mpMr;
        // auto& r_damping_reduced = *mpSr;   // Reduced damping matrix

        auto& r_force_vector_reduced = this->GetRHSr();
        auto& r_stiffness_matrix_reduced = this->GetKr();
        auto& r_mass_matrix_reduced = this->GetMr();
        auto& r_damping_reduced = this->GetDr();   // Reduced damping matrix

        // r_force_vector_reduced = (prod( rRHS, r_basis ));
        r_force_vector_reduced = prod( trans( r_basis ), rRHS );

        // KRATOS_WATCH( prod( matrix<double>(prod(trans(r_basis_r),r_stiffness_matrix)), r_basis_r))
        r_stiffness_matrix_reduced = prod( matrix< double >( prod( trans( r_basis ),rA ) ), r_basis );
        r_mass_matrix_reduced = prod( matrix< double >( prod( trans( r_basis ),rM ) ), r_basis );

        r_damping_reduced = prod( matrix< double >( prod( trans( r_basis ),rS ) ), r_basis );
        // KRATOS_WATCH(r_force_vector_reduced)
        KRATOS_WATCH(r_stiffness_matrix_reduced)
        KRATOS_WATCH(r_mass_matrix_reduced)
        // KRATOS_WATCH(r_basis)
        // KRATOS_WATCH(r_basis)

        std::cout << "MOR offline solve finished" << std::endl;
		//std::cout << "Compiled succesfully finally after long time"<<std::endl;
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
        return "MorSecondOrderKrylovStrategy";
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

    MorSecondOrderKrylovStrategy(const MorSecondOrderKrylovStrategy &Other){};

    ///@}

}; /* Class MorSecondOrderKrylovStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_SECOND_ORDER_KRYLOV_STRATEGY  defined */
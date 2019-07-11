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

#if !defined(MOR_OFFLINE_STRATEGY)
#define MOR_OFFLINE_STRATEGY

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
 * @class MorOfflineStrategy
 * @ingroup KratosCore
 * @brief This is the linear MOR matrix output strategy
 * @details This strategy builds the K and M matrices and outputs them
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class MorOfflineStrategy
    : public MorOfflineSecondOrderStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorOfflineStrategy);

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
    MorOfflineStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        vector< double > samplingPoints,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pNewLinearSolver, MoveMeshFlag)
    {
        KRATOS_TRY;

        // Saving the scheme
        mpScheme = pScheme;


        // Setting up the default builder and solver
         mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(
            new TBuilderAndSolverType(pNewLinearSolver)); 

        // Saving the linear solver
        mpLinearSolver = pNewLinearSolver;            

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(false);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        // GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(0);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpS = TSparseSpace::CreateEmptyMatrixPointer();
        mpM = TSparseSpace::CreateEmptyMatrixPointer();
        mpRHS = TSparseSpace::CreateEmptyVectorPointer();

        mpAr = TSparseSpace::CreateEmptyMatrixPointer();
        mpMr = TSparseSpace::CreateEmptyMatrixPointer();
        mpRHSr = TSparseSpace::CreateEmptyVectorPointer();
        // mpBasis = TDenseSpace::CreateEmptyMatrixPointer();
        mpBasis = TSparseSpace::CreateEmptyMatrixPointer();
        mpSr = TSparseSpace::CreateEmptyMatrixPointer();

        mSamplingPoints = samplingPoints;

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorOfflineStrategy() override
    {
        Clear();
    }

    /**
     * @brief Set method for the time scheme
     * @param pScheme The pointer to the time scheme considered
     */
    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        mpScheme = pScheme;
    };

    /**
     * @brief Get method for the time scheme
     * @return mpScheme: The pointer to the time scheme considered
     */
    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    /**
     * @brief Set method for the builder and solver
     * @param pNewBuilderAndSolver The pointer to the builder and solver considered
     */
    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    /**
     * @brief Get method for the builder and solver
     * @return mpBuilderAndSolver: The pointer to the builder and solver considered
     */
    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    /**
     * @brief This method sets the flag mInitializeWasPerformed
     * @param InitializePerformedFlag The flag that tells if the initialize has been computed
     */
    void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
    {
        mInitializeWasPerformed = InitializePerformedFlag;
    }

    /**
     * @brief This method gets the flag mInitializeWasPerformed
     * @return mInitializeWasPerformed: The flag that tells if the initialize has been computed
     */
    bool GetInitializePerformedFlag()
    {
        return mInitializeWasPerformed;
    }

    /**
     * @brief This method sets the flag mReformDofSetAtEachStep
     * @param Flag The flag that tells if each time step the system is rebuilt
     */
    void SetReformDofSetAtEachStepFlag(bool Flag)
    {
        mReformDofSetAtEachStep = Flag;
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }

    /**
     * @brief This method returns the flag mReformDofSetAtEachStep
     * @return The flag that tells if each time step the system is rebuilt
     */
    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic informations
     * - 2: Printing linear solver data
     * - 3: Print of debug informations: Echo of stiffness matrix, Dx, b...
     */

    void SetEchoLevel(int Level) override
    {
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * @brief Operation to predict the solution ... if it is not called a trivial predictor is used in which the
    values of the solution step of interest are assumed equal to the old values
     */
    // void Predict() override
    // {
    //     KRATOS_TRY
    //     //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
    //     //if the operations needed were already performed this does nothing
    //     if (mInitializeWasPerformed == false)
    //         Initialize();

    //     //initialize solution step
    //     if (mSolutionStepIsInitialized == false)
    //         InitializeSolutionStep();

    //     TSystemMatrixType& rA  = *mpA;
    //     TSystemMatrixType& rM = *mpM;
    //     TSystemVectorType& rRHS  = *mpRHS;

    //     DofsArrayType& r_dof_set = GetBuilderAndSolver()->GetDofSet();

    //     GetScheme()->Predict(BaseType::GetModelPart(), r_dof_set, rA, rRHS, rRHS);

    //     GetScheme()->Predict(BaseType::GetModelPart(), r_dof_set, rM, rRHS, rRHS);

    //     //move the mesh if needed
    //     if (this->MoveMeshFlag() == true)
    //         BaseType::MoveMesh();

    //     KRATOS_CATCH("")
    // }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY;

        if (mInitializeWasPerformed == false)
        {
            //pointers needed in the solution
            typename TSchemeType::Pointer p_scheme = GetScheme();

            //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (p_scheme->SchemeIsInitialized() == false)
                p_scheme->Initialize(BaseType::GetModelPart());

            //Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ElementsAreInitialized() == false)
                p_scheme->InitializeElements(BaseType::GetModelPart());

            //Initialize The Conditions - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ConditionsAreInitialized() == false)
                p_scheme->InitializeConditions(BaseType::GetModelPart());

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY;

        // if the preconditioner is saved between solves, it
        // should be cleared here.
        GetBuilderAndSolver()->GetLinearSystemSolver()->Clear();

        if (mpA != nullptr)
            SparseSpaceType::Clear(mpA);
        if (mpM != nullptr)
            SparseSpaceType::Clear(mpM);
        if (mpRHS != nullptr)
            SparseSpaceType::Clear(mpRHS);
        if (mpS != nullptr)
            SparseSpaceType::Clear(mpS);

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();
        GetScheme()->Clear();

        mInitializeWasPerformed = false;

        KRATOS_CATCH("");
    }


    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (mSolutionStepIsInitialized == false)
        {
            // std::cout << "system matrices before initialize solution step" << std::endl;
            // KRATOS_WATCH(*mpA)
            // KRATOS_WATCH(*mpM)
            // KRATOS_WATCH(*mpRHS)
            //pointers needed in the solution
            typename TSchemeType::Pointer p_scheme = GetScheme();
            typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

            const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

            //set up the system, operation performed just once unless it is required
            //to reform the dof set at each iteration
            BuiltinTimer system_construction_time;
            if (p_builder_and_solver->GetDofSetIsInitializedFlag() == false ||
                mReformDofSetAtEachStep == true)
            {
                //setting up the list of the DOFs to be solved
                BuiltinTimer setup_dofs_time;
                p_builder_and_solver->SetUpDofSet(p_scheme, BaseType::GetModelPart());
                KRATOS_INFO_IF("Setup Dofs Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_dofs_time.ElapsedSeconds() << std::endl;

                //shaping correctly the system
                BuiltinTimer setup_system_time;
                p_builder_and_solver->SetUpSystem(BaseType::GetModelPart());
                KRATOS_INFO_IF("Setup System Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << setup_system_time.ElapsedSeconds() << std::endl;

                //setting up the Vectors involved to the correct size
                BuiltinTimer system_matrix_resize_time;
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpA, mpRHS, mpRHS,
                                                                 BaseType::GetModelPart());
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpM, mpRHS, mpRHS,
                                                                 BaseType::GetModelPart()); 
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpS, mpRHS, mpRHS,
                                                                 BaseType::GetModelPart()); 

                KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << system_matrix_resize_time.ElapsedSeconds() << std::endl;
            }

            KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << system_construction_time.ElapsedSeconds() << std::endl;

            TSystemMatrixType& rA  = *mpA;
            TSystemMatrixType& rM = *mpM;
            TSystemVectorType& rRHS  = *mpRHS;

            //initial operations ... things that are constant over the Solution Step
            p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(), rA, rRHS, rRHS);
            p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(), rM, rRHS, rRHS);

            //initial operations ... things that are constant over the Solution Step
            p_scheme->InitializeSolutionStep(BaseType::GetModelPart(), rA, rRHS, rRHS);
            p_scheme->InitializeSolutionStep(BaseType::GetModelPart(), rM, rRHS, rRHS);

            mSolutionStepIsInitialized = true;
            // std::cout << "system matrices after initialize solution step" << std::endl;
            // KRATOS_WATCH(*mpA)
            // KRATOS_WATCH(*mpM)
            // KRATOS_WATCH(*mpRHS)
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

        TSystemMatrixType& rA  = *mpA;
        TSystemMatrixType& rM = *mpM;
        TSystemVectorType& rRHS  = *mpRHS;
        TSystemMatrixType& rS  = *mpS;

        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation

        p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), rA, rRHS, rRHS);
        p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), rM, rRHS, rRHS);
        p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), rA, rRHS, rRHS);
        p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), rM, rRHS, rRHS);

        p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), rS, rRHS, rRHS);
        p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), rS, rRHS, rRHS);

        //Cleaning memory after the solution
        p_scheme->Clean();

        //reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            SparseSpaceType::Clear(mpA);
            SparseSpaceType::Clear(mpRHS);
            SparseSpaceType::Clear(mpM);
            SparseSpaceType::Clear(mpS);

            this->Clear();
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY;
        std::cout << "hello! this is where the MOR magic happens" << std::endl;
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        TSystemMatrixType& rA  = *mpA;
        TSystemMatrixType& rM = *mpM;
        TSystemMatrixType& rS  = *mpS;
        TSystemVectorType& rRHS  = *mpRHS;
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
        // mpBasis = SparseSpaceType::CreateEmptyMatrixPointer();
        // auto r_basis = *mpBasis;
        // SparseSpaceType::Resize(r_basis, system_size, reduced_system_size);
        // DenseSpaceType::Resize(r_basis, system_size, reduced_system_size);

        vector< double > aux;
        vector <double > Q_vector;
        
        for( size_t i = 0; i < n_sampling_points; ++i )
        {
            KRATOS_WATCH( mSamplingPoints(i) )
            // KRATOS_WATCH( std::pow( mSamplingPoints(i), 2.0) )
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
        auto& r_basis = *mpBasis;
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
        auto& r_force_vector_reduced = *mpRHSr;
        // mpForceVectorReduced = ZeroMatrix( 3*n_sampling_points );
        r_force_vector_reduced = (prod( rRHS, r_basis ));
        auto& r_stiffness_matrix_reduced = *mpAr;
        auto& r_mass_matrix_reduced = *mpMr;
        auto& r_damping_reduced = *mpSr;   // Reduced damping matrix

        // KRATOS_WATCH( prod( matrix<double>(prod(trans(r_basis_r),r_stiffness_matrix)), r_basis_r))
        r_stiffness_matrix_reduced = prod( matrix< double >( prod( trans( r_basis ),rA ) ), r_basis );
        r_mass_matrix_reduced = prod( matrix< double >( prod( trans( r_basis ),rM ) ), r_basis );

        r_damping_reduced = prod( matrix< double >( prod( trans( r_basis ),rS ) ), r_basis );
        // KRATOS_WATCH(r_force_vector_reduced)
        // KRATOS_WATCH(r_stiffness_matrix_reduced)
        // KRATOS_WATCH(r_mass_matrix_reduced)
        // KRATOS_WATCH(r_basis)
        std::cout << "MOR offline solve finished" << std::endl;
		//std::cout << "Compiled succesfully finally after long time"<<std::endl;
		return true;
        KRATOS_CATCH("");
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        GetBuilderAndSolver()->Check(BaseType::GetModelPart());

        GetScheme()->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH("")
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

    /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    TSystemMatrixType &GetSystemMatrix()
    {
        TSystemMatrixType &mA = *mpA;

        return mA;
    }

    TSystemMatrixType &GetMassMatrix()
    {
        TSystemMatrixType &mM = *mpM;

        return mM;
    }

    TSystemMatrixType &GetDampingMatrix()
    {
        TSystemMatrixType &mS = *mpS;

        return mS;
    }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    TSystemVectorType& GetSystemVector()
    {
        TSystemVectorType& mb = *mpRHS;

        return mb;
    }

    TSystemMatrixType& GetMr()
    {
        TSystemMatrixType &mMr = *mpMr;

        return mMr;
    }

    TSystemMatrixPointerType& pGetMr()
    {
        return mpMr;
    }

    TSystemMatrixType& GetKr()
    {
        TSystemMatrixType &mAr = *mpAr;

        return mAr;
    }

    TSystemMatrixPointerType& pGetKr()
    {
        return mpAr;
    }

    TSystemMatrixType& GetDr()
    {
        TSystemMatrixType &mSr = *mpSr;

        return mSr;
    };

    TSystemMatrixPointerType& pGetDr()
    {
        return mpSr;
    }

    TSystemVectorType& GetRHSr()
    {
        
        TSystemVectorType& mb = *mpRHSr;

        return mb;
    };

    TSystemVectorPointerType& pGetRHSr()
    {
        return mpRHSr;
    }

    TSystemMatrixType& GetBasis()
    {
        TSystemMatrixType &mBasis = *mpBasis;

        return mBasis;
    }

    TSystemMatrixPointerType& pGetBasis()
    {
        return mpBasis;
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
    typename TLinearSolver::Pointer mpLinearSolver; /// The pointer to the linear solver considered
    typename TSchemeType::Pointer mpScheme; /// The pointer to the time scheme employed
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver; /// The pointer to the builder and solver employe

    TSystemVectorPointerType mpRHS; /// The RHS vector of the system of equations
    TSystemMatrixPointerType mpA; /// The Stiffness matrix of the system of equations
    TSystemMatrixPointerType mpM; /// The Mass matrix of the system of equations
    TSystemMatrixPointerType mpS; /// The Damping matrix of the system of equations


    // reduced matrices
    TSystemVectorPointerType mpRHSr; //reduced RHS
    TSystemMatrixPointerType mpAr;
    TSystemMatrixPointerType mpMr;
    // TDenseMatrixPointerType mpBasis;
    TSystemMatrixPointerType mpBasis;
    TSystemMatrixPointerType mpSr;

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

    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     * @param IterationNumber The non linear iteration in the solution loop
     */
    virtual void EchoInfo(const unsigned int IterationNumber)
    {
        TSystemMatrixType& rA  = *mpA;
        TSystemMatrixType& rM = *mpM;
        TSystemVectorType& rRHS  = *mpRHS;
        TSystemMatrixType& rS  = *mpS;

        if (this->GetEchoLevel() == 2) //if it is needed to print the debug info
        {
            KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        }
        else if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
        {
            KRATOS_INFO("LHS") << "SystemMatrix = " << rA << std::endl;
            KRATOS_INFO("Dx")  << "Mass Matrix = " << mpM << std::endl;
            KRATOS_INFO("Sx")  << "Damping Matrix = " << mpS << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        }
        std::stringstream matrix_market_name;
        matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), rA, false);

        std::stringstream matrix_market_mass_name;
        matrix_market_mass_name << "M_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_mass_name.str()).c_str(), rM, false);   

        std::stringstream matrix_market_damping_name;
        matrix_market_name << "S_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_damping_name.str()).c_str(), rS, false);        

        std::stringstream matrix_market_vectname;
        matrix_market_vectname << "RHS_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm.rhs";
        TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), rRHS);
    }

    /**
     * @brief This method prints information after reach the max number of iterations
     * @todo Replace by logger
     */

    virtual void MaxIterationsExceeded()
    {
        if (this->GetEchoLevel() != 0 && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
        {
            std::cout << "***************************************************" << std::endl;
            std::cout << "******* ATTENTION: max iterations exceeded ********" << std::endl;
            std::cout << "***************************************************" << std::endl;
        }
    }

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

    MorOfflineStrategy(const MorOfflineStrategy &Other){};

    ///@}

}; /* Class MorOfflineStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_OFFLINE_STRATEGY  defined */
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

#if !defined(FREQUENCY_RESPONSE_ANALYSIS_STRATEGY)
#define FREQUENCY_RESPONSE_ANALYSIS_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"
#include "utilities/qr_utility.h"

//default builder and solver
#include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"

// Linear solvers
#include "linear_solvers/linear_solver.h"

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
 * @class FrequencyResponseAnalysisStrategy
 * @ingroup KratosCore
 * @brief This is the linear MOR matrix output strategy
 * @details This strategy builds the K and M matrices and outputs them
 * @author Srikkanth Varadharajan
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
        //   class TSolutionSpace
        //   bool TUseDamping
          >
class FrequencyResponseAnalysisStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(FrequencyResponseAnalysisStrategy);

    typedef TUblasSparseSpace<complex> ComplexSparseSpaceType;
    typedef TUblasDenseSpace<complex> ComplexDenseSpaceType;

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef SystemMatrixBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef TDenseSpace DenseSpaceType;

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

    typedef typename ComplexSparseSpaceType::MatrixType TSolutionMatrixType;

    typedef typename ComplexSparseSpaceType::MatrixPointerType TSolutionMatrixPointerType;

    typedef typename ComplexSparseSpaceType::VectorType TSolutionVectorType;

    typedef typename ComplexSparseSpaceType::VectorPointerType TSolutionVectorPointerType;

    typedef ComplexSparseSpaceType TSolutionSpace;

    typedef std::complex<double> ComplexType;


    // typedef UblasSpace<ComplexType, CompressedMatrix, boost::numeric::ublas::vector<std::complex<double>>> TComplexSparseSpaceType;

    // typedef typename TComplexSparseSpaceType::MatrixType TComplexSystemMatrixType;
    typedef LinearSolver<ComplexSparseSpaceType, ComplexDenseSpaceType> ComplexLinearSolverType;
    // typedef TLinearSolverType<std::complex<double>> ComplexLinearSolverType;

    // typedef TLinearSolver<ComplexSpaceType, ComplexLocalSpaceType> ComplexLinearSolverType;
    // typedef ComplexLinearSolverType::Pointer ComplexLinearSolverPointer;
    // typedef TComplexSparseSpace ComplexSparseSpaceType;

    // typename typedef TSparseSpace::MatrixType<ComplexType> ComplexMatrixType;

    // typedef DenseVector<ComplexType> ComplexVectorType;

    // typedef DenseMatrix<ComplexType> ComplexMatrixType;

    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Default constructor -> damped
     * @param rModelPart The model part of the problem
     * @param pScheme The integration schemed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    FrequencyResponseAnalysisStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        // typename TLinearSolver::Pointer pNewLinearSolver,
        typename ComplexLinearSolverType::Pointer pNewcomplexLinearSolver,
        bool MoveMeshFlag = false)
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag)
    {
        KRATOS_TRY;
        // typename ComplexLinearSolverType::Pointer pNewComplexLinearSolver = PastixComplexSolver<ComplexSpaceType,ComplexLocalSpaceType>())

        // Saving the scheme
        mpScheme = pScheme;

        // Setting default linear solver
        mpLinearSolver = typename TLinearSolver::Pointer(new TLinearSolver);

        // Setting up the default builder and solver
         mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(
            new TBuilderAndSolverType(mpLinearSolver)); 

        // Saving the linear solver
        // mpLinearSolver = pNewLinearSolver;            

        // Set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;
        mReformDofSetAtEachStep = false;

        // Tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(false);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(0);

        // Damping is included and the results will be complex
        mUseDamping = true;
        // KRATOS_WATCH(TSolutionSpace)
        mpComplexLinearSolver = pNewcomplexLinearSolver;
        
        mpA = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
        mpK = TSparseSpace::CreateEmptyMatrixPointer();
        mpC = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
        mpM = TSparseSpace::CreateEmptyMatrixPointer();
        mpRHS = ComplexSparseSpaceType::CreateEmptyVectorPointer();
        mpDx = ComplexSparseSpaceType::CreateEmptyVectorPointer();

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~FrequencyResponseAnalysisStrategy() override
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
            TSolutionSpace::Clear(mpA);
        if (mpK != nullptr)
            SparseSpaceType::Clear(mpK);
        if (mpM != nullptr)
            SparseSpaceType::Clear(mpM);
        if (mpC != nullptr)
            TSolutionSpace::Clear(mpC);
        if (mpRHS != nullptr)
            TSolutionSpace::Clear(mpRHS);
        if (mpDx != nullptr)
            TSolutionSpace::Clear(mpDx);

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
            typename TSchemeType::Pointer p_scheme = GetScheme();
            typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

            const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();
            //set up the system, operation performed just once unless it is required
            //to reform the dof set at each iteration
            TSolutionMatrixType& rA  = *mpA;
            TSystemMatrixType& rK  = *mpK;
            TSystemMatrixType& rM = *mpM;
            TSolutionMatrixType& rC  = *mpC;
            TSolutionVectorType& rRHS  = *mpRHS;
            // TSolutionVectorType& rDx  = *mpDx;

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
                TSystemVectorPointerType tmp_RHS = TSparseSpace::CreateEmptyVectorPointer();
                TSystemVectorType& r_tmp_RHS = *tmp_RHS;
                TSystemMatrixPointerType tmpC = TSparseSpace::CreateEmptyMatrixPointer();
                TSystemMatrixType& r_tmp_C = *tmpC;
                BuiltinTimer system_matrix_resize_time;
                // p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpA, mpDx, mpRHS,
                //                                                  BaseType::GetModelPart());
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpK, tmp_RHS, tmp_RHS,
                                                                 BaseType::GetModelPart());
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpM, tmp_RHS, tmp_RHS,
                                                                 BaseType::GetModelPart()); 
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, tmpC, tmp_RHS, tmp_RHS,
                                                                 BaseType::GetModelPart());
                if (mpDx->size() != p_builder_and_solver->GetEquationSystemSize())
                    mpDx->resize(p_builder_and_solver->GetEquationSystemSize(), false);


                KRATOS_INFO_IF("System Matrix Resize Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                    << system_matrix_resize_time.ElapsedSeconds() << std::endl;

                //set up system matrices
                TSystemVectorType tmp(rK.size1(), 0.0);

                p_builder_and_solver->BuildStiffnessMatrix(p_scheme, BaseType::GetModelPart(), rK, tmp);
                p_builder_and_solver->ApplyDirichletConditions(p_scheme, BaseType::GetModelPart(), rK, tmp, r_tmp_RHS);

                p_builder_and_solver->BuildMassMatrix(p_scheme, BaseType::GetModelPart(), rM, tmp);
                p_builder_and_solver->ApplyDirichletConditionsForMassMatrix(p_scheme, BaseType::GetModelPart(), rM);

                p_builder_and_solver->BuildDampingMatrix(p_scheme, BaseType::GetModelPart(), r_tmp_C, tmp);
                p_builder_and_solver->ApplyDirichletConditionsForDampingMatrix(p_scheme, BaseType::GetModelPart(), r_tmp_C);

                p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), r_tmp_RHS);

                rC = TSolutionMatrixType(r_tmp_C);
                rC *= complex(0,1);

                rRHS = TSolutionVectorType(r_tmp_RHS);
                KRATOS_WATCH(rC)
                // mpC *= complex(0,1);
                std::cout << "rK vorher" << std::endl;
                KRATOS_WATCH(rK)
                // rA = TSolutionMatrixType(mpK->size1(), mpK->size2(), mpK->nnz());
                rA = TSolutionMatrixType(rK);
                std::cout << "rK nachher" << std::endl;
                KRATOS_WATCH(rK)
                KRATOS_WATCH(rA)
            }

            KRATOS_INFO_IF("System Construction Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << system_construction_time.ElapsedSeconds() << std::endl;

            //initial operations ... things that are constant over the Solution Step
            // p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(), rA, rDx, rRHS);
            // p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(), rK, rDx, rRHS);
            // p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(), rM, rDx, rRHS);
            // p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(), rC, rDx, rRHS);

            //initial operations ... things that are constant over the Solution Step
            // p_scheme->InitializeSolutionStep(BaseType::GetModelPart(), rA, rDx, rRHS);
            // p_scheme->InitializeSolutionStep(BaseType::GetModelPart(), rK, rDx, rRHS);
            // p_scheme->InitializeSolutionStep(BaseType::GetModelPart(), rM, rDx, rRHS);
            // p_scheme->InitializeSolutionStep(BaseType::GetModelPart(), rC, rDx, rRHS);

            mSolutionStepIsInitialized = true;
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

        // TSolutionMatrixType& rA  = *mpA;
        // TSystemMatrixType& rK  = *mpK;
        // TSystemMatrixType& rM = *mpM;
        // TSystemMatrixType& rC  = *mpC;
        // TSystemVectorType& rRHS  = *mpRHS;
        // TSolutionVectorType& rDx  = *mpDx;

        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation
        //TODO!
        // p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), rA, rDx, rRHS);
        // p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), rK, rDx, rRHS);
        // p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), rM, rDx, rRHS);
        // p_scheme->FinalizeSolutionStep(BaseType::GetModelPart(), rC, rDx, rRHS);
        // p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), rA, rDx, rRHS);
        // p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), rK, rDx, rRHS);
        // p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), rM, rDx, rRHS);
        // p_builder_and_solver->FinalizeSolutionStep(BaseType::GetModelPart(), rC, rDx, rRHS);

        //Cleaning memory after the solution
        p_scheme->Clean();

        //reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            // TComplexSparseSpaceType::Clear(mpA);
            TSolutionSpace::Clear(mpA);
            SparseSpaceType::Clear(mpK);
            SparseSpaceType::Clear(mpM);
            TSolutionSpace::Clear(mpC);
            TSolutionSpace::Clear(mpRHS);
            TSolutionSpace::Clear(mpDx);

            this->Clear();
        }

        KRATOS_CATCH("");
    }

    void AssignVariables(TSystemVectorType& rDisplacement, int step=0)
    {
        auto& r_model_part = BaseType::GetModelPart();
        for( auto& node : r_model_part.Nodes() )
        {
            ModelPart::NodeType::DofsContainerType& rNodeDofs = node.GetDofs();
            
            for( auto it_dof = std::begin(rNodeDofs); it_dof != std::end(rNodeDofs); it_dof++ )
            {
                if( !it_dof->IsFixed() )
                {
                    it_dof->GetSolutionStepValue(step) = rDisplacement(it_dof->EquationId());
                }
                else
                {
                    it_dof->GetSolutionStepValue(step) = 0.0;
                }
            }
        }
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY;
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        TSolutionMatrixType& rA = *mpA;
        // TComplexSystemMatrixType& rA = *mpA;
        TSystemMatrixType& rK = *mpK;
        TSystemMatrixType& rM = *mpM;
        TSolutionMatrixType& rC  = *mpC;
        TSolutionVectorType& rRHS = *mpRHS;
        TSolutionVectorType& rDx = *mpDx;
        
        auto& r_process_info = BaseType::GetModelPart().GetProcessInfo();
        double excitation_frequency = r_process_info[FREQUENCY];

        //Building the dynamic stiffnes matrix
        // rA = rK + ComplexType(0, excitation_frequency) * rC - ( std::pow(excitation_frequency, 2.0 ) * rM );
        // complex i(0,1);

        // rA = rK + excitation_frequency*rC - std::pow(excitation_frequency, 2.0)*rM;
        // TSystemMatrixType tmp = excitation_frequency * rM;
        rA = - (std::pow(excitation_frequency, 2.0) * rM);
        rA += excitation_frequency*rC;
        rA += rK;
        // rA -= tmp;

        // KRATOS_WATCH(rK)
        // KRATOS_WATCH(i)
        // rA = TSolutionMatrixType(complex(rK,rK));
        // rA = prod(rK,)
        // TSolutionMatrixType tmp = rK;
        // rA = rA+tmp;
        // boost::numeric::ublas::compressed_matrix<std::complex<double>> tmp(rK);
        // tmp = tmp * complex(0,1);
        // KRATOS_WATCH(&tmp.index1_data())
        // KRATOS_WATCH(&tmp.index2_data())
        // KRATOS_WATCH(&rK.index1_data())
        // KRATOS_WATCH(&rK.index2_data())

        // KRATOS_WATCH(tmp.nnz())
        // KRATOS_WATCH(rK.nnz())


        // KRATOS_WATCH(tmp)
        // KRATOS_WATCH(rK)
        // rA = complex(0,rK);
        // rA = rK * 2;
        // rA = conj(rA);
        // rA = rA*i;
        // rA = rK - complex(0,1) * rC + ( std::pow(excitation_frequency, 2.0 ) * rM );
        std::cout << "solve" << std::endl;
        // KRATOS_WATCH(rK)
        // KRATOS_WATCH(rA)
        // p_builder_and_solver->GetLinearSystemSolver()->Solve( rA, rDx, rRHS );
        // TSolutionVectorType tmpRHS = TSolutionVectorType(rRHS);
        // KRATOS_WATCH(rA.index1_data_())
        // std::cout << *rA.index1_data() << std::endl;
        KRATOS_WATCH(rDx)
        KRATOS_WATCH(rRHS)
        mpComplexLinearSolver->Initialize( rA, rDx, rRHS);
        mpComplexLinearSolver->Solve( rA, rDx, rRHS);

        //assigning the computed displacement
        // AssignVariables(rDx);

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
        TSystemMatrixType &mA = *mpK;

        return mA;
    }

    TSystemMatrixType &GetMassMatrix()
    {
        TSystemMatrixType &mM = *mpM;

        return mM;
    }

    // TSystemMatrixType &GetDampingMatrix()
    // {
    //     TSystemMatrixType &mS = *mpC;

    //     return mS;
    // }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    // TSystemVectorType& GetSystemVector()
    // {
    //     TSystemVectorType& mb = *mpRHS;

    //     return mb;
    // }

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
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver; /// The pointer to the builder and solver
    ComplexLinearSolverType::Pointer mpComplexLinearSolver;

    TSolutionVectorPointerType mpRHS; /// The RHS vector
    // TSystemVectorPointerType mpDx; /// The solution vector
    // TSystemMatrixPointerType mpA; /// The system matrix
    TSolutionVectorPointerType mpDx;
    TSolutionMatrixPointerType mpA;
    TSystemMatrixPointerType mpK; /// The stiffness matrix
    TSystemMatrixPointerType mpM; /// The mass matrix
    TSolutionMatrixPointerType mpC; /// The damping matrix

    bool mReformDofSetAtEachStep;

    bool mSolutionStepIsInitialized; /// Flag to set as initialized the solution step

    bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

    bool mUseDamping;

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     * @param IterationNumber The non linear iteration in the solution loop
     */
    virtual void EchoInfo(const unsigned int IterationNumber)
    {
        TSystemMatrixType& rK  = *mpK;
        TSystemMatrixType& rM = *mpM;
        TSolutionVectorType& rRHS  = *mpRHS;
        TSolutionMatrixType& rC  = *mpC;

        if (this->GetEchoLevel() == 2) //if it is needed to print the debug info
        {
            KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        }
        else if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
        {
            KRATOS_INFO("K") << "SystemMatrix = " << rK << std::endl;
            KRATOS_INFO("M")  << "Mass Matrix = " << mpM << std::endl;
            KRATOS_INFO("C")  << "Damping Matrix = " << mpC << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << rRHS << std::endl;
        }
        std::stringstream matrix_market_name;
        matrix_market_name << "K_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), rK, false);

        std::stringstream matrix_market_mass_name;
        matrix_market_mass_name << "M_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_mass_name.str()).c_str(), rM, false);   

        std::stringstream matrix_market_damping_name;
        matrix_market_name << "C_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_damping_name.str()).c_str(), rC, false);        

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

    FrequencyResponseAnalysisStrategy(const FrequencyResponseAnalysisStrategy &Other){};

    ///@}

}; /* Class MorOfflineStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_OFFLINE_STRATEGY  defined */
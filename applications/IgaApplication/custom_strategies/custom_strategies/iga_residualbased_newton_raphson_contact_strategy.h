// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ /
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_IGA_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_STRATEGY )
#define  KRATOS_IGA_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"

/* Strategies */
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"


// edxtra for couple contact geometries
// #include "geometries/nurbs_coupling_geometry_2d.h"

#include "custom_conditions/support_contact_2D_condition.h"

namespace Kratos {

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
 * @class IgaResidualBasedNewtonRaphsonContactStrategy
 * @ingroup ContactStructuralMechanicsApplication
 * @brief  Contact Newton Raphson class
 * @details This class is a specialization of the Newton Raphson strategy with some custom modifications for contact problems
 * @author Vicente Mataix Ferrandiz
*/
template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class IgaResidualBasedNewtonRaphsonContactStrategy :
    public ResidualBasedNewtonRaphsonStrategy< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    /// Definizione dei puntatori di classe
    KRATOS_CLASS_POINTER_DEFINITION(IgaResidualBasedNewtonRaphsonContactStrategy);

    /// Alias per la classe base
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef SolvingStrategy<TSparseSpace, TDenseSpace> SolvingStrategyType;

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    //------------------------
    // EVENTUALLY REMOVE, GEOEMTRY STUFF
    typedef Node                                             NodeType;
    typedef Geometry<NodeType>                                  GeometryType;
    typedef GeometryType::Pointer                               GeometryPointerType;
    typedef typename GeometryType::GeometriesArrayType          GeometriesArrayType;
    typedef typename GeometryType::CoordinatesArrayType         CoordinatesArrayType;
    typedef typename GeometryType::IntegrationPointsArrayType   IntegrationPointsArrayType;

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    typedef typename Properties::Pointer PropertiesPointerType;

    typedef typename ModelPart::ElementsContainerType ElementsContainerType;
    typedef typename ModelPart::ConditionsContainerType ConditionsContainerType;

    using PointType = Node;
    

    /// Costruttore di default che richiama il costruttore della classe base
    explicit IgaResidualBasedNewtonRaphsonContactStrategy(ModelPart& rModelPart)
        : BaseType(rModelPart)
    {
    }

    /// Costruttore di default con parametri che richiama il costruttore della classe base
    explicit IgaResidualBasedNewtonRaphsonContactStrategy(ModelPart& rModelPart, Parameters ThisParameters)
        : BaseType(rModelPart, ThisParameters)
    {
    }

    /**
     * Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit IgaResidualBasedNewtonRaphsonContactStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, MaxIterations, CalculateReactions,
                   ReformDofSetAtEachStep, MoveMeshFlag)
    {
    }

    /**
     * @brief Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit IgaResidualBasedNewtonRaphsonContactStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions,
                   ReformDofSetAtEachStep, MoveMeshFlag)
    {
    }

    /**
     * @brief Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit IgaResidualBasedNewtonRaphsonContactStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters igaContactParameters,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions,
                   ReformDofSetAtEachStep, MoveMeshFlag)
    {
        mIgaContactParameters = igaContactParameters;
    }

    /**
     * @brief Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    KRATOS_DEPRECATED_MESSAGE("Constructor deprecated, please use the constructor without linear solver")
    explicit IgaResidualBasedNewtonRaphsonContactStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions,
                   ReformDofSetAtEachStep, MoveMeshFlag)
    {
    }

    /**
     * Constructor with Parameters
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param Settings Settings used in the strategy
     */
    IgaResidualBasedNewtonRaphsonContactStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        Parameters Settings)
        : BaseType(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, Settings)
    {
    }

    /**
     * @brief Constructor specifying the builder and solver and using Parameters
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param Settings Settings used in the strategy
     */
    // IgaResidualBasedNewtonRaphsonContactStrategy(
    //     ModelPart& rModelPart,
    //     typename TSchemeType::Pointer pScheme,
    //     typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
    //     typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
    //     Parameters Settings)
    //     : BaseType(rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, Settings)
    // {
    // }

    /**
     * @brief Constructor specifying the builder and solver and using Parameters
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param Parameters Settings used in the strategy
     */
    KRATOS_DEPRECATED_MESSAGE("Constructor deprecated, please use the constructor without linear solver")
    IgaResidualBasedNewtonRaphsonContactStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters Settings)
        : BaseType(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, pNewBuilderAndSolver, Settings)
    {
    }

    /// Distruttore di default
    ~IgaResidualBasedNewtonRaphsonContactStrategy() override = default;

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        // Pointers needed in the solution

        typename TSchemeType::Pointer p_scheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        ModelPart& r_model_part = BaseType::GetModelPart();

        // Set up the system, operation performed just once unless it is required
        // to reform the dof set at each iteration
        BuiltinTimer system_construction_time;
        if (!p_builder_and_solver->GetDofSetIsInitializedFlag() || BaseType::mReformDofSetAtEachStep)
        {
            // Setting up the list of the DOFs to be solved
            BuiltinTimer setup_dofs_time;
            p_builder_and_solver->SetUpDofSet(p_scheme, r_model_part);
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "Setup Dofs Time: "
                << setup_dofs_time << std::endl;
            

            // Shaping correctly the system
            BuiltinTimer setup_system_time;
            p_builder_and_solver->SetUpSystem(r_model_part);
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "Setup System Time: "
                << setup_system_time << std::endl;

            // Setting up the Vectors involved to the correct size
            BuiltinTimer system_matrix_resize_time;
            p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, BaseType::mpA, BaseType::mpDx, BaseType::mpb,
                                                                r_model_part);

            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "System Matrix Resize Time: "
                << system_matrix_resize_time << std::endl;

        }

        KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", BaseType::GetEchoLevel() > 0) << "System Construction Time: "
            << system_construction_time << std::endl;


        TSystemMatrixType& rA  = *BaseType::mpA;
        TSystemVectorType& rDx = *BaseType::mpDx;
        TSystemVectorType& rb  = *BaseType::mpb;

        // Initial operations ... things that are constant over the Solution Step
        p_builder_and_solver->InitializeSolutionStep(r_model_part, rA, rDx, rb);


        // Initial operations ... things that are constant over the Solution Step
        p_scheme->InitializeSolutionStep(r_model_part, rA, rDx, rb);
        
        // Initialisation of the convergence criteria
        if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag())
        {
            TSparseSpace::SetToZero(rb);
            p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
        }

        BaseType::mpConvergenceCriteria->InitializeSolutionStep(r_model_part, p_builder_and_solver->GetDofSet(), rA, rDx, rb);

        if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag()) {
            TSparseSpace::SetToZero(rb);
        }

        KRATOS_CATCH("");
    }


    bool SolveSolutionStep() override
    {
        // Pointers needed in the solution
        ModelPart& r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer p_scheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        auto& r_dof_set = p_builder_and_solver->GetDofSet();

        TSystemMatrixType& rA  = *BaseType::mpA;
        TSystemVectorType& rDx = *BaseType::mpDx;
        TSystemVectorType& rb  = *BaseType::mpb;

        //initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number = 1;

        // UpdateContactPairs(&r_model_part);
        // BaseType::Initialize();
        // InitializeSolutionStep();

        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool residual_is_updated = false;
        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        BaseType::mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);


        bool is_converged = BaseType::mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);



        // Function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            if (BaseType::mUseOldStiffnessInFirstIteration){
                p_builder_and_solver->BuildAndSolveLinearizedOnPreviousIteration(p_scheme, r_model_part, rA, rDx, rb,BaseType::MoveMeshFlag());
            } else {
                p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
            }
        } else {
            TSparseSpace::SetToZero(rDx);  // Dx = 0.00;
            TSparseSpace::SetToZero(rb);

            p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
        }

        // Debugging info
        BaseType::EchoInfo(iteration_number);

        // Updating the results stored in the database
        BaseType::UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

        p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
        BaseType::mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

        if (is_converged) {
            if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag()) {
                TSparseSpace::SetToZero(rb);

                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
            }

            is_converged = BaseType::mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
        }

        //Iteration Cycle... performed only for NonLinearProblems
        while (is_converged == false &&
               iteration_number++ < BaseType::mMaxIterationNumber)
        {
            //setting the number of iteration
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
            BaseType::mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            // InitializeSolutionStep();

            is_converged = BaseType::mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

            //call the linear system solver to find the correction mDx for the
            //it is not called if there is no system to solve
            if (SparseSpaceType::Size(rDx) != 0)
            {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
                {
                    if (BaseType::GetKeepSystemConstantDuringIterations() == false)
                    {
                        //A = 0.00;
                        TSparseSpace::SetToZero(rA);
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    }
                    else
                    {
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    }
                }
                else
                {
                    TSparseSpace::SetToZero(rDx);
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                }
            }
            else
            {
                KRATOS_WARNING("NO DOFS") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Debugging info
            BaseType::EchoInfo(iteration_number);

            // Updating the results stored in the database
            BaseType::UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
            BaseType::mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            residual_is_updated = false;

            if (is_converged == true)
            {
                if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                    residual_is_updated = true;
                }

                is_converged = BaseType::mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
            }
        }

        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= BaseType::mMaxIterationNumber) {
            BaseType::MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / "
                << BaseType::mMaxIterationNumber << " iterations" << std::endl;
        }

        //recalculate residual if needed
        //(note that some convergence criteria need it to be recalculated)
        if (residual_is_updated == false)
        {
            // NOTE:
            // The following part will be commented because it is time consuming
            // and there is no obvious reason to be here. If someone need this
            // part please notify the community via mailing list before uncommenting it.
            // Pooyan.

            //    TSparseSpace::SetToZero(mb);
            //    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, mb);
        }

        //calculate reactions if required
        if (BaseType::mCalculateReactionsFlag == true)
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);

        return is_converged;
    }

    void UpdateContactPairs(ModelPart* r_model_part) 
    {
        std::string contact_model_part_name = "ContactInterface";

        // Model& mainModel = r_model_part->GetModel();

        ModelPart& contact_model_part = r_model_part->HasSubModelPart(contact_model_part_name)
                                      ? r_model_part->GetSubModelPart(contact_model_part_name)
                                      : r_model_part->CreateSubModelPart(contact_model_part_name);


        contact_model_part.RemoveConditionsFromAllLevels();

        //-----------------------------------------------------------------------------------------------
        for (IndexType contact_scenario_id = 0; contact_scenario_id < mIgaContactParameters.size(); contact_scenario_id++) {
            // Obtain SLAVE interface b_reps
            const std::string slave_model_part_name = mIgaContactParameters[contact_scenario_id]["contact_parameters"]["slave_model_part"]
                                                                            ["sub_model_part_name"].GetString();


            KRATOS_ERROR_IF_NOT(r_model_part->HasSubModelPart(slave_model_part_name)) << "ERROR: SLAVE MODEL PART " 
                                                     << slave_model_part_name << "NOT CREATED" << std::endl; 

            ModelPart& slave_model_part = r_model_part->GetSubModelPart(slave_model_part_name);
            GeometriesArrayType geometry_list_slave;
            GetCadGeometryList(geometry_list_slave, slave_model_part, mIgaContactParameters[contact_scenario_id]["contact_parameters"]
                                                                                                                ["slave_model_part"]);

            const IndexType slave_property_id = mIgaContactParameters[contact_scenario_id]["contact_parameters"]["slave_model_part"]
                                                                                ["property_id"].GetInt();

            Properties::Pointer p_prop_slave = slave_model_part.pGetProperties(slave_property_id);

            // Obtain MASTER interface b_reps
            const std::string master_model_part_name = mIgaContactParameters[contact_scenario_id]["contact_parameters"]["master_model_part"]
                                                                                                 ["sub_model_part_name"].GetString();

            KRATOS_ERROR_IF_NOT(r_model_part->HasSubModelPart(master_model_part_name)) << "ERROR: MASTER MODEL PART " 
                                                    << master_model_part_name << "NOT CREATED" << std::endl; 

            ModelPart& master_model_part = r_model_part->GetSubModelPart(master_model_part_name);
            GeometriesArrayType geometry_list_master;
            GetCadGeometryList(geometry_list_master, master_model_part, mIgaContactParameters[contact_scenario_id]["contact_parameters"]
                                                                                                                  ["master_model_part"]);

            const IndexType master_property_id = mIgaContactParameters[contact_scenario_id]["contact_parameters"]["master_model_part"]
                                                                                                                 ["property_id"].GetInt();
            Properties::Pointer p_prop_master = master_model_part.pGetProperties(master_property_id);

            auto couplingGeometry = Kratos::make_shared<NurbsCouplingGeometry2D<PointType, PointerVector<NodeType>>>(geometry_list_slave, geometry_list_master);
            
            
            //---------------------------------------------------------------
            // ÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇ
            GeometriesArrayType geometries;
            SizeType shape_function_derivatives_order = 1;
            if (mIgaContactParameters[contact_scenario_id].Has("shape_function_derivatives_order")) {
                shape_function_derivatives_order = mIgaContactParameters[contact_scenario_id]["shape_function_derivatives_order"].GetInt();
            }
            else {
                // KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 4)
                //     << "shape_function_derivatives_order is not provided and thus being considered as 1. " << std::endl;
            }
            std::string quadrature_method = mIgaContactParameters[contact_scenario_id].Has("quadrature_method")
            ? mIgaContactParameters[contact_scenario_id]["integration_rule"].GetString()
            : "GAUSS";

            IntegrationInfo integration_info = couplingGeometry->GetDefaultIntegrationInfo();
            for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
                if (quadrature_method == "GAUSS") {
                    integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
                }
                else if (quadrature_method == "EXTENDED_GAUSS") {
                    integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS);
                }
                else if (quadrature_method == "GRID") {
                    integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GRID);
                }
                // EVENTUALLY ADD THE PIECEWISE INTEGRATION
                else {
                    KRATOS_INFO("CreateQuadraturePointGeometries") << "Quadrature method: " << quadrature_method
                        << " is not available. Available options are \"GAUSS\" and \"GRID\". Default quadrature method is being considered." << std::endl;
                }
            }

            couplingGeometry->CreateQuadraturePointGeometries(geometries, shape_function_derivatives_order, integration_info);

            // KRATOS_WATCH("NOTHING HERE")

            SizeType id = 1;
            // if (contact_model_part.GetRootModelPart().Conditions().size() > 0)
            //     id = contact_model_part.GetRootModelPart().Conditions().back().Id() + 1;
            if (contact_model_part.GetRootModelPart().Conditions().size() > 0)
                id = contact_model_part.GetRootModelPart().Conditions().back().Id() + 1;
            
            KRATOS_ERROR_IF_NOT(mIgaContactParameters[contact_scenario_id].Has("name"))
                << "\"name\" need to be specified." << std::endl;
            std::string name = mIgaContactParameters[contact_scenario_id]["name"].GetString();
            
            this->CreateConditions(
                        geometries.ptr_begin(), geometries.ptr_end(),
                        contact_model_part, name, id, p_prop_master, p_prop_slave);
        }
    }



     void CreateConditions(
        typename GeometriesArrayType::ptr_iterator rGeometriesBegin,
        typename GeometriesArrayType::ptr_iterator rGeometriesEnd,
        ModelPart& rModelPart,
        std::string& rConditionName,
        SizeType& rIdCounter,
        PropertiesPointerType pPropMaster,
        PropertiesPointerType pPropSlave) const
    {
        // const SupportContact2DCondition rReferenceCondition = SupportContact2DCondition();
        // const Condition& rReferenceCondition = KratosComponents<Condition>::Get(rConditionName);

        ModelPart::ConditionsContainerType new_condition_list;

        // KRATOS_INFO_IF("CreateConditions", mEchoLevel > 2)
        //     << "Creating conditions of type " << rConditionName
        //     << " in " << rModelPart.Name() << "-SubModelPart." << std::endl;

        for (auto it = rGeometriesBegin; it != rGeometriesEnd; ++it)
        {
            new_condition_list.push_back(
                Kratos::make_intrusive<SupportContact2DCondition>(
                rIdCounter, (*it), pPropMaster, pPropSlave));
            for (SizeType i = 0; i < (*it)->size(); ++i) {
                // rModelPart.AddNode((*it)->pGetPoint(i));
                rModelPart.Nodes().push_back((*it)->pGetPoint(i));
            }
            rIdCounter++;
        }

        rModelPart.AddConditions(new_condition_list.begin(), new_condition_list.end());

        EntitiesUtilities::InitializeEntities<Condition>(rModelPart);

    }

    ///@name CAD functionalities
    ///@{

    void GetCadGeometryList(
        GeometriesArrayType& rGeometryList,
        ModelPart& rModelPart,
        const Parameters rParameters) const
    {

        static int starting_brep_ids;

        // get surrogate model part name
        // std::string surrogate_model_part_name;
        // if (!mParameters.Has("surrogate_model_part_name")) surrogate_model_part_name = "surrogate_model_part";
        // else {
        //     surrogate_model_part_name = mParameters["surrogate_model_part_name"].GetString();
        // }

        if (rParameters.Has("brep_id")) {
            rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_id"].GetInt()));
        }
        if (rParameters.Has("brep_ids")) {
            for (SizeType i = 0; i < rParameters["brep_ids"].size(); ++i) {
                rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_ids"][i].GetInt()));
            }
            // int lastIndex
            starting_brep_ids = rParameters["brep_ids"][rParameters["brep_ids"].size()-1].GetInt() + 1;

            // SBM THINGS
            // // OUTER
            // std::string conditionName = rParameters["iga_model_part"].GetString();
            // if (conditionName.rfind("SBM", 0) == 0) { 
            //     ModelPart& surrogateModelPart_outer = mpModel->GetModelPart(surrogate_model_part_name +"_outer");
            //     if (surrogateModelPart_outer.Conditions().size() > 0) {
            //         // 2D
            //         if (surrogateModelPart_outer.GetCondition(1).GetGeometry().size() == 2) {
            //             int sizeSurrogateLoop_outer = surrogateModelPart_outer.Nodes().size();
            //             for (int j = 0; j < (sizeSurrogateLoop_outer-1); ++j) {
            //                 // Add the brep_ids of the internal boundary for SBMLaplacianCondition
            //                 rGeometryList.push_back(rModelPart.pGetGeometry(starting_brep_ids));
            //                 starting_brep_ids++;
            //             }
                        
            //         } else { // 3D
            //             int sizeSurrogateLoop = surrogateModelPart_outer.Conditions().size(); //lastSurrogateNode.Id() - firstSurrogateNode.Id() + 1 ;
            //             for (SizeType j = 0; j < sizeSurrogateLoop-1; ++j) {
            //                 // Add the brep_ids of the internal boundary for SBMLaplacianCondition
            //                 rGeometryList.push_back(rModelPart.pGetGeometry(starting_brep_ids));
            //                 starting_brep_ids++;
            //             }
            //         }     
            //     }
            // }
        }
        // else { // SBM THINGS
        //     // INNER   
        //     ModelPart& surrogateModelPart_inner = mpModel->GetModelPart(surrogate_model_part_name + "_inner");
        //     if (surrogateModelPart_inner.Elements().size() > 0) {

        //         int sizeSurrogateLoop = surrogateModelPart_inner.Nodes().size();

        //         for (int iel = 1; iel < surrogateModelPart_inner.Elements().size()+1; iel++) {
        //             // Each element in the surrogate_model_part represents a surrogate boundary loop. First "node" is the initial ID of the first surrogate node and
        //             // the second "node" is the last surrogate node of that loop. (We have done this in the case we have multiple surrogate boundaries and 1 model part)
        //             Node& firstSurrogateNode = surrogateModelPart_inner.pGetElement(iel)->GetGeometry()[0]; // Element 1 because is the only surrogate loop
        //             Node& lastSurrogateNode = surrogateModelPart_inner.pGetElement(iel)->GetGeometry()[1];  // Element 1 because is the only surrogate loop
        //             int sizeSurrogateLoop = lastSurrogateNode.Id() - firstSurrogateNode.Id() + 1 ;

        //             for (SizeType j = 0; j < sizeSurrogateLoop; ++j) {
        //                 // Add the brep_ids of the internal boundary for SBMLaplacianCondition
        //                 rGeometryList.push_back(rModelPart.pGetGeometry(starting_brep_ids));
        //                 starting_brep_ids++;
        //             }
        //         }
        //     } else { // 3D
        //         int sizeSurrogateLoop = surrogateModelPart_inner.Conditions().size(); //lastSurrogateNode.Id() - firstSurrogateNode.Id() + 1 ;
        //         for (SizeType j = 0; j < sizeSurrogateLoop; ++j) {
        //             // Add the brep_ids of the internal boundary for SBMLaplacianCondition
        //             rGeometryList.push_back(rModelPart.pGetGeometry(starting_brep_ids));
        //             starting_brep_ids++;
        //         }
        //     }
        // }
        if (rParameters.Has("brep_name")) {
            rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_name"].GetString()));
        }
        if (rParameters.Has("brep_names")) {
            for (SizeType i = 0; i < rParameters["brep_names"].size(); ++i) {
                rGeometryList.push_back(rModelPart.pGetGeometry(rParameters["brep_names"][i].GetString()));
            }
        }
        KRATOS_ERROR_IF(rGeometryList.size() == 0)
            << "Empty geometry list. Either \"brep_id\", \"brep_ids\", \"brep_name\" or \"brep_names\" are the possible options." << std::endl;
    }



    /// Turn back information as a string.
    std::string Info() const override
    {
        return "IgaResidualBasedNewtonRaphsonContactStrategy";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    protected: 
    /**
     * Copy constructor.
     */

    IgaResidualBasedNewtonRaphsonContactStrategy(const IgaResidualBasedNewtonRaphsonContactStrategy &Other){};

    private:
    
    Parameters mIgaContactParameters;

    ///@}

}; /* Class IgaResidualBasedNewtonRaphsonContactStrategy */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  // namespace Kratos


#endif /* KRATOS_IGA_RESIDUALBASED_NEWTON_RAPHSON_CONTACT_STRATEGY  defined */
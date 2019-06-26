//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


#if !defined(KRATOS_MPM_STRATEGY )
#define  KRATOS_MPM_STRATEGY

/* System includes */
#include <set>

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/kratos_flags.h"
#include "geometries/geometry.h"
#include "includes/element.h"
#include "boost/numeric/ublas/matrix.hpp"

// Application includes
#include "particle_mechanics_application_variables.h"

// Geometry utilities
#include "utilities/geometry_utilities.h"

// Custom includes
#include "custom_strategies/schemes/MPM_residual_based_bossak_scheme.hpp"
#include "custom_strategies/strategies/MPM_residual_based_newton_raphson_strategy.hpp"
#include "custom_utilities/mpm_search_element_utility.h"
#include "custom_utilities/mpm_particle_generator_utility.h"

// Core includes
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */
/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.
Detail class definition.

 */
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver, //= LinearSolver<TSparseSpace,TDenseSpace>
         unsigned int TDim>
class MPMStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    /**@name Type Definitions */
    /*@{ */
    //		typedef std::set<Dof::Pointer,ComparePDof> DofSetType;

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef typename TSparseSpace::DataType TDataType;
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TSparseSpace::MatrixPointerType TSystemMatrixPointerType;
    typedef typename TSparseSpace::VectorPointerType TSystemVectorPointerType;

    typedef Node < 3 > NodeType;
    typedef Geometry<NodeType> GeometryType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    typedef Scheme<TSparseSpace, TDenseSpace> TSchemeType;
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> TBuilderAndSolverType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > SolvingStrategyType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MPMStrategy);

    typedef typename ModelPart::DofType TDofType;
    typedef typename ModelPart::DofsArrayType DofsArrayType;

    typedef typename PointerVectorSet<TDofType, IndexedObject>::iterator DofIterator;
    typedef typename PointerVectorSet<TDofType, IndexedObject>::const_iterator DofConstantIterator;


    /**
     * @brief Default constructor of MPM Strategy.
     * @details
     * The grid model part contains all the information about ID, geometry and initial/boundary conditions
     * of the computational mesh.
     * In the costructor of time scheme the model part of material points is defined as:
     *
     * STEP 1:
     * The nodes, properties and process info of grid_model_part are assigned to the material points' mdpa.
     *
     * STEP 2:
     * loop over grid elements to evaluate:
     * - the MP integration weight and MP mass
     * - rGeo : connectivity of the grid element
     * - shape function values of the integration points of the grid element
     *
     * STEP 3:
     * loop over the integration points of a grid element to
     * - create a new MP element
     * - evaluate xg which is the coordinate of the integration point
     * - to assign all MP variables
     *
    */

    /*@{ */

    MPMStrategy(ModelPart& grid_model_part, ModelPart& initial_model_part, ModelPart& mpm_model_part, typename TLinearSolver::Pointer plinear_solver,
        std::string SolutionType = "static", int MaxIteration = 10, bool ComputeReaction = false, bool BlockBuilder = false,
        bool IsMixedFormulation = false, bool MoveMeshFlag = false)
        : SolvingStrategyType(grid_model_part, MoveMeshFlag), mr_grid_model_part(grid_model_part), mr_initial_model_part(initial_model_part),
        mr_mpm_model_part(mpm_model_part)
    {

        // Prepare Dimension and Block Size
        unsigned int TBlock = TDim;
        if (IsMixedFormulation) TBlock ++;

        KRATOS_INFO("MPM_Strategy") << "Dimension Size = " << TDim << " and Block Size = " << TBlock << std::endl;

        // Define a standard static strategy to be used in the calculation
        if(SolutionType == "static" || SolutionType == "Static")
        {
            typename TSchemeType::Pointer pscheme = typename TSchemeType::Pointer( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace >() );

            typename TBuilderAndSolverType::Pointer pBuilderAndSolver;
            if(BlockBuilder == true){
                KRATOS_INFO("MPM_Strategy") << "Block Builder is used" << std::endl;
                pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            }
            else{
                KRATOS_INFO("MPM_Strategy") << "Block Builder is not used" << std::endl;
                pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            }

            const double ratio_tolerance = 1e-04;
            const double always_converged_norm = 1e-09;
            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = typename TConvergenceCriteriaType::Pointer(new ResidualCriteria< TSparseSpace, TDenseSpace >(ratio_tolerance,always_converged_norm));

            bool reform_DOF_at_each_iteration = false;

            mp_solving_strategy = typename SolvingStrategyType::Pointer( new ResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(mr_mpm_model_part,pscheme,plinear_solver,pConvergenceCriteria,pBuilderAndSolver,MaxIteration,ComputeReaction,reform_DOF_at_each_iteration,MoveMeshFlag) );
        }

        // Define a quasi-static strategy to be used in the calculation
        else if(SolutionType == "quasi_static" || SolutionType == "Quasi-static")
        {
            double alpha_M;
            double dynamic;
            typename TSchemeType::Pointer pscheme = typename TSchemeType::Pointer( new MPMResidualBasedBossakScheme< TSparseSpace,TDenseSpace >(mr_grid_model_part, TDim, TBlock, alpha_M = 0.00, dynamic=0) );

            typename TBuilderAndSolverType::Pointer pBuilderAndSolver;
            if(BlockBuilder == true){
                KRATOS_INFO("MPM_Strategy") << "Block Builder is used" << std::endl;
                pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            }
            else{
                KRATOS_INFO("MPM_Strategy") << "Block Builder is not used" << std::endl;
                pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            }

            const double ratio_tolerance = 0.0001;
            const double always_converged_norm = 1e-09;
            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = typename TConvergenceCriteriaType::Pointer(new ResidualCriteria< TSparseSpace, TDenseSpace >(ratio_tolerance,always_converged_norm));

            bool reform_DOF_at_each_iteration = false;

            mp_solving_strategy = typename SolvingStrategyType::Pointer( new ResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(mr_mpm_model_part,pscheme,plinear_solver,pConvergenceCriteria,pBuilderAndSolver,MaxIteration,ComputeReaction,reform_DOF_at_each_iteration,MoveMeshFlag) );
        }

        // Define a dynamic strategy to be used in the calculation
        else if(SolutionType == "dynamic" || SolutionType == "Dynamic")
        {
            double alpha_M;
            double dynamic;
            typename TSchemeType::Pointer pscheme = typename TSchemeType::Pointer( new MPMResidualBasedBossakScheme< TSparseSpace,TDenseSpace >(mr_grid_model_part, TDim, TBlock, alpha_M = 0.0, dynamic=1) );

            typename TBuilderAndSolverType::Pointer pBuilderAndSolver;
            if(BlockBuilder == true){
                KRATOS_INFO("MPM_Strategy") << "Block Builder is used" << std::endl;
                pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            }
            else{
                KRATOS_INFO("MPM_Strategy") << "Block Builder is not used" << std::endl;
                pBuilderAndSolver = typename TBuilderAndSolverType::Pointer(new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(plinear_solver) );
            }

            const double ratio_tolerance = 0.00005;
            const double always_converged_norm = 1e-09;

            typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = typename TConvergenceCriteriaType::Pointer(new ResidualCriteria< TSparseSpace, TDenseSpace >(ratio_tolerance,always_converged_norm));
            bool reform_DOF_at_each_iteration = false;

            mp_solving_strategy = typename SolvingStrategyType::Pointer( new ResidualBasedNewtonRaphsonStrategy<TSparseSpace,TDenseSpace,TLinearSolver >(mr_mpm_model_part,pscheme,plinear_solver,pConvergenceCriteria,pBuilderAndSolver,MaxIteration,ComputeReaction,reform_DOF_at_each_iteration,MoveMeshFlag) );
        }

    }

    /*@} */

    /** Destructor.
     */

    /*@{ */
    virtual ~MPMStrategy()
    {
    }
    /*@} */


    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT:*/

    /**
     * @brief This sets the level of echo for the solution strategy
     * @param Level of echo for the solution strategy
     * @details
     * {
     * 0->Mute... no echo at all
     * 1->Printing time and basic informations
     * 2->Printing linear solver data
     * 3->Print of debug informations: Echo of stiffness matrix, Dx, b...
     * }
     */
    void SetEchoLevel(const int Level) override
    {
        BaseType::mEchoLevel = Level;
        mp_solving_strategy->SetEchoLevel(Level);
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        // Initialize solving strategy: only to be done at the beginning of time step
        mp_solving_strategy->Initialize();
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     */
    void InitializeSolutionStep() override
    {
        // The nodal initial conditions are computed
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - InitializeSolutionStep" <<std::endl;
        mp_solving_strategy->InitializeSolutionStep();
    }

    /**
     * @brief Operation to predict the solution
     */
    void Predict() override
    {
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - Predict" <<std::endl;
        mp_solving_strategy->Predict();
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        // Do solution iterations
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - SolveSolutionStep" <<std::endl;
        return mp_solving_strategy->SolveSolutionStep();
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     */
    void FinalizeSolutionStep() override
    {
        // The nodal solution are mapped from mesh to MP
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - FinalizeSolutionStep" <<std::endl;
        mp_solving_strategy->FinalizeSolutionStep();
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_INFO_IF("MPM_Strategy", this->GetEchoLevel() > 1) << "Main Solve - Clear" <<std::endl;
        mp_solving_strategy->Clear();
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        const auto& mpm_process_info  = mr_mpm_model_part.GetProcessInfo();
        for (auto& r_element : mr_mpm_model_part.Elements())
            r_element.Check(mpm_process_info);

        const auto& grid_process_info = mr_grid_model_part.GetProcessInfo();
        for (auto& r_condition : mr_grid_model_part.Conditions())
            r_condition.Check(grid_process_info);

        return 0;
        KRATOS_CATCH("")
    }

protected:
    /**@name Protected static Member Variables */
    /*@{ */

    ModelPart& mr_grid_model_part;
    ModelPart& mr_initial_model_part;
    ModelPart& mr_mpm_model_part;

    SolvingStrategyType::Pointer mp_solving_strategy;

    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */


private:

    /*@} */
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /** Copy constructor.
     */

    /*@} */

}; /* Class NewSolvingStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_MPM_STRATEGY  defined */


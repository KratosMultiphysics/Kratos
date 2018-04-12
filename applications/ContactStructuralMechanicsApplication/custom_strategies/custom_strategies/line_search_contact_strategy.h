// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_LINE_SEARCH_CONTACT_STRATEGY)
#define KRATOS_LINE_SEARCH_CONTACT_STRATEGY

/* System Includes */

/* External Includes */

/* Project includes */
#include "includes/kratos_parameters.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/line_search_strategy.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

// Default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

// TODO: Extend the descriptions

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
    
/** \brief  Short class definition.
This class 
*/

template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
         
class LineSearchContactStrategy :
    public LineSearchStrategy< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    
    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( LineSearchContactStrategy );

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>              StrategyBaseType;
    
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> NRBaseType;
    
    typedef LineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver>                   BaseType;
    
    typedef typename BaseType::TBuilderAndSolverType                          TBuilderAndSolverType;
 
    typedef typename BaseType::TDataType                                                  TDataType;

    typedef TSparseSpace                                                            SparseSpaceType;

    typedef typename BaseType::TSchemeType                                              TSchemeType;

    typedef typename BaseType::DofsArrayType                                          DofsArrayType;

    typedef typename BaseType::TSystemMatrixType                                  TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                                  TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType                          LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType                          LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType                    TSystemMatrixPointerType;
    
    typedef typename BaseType::TSystemVectorPointerType                    TSystemVectorPointerType;
    
    typedef ModelPart::NodesContainerType                                            NodesArrayType;
    
    typedef ModelPart::ConditionsContainerType                                  ConditionsArrayType;
    
    typedef std::size_t                                                                   IndexType;
    
    /**
     * Default constructor 
     * @param rModelPart: The model part of the problem
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pNewConvergenceCriteria: The convergence criteria employed
     * @param MaxIterationNumber: The maximum number of iterations
     * @param CalculateReactions: The flag for the reaction calculation
     * @param ReformDofSetAtEachStep: The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag: The flag that allows to move the mesh
     */
    
    LineSearchContactStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        IndexType MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        Parameters ThisParameters =  Parameters(R"({})")
    )
        : LineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        KRATOS_TRY;

        Parameters DefaultParameters = Parameters(R"(
        {
        })" );

        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);

        KRATOS_CATCH("");
    }

    /**
     * Default constructor 
     * @param rModelPart: The model part of the problem
     * @param pScheme: The integration scheme
     * @param pNewLinearSolver: The linear solver employed
     * @param pNewConvergenceCriteria: The convergence criteria employed
     * @param MaxIterationNumber: The maximum number of iterations
     * @param CalculateReactions: The flag for the reaction calculation
     * @param ReformDofSetAtEachStep: The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag: The flag that allows to move the mesh
     */
    
    LineSearchContactStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        IndexType MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        Parameters ThisParameters =  Parameters(R"({})")
    )
        : LineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag )
    {
        KRATOS_TRY;

        Parameters DefaultParameters = Parameters(R"(
        {
        })" );

        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);

        KRATOS_CATCH("");
    }

    /** 
     * Destructor.
     */
    
    ~LineSearchContactStrategy() override
    = default;
        
    ///@}
    ///@name Access
    ///@{
    
    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    
    bool mRecalculateFactor;         // To check if we recalculate or not the scale factor

    ///@}
    ///@name Protected Operators
    ///@{
    
    /**
     * Performs all the required operations that should be done (for each step) 
     * before solving the solution step.
     * A member variable should be used as a flag to make sure this function is called only once per step.
     */
        
    void InitializeSolutionStep() override
    {
        BaseType::InitializeSolutionStep();
        
        // TODO: Add something if necessary
    }
    
    /**
     * Here the database is updated
     */
     
    void UpdateDatabase(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        const bool MoveMesh
    ) override
    {
        typename TSchemeType::Pointer pScheme = this->GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = this->GetBuilderAndSolver(); // FIXME: Separate in the parts of LM and displacement

        TSystemVectorType aux(b.size()); //TODO: do it by using the space
        TSparseSpace::Assign(aux, 0.5, Dx);

        TSystemVectorType DxDisp(b.size()); 
        TSystemVectorType DxLM(b.size()); 
        ComputeSplitDx(Dx, DxDisp, DxLM);
        
        // Compute residual without update
        TSparseSpace::SetToZero(b);
        pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
        double roDisp;
        double roLM;
        ComputeMixedResidual(b, roDisp, roLM);
        
        // Compute half step residual
        NRBaseType::UpdateDatabase(A,aux,b,MoveMesh);
        TSparseSpace::SetToZero(b);
        pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
        double rhDisp;
        double rhLM;
        ComputeMixedResidual(b, rhDisp, rhLM);

        // Compute full step residual (add another half Dx to the previous half)
        NRBaseType::UpdateDatabase(A,aux,b,MoveMesh);
        TSparseSpace::SetToZero(b);
        pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
        double rfDisp;
        double rfLM;
        ComputeMixedResidual(b, rfDisp, rfLM);

        // We compute the parabola        
        double XminDisp = 1e-3;
        double XmaxDisp = 1.0;
        double XminLM = 1e-3;
        double XmaxLM = 1.0;
        
        ComputeParabola(XminDisp, XmaxDisp, rfDisp, roDisp, rhDisp);
        ComputeParabola(XminLM, XmaxLM, rfLM, roLM, rhLM);
        
        // Perform final update
        TSparseSpace::Assign(aux,-(1.0 - XmaxDisp), DxDisp);
        TSparseSpace::UnaliasedAdd(aux,-(1.0 - XmaxLM), DxLM);
        NRBaseType::UpdateDatabase(A,aux,b,MoveMesh);
    }

    /**
     * This method split the vector of increment of DoF in displacement and LM
     * @param Dx The increment of displacements and LM
     * @param DxDisp The increment of displacements
     * @param DxLM The increment of LM
     */
        
    void ComputeSplitDx(
        TSystemVectorType& Dx,
        TSystemVectorType& DxDisp,
        TSystemVectorType& DxLM
        )
    {        
        // Now we iterate over all the nodes
        NodesArrayType& nodes_array = StrategyBaseType::GetModelPart().Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size()); 
        
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i) 
        {
            auto it_node = nodes_array.begin() + i;
    
            for(auto itDoF = it_node->GetDofs().begin() ; itDoF != it_node->GetDofs().end() ; itDoF++)
            {
                const int j = (itDoF)->EquationId();
                std::size_t CurrVar = (itDoF)->GetVariable().Key();
                
                if ((CurrVar == DISPLACEMENT_X) || (CurrVar == DISPLACEMENT_Y) || (CurrVar == DISPLACEMENT_Z))
                {          
                    DxDisp[j] = Dx[j];
                    DxLM[j] = 0.0;
                }
                else // Corresponding with contact
                {
                    DxDisp[j] = 0.0;
                    DxLM[j] = Dx[j];
                }
            }
        }
    }
    
    /**
     * This method calculates the norm considering one norm for the displacement and other norm for the LM
     * @param b The residual vector
     * @param normDisp normDisp: The norm of the displacement
     * @param normLM The norm of the LM
     */
        
    void ComputeMixedResidual(
        TSystemVectorType& b,
        double& normDisp, 
        double& normLM
        )
    {        
        // Now we iterate over all the nodes
        NodesArrayType& nodes_array = StrategyBaseType::GetModelPart().Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size()); 
        
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i) 
        {
            auto it_node = nodes_array.begin() + i;
    
            for(auto itDoF = it_node->GetDofs().begin() ; itDoF != it_node->GetDofs().end() ; itDoF++)
            {
                const int j = (itDoF)->EquationId();
                std::size_t CurrVar = (itDoF)->GetVariable().Key();
                
                if ((CurrVar == DISPLACEMENT_X) || (CurrVar == DISPLACEMENT_Y) || (CurrVar == DISPLACEMENT_Z))
                {          
                    #pragma omp atomic
                    normDisp += b[j] * b[j];
                }
                else // Corresponding with contact
                {
                    #pragma omp atomic
                    normLM += b[j] * b[j];
                }
            }
        }
        
        normDisp = std::sqrt(normDisp);
        normLM = std::sqrt(normLM);
    }
    
    /**
     * This method computes the parabola necessary for the line search
     * @param Xmax The maximal abscissa
     * @param Xmin The norm of the LM
     * @param rf The residual norm of the full step
     * @param ro The residual norm without step
     * @param rh The residual norm of the half step
     */
        
    void ComputeParabola(
        double& Xmax,
        double& Xmin,
        const double rf,
        const double ro,
        const double rh
        )
    {   
        // Compute optimal (limited to the range 0-1)
        // Parabola is y = a*x^2 + b*x + c -> min/max for
        // x=0   --> r=ro
        // x=1/2 --> r=rh
        // x=1   --> r =
        // c= ro,     b= 4*rh -rf -3*ro,  a= 2*rf - 4*rh + 2*ro
        // max found if a>0 at the position  Xmax = (rf/4 - rh)/(rf - 2*rh);
        
        const double parabole_a = 2 * rf + 2 * ro - 4 * rh;
        const double parabole_b = 4 * rh - rf - 3 * ro;
        
        if( parabole_a > 0.0) //  If parabola has a local minima
        {
            Xmax = -0.5 * parabole_b/parabole_a; // -b / 2a
            if( Xmax > 1.0)
                Xmax = 1.0;
            else if(Xmax < -1.0)
                Xmax = -1.0;
        }
        else // Parabola degenerates to either a line or to have a local max. best solution on either extreme
        {
            if(rf < ro)
                Xmax = 1.0;
            else
                Xmax = Xmin; // Should be zero, but otherwise it will stagnate
        }
    }
    
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
    ///@{

    /** 
     * Copy constructor.
     */
    
    LineSearchContactStrategy(const LineSearchContactStrategy& Other)
    {
    };

private:

    ///@name Static Member Variables
    ///@{
    
    ///@}
    ///@name Member Variables
    ///@{

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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; /* Class LineSearchContactStrategy */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  // namespace Kratos

#endif /* KRATOS_LINE_SEARCH_CONTACT_STRATEGY */

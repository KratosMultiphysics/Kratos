// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_LINE_SEARCH_CONTACT_STRATEGY)
#define KRATOS_LINE_SEARCH_CONTACT_STRATEGY

/* System Includes */

/* External Includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/kratos_parameters.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "custom_utilities/contact_utilities.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/line_search_strategy.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

// Default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

// TODO: Extend the descriptions

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

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>            StrategyBaseType;
    
    typedef LineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver>                 BaseType;
    
    typedef typename BaseType::TBuilderAndSolverType                        TBuilderAndSolverType;

    typedef typename BaseType::TDataType                                                TDataType;

    typedef TSparseSpace                                                          SparseSpaceType;

    typedef typename BaseType::TSchemeType                                            TSchemeType;

    typedef typename BaseType::DofsArrayType                                        DofsArrayType;

    typedef typename BaseType::TSystemMatrixType                                TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                                TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType                        LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType                        LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType                  TSystemMatrixPointerType;
    
    typedef typename BaseType::TSystemVectorPointerType                  TSystemVectorPointerType;
    
    typedef ModelPart::NodesContainerType                                          NodesArrayType;
    
    typedef ModelPart::ConditionsContainerType                                ConditionsArrayType;
    
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
        unsigned int MaxIterations = 30,
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
            "rescale_factor"                   : false
        })" );

        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        
        mRecalculateFactor = ThisParameters["rescale_factor"].GetBool();

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
        unsigned int MaxIterations = 30,
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
            "rescale_factor"                   : false
        })" );

        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        
        mRecalculateFactor = ThisParameters["rescale_factor"].GetBool();

        KRATOS_CATCH("");
    }

    /** 
     * Destructor.
     */
    
    virtual ~LineSearchContactStrategy()
    {
    }
        
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
        
        // Now we rescale the scale factor
        if (mRecalculateFactor == true && StrategyBaseType::GetModelPart().GetProcessInfo()[TIME_STEPS] == 1)
        {
            RescaleFactor();
        }
    }
    
    /**
     * We rescale the scale factor in function of the norm of the RHS
     */
    
    void RescaleFactor()
    {
        // We get the scale factor
        double& ScaleFactor = StrategyBaseType::GetModelPart().GetProcessInfo()[SCALE_FACTOR]; 
        if (ScaleFactor == 0.0)
        {
            KRATOS_ERROR << "You don't have any value assigned to SCALE_FACTOR" << std::endl;
        }
        
        // Pointers needed in the solution
        typename TSchemeType::Pointer pScheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
        
        // We recalculate the RHS
        TSystemVectorType& b = *BaseType::mpb;
        TSparseSpace::SetToZero(b);
        pBuilderAndSolver->BuildRHS(pScheme, StrategyBaseType::GetModelPart(), b);
        
        // We initialize the values of the contact and non contact norm
        double AuxContact = 0.0;
        double AuxNonContact = 0.0;
        
        // Now we iterate over all the nodes
        NodesArrayType& pNode = StrategyBaseType::GetModelPart().GetSubModelPart("Contact").Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
    
            for(auto itDoF = itNode->GetDofs().begin() ; itDoF != itNode->GetDofs().end() ; itDoF++)
            {
                const int j = (itDoF)->EquationId();
                
                if (((itDoF)->GetVariable().Name()).find("VECTOR_LAGRANGE") != std::string::npos || ((itDoF)->GetVariable().Name()).find("NORMAL_CONTACT_STRESS") == std::string::npos) // Corresponding with contact
                {          
                    #pragma omp atomic
                    AuxContact += b[j] * b[j];
                }
                else 
                {
                    #pragma omp atomic
                    AuxNonContact += b[j] * b[j];
                }
            }
        }
        
        ScaleFactor *= std::sqrt(AuxContact/AuxNonContact); // NOTE: The inverse¿?
//         ScaleFactor *= std::sqrt(AuxNonContact/AuxContact);
        
        if (StrategyBaseType::mEchoLevel > 0)
        {
            std::cout << "The new SCALE_FACTOR is: " << ScaleFactor << std::endl;
        }
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
        BaseType::UpdateDatabase(A,Dx,b,MoveMesh);
        
        CalculateContactReactions(b);
    }
    
    /**
     * We activate all yhe nodes in the active conditions
     */
    
    void ReactivateNodes()
    {
        // Now we iterate over all the conditions and we set all the nodes as ACTIVE in the active conditions 
        ConditionsArrayType& pCond = StrategyBaseType::GetModelPart().GetSubModelPart("Contact").Conditions();
        auto numConditions = pCond.end() - pCond.begin();
        
        #pragma omp parallel for
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pCond.begin() + i;
            
            if (itCond->Is(ACTIVE) == true)
            {
                for (unsigned int iNode = 0; iNode < itCond->GetGeometry().size(); iNode++)
                {
                    itCond->GetGeometry()[iNode].SetLock();
                    itCond->GetGeometry()[iNode].Set(ACTIVE, true);
                    itCond->GetGeometry()[iNode].UnSetLock();
                }
            }
        }
    }
    
    /**
     * This method calculates the reactions concerning the contact (residual of the contact)
     */
    
    void CalculateContactReactions(TSystemVectorType& b)
    {       
        // First we activate all the nodes in the active conditions
        ReactivateNodes();
        
        // Pointers needed in the solution
        typename TSchemeType::Pointer pScheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = BaseType::GetBuilderAndSolver();
        
        const double ScaleFactor = (StrategyBaseType::GetModelPart().GetProcessInfo()[SCALE_FACTOR] > 0.0) ? StrategyBaseType::GetModelPart().GetProcessInfo()[SCALE_FACTOR] : 1.0;
        
        // We recalculate the RHS (NOTE: EXPENSIVE)
        TSparseSpace::SetToZero(b);
        pBuilderAndSolver->BuildRHS(pScheme, StrategyBaseType::GetModelPart(), b);
        
        // Now we iterate over all the nodes
        NodesArrayType& pNode = StrategyBaseType::GetModelPart().GetSubModelPart("Contact").Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            if (itNode->Is(ACTIVE) == true)
            {                
                for(auto itDoF = itNode->GetDofs().begin() ; itDoF != itNode->GetDofs().end() ; itDoF++)
                {
                    const int j = (itDoF)->EquationId();
                    
                    if (((itDoF)->GetReaction().Name()).find("WEIGHTED") != std::string::npos) // Corresponding with contact
                    {                        
                        (itDoF)->GetSolutionStepReactionValue() = b[j]/ScaleFactor;
                    }
                    else if ((itDoF)->GetReaction().Name() != "NONE") // The others
                    {                        
                        (itDoF)->GetSolutionStepReactionValue() = -b[j];
                    }
                }
            }
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

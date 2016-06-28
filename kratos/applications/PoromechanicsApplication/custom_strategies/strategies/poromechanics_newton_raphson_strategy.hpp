//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_POROMECHANICS_NEWTON_RAPHSON_STRATEGY)
#define KRATOS_POROMECHANICS_NEWTON_RAPHSON_STRATEGY

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace,class TDenseSpace,class TLinearSolver>

class PoromechanicsNewtonRaphsonStrategy : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(PoromechanicsNewtonRaphsonStrategy);
    
    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> MotherType;
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
    typedef typename BaseType::TSchemeType TSchemeType;
    using MotherType::mInitializeWasPerformed;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    PoromechanicsNewtonRaphsonStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme, pNewLinearSolver,
                pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
        {
            
        }

    //------------------------------------------------------------------------------------

    ///Destructor
    virtual ~PoromechanicsNewtonRaphsonStrategy() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize()
    {
        KRATOS_TRY

        if (mInitializeWasPerformed == false)
		{
            MotherType::Initialize();
            
            //Initialize ProcessInfo variables
            BaseType::GetModelPart().GetProcessInfo()[IS_CONVERGED] = true;
        }

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check()
    {
        KRATOS_TRY
        
        int ierr = MotherType::Check();
        if(ierr != 0) return ierr;
        
        if(IS_CONVERGED.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"IS_CONVERGED has Key zero! (check if the application is correctly registered", "" )
        
        return ierr;

        KRATOS_CATCH( "" )
    }

}; // Class PoromechanicsNewtonRaphsonStrategy

} // namespace Kratos

#endif // KRATOS_POROMECHANICS_NEWTON_RAPHSON_STRATEGY  defined

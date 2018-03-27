//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_COMPONENT_WISE_NEWTON_RAPHSON_STRATEGY_H_INCLUDED )
#define  KRATOS_COMPONENT_WISE_NEWTON_RAPHSON_STRATEGY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_strategies/strategies/newton_raphson_strategy.hpp"

//default builder and solver
#include "custom_strategies/builders_and_solvers/component_wise_builder_and_solver.hpp"

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
 * @class ComponentWise NewtonRaphsonStrategy
 * @brief This is the base Newton Raphson strategy
 * @details This strategy iterates until the convergence is achieved (or the maximum number of iterations is surpassed) using a Newton Raphson algorithm
 */

template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ComponentWiseNewtonRaphsonStrategy : public NewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
 public:
  ///@name Type Definitions
  ///@{
  
  // Counted pointer of ClassName
  KRATOS_CLASS_POINTER_DEFINITION( ComponentWiseNewtonRaphsonStrategy );

  typedef NewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>       BaseType;

  typedef ConvergenceCriteria<TSparseSpace, TDenseSpace>        ConvergenceCriterionType;

  typedef typename BaseType::BuilderAndSolverType                   BuilderAndSolverType;

  typedef typename BaseType::SchemeType                                       SchemeType;

  typedef TLinearSolver                                                 LinearSolverType;
  
  ///@}
  ///@name Life Cycle

  ///@{

  /// Constructor.
  ComponentWiseNewtonRaphsonStrategy(ModelPart& rModelPart,
                                     typename SchemeType::Pointer pScheme,
                                     typename BuilderAndSolverType::Pointer pBuilderAndSolver,
                                     typename ConvergenceCriterionType::Pointer pConvergenceCriterion,
                                     Flags& rOptions,
                                     unsigned int MaxIterations = 30)
      : NewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pBuilderAndSolver, pConvergenceCriterion, rOptions, MaxIterations)
  {
    
    KRATOS_TRY
            
    //component-wise options
    typename BuilderAndSolverType::GlobalSystemComponents& rGlobalSystem = this->mpBuilderAndSolver->GetGlobalSystemComponents();
	
    rGlobalSystem.SetRHS_Element_Components( this->mpConvergenceCriteria->GetRHS_Element_Components() );
    rGlobalSystem.SetRHS_Element_Variables( this->mpConvergenceCriteria->GetRHS_Element_Variables() );

    rGlobalSystem.SetRHS_Condition_Components( this->mpConvergenceCriteria->GetRHS_Condition_Components() );
    rGlobalSystem.SetRHS_Condition_Variables( this->mpConvergenceCriteria->GetRHS_Condition_Variables() );
    //component-wise options

    KRATOS_CATCH( "" )
  }


  /// Constructor.
  ComponentWiseNewtonRaphsonStrategy(ModelPart& rModelPart,
                                     typename SchemeType::Pointer pScheme,
                                     typename LinearSolverType::Pointer pLinearSolver,
                                     typename ConvergenceCriterionType::Pointer pConvergenceCriterion,
                                     Flags& rOptions,
                                     unsigned int MaxIterations = 30)
      : ComponentWiseNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, typename BuilderAndSolverType::Pointer(new ComponentWiseBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pLinearSolver)), pConvergenceCriterion, rOptions, MaxIterations)
  {
    KRATOS_TRY
        
    KRATOS_CATCH( "" )
  }



  /// Destructor.
  virtual ~ComponentWiseNewtonRaphsonStrategy()
  {
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
  ///@name Inquiry
  ///@{
  ///@}
  ///@name Friends
  ///@{
  ///@}

 protected:
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
  ///@name Private Inquiry
  ///@{
  ///@}
  ///@name Un accessible methods
  ///@{

  /// Copy constructor.
  ComponentWiseNewtonRaphsonStrategy(const ComponentWiseNewtonRaphsonStrategy& Other){};

  ///@}

}; /// Class ComponentWisNewtonRaphsonStrategy

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_COMPONENT_WISE_NEWTON_RAPHSON_STRATEGY_H_INCLUDED  defined


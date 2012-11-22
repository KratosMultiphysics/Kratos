//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 2008-10-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//

#if !defined(KRATOS_MPI_EXPLICIT_SOLVER_STRATEGY)
#define  KRATOS_MPI_EXPLICIT_SOLVER_STRATEGY

// /* External includes */
// #include "boost/smart_ptr.hpp"

// System includes

// Project includes
#include "utilities/timer.h"
#include "custom_utilities/neighbours_calculator.h"
#include "custom_utilities/mpi_neighbours_calculator.h"
#include "custom_utilities/create_and_destroy.h"

#include "custom_elements/spheric_particle.h" //M: le afegit jo.. no hi era. cal que hi sigui oi???
#include "includes/variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
#include <iostream>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/schemes/scheme.h"

// #include "custom_utilities/neighbours_calculator.h"
#include "custom_utilities/mpi_neighbours_calculator.h"
#include "custom_strategies/schemes/integration_scheme.h"
#include "custom_strategies/strategies/explicit_solver_strategy.h"

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

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
  
  /// Short class definition.
  /** Detail class definition.
  */
  template<
  class TSparseSpace,
  class TDenseSpace, 
  class TLinearSolver> 
  class MpiExplicitSolverStrategy : public ExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
     {
      public:
      ///@name Type Definitions
      ///@{
 
      typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>   BaseType;
      typedef typename BaseType::TDataType                              TDataType;
      typedef typename BaseType::TBuilderAndSolverType                  TBuilderAndSolverType;
      typedef typename BaseType::TSchemeType                            TSchemeType;
      typedef typename BaseType::DofsArrayType                          DofsArrayType;
      typedef typename Element::DofsVectorType                          DofsVectorType;
      
      typedef ModelPart::NodesContainerType                             NodesArrayType;
      typedef ModelPart::ElementsContainerType                          ElementsArrayType;
      typedef ModelPart::ConditionsContainerType                        ConditionsArrayType;     
      
      typedef ModelPart::NodesContainerType::ContainerType              NodesContainerType;
      typedef ModelPart::ElementsContainerType::ContainerType           ElementsContainerType;
      typedef ModelPart::ConditionsContainerType::ContainerType         ConditionsContainerType;
      
      typedef DiscreteElement                                           ParticleType;
      
      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(MpiExplicitSolverStrategy);
 
      
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      MpiExplicitSolverStrategy(){}
      
      MpiExplicitSolverStrategy(ModelPart& model_part,
                             ModelPart& contacts_model_part, 
                             const int dimension,
                             const double enlargement_factor,
                             const double damping_ratio,
                             const double fraction_delta_time,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const bool MoveMeshFlag,
                             const bool delta_option,
                             const bool continuum_simulating_option,
                             typename IntegrationScheme::Pointer pScheme
      ) : ExplicitSolverStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(
                             model_part,
                             contacts_model_part,
                             dimension,
                             enlargement_factor,
                             damping_ratio,
                             fraction_delta_time,
                             max_delta_time,
                             n_step_search,
                             safety_factor,
                             MoveMeshFlag,
                             delta_option,
                             continuum_simulating_option,
                             pScheme )
      {
      }

      /// Destructor.
      virtual ~MpiExplicitSolverStrategy(){}
          
    protected:

    private:
 
    typename IntegrationScheme::Pointer mpScheme;
    
    virtual void Synchronize(ModelPart& r_model_part)
    {
        r_model_part.GetCommunicator().SynchronizeNodalSolutionStepsData();
        r_model_part.GetCommunicator().SynchronizeDofs();
    }
    
    virtual void Repart(ModelPart& r_model_part)
    {
        typedef Mpi_Neighbours_Calculator<ParticleType> NeighboursCalculatorType;
        
        NeighboursCalculatorType::Parallel_partitioning(r_model_part,true);
    }
    
    virtual ElementsArrayType& GetElements(ModelPart& r_model_part)
    {
        return r_model_part.GetCommunicator().LocalMesh().Elements();
    }

    virtual void SearchIniNeighbours(ModelPart& r_model_part,bool extension_option)
    { 
        typedef Mpi_Neighbours_Calculator<ParticleType> NeighboursCalculatorType;
              
        NeighboursCalculatorType neighbourCalc;
        neighbourCalc.Search_Ini_Neighbours(r_model_part, extension_option);
    }//SearchIniNeighbours


    virtual void SearchNeighbours(ModelPart& r_model_part,bool extension_option)
    {
        typedef Mpi_Neighbours_Calculator<ParticleType> NeighboursCalculatorType;
              
        NeighboursCalculatorType neighbourCalc;
        neighbourCalc.Search_Neighbours(r_model_part, extension_option);
    }//SearchNeighbours

  
  }; // Class MpiExplicitSolverStrategy  


        
 /*
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
                    MpiExplicitSolverStrategy& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
                    const MpiExplicitSolverStrategy& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
    */
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 





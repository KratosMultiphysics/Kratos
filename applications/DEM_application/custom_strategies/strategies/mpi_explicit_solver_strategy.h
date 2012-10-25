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

#include "custom_utilities/neighbours_calculator.h"
#include "custom_strategies/schemes/integration_scheme.h"

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
  class MpiExplicitSolverStrategy : MpiExplicitSolverStrategy
     {
      public:
      ///@name Type Definitions
      ///@{
 
      typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;
      typedef typename BaseType::TDataType TDataType;
      typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
      typedef typename BaseType::TSchemeType TSchemeType;
      typedef typename BaseType::DofsArrayType DofsArrayType;
      typedef typename Element::DofsVectorType DofsVectorType;
      
      typedef ModelPart::NodesContainerType NodesArrayType;
      typedef ModelPart::ElementsContainerType ElementsArrayType;
      typedef ModelPart::ConditionsContainerType ConditionsArrayType;     
      
      typedef ModelPart::NodesContainerType::ContainerType      NodesContainerType;
      typedef ModelPart::ElementsContainerType::ContainerType   ElementsContainerType;
      typedef ModelPart::ConditionsContainerType::ContainerType ConditionsContainerType;
      
      typedef DiscreteElement                         ParticleType;
      typedef Mpi_Neighbours_Calculator<ParticleType> NeighboursCalculatorType;
      
      /// Pointer definition of ExplicitSolverStrategy
      KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolverStrategy);
 
      
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      ExplicitSolverStrategy(){}
      
      ExplicitSolverStrategy(ModelPart& model_part,   
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
      ) : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag),mdimension(dimension) //inicialitzacio de variables const. no poden inicialitzarse a l'esquerra d'un igual.
      {
          mdelta_option                = delta_option;
          mcontinuum_simulating_option = continuum_simulating_option;
          mvirtual_mass                = false;  //M: it has to be implemented.
          mElementsAreInitialized      = false;
          mConditionsAreInitialized    = false;
          mCalculateOldTime            = false;
          mSolutionStepIsInitialized   = false;
          mInitializeWasPerformed      = false;
          mComputeTime                 = false;
          mInitialConditions           = false;
          mEnlargementFactor           = enlargement_factor;
          mdamping_ratio               = damping_ratio; //not in use
          mfraction_delta_time         = fraction_delta_time;
          mmax_delta_time              = max_delta_time;
          molddelta_time               = 0.00;
          mtimestep                    = 0.00;
          mpScheme                     = pScheme;
          
          mtimestep                    = max_delta_time;
          mnstepsearch                 = n_step_search;
          msafety_factor               = safety_factor;
      }

      /// Destructor.
      virtual ~ExplicitSolverStrategy(){}
      
      void Initialized()
      {
          KRATOS_TRY

          //M: faig una primera búsqueda abans de inicialitzar elements pk allí guardaré veins inicials i altres coses.
          ModelPart& r_model_part           = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo  = r_model_part.GetProcessInfo();

          //1. Search Neighbours with tolerance (Not in mpi.)
          bool extension_option = false;

          if(rCurrentProcessInfo[CASE_OPTION] != 0)
          {
              extension_option = true;
          }

          //2. Initializing elements and perform the 1st repartition
          if(mElementsAreInitialized == false)
          {
//               #ifdef _OPENMPI
//               Neighbours_Calculator<2, DiscreteElement>::parallel_partitioning(r_model_part,true);
//               for(int i = 0; i < 30; i++)
//                   Neighbours_Calculator<2, DiscreteElement>::parallel_partitioning(r_model_part,true);
//               #endif
              InitializeElements();
          }
          
          //3. Search Neighbours with tolerance (afther first repartition process)
          SearchIniNeighbours(r_model_part,extension_option);  
          mInitializeWasPerformed   = true;

          // 3. Set Initial Contacts
          if(mdelta_option || mcontinuum_simulating_option)
          {
              Set_Initial_Contacts(mdelta_option, mcontinuum_simulating_option);  //delta option no fa falta i fer el continuu
          }

          //4.Final operations
          FinalizeSolutionStep();

          KRATOS_CATCH("")
      }
       
      double Solve()
      {
          KRATOS_TRY

          std::cout<<std::fixed<<std::setw(15)<<std::scientific<<std::setprecision(5);
          
          ModelPart& r_model_part          = BaseType::GetModelPart();
          ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();
          
          bool extension_option = false;

          if(rCurrentProcessInfo[CASE_OPTION] )
          {
              extension_option = true;
          }

          int time_step = rCurrentProcessInfo[TIME_STEPS];

//        std::cout<<"------------------------------------------------------------------------"<<std::endl;
//        std::cout<<"                 KRATOS DEM APPLICATION. TIME STEPS = " <<  time_step    <<std::endl;
//        std::cout<<"------------------------------------------------------------------------"<<std::endl;

          //STRATEGY:
          //0.0
          Synchronize(r_model_part);
          
          //1.0
          InitializeSolutionStep();

          //1. Get and Calculate the forces
          GetForce();
  
          //1.1. Calculate Local Dampings
          int rota_damp_id            = rCurrentProcessInfo[ROTA_DAMP_TYPE];
          int rotation_OPTION         = rCurrentProcessInfo[ROTATION_OPTION];
      
          if ( (rotation_OPTION == 1) && (rota_damp_id == 1) ) 
          {
              ApplyRotationalDampings();
          }
          
          //2. Motion Integration
          ComputeIntermedialVelocityAndNewDisplacement(); //llama al scheme, i aquesta ja fa el calcul dels despaçaments i tot
          
          //3. Neighbouring search. Every N times. +bounding box destruction
          if( time_step == 1)
          {
              mParticle_Creator_Destructor.CalculateSurroundingBoundingBox(r_model_part, mEnlargementFactor);
          }
          
          if ( (time_step + 1)%mnstepsearch == 0 )
          {
              if ( (time_step + 1)%(mnstepsearch*10) == 0 )
              {
                  Repart()
              }
              if(rCurrentProcessInfo[BOUNDING_BOX_OPTION]==1)
              {
                  BoundingBoxUtility(mEnlargementFactor);
              }
              
              SearchNeighbours(r_model_part,extension_option); //extension option false;
          }
          
          //4.Final operations
          FinalizeSolutionStep();
          
          return 0.00;
          KRATOS_CATCH("")
      }
          
    protected:

    private:
 
    typename IntegrationScheme::Pointer mpScheme;
    
    virtual void Synchronize(ModelPart& r_model_part)
    {
        r_model_part.GetCommunicator().SynchronizeNodalSolutionStepsData();
        r_model_part.GetCommunicator().SynchronizeDofs();
    }
    
    virtual void Repart()
    {
        NeighboursCalculatorType::Parallel_partitioning(r_model_part,true);
        InitializeElements(); //TODO: REMOVE
    }
    
    virtual ElementsArrayType& GetElements(ModelPart& r_model_part)
    {
        return r_model_part.GetCommunicator().LocalMesh().Elements();
    }

    virtual void SearchIniNeighbours(ModelPart r_model_part,bool extension_option)
    {        
        NeighboursCalculatorType neighbourCalc;
        neighbourCalc.Search_Ini_Neighbours(r_model_part, extension_option);

    }//SearchIniNeighbours


    virtual void SearchNeighbours(ModelPart r_model_part,bool extension_option)
    {
        NeighboursCalculatorType neighbourCalc;
        neighbourCalc.Search_Neighbours(r_model_part, extension_option);

    }//SearchNeighbours

  
  }; // Class ExplicitSolverStrategy  


        
 /*
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
                    ExplicitSolverStrategy& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
                    const ExplicitSolverStrategy& rThis)
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





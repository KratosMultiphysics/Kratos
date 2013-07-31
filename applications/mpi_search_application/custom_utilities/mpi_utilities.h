//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Carlos A. Roig $
//   Date:                $Date: 2008-10-13 08:56:42 $
//   Revision:            $Revision: 1 $
//
//

#if !defined(KRATOS_MPI_UTILITIES)
#define  KRATOS_MPI_UTILITIES


// System includes

// Project includes
#include "utilities/timer.h"

/* System includes */
#include <limits>
#include <iomanip>
#include <iostream>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/line_3d_2.h"

/* Search */
#include "custom_utilities/lloyd_parallel_partitioner.h"

#define CUSTOMTIMER 1

/* Timer defines */
#ifdef CUSTOMTIMER
#define KRATOS_TIMER_START(t) Timer::Start(t);
#define KRATOS_TIMER_STOP(t) Timer::Stop(t);
#else
#define KRATOS_TIMER_START(t)
#define KRATOS_TIMER_STOP(t)
#endif

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
  
  class MpiUtilities
  {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of MpiUtilities
      KRATOS_CLASS_POINTER_DEFINITION(MpiUtilities);
      
      typedef SpatialSearch                                         SearchType;
      typedef MpiDiscreteParticleConfigure<3>                       Configure;

      typedef SearchType::ElementsContainerType::ContainerType      ElementsContainerType;
      typedef SearchType::NodesContainerType::ContainerType         NodesContainerType;
      typedef SearchType::ConditionsContainerType::ContainerType    ConditionsContainerType;
      
//       typedef ElementsContainerType::value_type                     ElementPointerType;
//       typedef SearchType::IteratorType                              IteratorType;

      
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      MpiUtilities(){}
      
      MpiUtilities(ModelPart& model_part,
                   ModelPart& contacts_model_part
      ) 
      {
      }

      /// Destructor.
      virtual ~MpiUtilities(){}
      
      void TransferModelElements(ModelPart& rModelPart)
      {
          KRATOS_TRY
          
          rModelPart.GetCommunicator().TransferModelElements(rModelPart);
          
          KRATOS_CATCH("")
      }
      
      void TransferModelNodes(ModelPart& rModelPart)
      {
          KRATOS_TRY
          
          rModelPart.GetCommunicator().TransferModelNodes(rModelPart);
          
          KRATOS_CATCH("")
      }
    
      void ParallelPartitioning(ModelPart& rModelPart, bool extension_option, int CalculateBoundry)
      {
          KRATOS_TRY
          
          KRATOS_TIMER_START("PART")
                    
          ElementsContainerType& pLocalElements = rModelPart.GetCommunicator().LocalMesh().ElementsArray();

          ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
          
          double radius_extend = 0.0;
          if (extension_option) radius_extend = rCurrentProcessInfo[SEARCH_RADIUS_EXTENSION];
          
          static double MaxNodeRadius = 0.0f;
          if(MaxNodeRadius == 0.0f) //TODO
              for (ElementsContainerType::iterator particle_pointer_it = pLocalElements.begin(); particle_pointer_it != pLocalElements.end(); ++particle_pointer_it)
              {
                  double NodeRaidus = (1.0 + radius_extend) * (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                  MaxNodeRadius = NodeRaidus > MaxNodeRadius ? NodeRaidus : MaxNodeRadius;
              }
          
          LloydParallelPartitioner<Configure> partitioner;
          partitioner.LloydsBasedParitioner(rModelPart,MaxNodeRadius,CalculateBoundry);
          
          KRATOS_TIMER_STOP("PART")
          
          Timer::SetOuputFile("TimesPartitioner");
          Timer::PrintTimingInformation();
          
          KRATOS_CATCH("")
      }
      
      void CalculateModelNewIds(ModelPart& mModelPart, int offset)
      {
          CalculateElementsNewId(mModelPart,offset);
          CalculateNodesNewId(mModelPart,offset);
          CalculateConditionsNewId(mModelPart,offset);
      }
    
      void CalculateElementsNewId(ModelPart& mModelPart, int offset)
      {
          KRATOS_TRY
        
          int element_size = mModelPart.GetCommunicator().LocalMesh().Elements().size();
          int iteratorId = -1;
        
          Configure::ReduceIds(element_size,offset,iteratorId);
          
          if(iteratorId == -1)
              std::cout << "Invalid starting Id" << std::endl;
          
          for (ElementsContainerType::iterator it = mModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin(); it != mModelPart.GetCommunicator().LocalMesh().Elements().ptr_end(); ++it)
              (*it)->SetId(iteratorId++);
          
          KRATOS_CATCH("")
      }
      
      void CalculateNodesNewId(ModelPart& mModelPart, int offset)
      {
          KRATOS_TRY
        
          int nodes_size = mModelPart.GetCommunicator().LocalMesh().Nodes().size();
          int iteratorId = -1;
          
          Configure::ReduceIds(nodes_size,offset,iteratorId);
          
          if(iteratorId == -1)
              std::cout << "Invalid starting Id" << std::endl;
          
          std::cout << "STARTING NODE ID: " << iteratorId << std::endl;
          std::cout << "ENDING   NODE ID: " << iteratorId + nodes_size << std::endl;

          for (NodesContainerType::iterator it = mModelPart.GetCommunicator().LocalMesh().Nodes().ptr_begin(); it != mModelPart.GetCommunicator().LocalMesh().Nodes().ptr_end(); ++it)
          {
//             std::cout << (*it)->Id() << "\t-->\t" << iteratorId << std::endl;
              (*it)->SetId(iteratorId);
              iteratorId++;
          }
          
          std::cout << "New ids assigned" << std::endl;
          
          KRATOS_CATCH("")
      }
      
      void CalculateConditionsNewId(ModelPart& mModelPart, int offset)
      {
          KRATOS_TRY
        
          int conditions_size = mModelPart.GetCommunicator().LocalMesh().Conditions().size();
          int iteratorId = -1;
        
          Configure::ReduceIds(conditions_size,offset,iteratorId);
          
          if(iteratorId == -1)
              std::cout << "Invalid starting Id" << std::endl;
          
          for (ConditionsContainerType::iterator it = mModelPart.GetCommunicator().LocalMesh().Conditions().ptr_begin(); it != mModelPart.GetCommunicator().LocalMesh().Conditions().ptr_end(); ++it)
              (*it)->SetId(iteratorId++);
          
          KRATOS_CATCH("")
      }
          
    protected:

    private:
  
  }; // Class MpiUtilities  


        
 /*
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
                    MpiUtilities& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
                    const MpiUtilities& rThis)
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





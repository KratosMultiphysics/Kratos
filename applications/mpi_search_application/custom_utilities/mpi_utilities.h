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
      
      typedef MpiDiscreteParticleConfigure < 3 >                Configure;

      typedef Configure::ElementsContainerType::ContainerType   ElementsContainerType;
      typedef Configure::PointerType                            PointerType;
      typedef Configure::IteratorType                           IteratorType;
      typedef Configure::ResultIteratorType                     ResultIteratorType;

      
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
      
      void TransferObjects(ModelPart& rModelPart)
      {
          KRATOS_TRY
          
          rModelPart.GetCommunicator().TransferObjects(rModelPart);
          
          KRATOS_CATCH("")
      }
    
      void ParallelPartitioning(ModelPart& rModelPart, bool extension_option, int CalculateBoundry)
      {
          KRATOS_TRY
                    
          ElementsContainerType& pLocalElements = rModelPart.GetCommunicator().LocalMesh().ElementsArray();

          ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
          
          double radius_extend = 0.0;
          if (extension_option) radius_extend = rCurrentProcessInfo[SEARCH_RADIUS_EXTENSION];
          
          static double MaxNodeRadius = 0.0f;
          if(MaxNodeRadius == 0.0f) //TODO
              for (IteratorType particle_pointer_it = pLocalElements.begin(); particle_pointer_it != pLocalElements.end(); ++particle_pointer_it)
              {
                  double NodeRaidus = (1.0 + radius_extend) * (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                  MaxNodeRadius = NodeRaidus > MaxNodeRadius ? NodeRaidus : MaxNodeRadius;
              }
          
          LloydParallelPartitioner<Configure> partitioner;
          partitioner.LloydsBasedParitioner(rModelPart,MaxNodeRadius,CalculateBoundry);
          
          KRATOS_CATCH("")
      }
    
      void CompactIds(ModelPart& mContactsModelPart)
      {
          KRATOS_TRY
        
          int contacts_model_part_size = mContactsModelPart.GetCommunicator().LocalMesh().Elements().size();
          int iteratorId = -1;
        
          Configure::ReduceIds(contacts_model_part_size,iteratorId);
          
          if(iteratorId == -1)
              std::cout << "Something went wrong :(" << std::endl;
          
          for (IteratorType it = mContactsModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin(); it != mContactsModelPart.GetCommunicator().LocalMesh().Elements().ptr_end(); ++it)
          {
              (*it)->SetId(iteratorId++);
          }
          
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





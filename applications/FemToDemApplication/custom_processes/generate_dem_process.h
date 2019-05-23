//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_GENERATE_DEM_PROCESS)
#define KRATOS_GENERATE_DEM_PROCESS


#include "includes/model_part.h"
#include "processes/process.h"
#include "fem_to_dem_application_variables.h"

namespace Kratos {

typedef std::size_t SizeType;

class GenerateDemProcess : public Process 
{
 public:

  /// Pointer definition of ApplyMultipointConstraintsProcess
  KRATOS_CLASS_POINTER_DEFINITION(GenerateDemProcess);

  typedef std::size_t IndexType;

  // Constructor
  GenerateDemProcess(ModelPart& rModelPart, ModelPart& rDemModelPart);

  // Destructor
  ~GenerateDemProcess() override = default;

  void operator()() { Execute(); }

  void Execute() override;

protected:

  // Member Variables
  ModelPart& mrModelPart;
  ModelPart& mrDEMModelPart;

};  // Class

}  // namespace Kratos
#endif /* KRATOS_EXTEND_PRESSURE_PROCESS defined */
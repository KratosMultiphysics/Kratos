//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license:
//kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_ASSIGN_PRESSURE_ID_PROCESS)
#define KRATOS_ASSIGN_PRESSURE_ID_PROCESS


#include "includes/model_part.h"
#include "processes/process.h"
#include "fem_to_dem_application_variables.h"

namespace Kratos {

typedef std::size_t SizeType;

class ComputeInitialVolumeProcess : public Process 
{
 public:

  /// Pointer definition of ApplyMultipointConstraintsProcess
  KRATOS_CLASS_POINTER_DEFINITION(ComputeInitialVolumeProcess);

  typedef std::size_t IndexType;

  // Constructor
  ComputeInitialVolumeProcess(ModelPart& rModelPart);

  // Destructor
  ~ComputeInitialVolumeProcess() override = default;

  void operator()() { Execute(); }

  void Execute() override;

  void AssignPressureIdToNodes(std::string rSubModelPartName, const int PressureId);

protected:

  // Member Variables
  ModelPart& mrModelPart;
  int mDimension;

};  // Class

}  // namespace Kratos
#endif /* KRATOS_EXTEND_PRESSURE_PROCESS defined */
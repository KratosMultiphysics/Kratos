//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//    Kratos default license:
//  kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_UPDATE_FLAG_NO_REMESH_FEMDEM_BOUNDARY_PROCESS)
#define KRATOS_UPDATE_FLAG_NO_REMESH_FEMDEM_BOUNDARY_PROCESS


#include "processes/process.h"
#include "fem_to_dem_application_variables.h"

namespace Kratos {

typedef std::size_t SizeType;
typedef Node<3> NodeType;
typedef ModelPart::NodesContainerType::iterator NodeIteratorType;

class UpdateFlagNoRemeshFemDemBoundaryProcess : public Process 
{
 public:

  /// Pointer definition of ApplyMultipointConstraintsProcess
  KRATOS_CLASS_POINTER_DEFINITION(UpdateFlagNoRemeshFemDemBoundaryProcess);

  typedef std::size_t IndexType;

  // Constructor
  UpdateFlagNoRemeshFemDemBoundaryProcess(ModelPart& rModelPart);

  // Destructor
  ~UpdateFlagNoRemeshFemDemBoundaryProcess() override = default;

  void operator()() { Execute(); }

  void Execute() override;

protected:

  // Member Variables
  ModelPart& mrModelPart;

};  // Class

}  // namespace Kratos
#endif /* KRATOS_UPDATE_FLAG_NO_REMESH_FEMDEM_BOUNDARY_PROCESS defined */
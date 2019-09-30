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

#if !defined(KRATOS_UPDATE_DEM_KINEMATICS_PROCESS)
#define KRATOS_UPDATE_DEM_KINEMATICS_PROCESS


#include "includes/model_part.h"
#include "processes/process.h"
#include "fem_to_dem_application_variables.h"

namespace Kratos {

typedef std::size_t SizeType;
typedef Node<3> NodeType;
typedef ModelPart::NodesContainerType::iterator NodeIteratorType;

class UpdateDemKinematicsProcess : public Process 
{
 public:

  /// Pointer definition of ApplyMultipointConstraintsProcess
  KRATOS_CLASS_POINTER_DEFINITION(UpdateDemKinematicsProcess);

  typedef std::size_t IndexType;

  // Constructor
  UpdateDemKinematicsProcess(ModelPart& rModelPart, ModelPart& rDemModelPart);

  // Destructor
  ~UpdateDemKinematicsProcess() override = default;

  void operator()() { Execute(); }

  void Execute() override;
  void UpdateKinematics(const NodeIteratorType& rNode, NodeType& rDEMNode);

protected:

  // Member Variables
  ModelPart& mrModelPart;
  ModelPart& mrDEMModelPart;

};  // Class

}  // namespace Kratos
#endif /* KRATOS_GENERATE_DEM_PROCESS defined */
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

#if !defined(KRATOS_TRANSFER_NODAL_FORCES_TO_FEM_PROCESS)
#define KRATOS_TRANSFER_NODAL_FORCES_TO_FEM_PROCESS


#include "includes/model_part.h"
#include "processes/process.h"
#include "fem_to_dem_application_variables.h"
#include "solid_mechanics_application_variables.h"
#include "DEM_application_variables.h"
#include "includes/variables.h"

namespace Kratos {

typedef std::size_t SizeType;
typedef Node<3> NodeType;
typedef ModelPart::NodesContainerType::iterator NodeIteratorType;

class TransferNodalForcesToFem : public Process 
{
 public:

  /// Pointer definition of ApplyMultipointConstraintsProcess
  KRATOS_CLASS_POINTER_DEFINITION(TransferNodalForcesToFem);

  typedef std::size_t IndexType;

  // Constructor
  TransferNodalForcesToFem(ModelPart& rModelPart, ModelPart& rDemModelPart);

  // Destructor
  ~TransferNodalForcesToFem() override = default;

  void operator()() { Execute(); }

  void Execute() override;

protected:

  // Member Variables
  ModelPart& mrModelPart;
  ModelPart& mrDEMModelPart;

};  // Class

}  // namespace Kratos
#endif /* KRATOS_TRANSFER_NODAL_FORCES_TO_FEM_PROCESS defined */
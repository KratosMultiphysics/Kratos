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

#if !defined(KRATOS_REMOVE_ALONE_DEM_ELEMENTS_PROCESS)
#define KRATOS_REMOVE_ALONE_DEM_ELEMENTS_PROCESS

#include "processes/process.h"
#include "fem_to_dem_application_variables.h"
#include "solid_mechanics_application_variables.h"
#include "includes/variables.h"

namespace Kratos {

typedef std::size_t SizeType;
typedef Node<3> NodeType;
typedef ModelPart::NodesContainerType::iterator NodeIteratorType;

class RemoveAloneDEMElementsProcess : public Process 
{
 public:

    /// Pointer definition of ApplyMultipointConstraintsProcess
  KRATOS_CLASS_POINTER_DEFINITION(RemoveAloneDEMElementsProcess);

  typedef std::size_t IndexType;

  // Constructor
  RemoveAloneDEMElementsProcess(ModelPart& rModelPart, ModelPart& rDemModelPart);

  // Destructor
  ~RemoveAloneDEMElementsProcess() override = default;

  void operator()() { Execute(); }

  void Execute() override;

protected:

  // Member Variables
  ModelPart& mrModelPart;
  ModelPart& mrDEMModelPart;

};  // Class

}  // namespace Kratos
#endif /* KRATOS_REMOVE_ALONE_DEM_ELEMENTS_PROCESS defined */
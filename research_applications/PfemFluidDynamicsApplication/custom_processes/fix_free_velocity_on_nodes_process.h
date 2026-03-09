//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics PfemFluidDynamics Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Carlos Eulogio Flores
//

#if !defined(KRATOS_PFEM_FIX_FREE_VELOCITY_ON_NODES_PROCESS)
#define KRATOS_PFEM_FIX_FREE_VELOCITY_ON_NODES_PROCESS

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos {

class KRATOS_API(PFEM_FLUID_DYNAMICS_APPLICATION) PFEMFixFreeVelocityOnNodesProcess : public Process
{
public:

  /// Pointer definition of PFEMFixFreeVelocityOnNodesProcess
  KRATOS_CLASS_POINTER_DEFINITION(PFEMFixFreeVelocityOnNodesProcess);

  // Constructor
  PFEMFixFreeVelocityOnNodesProcess(ModelPart& rModelPart, const bool rFreeOrFix);

  void operator()() { Execute(); }

  void Execute() override;

protected:

  // Member Variables
  ModelPart& mrModelPart;
  bool mFreeOrFix;

};

}  // namespace Kratos

#endif // KRATOS_FIX_FREE_VELOCITY_ON_NODES_PROCESS defined

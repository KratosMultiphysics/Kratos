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

#if !defined(KRATOS_COMPUTE_SAND_PRODUCTION_PROCESS)
#define KRATOS_COMPUTE_SAND_PRODUCTION_PROCESS


#include "includes/model_part.h"
#include "processes/process.h"
#include "fem_to_dem_application_variables.h"

namespace Kratos {

typedef std::size_t SizeType;

class ComputeSandProduction : public Process 
{
 public:

  /// Pointer definition of ApplyMultipointConstraintsProcess
  KRATOS_CLASS_POINTER_DEFINITION(ComputeSandProduction);

  typedef std::size_t IndexType;

  // Constructor
  ComputeSandProduction(ModelPart& rModelPart);

  // Destructor
  ~ComputeSandProduction() override = default;

  void operator()() { Execute(); }

  void Execute() override;


protected:

  // Member Variables
  ModelPart& mrModelPart;

};  // Class ComputeSandProduction

}  // namespace Kratos
#endif /* KRATOS_EXTEND_PRESSURE_PROCESS defined */
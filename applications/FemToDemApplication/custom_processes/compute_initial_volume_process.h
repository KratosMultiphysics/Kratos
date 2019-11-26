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

#if !defined(KRATOS_COMPUTE_INITIAL_VOLUME_PROCESS)
#define KRATOS_COMPUTE_INITIAL_VOLUME_PROCESS

#include "processes/process.h"
#include "fem_to_dem_application_variables.h"

namespace Kratos {

typedef std::size_t SizeType;

/** 
 * @class ComputeInitialVolumeProcess
 * @ingroup FemToDemApplication 
 * @brief Computes the initial volume of the blast circle and
 * assigns it to the nodes
 * @author Alejandro Cornejo
 */
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

  /**
   * @brief Computes the initial volume of the blast circle and assigns it to the nodes
   */
  void Execute() override;

  /**
   * @brief Gets the pressure id of a SubModelPart
   */
  int GetPressureIdSubModel(const std::string &rSubModelName);

  /**
   * @brief Computes the initial volume of the blast circle and assigns it to the nodes
   */
  double ComputeInitialVolumeSubModel(const ModelPart &rSubModel);

  /**
   * @brief assigns the initial volume to the nodes -> SetValue(--)
   */
  void AssignInitialVolumeToNodes(const ModelPart &rSubModel, const double InitialVolume);

protected:

  // Member Variables
  ModelPart& mrModelPart;
  std::size_t mDimension;
  std::string mPressureName;

};  // Class

}  // namespace Kratos
#endif /* KRATOS_COMPUTE_INITIAL_VOLUME_PROCESS defined */
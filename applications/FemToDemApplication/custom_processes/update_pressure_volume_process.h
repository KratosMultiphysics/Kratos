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

#if !defined(KRATOS_UPDATE_PRESSURE_VOLUME_PROCESS)
#define KRATOS_UPDATE_PRESSURE_VOLUME_PROCESS

#include "processes/process.h"
#include "fem_to_dem_application_variables.h"

namespace Kratos {

typedef std::size_t SizeType;

/** 
 * @class UpdatePressureVolumeProcess
 * @ingroup FemToDemApplication 
 * @brief updates the volume according to a pressure load
 * @author Alejandro Cornejo
 */
class UpdatePressureVolumeProcess : public Process 
{
 public:

  /// Pointer definition of ApplyMultipointConstraintsProcess
  KRATOS_CLASS_POINTER_DEFINITION(UpdatePressureVolumeProcess);

  typedef std::size_t IndexType;

  // Constructor
  UpdatePressureVolumeProcess(ModelPart& rModelPart);

  // Destructor
  ~UpdatePressureVolumeProcess() override = default;

  void operator()() { Execute(); }

  /**
   * @brief updates the volume according to a pressure load
   */
  void Execute() override;

  /**
   * @brief returns the pressure Id of the element
   */
  int GetPressureId(ModelPart::ElementIterator itElem);

protected:

  // Member Variables
  ModelPart& mrModelPart;
  std::size_t mDimension;
  std::string mPressureName;

};  // Class

}  // namespace Kratos
#endif /* KRATOS_UPDATE_PRESSURE_VOLUME_PROCESS defined */
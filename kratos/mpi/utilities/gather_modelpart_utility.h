//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

#if !defined(KRATOS_GATHER_MODELPART_UTILITY)
#define  KRATOS_GATHER_MODELPART_UTILITY

// System includes

/* Project includes */
#include "containers/variable.h"
#include "containers/array_1d.h"
#include "includes/model_part.h"

namespace Kratos
{
class KRATOS_API(KRATOS_MPI_CORE) GatherModelPartUtility
{
public:

  KRATOS_CLASS_POINTER_DEFINITION(GatherModelPartUtility);

  typedef ModelPart::NodesContainerType NodesContainerType;

  typedef ModelPart::ElementsContainerType ElementsContainerType;

  typedef ModelPart::ConditionsContainerType ConditionsContainerType;

  ///This function is designed to obtain data from "origin_model_part.GetMesh(mesh_id)", copy it to a new model part
  ///and make rank "gather_rank" to have a copy of it. Transferred nodes will be treated as ghost on the gather_rank
  ///@param gather_rank --> mpi rank to which the model part is gathered
  ///@param origin_model_part --> model part on which the origin mesh is contained
  ///@param mesh_id --> id of the mesh which contains the data
  ///@param destination_model_part --> model part to which we gather the data
  GatherModelPartUtility(int gather_rank,
                         ModelPart& origin_model_part,
                         int mesh_id,
                         ModelPart& destination_model_part);

  void GatherOnMaster();

  template <class TDataType>
  void GatherOnMaster(const Variable<TDataType>& ThisVariable);

  template <class TDataType>
  void ScatterFromMaster(const Variable<TDataType>& ThisVariable);

  private:
  ModelPart& mr_model_part;
  int mgather_rank;

};

extern template void GatherModelPartUtility::GatherOnMaster(const Variable<double>&);
extern template void GatherModelPartUtility::GatherOnMaster(const Variable<array_1d<double, 3>>&);
extern template void GatherModelPartUtility::ScatterFromMaster(const Variable<double>&);
extern template void GatherModelPartUtility::ScatterFromMaster(const Variable<array_1d<double, 3>>&);

} // namespace Kratos.

#endif // KRATOS_GATHER_MODELPART_UTILITY  defined

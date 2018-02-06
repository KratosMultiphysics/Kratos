//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_HDF5_PARTITIONED_MODEL_PART_IO_H_INCLUDED)
#define KRATOS_HDF5_PARTITIONED_MODEL_PART_IO_H_INCLUDED

// System includes
#include <tuple>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/io.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_model_part_io_base.h"

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

/// A class for partitioned IO of a model part in HDF5.
class PartitionedModelPartIO : public ModelPartIOBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PartitionedModelPartIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    PartitionedModelPartIO(Parameters Settings, File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{
    bool ReadNodes(NodesContainerType& rNodes) override;

    void WriteNodes(NodesContainerType const& rNodes) override;

    void ReadElements(NodesContainerType& rNodes,
                      PropertiesContainerType& rProperties,
                      ElementsContainerType& rElements) override;

    void WriteElements(ElementsContainerType const& rElements) override;

    void ReadConditions(NodesContainerType& rNodes,
                        PropertiesContainerType& rProperties,
                        ConditionsContainerType& rConditions) override;

    void WriteConditions(ConditionsContainerType const& rConditions) override;

    void ReadModelPart(ModelPart& rModelPart) override;

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void Check();

    ///@}

private:
    ///@name Member Variables
    ///@{

    ///@}

    ///@name Private Operations
    ///@{
    std::tuple<unsigned, unsigned> GetPartitionStartIndexAndBlockSize(std::string Path) const;

    void WritePartitionIndex(std::string Path, NodesContainerType const& rGhostNodes);

    void ReadAndAssignPartitionIndex(std::string Path, ModelPart& rModelPart) const;
    ///@}
};

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_PARTITIONED_HDF5_MODEL_PART_IO_H_INCLUDED defined

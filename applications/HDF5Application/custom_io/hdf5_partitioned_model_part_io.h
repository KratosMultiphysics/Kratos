//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//                  Suneth Warnakulasuriya
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_model_part_io.h"

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

/// A class for partitioned IO of a model part in HDF5.
class PartitionedModelPartIO : public ModelPartIO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PartitionedModelPartIO);

    typedef ModelPartIO BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    PartitionedModelPartIO(
        Parameters Settings,
        File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    bool ReadNodes(NodesContainerType& rNodes) override;

    void WriteNodes(NodesContainerType const& rNodes) override;

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void Check();

    void ReadParitionIndices(ModelPart& rModelPart) override;

    void SetCommunicator(ModelPart& rModelPart) const override;

    ///@}

private:
    ///@name Private Operations
    ///@{

    void WritePartitionIndex(const std::string& rPath, NodesContainerType const& rGhostNodes);

    ///@}
};

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.
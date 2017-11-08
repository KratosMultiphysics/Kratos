//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_HDF5_PARTITIONED_MODEL_PART_IO_H_INCLUDED)
#define KRATOS_HDF5_PARTITIONED_MODEL_PART_IO_H_INCLUDED

// System includes
#include <vector>
#include <string>
#include <tuple>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_file.h"

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

/// A class for partitioned IO of a model part in HDF5.
class PartitionedModelPartIO : public IO
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
    PartitionedModelPartIO(Parameters& rParams, File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{
    bool ReadNodes(NodesContainerType& rNodes) override;

    std::size_t ReadNodesNumber() override;

    void WriteNodes(NodesContainerType const& rNodes) override;

    void ReadProperties(PropertiesContainerType& rProperties) override;
    
    void WriteProperties(Properties const& rProperties) override;
    
    void WriteProperties(PropertiesContainerType const& rProperties) override;

    void ReadElements(NodesContainerType& rNodes,
                      PropertiesContainerType& rProperties,
                      ElementsContainerType& rElements) override;

    void WriteElements(ElementsContainerType const& rElements) override;

    void ReadConditions(NodesContainerType& rNodes,
                        PropertiesContainerType& rProperties,
                        ConditionsContainerType& rConditions) override;

    void WriteConditions(ConditionsContainerType const& rConditions) override;

    void ReadModelPart(ModelPart& rModelPart) override;

    void WriteModelPart(ModelPart& rModelPart) override;

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    File& GetFile() const;

    void Check();

    ///@}

private:
    ///@name Member Variables
    ///@{

    File::Pointer mpFile;
    std::string mPrefix;
    std::vector<std::string> mElementNames;
    std::vector<const Element*> mElementPointers;
    std::vector<std::string> mConditionNames;
    std::vector<const Condition*> mConditionPointers;

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

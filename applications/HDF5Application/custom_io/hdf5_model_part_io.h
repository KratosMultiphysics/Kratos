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
#include <string>
#include <tuple>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "includes/kratos_parameters.h"

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

/// A class for serial IO of a model part in HDF5.
class KRATOS_API(HDF5_APPLICATION) ModelPartIO : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ModelPartIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ModelPartIO(
        Parameters Settings,
        File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{

    bool ReadNodes(NodesContainerType& rNodes) override;

    std::size_t ReadNodesNumber() override;

    void WriteNodes(NodesContainerType const& rNodes) override;

    void ReadProperties(PropertiesContainerType& rThisProperties) override;

    void WriteProperties(PropertiesContainerType const& rThisProperties) override;

    void ReadElements(
        NodesContainerType& rNodes,
        PropertiesContainerType& rProperties,
        ElementsContainerType& rElements) override;

    void WriteElements(ElementsContainerType const& rElements) override;

    void ReadConditions(
        NodesContainerType& rNodes,
        PropertiesContainerType& rProperties,
        ConditionsContainerType& rConditions) override;

    void WriteConditions(ConditionsContainerType const& rConditions) override;

    void WriteModelPart(ModelPart& rModelPart) override;

    void ReadModelPart(ModelPart& rModelPart) override;

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    virtual void ReadParitionIndices(ModelPart& rModelPart);

    /**
     * @brief Set the Communicator for the given @p rModelPart
     * @details This method sets the communicator of the given @p rModelPart,
     *          sets the local mesh entities correctly. This is done recursively for
     *          all the sub model parts as well.
     *
     * @param rModelPart Model part to set the communicator.
     */
    virtual void SetCommunicator(ModelPart& rModelPart) const;

    ///@}
    ///@name Member Variables
    ///@{

    File::Pointer mpFile;

    std::string mPrefix;

    ///@}

private:
    ///@name Private member variables
    ///@{

    bool mWriteEntityPropertyIds;

    Parameters mCustomAttributes;

    ///@}
    ///@name Private life cycle
    ///@{

    ModelPartIO(
        const std::string& rPrefix,
        const bool WriteEntityPropertyIds,
        File::Pointer pFile);

    ///@}
    ///@name Private Operations
    ///@{

    std::vector<std::size_t> ReadContainerIds(std::string const& rPath) const;

    void WriteNodeIds(
        const std::string& rPath,
        const NodesContainerType& rNodes);

    void WriteSubModelParts(
        const ModelPart::SubModelPartsContainerType& rSubModelPartsContainer,
        const std::string& GroupName);

    void ReadSubModelParts(
        ModelPart& rModelPart,
        const std::string& rPath);

    ///@}
};

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.
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
//                  Suneth Warnakulasuriya, https://github.com/sunethwarna
//

#pragma once

// System includes
#include <string>
#include <vector>
#include <map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "custom_io/hdf5_file.h"

namespace Kratos
{
class Parameters;

namespace HDF5
{
///@addtogroup HDF5Application
///@{

///@name Kratos Classes
///@{

/// A class for IO of element data in HDF5.
template <class TContainerType, class TContainerIO, class... TComponents>
class KRATOS_API(HDF5_APPLICATION) ContainerComponentIO
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ContainerComponentIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /**
     * @brief Construct a new Container Component IO to read and write from kratos containers using variables
     * @details This class is the inner most class which can be used to read from or write to
     *          different containers in Kratos. followings are the supported containers.
     *
     *               1. Nodal historical data value container
     *               2. Nodal non-historical data value container
     *               3. Condition data value container
     *               4. Element data value container
     *               5. Condition container with gauss point data values
     *               6. Element container gauss point data values
     *               7. Vertex container
     *               8. PropertiesContainer
     *
     *          This supports following variable types:
     *              - int
     *              - double
     *              - array_1d<double, 3>
     *              - array_1d<double, 4>
     *              - array_1d<double, 6>
     *              - array_1d<double, 9>
     *              - Vector
     *              - Matrix
     *
     *          This also allows writing custom attributes to the HDF5 dataset. These can be anything,
     *          and should be passed as a @ref Kratos::Parameters object. The attribute keys starting
     *          with "__" are reserved for the HDF5Application and must not be used as custom attribute
     *          keys. Followings are list of reserved keys:
     *              - "__data_dimension"
     *              - "__data_shape"
     *              - "__data_availability"
     *              - "__container_type"
     *              - "__data_name"
     *              - "__data_location"
     *              - "__mesh_location"
     *
     * @throws Runtime error when any one of the given components are not found in at least one entity
     *         in at least on rank.
     *
     * @todo Remove the rLegacySuffix once the python side is properly fixed.
     *
     * @param Settings          Settings to be used with this Component IO.
     * @param pFile             HDF5 file to be used for reading and writing components.
     * @param rLegacySuffix     Suffix used in legacy HDF5 writing. Planned to be removed once python side is fixed.
     */
    ContainerComponentIO(
        Parameters Settings,
        File::Pointer pFile,
        const std::string& rLegacySuffix = "");

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Write data from the local container having the type TContainerType to the HDF5 file for the components specified in Settings.
     * @details This method can be used in either shared memory/distributed memory parallelized architectures.
     *
     * @param rModelPart            Model part to get the local container of type TContainerType to get the data.
     * @param rContainerDataIO      Reading mechanism from where the data will be read from the TContainerType.
     * @param Attributes            Custom attributes to be written in to the HDF5 dataset.
     */
    void Write(
        const ModelPart& rModelPart,
        const TContainerIO& rContainerDataIO,
        const Parameters Attributes);

    /**
     * @brief Read an hdf5 dataset and populate data to the container with type TContainerType for the components specified in Settings.
     * @details This method can be used in either shared memory/distributed memory parallelized architectures.
     *          It will synchronize the written data to populate the ghost meshes correctly.
     *          The returning map will have the component name as the key, and its attributes as values.
     *
     * @param rModelPart                            Model part to read the data into from the HDF5 dataset.
     * @param rContainerDataIO                      Writing mechanism to where the data will be written in rModelPart from the HDF5 dataset.
     * @return std::map<std::string, Parameters>    A map containing component name as the key and each components' attributes as values.
     */
    std::map<std::string, Parameters> Read(
        ModelPart& rModelPart,
        const TContainerIO& rContainerDataIO);

    /**
     * @brief Write data from the passed local data container to the HDF5 dataset for the components specified in Settings.
     *
     * @param rLocalContainer       Local container with the data.
     * @param rContainerDataIO      Reading mechanism specifying where to read the data from the local container.
     * @param Attributes            Custom attributes to be written in to the HDF5 dataset.
     */
    void Write(
        const TContainerType& rLocalContainer,
        const TContainerIO& rContainerDataIO,
        const Parameters Attributes);

    /**
     * @brief Read an hdf5 dataset and populate data to the provided local container for the components specified in Settings.
     * @details This method can be used in either shared memory/distributed memory parallelized architectures.
     *          It will synchronize the written data to populate the ghost meshes correctly.
     *          The returning map will have the component name as the key, and its attributes as values.
     *
     * @param rLocalContainer                       Local container with the data.
     * @param rContainerDataIO                      Writing mechanism specifying where to write the data in the local container.
     * @param rCommunicator                         Communicator to synchronize and populate ghost meshes.
     * @return std::map<std::string, Parameters>    A map containing component name as the key and each components' attributes as values.
     */
    std::map<std::string, Parameters> Read(
        TContainerType& rLocalContainer,
        const TContainerIO& rContainerDataIO,
        Communicator& rCommunicator);

    /**
     * @brief Reads only the attributes for components specified in the Settings.
     * @details This method only reads the attributes for the specified components. Does not read the dataset.
     *
     * @return std::map<std::string, Parameters>    A map containing component name as the key and each components' attributes as values.
     */
    std::map<std::string, Parameters> ReadAttributes();

    ///@}

protected:
    ///@name Protected member variables
    ///@{

    File::Pointer mpFile;

    std::vector<std::string> mComponentNames;

    std::string mComponentPrefix;

    ///@}

private:
    ///@name Private static member variables
    ///@{

    const static inline std::vector<std::string> ReservedAttributeKeys = {
                                            "__data_dimension",
                                            "__data_shape",
                                            "__data_availability",
                                            "__container_type",
                                            "__data_name",
                                            "__data_location",
                                            "__mesh_location"
                                        };

    ///}
    ///@name Private static operations
    ///@{

    static void CheckReservedAttributes(const Parameters Attributes);

    static void RemoveReservedAttributes(Parameters Attributes);

    ///@}
    ///@name Private Operations
    ///@{

    template<class TComponentType>
    bool WriteComponents(
        const std::string& rComponentName,
        const TContainerIO& rContainerDataIO,
        const TContainerType& rLocalContainer,
        Parameters Attributes,
        WriteInfo& rInfo);

    template<class TComponentType>
    bool ReadComponents(
        const std::string& rComponentName,
        const TContainerIO& rContainerDataIO,
        TContainerType& rLocalContainer,
        Communicator& rCommunicator,
        std::map<std::string, Parameters>& rAttributesMap,
        const IndexType StartIndex,
        const IndexType BlockSize);

    ///@}

}; // class ContainerComponentIO.

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

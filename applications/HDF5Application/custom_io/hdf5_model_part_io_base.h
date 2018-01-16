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

#if !defined(KRATOS_HDF5_MODEL_PART_IO_BASE_H_INCLUDED)
#define KRATOS_HDF5_MODEL_PART_IO_BASE_H_INCLUDED

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "includes/kratos_parameters.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_file.h"
#include "custom_utilities/hdf5_connectivities_io.h"

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

class ModelPartIOBase : public IO
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ModelPartIOBase(Parameters Settings, File::Pointer pFile);

    ///@}
    ///@name Operations
    ///@{
    std::size_t ReadNodesNumber() override;

    void ReadProperties(PropertiesContainerType& rThisProperties) override;

    void WriteProperties(Properties const& rThisProperties) override;

    void WriteProperties(PropertiesContainerType const& rThisProperties) override;

    void WriteModelPart(ModelPart& rModelPart) override;

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    ~ModelPartIOBase() = default;

    File& GetFile() const;

    ///@}
    ///@name Member Variables
    ///@{

    File::Pointer mpFile;
    std::string mPrefix;
    std::vector<Internals::FileItem> mElementIO;
    std::vector<Internals::FileItem> mConditionIO;

    ///@}

    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_MODEL_PART_IO_BASE_H_INCLUDED defined

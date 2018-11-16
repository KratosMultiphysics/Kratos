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

#if !defined(KRATOS_HDF5_XDMF_CONNECTIVITIES_WRITER_PROCESS_H_INCLUDED)
#define KRATOS_HDF5_XDMF_CONNECTIVITIES_WRITER_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"

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

/// Writes Xdmf connectivities.
/**
 * In partitioned simulations, the node ids are not sorted globally. Thus
 * the Kratos element connectivities cannot be used directly to index the
 * node arrays in the HDF5 file. This utility transforms Kratos element
 * connectivities to Xdmf connectivities, which can be used directly to
 * index the node arrays.
 */
class XdmfConnectivitiesWriterProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(XdmfConnectivitiesWriterProcess);
    typedef std::unordered_map<int, int> IdMapType;
    ///@}
    ///@name Life Cycle
    ///@{
    XdmfConnectivitiesWriterProcess(const std::string& rFileName, const std::string& rPrefix);

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    ///@}
private:
    ///@name Member Variables
    ///@{

    File::Pointer mpFile;
    std::string mPrefix;
    IdMapType mKratosToXdmfIdMap;

    ///@}
    ///@name Private Operations
    ///@{

    void CreateXdmfConnectivities(const std::string& rKratosConnectivitiesPath, const std::string& rXdmfConnectivitiesPath) const;

    ///@}

}; // class XdmfConnectivitiesWriterProcess

///@} // Kratos Classes
///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_XDMF_CONNECTIVITIES_WRITER_PROCESS_H_INCLUDED defined

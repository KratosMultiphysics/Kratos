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

#if !defined(KRATOS_HDF5_POINTS_DATA_H_INCLUDED)
#define KRATOS_HDF5_POINTS_DATA_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"

// Application includes
#include "hdf5_application_define.h"

namespace Kratos
{
namespace HDF5
{

class File;
struct WriteInfo;

namespace Internals
{
///@addtogroup HDF5Application
///@{
///@name Kratos Classes
///@{

/// A class representing points in a mesh.
class PointsData
{
public:
    ///@name Type Definitions
    ///@{
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PointsData);

    ///@}
    ///@name Life Cycle
    ///@{
    ///@}
    ///@name Operations
    ///@{

    inline unsigned size() const
    {
        return mIds.size();
    }

    void ReadData(File& rFile, const std::string& rPath, unsigned StartIndex, unsigned BlockSize);

    void WriteData(File& rFile, const std::string& rPath, WriteInfo& rInfo);

    void CreateNodes(NodesContainerType& rNodes);

    // Fill data from nodes.
    void SetData(NodesContainerType const& rNodes);

    void Clear();
    ///@}
private:
    ///@name Member Variables
    ///@{
    Vector<int> mIds;
    Vector<array_1d<double, 3>> mCoords;
    ///@}
}; // class PointsData

///@} // Kratos Classes
///@} addtogroup
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_POINTS_DATA_H_INCLUDED defined
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

#if !defined(KRATOS_HDF5_POINTS_DATA_H_INCLUDED)
#define KRATOS_HDF5_POINTS_DATA_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "includes/model_part.h"

// Application includes
#include "custom_io/hdf5_file.h"
#include "custom_utilities/hdf5_utils.h"

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{

class PointsData
{
public:
    ///@name Type Definitions
    ///@{
    /// Pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(PointsData);

    typedef ModelPart::NodesContainerType NodesContainerType;

    typedef ModelPart::NodeType NodeType;

    ///@}
    ///@name Life Cycle
    ///@{
    ///@}
    ///@name Operations
    ///@{
    inline Vector<int> const& GetIds() const
    {
        return mIds;
    }

    inline Vector<array_1d<double, 3>> const& GetCoordinates() const
    {
        return mCoords;
    }

    inline unsigned size() const
    {
        return mIds.size();
    }

    void ReadData(File& rFile, std::string Path, unsigned StartIndex, unsigned BlockSize);

    void WriteData(File& rFile, std::string Path);

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

///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_POINTS_DATA_H_INCLUDED defined
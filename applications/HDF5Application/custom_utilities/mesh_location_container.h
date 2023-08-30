//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>
#include <tuple>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{
namespace HDF5
{

class MeshLocationContainer
{
public:
    ///@name Type defintions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(MeshLocationContainer);

    ///@}
    ///@name Public operations
    ///@{

    void Set(
        const int HDF5RankId,
        const int HDF5ProcessId,
        const std::string& rMeshLocation);

    bool Has(
        const int HDF5RankId,
        const int HDF5ProcessId) const;

    std::string Get(
        const int HDF5RankId,
        const int HDF5ProcessId) const;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    std::vector<std::pair<int, int>> mProcessIds;

    std::vector<std::string> mMeshLocations;

    ///@}
    ///@name Friends
    ///@{

    friend class Kratos::Serializer;

    void save(Kratos::Serializer& rSerializer) const;

    void load(Kratos::Serializer& rSerializer);

    ///@}
};

void SetMeshLocationContainer(
    ModelPart& rModelPart,
    MeshLocationContainer::Pointer pMeshLocationContainer);

bool HasMeshLocationContainer(const ModelPart& rModelPart);

MeshLocationContainer& GetMeshLocationContainer(ModelPart& rModelPart);

const MeshLocationContainer& GetMeshLocationContainer(const ModelPart& rModelPart);

void AddProcessId(
    Parameters Settings,
    const int HDF5RankId,
    const int HDF5ProcessId);

bool HasProcessId(const Parameters Settings);

std::pair<int, int> GetProcessId(const Parameters Settings);

std::ostream& operator << (
    std::ostream& rOStream,
    const MeshLocationContainer& rThis);

} // namespace HDF5

} // namespace Kratos
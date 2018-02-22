//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_HDF5_TEST_UTILS_H_INCLUDED)
#define KRATOS_HDF5_TEST_UTILS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_io/hdf5_file.h"

namespace Kratos
{
namespace Testing
{

template<typename T>
HDF5::File::Vector<T> TestVector(std::size_t n=3)
{
    HDF5::File::Vector<T> vec(n);
    for (std::size_t i = 0; i < n; ++i)
        vec(i) = i;
    return vec;
}

template <>
HDF5::File::Vector<array_1d<double, 3>> TestVector(std::size_t n);

template<typename T>
HDF5::File::Matrix<T> TestMatrix(std::size_t m=3, std::size_t n=3)
{
    HDF5::File::Matrix<T> mat(m, n);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < n; ++j)
            mat(i, j) = i * n + j;
    return mat;
}

void CreateTestMesh(HDF5::NodesContainerType& rNodes,
                    HDF5::PropertiesContainerType& rProperties,
                    HDF5::ElementsContainerType& rElements,
                    HDF5::Vector<int>& rElementIds,
                    HDF5::Vector<int>& rPropertiesIds,
                    HDF5::Matrix<int>& rConnectivities);

void CreateTestMesh(HDF5::NodesContainerType& rNodes,
                    HDF5::PropertiesContainerType& rProperties,
                    HDF5::ConditionsContainerType& rConditions,
                    HDF5::Vector<int>& rConditionIds,
                    HDF5::Vector<int>& rPropertiesIds,
                    HDF5::Matrix<int>& rConnectivities);

} // namespace Testing
} // namespace Kratos.

#endif

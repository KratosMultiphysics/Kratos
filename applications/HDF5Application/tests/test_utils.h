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
#include <string>

// External includes

// Project includes
#include "includes/model_part.h"
#include "containers/data_value_container.h"

// Application includes
#include "custom_io/hdf5_file.h"
#include "custom_io/hdf5_file_serial.h"

namespace Kratos
{

class ModelPart;

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

class TestModelPartFactory
{
    public:
        static void CreateModelPart(
            ModelPart& rTestModelPart,
            std::vector<std::string> const& rElements = {},
            std::vector<std::string> const& rConditions = {},
            std::vector<std::string> const& rNodalVariables = {});

        static void AssignNonHistoricalNodalTestData(ModelPart& rTestModelPart,
                                                     std::vector<std::string> const& rNodalVariables = {});

        static void AssignDataValueContainer(DataValueContainer& rData,
                                             Flags& rFlags,
                                             std::vector<std::string> const& rVariables = {});

    private:
        explicit TestModelPartFactory(ModelPart& rTestModelPart);

        void AddNodalVariables(std::vector<std::string> const& rNodalVariables);

        std::size_t AddNodes(std::size_t NumNodes);

        void SetBufferSize(std::size_t BufferSize);

        void AssignNodalTestData(std::vector<std::string> const& rNodalVariables);

        std::size_t AddElements(std::string const& rElement, std::size_t NumElems);

        std::size_t AddConditions(std::string const& rCondition, std::size_t NumConds);

        void AddSubModelParts();

        void AddEmptySubModelPart();

        void AddElementsSubModelPart();

        void AddConditionsSubModelPart();

        void AddElementsAndConditionsSubModelPart();

        ModelPart& mrTestModelPart;
};

void CompareNodes(HDF5::NodesContainerType& rNodes1, HDF5::NodesContainerType& rNodes2);

void CompareElements(HDF5::ElementsContainerType& rElements1, HDF5::ElementsContainerType& rElements2);

void CompareConditions(HDF5::ConditionsContainerType& rConditions1, HDF5::ConditionsContainerType& rConditions2);

void CompareModelParts(ModelPart& rModelPart1, ModelPart& rModelPart2);

void CompareDataValueContainers(DataValueContainer const& rData1, Flags const& rFlags1, DataValueContainer const& rData2, Flags const& rFlags2);

void CompareNonHistoricalNodalData(HDF5::NodesContainerType& rNodes1,
                                   HDF5::NodesContainerType& rNodes2);

HDF5::File::Pointer pGetTestSerialFile();

HDF5::File GetTestFile();

HDF5::FileSerial GetTestSerialFile();

/// Silences HDF5 stderr messages for duration of local scope.
class H5_stderr_muter
{
    H5E_auto2_t old_func;
    void* old_client_data;

public:
    H5_stderr_muter()
    {
        H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
        H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    }
    ~H5_stderr_muter()
    {
        H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);
    }
};

} // namespace Testing
} // namespace Kratos.

#endif

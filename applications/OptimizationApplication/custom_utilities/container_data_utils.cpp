//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "includes/define.h"

// Application includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Include base h
#include "container_data_utils.h"

namespace Kratos
{


double ContainerDataUtils::NormInf(const ContainerDataBase& rContainer)
{
    const auto& r_data = rContainer.GetData();
    return rContainer.GetModelPart().GetCommunicator().GetDataCommunicator().MaxAll(IndexPartition<IndexType>(r_data.size()).for_each<MaxReduction<double>>([&](const IndexType Index) {
        return r_data[Index];
    }));
}

double ContainerDataUtils::InnerProduct(
    const ContainerDataBase& rContainer1,
    const ContainerDataBase& rContainer2)
{
    const auto& r_data_1 = rContainer1.GetData();
    const auto& r_data_2 = rContainer2.GetData();

    KRATOS_ERROR_IF(r_data_1.size() != r_data_2.size())
        << "Data size mismatch in InnerProduct calculation. "
        << "Followings are the given containers: \n"
        << "   Container 1: " << rContainer1 << "\n"
        << "   Container 2: " << rContainer2 << "\n";

    KRATOS_ERROR_IF_NOT(rContainer1.IsCompatibleWithContainerData(rContainer2))
        << "Incompatible containers given for InnerProduct. "
        << "Followings are the given containers: \n"
        << "   Container 1: " << rContainer1 << "\n"
        << "   Container 2: " << rContainer2 << "\n";

    return rContainer1.GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(IndexPartition<IndexType>(r_data_1.size()).for_each<SumReduction<double>>([&](const IndexType Index) {
        return r_data_1[Index] * r_data_2[Index];
    }));
}

}
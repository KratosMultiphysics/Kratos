//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi

#include <utility>
#include <iostream>
#include "testing/testing.h"
#include "utilities/parallel_utilities.h"

namespace Kratos {
namespace Testing {

// Basic Type
KRATOS_TEST_CASE_IN_SUITE(BlockPartitioner, KratosCoreFastSuite)
{
    int nsize = 1e3;
    std::vector<double> data_vector(nsize);
    for(auto& it : data_vector)
        it = 5.0;

    //here we raise every entry of a vector to the power 0.1
    //c++11 version
    BlockPartition<std::vector<double>>(data_vector).for_each(
                                         [](double& item)
    {
        item = std::pow(item, 0.1);
    });

    //c++17 version - no explicit template parameter in constructor and auto in lambda
    // BlockPartition(data_vector).for_each(
    //     [](auto& item){
    //         item = std::pow(item, 0.1);
    //     });

    //error check
    for(auto& item : data_vector)
    {
        KRATOS_CHECK_EQUAL(item, std::pow(5.0, 0.1));
    }

    //here we check for a reduction (computing the sum of all the entries)
    SumReduction<double> reducer;

    //compute a sum reduction - c++11 version
    BlockPartition<std::vector<double>>(data_vector).for_each(
        reducer,
        [](double& item, SumReduction<double>& reduction_helper)
        {
            reduction_helper.mvalue += item;
        }
    );

    //compute a sum reduction - c++17 version
    // BlockPartition(data_vector).for_each(
    //     reducer,
    //     [](auto& item, auto& reduction_helper){
    //         reduction_helper.mvalue += item;
    //         }
    //     );

    //get the value of the final sum
    double final_sum = reducer.mvalue;

    double expected_value = std::pow(5.0, 0.1)*nsize;
    KRATOS_CHECK_NEAR( std::abs(final_sum-expected_value)/std::abs(expected_value), 0.0, 1e-10  );
}

// Basic Type
KRATOS_TEST_CASE_IN_SUITE(IndexPartitioner, KratosCoreFastSuite)
{
    int nsize = 1e3;
    std::vector<double> data_vector(nsize), output(nsize);
    for(auto& it : data_vector)
        it = -1.0;

    //output = 2*data_vector (in parallel, and accessing by index)
    IndexPartition<unsigned int>(data_vector.size()).for_each(
        [&](unsigned int i){
            output[i] = 2.0*data_vector[i];
            }
        );

    for(unsigned int i=0; i<output.size(); ++i)
        KRATOS_CHECK_EQUAL(output[i], -2.0 );
}


} // namespace Testing
} // namespace Kratos

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

KRATOS_TEST_CASE_IN_SUITE(CustomReduction, KratosCoreFastSuite)
{
    int nsize = 1e3;
    std::vector<double> data_vector(nsize);
    for(int i=0; i<nsize; ++i)
        data_vector[i] = -i;

    //we want to find the maximum value and the maximum absolute value.
    //this is a "custom reduction"
    class CustomReducer{
        public:
            double max_value = -std::numeric_limits<double>::max();
            double max_abs = 0.0;
            void ThreadSafeMerge(CustomReducer& rOther){
                #pragma omp critical
                {
                this->max_value = std::max(this->max_value,rOther.max_value);
                this->max_abs   = std::max(this->max_abs,std::abs(rOther.max_abs));
                }
            }
    };

    CustomReducer Reducer;
    IndexPartition<unsigned int>(data_vector.size()).for_each(
        Reducer,
        [&](unsigned int i, CustomReducer& r){
                r.max_value = std::max(r.max_value, data_vector[i]);
                r.max_abs   = std::max(r.max_abs,std::abs(data_vector[i]));
            }
        );

    KRATOS_CHECK_EQUAL(Reducer.max_value, 0.0 );
    KRATOS_CHECK_EQUAL(Reducer.max_abs, nsize-1 );

    class CustomReducerReturnValueVersion{
        public:
            double max_value = std::numeric_limits<double>::lowest();
            double max_abs = 0.0;

            void LocalMerge(double function_return_value){
                this->max_value = std::max(this->max_value,function_return_value);
                this->max_abs   = std::max(this->max_abs,std::abs(function_return_value));
            }
            void ThreadSafeMerge(CustomReducerReturnValueVersion& rOther){
                #pragma omp critical
                {
                this->max_value = std::max(this->max_value,rOther.max_value);
                this->max_abs   = std::max(this->max_abs,std::abs(rOther.max_abs));
                }
            }
    };

    auto ReturnValueReducer = IndexPartition<unsigned int>(data_vector.size()).
        for_reduce<CustomReducerReturnValueVersion>(
            [&](unsigned int i)->double{
                return data_vector[i]; //note that here the lambda returns the values to be reduced
                }
            );
    KRATOS_CHECK_EQUAL(ReturnValueReducer.max_value, 0.0 );
    KRATOS_CHECK_EQUAL(ReturnValueReducer.max_abs, nsize-1 );

}


} // namespace Testing
} // namespace Kratos

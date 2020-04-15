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
    auto reducer = BlockPartition<std::vector<double>>(data_vector).for_reduce<SumReduction<double>>(
        [](double& item)
        {
            return item;
        }
    );

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

    double reference_max = std::numeric_limits<double>::lowest();
    double reference_min = std::numeric_limits<double>::max();
    double reference_sum = 0.0;
    double reference_sub = 0.0;
    for(auto item : data_vector)
    {
        reference_max = std::max(reference_max, item);
        reference_min = std::min(reference_min, item);
        reference_sum += item;
        reference_sub -= item;
    }
    class CustomReducer{
        public:
            double max_value = -std::numeric_limits<double>::max();
            double max_abs = 0.0;

            void LocalMerge(double function_return_value){
                this->max_value = std::max(this->max_value,function_return_value);
                this->max_abs   = std::max(this->max_abs,std::abs(function_return_value));
            }
            void ThreadSafeMerge(CustomReducer& rOther){
                #pragma omp critical
                {
                this->max_value = std::max(this->max_value,rOther.max_value);
                this->max_abs   = std::max(this->max_abs,std::abs(rOther.max_abs));
                }
            }
    };

    auto partition = IndexPartition<unsigned int>(data_vector.size());
    auto ReturnValueReducer = partition.for_reduce<CustomReducer>([&](unsigned int i){
            return data_vector[i]; //note that here the lambda returns the values to be reduced
        });

    KRATOS_CHECK_EQUAL(ReturnValueReducer.max_value, 0.0 );
    KRATOS_CHECK_EQUAL(ReturnValueReducer.max_abs, nsize-1 );

    typedef CombinedReduction< SumReduction<double>,
                               MinReduction<double>,
                               MaxReduction<double>,
                               SubReduction<double>
            > MultipleReduction;

    auto combined = IndexPartition<unsigned int>(data_vector.size()).
        for_reduce<MultipleReduction>(
            [&](unsigned int i){
                return std::tuple<double,double,double,double>{
                    data_vector[i],data_vector[i],data_vector[i],data_vector[i]}; //note that here the lambda returns the values to be reduced
                }
            );
    KRATOS_CHECK_EQUAL(std::get<0>(combined.GetValue()).mvalue, reference_sum );
    KRATOS_CHECK_EQUAL(std::get<1>(combined.GetValue()).mvalue, reference_min );
    KRATOS_CHECK_EQUAL(std::get<2>(combined.GetValue()).mvalue, reference_max );
    KRATOS_CHECK_EQUAL(std::get<3>(combined.GetValue()).mvalue, reference_sub );
    // KRATOS_CHECK_EQUAL(combined.GetValue().second, reference_min );
}


} // namespace Testing
} // namespace Kratos

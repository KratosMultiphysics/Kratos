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
#include "utilities/openmp_utils.h"

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
    BlockPartition<std::vector<double>>(data_vector).for_each(
                                         [](double& item)
    {
        item = std::pow(item, 0.1);
    });

    //error check
    for(auto& item : data_vector)
    {
        KRATOS_CHECK_EQUAL(item, std::pow(5.0, 0.1));
    }

    //shorter form
    block_for_each(data_vector, [](double& item){
            item = std::pow(5.0, 0.1);
    });

    //error check
    for(auto& item : data_vector)
    {
        KRATOS_CHECK_EQUAL(item, std::pow(5.0, 0.1));
    }

    //here we check for a reduction (computing the sum of all the entries)
    auto final_sum = BlockPartition<std::vector<double>>(data_vector).for_each<SumReduction<double>>(
        [](double& item)
        {
            return item;
        }
    );

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
            typedef std::tuple<double,double> value_type;
            double max_value = -std::numeric_limits<double>::max();
            double max_abs = 0.0;

            value_type GetValue()
            {
                value_type values;
                std::get<0>(values) = max_value;
                std::get<1>(values) = max_abs;
                return values;
            }

            void LocalReduce(double function_return_value){
                this->max_value = std::max(this->max_value,function_return_value);
                this->max_abs   = std::max(this->max_abs,std::abs(function_return_value));
            }
            void ThreadSafeReduce(CustomReducer& rOther){
                #pragma omp critical
                {
                this->max_value = std::max(this->max_value,rOther.max_value);
                this->max_abs   = std::max(this->max_abs,std::abs(rOther.max_abs));
                }
            }
    };

    auto partition = IndexPartition<unsigned int>(data_vector.size());

    double max_value,max_abs;
    std::tie(max_value,max_abs) = partition.for_each<CustomReducer>([&](unsigned int i){
            return data_vector[i]; //note that here the lambda returns the values to be reduced
        });

    KRATOS_CHECK_EQUAL(max_value, 0.0 );
    KRATOS_CHECK_EQUAL(max_abs, nsize-1 );

    //same but with short form with block version
    std::tie(max_value,max_abs) = block_for_each<CustomReducer>(data_vector,[&](double& item){
            return item; //note that here the lambda returns the values to be reduced
        });

    KRATOS_CHECK_EQUAL(max_value, 0.0 );
    KRATOS_CHECK_EQUAL(max_abs, nsize-1 );



    //******************************************************************************************
    //here we reduce at once, the sum,min,max and sub
    typedef CombinedReduction< SumReduction<double>,
                               MinReduction<double>,
                               MaxReduction<double>,
                               SubReduction<double>
            > MultipleReduction;

    //auto reduction_res
    double sum,min,max,sub;
    std::tie(sum,min,max,sub) = IndexPartition<unsigned int>(data_vector.size()).
        for_each<MultipleReduction>(
            [&](unsigned int i){
                    double to_sum = data_vector[i];
                    double to_max = data_vector[i];
                    double to_min = data_vector[i];
                    double to_sub = data_vector[i];
                    return std::make_tuple( to_sum, to_max, to_min, to_sub ); //note that these may have different types
                }
            );
    KRATOS_CHECK_EQUAL(sum, reference_sum );
    KRATOS_CHECK_EQUAL(min, reference_min );
    KRATOS_CHECK_EQUAL(max, reference_max );
    KRATOS_CHECK_EQUAL(sub, reference_sub );
}

KRATOS_TEST_CASE_IN_SUITE(OmpVsPureC11, KratosCoreFastSuite)
{
    int nsize = 1e7;
    std::vector<double> data_vector(nsize), output(nsize);
    for(int i=0; i<nsize; ++i)
        data_vector[i] = i;

    //check ability to handle exceptions in pure c++ - DELIBERATELY THROWING AN EXCEPTION!
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        IndexPartition<unsigned int>(data_vector.size()).for_pure_c11([&](unsigned int i){
                if(i==0)
                    KRATOS_ERROR << "test error on thread 0";
                }
            );
        ,
        "test error on thread 0"
        );

    //benchmark openmp vs pure c++11 impementation in a simple loop
    double start_omp = OpenMPUtils::GetCurrentTime();
    IndexPartition<unsigned int>(nsize).for_each(
            [&](unsigned int i){
                    output[i] = std::pow(data_vector[i],0.01);
                }
            );
    double stop_omp = OpenMPUtils::GetCurrentTime();
    std::cout << "omp time = " << stop_omp-start_omp << std::endl;

    double start_pure= OpenMPUtils::GetCurrentTime();
    IndexPartition<unsigned int>(nsize).for_pure_c11(
            [&](unsigned int i){
                    output[i] = std::pow(data_vector[i],0.01);
                }
            );
    double stop_pure = OpenMPUtils::GetCurrentTime();
    std::cout << "pure c++11 time = " << stop_pure-start_pure << std::endl;



}


} // namespace Testing
} // namespace Kratos

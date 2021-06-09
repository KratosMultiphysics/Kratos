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
//                   Philipp Bucher

// System includes
#include <utility>
#include <numeric>
#include <iostream>

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/parallel_utilities.h"
#include "utilities/builtin_timer.h"

namespace Kratos {
namespace Testing {

namespace { // internals used for testing

template< std::size_t TSize>
class RHSElement
{
    public:
    explicit RHSElement(const double Val) : mRHSVal(Val) {}
    void CalculateRHS(std::vector<double>& rVector)
    {
        if (rVector.size() != TSize) { rVector.resize(TSize); }
        std::fill(rVector.begin(), rVector.end(), mRHSVal);
    }
    double GetAccumRHSValue() {return mAccumRHSValue;}
    void SetAccumRHSValue(double Value) {mAccumRHSValue = Value;}

    private:
    double mRHSVal;
    double mAccumRHSValue = 0.0;
};

} // namespace Internals


// Basic Type
KRATOS_TEST_CASE_IN_SUITE(BlockPartitioner, KratosCoreFastSuite)
{
    int nsize = 1e3;
    std::vector<double> data_vector(nsize, 5.0);

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
KRATOS_TEST_CASE_IN_SUITE(BlockPartitionerConstContainer, KratosCoreFastSuite)
{
    int nsize = 1e3;
    const std::vector<double> data_vector(nsize, 5.0);

    //here we check for a reduction (computing the sum of all the entries)
    auto final_sum = BlockPartition<decltype(data_vector)>(data_vector).for_each<SumReduction<double>>(
        [](const double item)
        {
            return item;
        }
    );

    //here we check for a reduction (computing the sum of all the entries)
    auto final_sum_short = block_for_each<SumReduction<double>>(data_vector,
        [](const double item)
        {
            return item;
        }
    );

    const double expected_value = 5.0*nsize;
    KRATOS_CHECK_DOUBLE_EQUAL(final_sum, expected_value);
    KRATOS_CHECK_DOUBLE_EQUAL(final_sum_short, expected_value);
}

// Basic Type
KRATOS_TEST_CASE_IN_SUITE(IndexPartitioner, KratosCoreFastSuite)
{
    int nsize = 1e3;
    std::vector<double> data_vector(nsize, -1.0), output(nsize);

    //output = 2*data_vector (in parallel, and accessing by index)
    IndexPartition<unsigned int>(data_vector.size()).for_each(
        [&](unsigned int i){
            output[i] = 2.0*data_vector[i];
            }
        );

    for(unsigned int i=0; i<output.size(); ++i)
        KRATOS_CHECK_EQUAL(output[i], -2.0 );
}

KRATOS_TEST_CASE_IN_SUITE(BlockPartitionerThreadLocalStorage, KratosCoreFastSuite)
{
    constexpr std::size_t vec_size = 6;
    constexpr std::size_t n_elems = 10e2;
    constexpr double tol = 1e-8;

    using RHSElementType = RHSElement<vec_size>;

    std::vector<double> rhs_vals(n_elems);
    for (std::size_t i=0; i<n_elems; ++i) {
        rhs_vals[i] = (i%12) * 1.889;
    }

    std::vector<RHSElementType> elements;
    for (std::size_t i=0; i<rhs_vals.size(); ++i) {
        elements.push_back(RHSElementType(rhs_vals[i]));
    }

    const double exp_sum = (std::accumulate(rhs_vals.begin(), rhs_vals.end(), 0.0)) * vec_size;

    auto tls_lambda_manual_reduction = [](RHSElementType& rElem, std::vector<double>& rTLS)
    {
        rElem.CalculateRHS(rTLS);
        double rhs_sum = std::accumulate(rTLS.begin(), rTLS.end(), 0.0);
        rElem.SetAccumRHSValue(rhs_sum);
    };

    auto tls_lambda_reduction = [](RHSElementType& rElem, std::vector<double>& rTLS)
    {
        rElem.CalculateRHS(rTLS);
        return std::accumulate(rTLS.begin(), rTLS.end(), 0.0);
    };


    // Manual Reduction, long form
    // here the TLS is constructed on the fly. This is the "private" approach of OpenMP
    // the result is checked with a "manual reduction"
    BlockPartition<std::vector<RHSElementType>>(elements).for_each(std::vector<double>(), tls_lambda_manual_reduction);

    const double sum_elem_rhs_vals = std::accumulate(elements.begin(), elements.end(), 0.0, [](double acc, RHSElementType& rElem){
        return acc + rElem.GetAccumRHSValue();
    });

    KRATOS_CHECK_NEAR(sum_elem_rhs_vals, exp_sum, tol);


    // Manual Reduction, short form
    // here the TLS is constructed on the fly. This is the "private" approach of OpenMP
    // the result is checked with a "manual reduction"
    block_for_each(elements, std::vector<double>(), tls_lambda_manual_reduction);

    const double sum_elem_rhs_vals_short = std::accumulate(elements.begin(), elements.end(), 0.0, [](double acc, RHSElementType& rElem){
        return acc + rElem.GetAccumRHSValue();
    });

    KRATOS_CHECK_NEAR(sum_elem_rhs_vals_short, exp_sum, tol);


    // Reduction, long form
    // here the TLS is constructed beforehand. This is the "firstprivate" approach of OpenMP
    // checking the results using reduction
    std::vector<double> tls(6);
    const double final_sum = BlockPartition<std::vector<RHSElementType>>(elements).for_each<SumReduction<double>>(tls, tls_lambda_reduction);

    KRATOS_CHECK_NEAR(final_sum, exp_sum, tol);


    // Reduction, short form
    // here the TLS is constructed beforehand. This is the "firstprivate" approach of OpenMP
    // checking the results using reduction
    std::vector<double> tls_short(6);
    const double final_sum_short = block_for_each<SumReduction<double>>(elements, tls_short, tls_lambda_reduction);

    KRATOS_CHECK_NEAR(final_sum_short, exp_sum, tol);
}

KRATOS_TEST_CASE_IN_SUITE(IndexPartitionerThreadLocalStorage, KratosCoreFastSuite)
{
    constexpr std::size_t vec_size = 6;
    constexpr std::size_t n_elems = 10e2;
    constexpr double tol = 1e-8;

    using RHSElementType = RHSElement<vec_size>;

    std::vector<double> rhs_vals(n_elems);
    for (std::size_t i=0; i<n_elems; ++i) {
        rhs_vals[i] = (i%12) * 1.889;
    }

    std::vector<RHSElementType> elements;
    for (std::size_t i=0; i<rhs_vals.size(); ++i) {
        elements.push_back(RHSElementType(rhs_vals[i]));
    }

    const double exp_sum = (std::accumulate(rhs_vals.begin(), rhs_vals.end(), 0.0)) * vec_size;

    auto tls_lambda_manual_reduction = [&elements](std::size_t i, std::vector<double>& rTLS)
    {
        elements[i].CalculateRHS(rTLS);
        double rhs_sum = std::accumulate(rTLS.begin(), rTLS.end(), 0.0);
        elements[i].SetAccumRHSValue(rhs_sum);
    };

    auto tls_lambda_reduction = [&elements](std::size_t i, std::vector<double>& rTLS)
    {
        elements[i].CalculateRHS(rTLS);
        return std::accumulate(rTLS.begin(), rTLS.end(), 0.0);
    };


    // Manual Reduction, long form
    // here the TLS is constructed on the fly. This is the "private" approach of OpenMP
    // the result is checked with a "manual reduction"
    IndexPartition<std::size_t>(elements.size()).for_each(std::vector<double>(), tls_lambda_manual_reduction);

    const double sum_elem_rhs_vals = std::accumulate(elements.begin(), elements.end(), 0.0, [](double acc, RHSElementType& rElem){
        return acc + rElem.GetAccumRHSValue();
    });

    KRATOS_CHECK_NEAR(sum_elem_rhs_vals, exp_sum, tol);

    const double sum_elem_rhs_vals_short = std::accumulate(elements.begin(), elements.end(), 0.0, [](double acc, RHSElementType& rElem){
        return acc + rElem.GetAccumRHSValue();
    });

    KRATOS_CHECK_NEAR(sum_elem_rhs_vals_short, exp_sum, tol);


    // Reduction, long form
    // here the TLS is constructed beforehand. This is the "firstprivate" approach of OpenMP
    // checking the results using reduction
    std::vector<double> tls(6);
    const double final_sum = IndexPartition<std::size_t>(elements.size()).for_each<SumReduction<double>>(tls, tls_lambda_reduction);

    KRATOS_CHECK_NEAR(final_sum, exp_sum, tol);
}

KRATOS_TEST_CASE_IN_SUITE(ParallelUtilsContinue, KratosCoreFastSuite)
{
    int nsize = 1e3;
    std::vector<double> data_vector(nsize, -1.0), output(nsize, 3.3);

    IndexPartition<unsigned int>(data_vector.size()).for_each(
        [&](unsigned int i){
            if (i%4 == 0) return; // this is the equivalent of "continue" in a regular loop

            output[i] = 2.0*data_vector[i];
            }
        );

    for(unsigned int i=0; i<output.size(); ++i) {
        if (i%4 == 0) {
            KRATOS_CHECK_DOUBLE_EQUAL(output[i], 3.3);
        } else {
            KRATOS_CHECK_DOUBLE_EQUAL(output[i], -2.0);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(AccumReductionVector, KratosCoreFastSuite)
{
    int nsize = 1e3;
    std::vector<int> input_data_vector(nsize);
    std::vector<int> expct_data_vector(nsize);

    std::iota(input_data_vector.begin(), input_data_vector.end(), 0);
    std::iota(expct_data_vector.begin(), expct_data_vector.end(), 1);

    auto assembled_vector = block_for_each<AccumReduction<int>>(input_data_vector, [](int& rValue) {
        return rValue+1;
    });

    std::sort(assembled_vector.begin(), assembled_vector.end());

    KRATOS_CHECK_VECTOR_EQUAL(assembled_vector, expct_data_vector);
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
            typedef double value_type;
            typedef std::tuple<double,double> return_type;

            value_type max_value = -std::numeric_limits<double>::max();
            value_type max_abs = 0.0;

            return_type GetValue()
            {
                return_type values;
                std::get<0>(values) = max_value;
                std::get<1>(values) = max_abs;
                return values;
            }

            void LocalReduce(value_type function_return_value){
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
    const auto timer_omp = BuiltinTimer();
    IndexPartition<unsigned int>(nsize).for_each(
            [&](unsigned int i){
                    output[i] = std::pow(data_vector[i],0.01);
                }
            );
    std::cout << "OMP time = " << timer_omp.ElapsedSeconds() << std::endl;

    const auto timer_pure = BuiltinTimer();
    IndexPartition<unsigned int>(nsize).for_pure_c11(
            [&](unsigned int i){
                    output[i] = std::pow(data_vector[i],0.01);
                }
            );
    std::cout << "Pure c++11 time = " << timer_pure.ElapsedSeconds() << std::endl;



}


} // namespace Testing
} // namespace Kratos

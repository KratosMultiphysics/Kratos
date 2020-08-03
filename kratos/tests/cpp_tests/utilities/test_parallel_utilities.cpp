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

#include <utility>
#include <numeric>
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

KRATOS_TEST_CASE_IN_SUITE(IndexPartitionerThreadLocalStorage, KratosCoreFastSuite)
{
    // In this example we use a non thread safe search structure
    //In order to use in a parallel context we need to create a search structure
    //for each thread. This is done with the use of TLS
    struct LocalSearch
    {
      public:
        explicit LocalSearch(const std::vector<double>& rPoints) : mPoints(rPoints) {}

        // the TLS is created using the copy constructor
        // hence we implement the initialization of the
        // search structure in the copy constructor
        LocalSearch(const LocalSearch& rSearchPrototype) : mPoints(rSearchPrototype.mPoints)
        {
            // here we initalize the search structure
            // this is a potentially expensive operation!
            // e.g.
            // mBins = BinsDynamicSearch(...)
            // as this is a very basic example we don't do anything
        }

        double GetIndexOfClosestPoint(const double Coord)
        {
            // alternatively the initialization of the search structure
            // could be done the first time this function is called

            std::size_t closest_index = 0;
            double closest_distance = std::abs(mPoints[0] - Coord);

            for (std::size_t i=1; i<mPoints.size(); ++i) {
                double cur_dist = std::abs(mPoints[i] - Coord);
                if (std::abs(mPoints[i] - Coord) < closest_distance) {
                    closest_index = i;
                    closest_distance = cur_dist;
                }
            }
            return closest_index;
        }

      private:
        std::vector<double> mPoints;
    };

    std::vector<double> point_coords {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

    constexpr std::size_t n_seach_points = 10e2;

    std::vector<double> search_point_coords (n_seach_points);
    std::vector<std::size_t> exp_results (n_seach_points);

    for (std::size_t i=0; i<n_seach_points; ++i) {
        search_point_coords[i] = i%6 + 0.2;
        exp_results[i] = i%6;
    }

    std::vector<std::size_t> results(search_point_coords.size());

    IndexPartition<unsigned int>(search_point_coords.size()).for_each(LocalSearch(point_coords),
        [&search_point_coords, &results](unsigned int i, LocalSearch& rLocalSearch){
            results[i] = rLocalSearch.GetIndexOfClosestPoint(search_point_coords[i]);
        }
    );

    KRATOS_CHECK_EQUAL(results.size(), exp_results.size());
    for (std::size_t i=0; i<results.size(); ++i) {
        KRATOS_CHECK_EQUAL(results[i], exp_results[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(BlockPartitionerThreadLocalStorage, KratosCoreFastSuite)
{
    constexpr std::size_t vec_size = 6;
    constexpr std::size_t n_elems = 10e2;
    constexpr double tol = 1e-8;

    class RHSElement
    {
      public:
        explicit RHSElement(const double Val) : mRHSVal(Val) {}
        void CalculateRHS(std::vector<double>& rVector)
        {
            if (rVector.size() != vec_size) { rVector.resize(vec_size); }
            std::fill(rVector.begin(), rVector.end(), mRHSVal);
        }
        double GetAccumRHSValue() {return mAccumRHSValue;}
        void SetAccumRHSValue(double Value) {mAccumRHSValue = Value;}

      private:
        double mRHSVal;
        double mAccumRHSValue = 0.0;
    };

    std::vector<double> rhs_vals(n_elems);
    for (std::size_t i=0; i<n_elems; ++i) {
        rhs_vals[i] = (i%12) * 1.889;
    }

    std::vector<RHSElement> elements;
    for (std::size_t i=0; i<rhs_vals.size(); ++i) {
        elements.push_back(RHSElement(rhs_vals[i]));
    }

    const double exp_sum = (std::accumulate(rhs_vals.begin(), rhs_vals.end(), 0.0)) * vec_size;

    auto tls_lambda_manual_reduction = [](RHSElement& rElem, std::vector<double>& rTLS)
    {
        rElem.CalculateRHS(rTLS);
        double rhs_sum = std::accumulate(rTLS.begin(), rTLS.end(), 0.0);
        rElem.SetAccumRHSValue(rhs_sum);
    };

    auto tls_lambda_reduction = [](RHSElement& rElem, std::vector<double>& rTLS)
    {
        rElem.CalculateRHS(rTLS);
        return std::accumulate(rTLS.begin(), rTLS.end(), 0.0);
    };


    // Manual Reduction, long form
    // here the TLS is constructed on the fly. This is the "private" approach of OpenMP
    // the result is checked with a "manual reduction"
    BlockPartition<std::vector<RHSElement>>(elements).for_each(std::vector<double>(), tls_lambda_manual_reduction);

    const double sum_elem_rhs_vals = std::accumulate(elements.begin(), elements.end(), 0.0, [](double acc, RHSElement& rElem){
        return acc + rElem.GetAccumRHSValue();
    });

    KRATOS_CHECK_NEAR(sum_elem_rhs_vals, exp_sum, tol);


    // Manual Reduction, short form
    // here the TLS is constructed on the fly. This is the "private" approach of OpenMP
    // the result is checked with a "manual reduction"
    block_for_each(elements, std::vector<double>(), tls_lambda_manual_reduction);

    const double sum_elem_rhs_vals_short = std::accumulate(elements.begin(), elements.end(), 0.0, [](double acc, RHSElement& rElem){
        return acc + rElem.GetAccumRHSValue();
    });

    KRATOS_CHECK_NEAR(sum_elem_rhs_vals_short, exp_sum, tol);


    // Reduction, long form
    // here the TLS is constructed beforehand. This is the "firstprivate" approach of OpenMP
    // checking the results using reduction
    std::vector<double> tls(6);
    const double final_sum = BlockPartition<std::vector<RHSElement>>(elements).for_each<SumReduction<double>>(tls, tls_lambda_reduction);

    KRATOS_CHECK_NEAR(final_sum, exp_sum, tol);


    // Reduction, short form
    // here the TLS is constructed beforehand. This is the "firstprivate" approach of OpenMP
    // checking the results using reduction
    std::vector<double> tls_short(6);
    const double final_sum_short = block_for_each<SumReduction<double>>(elements, tls_short, tls_lambda_reduction);

    KRATOS_CHECK_NEAR(final_sum_short, exp_sum, tol);
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

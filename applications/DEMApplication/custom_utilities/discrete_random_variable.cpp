// Authors:
// Guillermo Casas gcasas@cimne-upc.edu

#include "discrete_random_variable.h"
#include "includes/checks.h"
#include <numeric>

namespace Kratos {
    DiscreteRandomVariable::DiscreteRandomVariable():
    mRandomNumberGenerator(std::random_device{}()){}

    DiscreteRandomVariable::DiscreteRandomVariable(const Parameters rParameters):
    DiscreteRandomVariable(rParameters, std::random_device{}()){};

    DiscreteRandomVariable::DiscreteRandomVariable(const Parameters rParameters, const int seed)
    : mRandomNumberGenerator(seed)
    {
        Parameters default_parameters = Parameters(R"(
        {
            "possible_values"              : [0, 1],
            "relative_frequencies"         : [1, 1],
            "do_use_seed"                  : false,
            "seed"                         : 1,
            "relative_closeness_tolerance" : 1e-6
        })" );

        default_parameters.ValidateAndAssignDefaults(rParameters);

        if (!rParameters.Has("relative_closeness_tolerance")){
            mRelativeClosenessTolerance = 1.0e-6;
        }
        else {
            mRelativeClosenessTolerance = rParameters["relative_closeness_tolerance"].GetDouble();
        }

        const auto possible_values = rParameters["possible_values"].GetVector();
        const auto values = rParameters["relative_frequencies"].GetVector();

        if (possible_values.size() != values.size()){
            KRATOS_ERROR << "Check failed because the possible_values list has a different size than the relative_frequencies list.";
        }

        const std::size_t n_points = possible_values.size();
        mPossibleValues.resize(n_points);
        mRelativeFrequencies.resize(n_points);

        for (std::size_t i = 0; i < n_points; ++i) {
            mPossibleValues[i] = possible_values[i];
            mRelativeFrequencies[i] = values[i];
        }

        mTrapezoidsDiscreteDistribution = std::discrete_distribution<int>(mRelativeFrequencies.begin(), mRelativeFrequencies.end());

        size_t low_index;
        size_t high_index;
        CalculateFirstAndLastIndicesWithNonzeroValue<double>(mRelativeFrequencies, low_index, high_index);
        SetSupport(mPossibleValues[low_index], mPossibleValues[high_index]);

        Check();

        Normalize();
    }

   void DiscreteRandomVariable::Check()
   {
       for (std::size_t i = 0; i < mRelativeFrequencies.size(); ++i){
           if (mRelativeFrequencies[i] < 0.0){
               KRATOS_ERROR << "Check failed because the " << i << "-th entry of the relative_frequencies list is negative.";
           }
       }

       const double support = mPossibleValues[mPossibleValues.size() - 1] - mPossibleValues[0];

       for (std::size_t i = 0; i < mPossibleValues.size() - 1; ++i){
           KRATOS_CHECK_GREATER(mPossibleValues[i+1], mPossibleValues[i]);
           if (std::abs(mPossibleValues[i+1] - mPossibleValues[i]) < mRelativeClosenessTolerance * support){
               KRATOS_ERROR << "Check failed because the " << i << "-th" << " and " << i+1 << "-th consecutive possible_values are too close (" \
                            << mPossibleValues[i] << " ~ " << mPossibleValues[i+1] << "). Increase the separation between them or pick a " \
                            << "smaller tolerance.";
           }
       }
   }

   double DiscreteRandomVariable::ProbabilityDensity(const double x)
   {
       const double& x_min = mPossibleValues[0];
       const double& x_max = mPossibleValues[mPossibleValues.size() - 1];

       if (x < x_min || x > x_max){
           return 0.0;
       }

       for (std::size_t i = 0; i < mPossibleValues.size() - 1; ++i){
           const double& x2 = mPossibleValues[i + 1];
           if (x <= x2 + mRelativeClosenessTolerance && x > x2 - mRelativeClosenessTolerance){
               return mRelativeFrequencies[i];
           }
       }

       return 0.0;
   }

    double DiscreteRandomVariable::Sample(){
        return mPossibleValues[mTrapezoidsDiscreteDistribution(mRandomNumberGenerator)];
    }

    void DiscreteRandomVariable::Normalize(){
        const double total_area = std::accumulate(mRelativeFrequencies.begin(), mRelativeFrequencies.end(), 0);
        // Normalization of probability function
        for (std::size_t i = 0; i < mRelativeFrequencies.size(); ++i) {
            mRelativeFrequencies[i] /= total_area;
        }
    }

    double DiscreteRandomVariable::GetMean(){
        if (!mMeanHasAlreadyBeenCalculated){
            const auto n = mRelativeFrequencies.size();
            mMean = std::accumulate(mRelativeFrequencies.begin(), mRelativeFrequencies.end(), 0.0) / n;
            mMeanHasAlreadyBeenCalculated = true;
        }

        return mMean;
    }

}

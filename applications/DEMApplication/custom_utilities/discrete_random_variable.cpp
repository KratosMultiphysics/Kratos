// Authors:
// Guillermo Casas gcasas@cimne-upc.edu

#include "discrete_random_variable.h"
#include "includes/checks.h"

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
            "possible_values" : [0, 1],
            "pdf_values"      : [1, 1],
            "do_use_seed"     : false,
            "seed"            : 1,
            "relative_closeness_tolerance" : 1e-6
        })" );

        default_parameters.ValidateAndAssignDefaults(rParameters);

        if (!rParameters.Has("relative_closeness_tolerance")){
            mRelativeClosenessTolerance = 1.0e-6;
        }
        else {
            mRelativeClosenessTolerance = rParameters["relative_closeness_tolerance"].GetDouble();
        }

        const auto breakpoints = rParameters["possible_values"].GetVector();
        const auto values = rParameters["pdf_values"].GetVector();

        if (breakpoints.size() != values.size()){
            KRATOS_ERROR << "Check failed because the breakpoints list has a different size than the pdf_values list.";
        }

        const std::size_t n_points = breakpoints.size();
        mPDFBreakpoints.resize(n_points);
        mPDFValues.resize(n_points);

        for (std::size_t i = 0; i < n_points; ++i) {
            mPDFBreakpoints[i] = breakpoints[i];
            mPDFValues[i] = values[i];
        }

        mTrapezoidsDiscreteDistribution = std::discrete_distribution<int>(mPDFValues.begin(), mPDFValues.end());

        Check();

        Normalize();
    }

   void DiscreteRandomVariable::Check()
   {
       for (std::size_t i = 0; i < mPDFValues.size(); ++i){
           if (mPDFValues[i] < 0.0){
               KRATOS_ERROR << "Check failed because the " << i << "-th entry of the pdf_values list is negative.";
           }
       }

       const double support = mPDFBreakpoints[mPDFBreakpoints.size() - 1] - mPDFBreakpoints[0];

       for (std::size_t i = 0; i < mPDFBreakpoints.size() - 1; ++i){
           KRATOS_CHECK_GREATER(mPDFBreakpoints[i+1], mPDFBreakpoints[i]);
           if (std::abs(mPDFBreakpoints[i+1] - mPDFBreakpoints[i]) < mRelativeClosenessTolerance * support){
               KRATOS_ERROR << "Check failed because the " << i << "-th" << " and " << i+1 << "-th consecutive breakpoints are too close (" \
                            << mPDFBreakpoints[i] << " ~ " << mPDFBreakpoints[i+1] << "). Increase the separation between them or pick a " \
                            << "smaller tolerance.";
           }
       }
   }

   double DiscreteRandomVariable::ProbabilityDensity(const double x)
   {
       const double& x_min = mPDFBreakpoints[0];
       const double& x_max = mPDFBreakpoints[mPDFBreakpoints.size() - 1];

       if (x < x_min || x > x_max){
           return 0.0;
       }

       for (std::size_t i = 0; i < mPDFBreakpoints.size() - 1; ++i){
           const double& x2 = mPDFBreakpoints[i + 1];
           if (x <= x2 + mRelativeClosenessTolerance && x > x2 - mRelativeClosenessTolerance){
               return mPDFValues[i];
           }
       }

       return 0.0;
   }

    double DiscreteRandomVariable::Sample(){
        return mPDFBreakpoints[mTrapezoidsDiscreteDistribution(mRandomNumberGenerator)];
    }

    void DiscreteRandomVariable::Normalize(){
        const double total_area = std::accumulate(mPDFValues.begin(), mPDFValues.end(), 0);
        // Normalization of probability function
        for (std::size_t i = 0; i < mPDFValues.size(); ++i) {
            mPDFValues[i] /= total_area;
        }
    }

    double DiscreteRandomVariable::GetMean(){
        if (!mMeanHasAlreadyBeenCalculated){
            const auto n = mPDFValues.size();
            mMean = std::accumulate(mPDFValues.begin(), mPDFValues.end(), 0.0) / n;
            mMeanHasAlreadyBeenCalculated = true;
        }

        return mMean;
    }

}

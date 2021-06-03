// Authors:
// Guillermo Casas gcasas@cimne-upc.edu

#include "piecewise_linear_random_variable.h"
#include "includes/checks.h"

namespace Kratos {
    PiecewiseLinearRandomVariable::PiecewiseLinearRandomVariable():
    mRandomNumberGenerator(std::random_device{}()){}

    PiecewiseLinearRandomVariable::PiecewiseLinearRandomVariable(const Parameters rParameters):
    PiecewiseLinearRandomVariable(rParameters, std::random_device{}()){};

    PiecewiseLinearRandomVariable::PiecewiseLinearRandomVariable(const Parameters rParameters, const int seed)
    : mRandomNumberGenerator(seed)
    {
        Parameters default_parameters = Parameters(R"(
        {
            "pdf_breakpoints" : [0, 1],
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

        const auto breakpoints = rParameters["pdf_breakpoints"].GetVector();
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

        Check();

        CalculateTrapezoidProbabilitiesAndNormalize();
    }

   void PiecewiseLinearRandomVariable::Check()
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
   };

   double PiecewiseLinearRandomVariable::ProbabilityDensity(const double x)
   {
       const double& x_min = mPDFBreakpoints[0];
       const double& x_max = mPDFBreakpoints[mPDFBreakpoints.size() - 1];

       if (x < x_min || x > x_max){
           return 0.0;
       }

       for (std::size_t i = 0; i < mPDFBreakpoints.size() - 1; ++i){
           const double& x2 = mPDFBreakpoints[i + 1];
           if (x <= x2){
               const double& x1 = mPDFBreakpoints[i];
               const double alpha =  (x - x1) / (x2 - x1);
               const double value1 = mPDFValues[i];
               const double value2 = mPDFValues[i + 1];
               return (1 - alpha) * value1 + alpha * value2;
           }
       }

       return 0.0;
   }

    double PiecewiseLinearRandomVariable::Sample(){
        const auto i_bin = SampleTrapezoidChoice();
        const double x0 = mPDFBreakpoints[i_bin];
        const double H = mPDFBreakpoints[i_bin + 1] - x0;
        const double B1 = mPDFValues[i_bin];
        const double B2 = mPDFValues[i_bin + 1];
        const double x_within = SampleWithinTrapezoid(H, B1, B2);
        return x0 + x_within;
    }

    void PiecewiseLinearRandomVariable::CalculateTrapezoidProbabilitiesAndNormalize(){
        double total_area = 0.0;
        std::vector<double> areas(mPDFBreakpoints.size() - 1);

        // Area under each straight bin
        for (std::size_t i = 0; i < areas.size(); ++i) {
            const double trapezoid_area = 0.5 * (mPDFBreakpoints[i + 1] - mPDFBreakpoints[i]) * (mPDFValues[i + 1] + mPDFValues[i]);
            total_area += trapezoid_area;
            areas[i] = trapezoid_area;
        }

        // Normalization to obtain probability of hitting each particular straight bin
        for (std::size_t i = 0; i < areas.size(); ++i) {
            areas[i] /= total_area;
        }

        // Normalization of probability function
        for (std::size_t i = 0; i < mPDFValues.size(); ++i) {
            mPDFValues[i] /= total_area;
        }

        mTrapezoidsDiscreteDistribution = std::discrete_distribution<int>(areas.begin(), areas.end());
    }

    std::size_t PiecewiseLinearRandomVariable::SampleTrapezoidChoice(){
        return mTrapezoidsDiscreteDistribution(mRandomNumberGenerator);
    }

    double PiecewiseLinearRandomVariable::GetMean(){

        if (!mMeanHasAlreadyBeenCalculated){
            std::vector<double> areas(mPDFBreakpoints.size() - 1);

            // Compute the mean by weighing the centroid coordinate of every trapezoid by its area
            mMean = 0.0;
            for (std::size_t i = 0; i < areas.size(); ++i) {
                const double& x1 = mPDFBreakpoints[i];
                const double& x2 = mPDFBreakpoints[i + 1];
                const double& v1 = mPDFValues[i];
                const double& v2 = mPDFValues[i+1];
                const double h = x2 - x1;
                const double trapezoid_area = 0.5 * (v1 + v2) * h;
                const double square_area = std::min(v1, v2) * h;
                const double delta_v = v2 - v1;
                const double triangle_area = 0.5 * std::abs(delta_v)  * h;
                const int sign_delta_v = (delta_v > 0) - (delta_v < 0);
                const double x_triangle = (0.5 + 1.0/6 * sign_delta_v) * h; // covers triangles sloping up or down
                const double x_square = 0.5 * h;
                const double x_centroid_within_trapezoid = (triangle_area * x_triangle + square_area * x_square) / trapezoid_area;
                const double x_centroid = x1 + x_centroid_within_trapezoid;
                mMean += trapezoid_area * x_centroid;
            }

            mMeanHasAlreadyBeenCalculated = true;
        }

        return mMean;
    }

    double PiecewiseLinearRandomVariable::SampleWithinTrapezoid(const double H, const double B1, const double B2){
        double x;
        if (B1 == 0){ //TODO: improve
            x = SamplePositiveSlopingStandardTriangle();
        }

        else {
            const double beta = B2/B1;
            const double b = 2.0 / (1 + beta);
            x = SampleWithinStandardTrapezoid(b);
        }
        return H * x;
    }

    double PiecewiseLinearRandomVariable::SampleWithinStandardTrapezoid(const double b){
        std::uniform_real_distribution<> uniform_distribution(0, 1);
        const double alpha = uniform_distribution(mRandomNumberGenerator);
        double x;
        if (alpha < b/2){
            x = SampleNegativeSlopingStandardTriangle();
        }

        else {
            x = SamplePositiveSlopingStandardTriangle();
        }

        return x;
    }

    double PiecewiseLinearRandomVariable::SamplePositiveSlopingStandardTriangle(){
        std::uniform_real_distribution<> uniform_distribution(0, 1);
        const double x = uniform_distribution(mRandomNumberGenerator);
        return std::sqrt(x);
    }

    double PiecewiseLinearRandomVariable::SampleNegativeSlopingStandardTriangle(){
        return 1.0 - SamplePositiveSlopingStandardTriangle();
    }

}

#include "piecewise_linear_random_variable.h"

namespace Kratos {
    PiecewiseLinearRandomVariable::PiecewiseLinearRandomVariable():
    mRandomNumberGenerator(std::random_device{}()){}

    PiecewiseLinearRandomVariable::PiecewiseLinearRandomVariable(const Parameters rParameters):
    PiecewiseLinearRandomVariable(rParameters, std::random_device{}()){};

    PiecewiseLinearRandomVariable::PiecewiseLinearRandomVariable(const Parameters rParameters, const int seed)
    : mRandomNumberGenerator(seed){

        const auto breakpoints = rParameters["pdf_breakpoints"].GetVector();
        const auto values = rParameters["pdf_values"].GetVector();
        const std::size_t n_points = breakpoints.size();
        mPDFBreakpoints.resize(n_points);
        mPDFValues.resize(n_points);

        for (std::size_t i = 0; i < n_points; ++i) {
            mPDFBreakpoints[i] = breakpoints[i];
            mPDFValues[i] = values[i];
        }

        CalculateTrapezoidProbabilitiesAndNormalize();
    }


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

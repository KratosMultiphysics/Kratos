#include "piecewise_linear_random_variable.h"

namespace Kratos {
    PiecewiseLinearRandomVariable::PiecewiseLinearRandomVariable(){
    }
    PiecewiseLinearRandomVariable::PiecewiseLinearRandomVariable(const Parameters rParameters){

        const auto breakpoints = rParameters["pdf_breakpoints"].GetVector();
        const auto values = rParameters["pdf_values"].GetVector();
        const std::size_t n_points = breakpoints.size();
        mPDFRange.resize(n_points);
        mPDFValues.resize(n_points);

        for (std::size_t i = 0; i < n_points; ++i) {
            mPDFRange[i] = breakpoints[i];
            mPDFValues[i] = values[i];
        }

        CalculateTrapezoidProbabilitiesAndNormalize();
    }

    double PiecewiseLinearRandomVariable::Sample(){
        const auto i_bin = SampleTrapezoidChoice();
        const double x0 = mPDFRange[i_bin];
        const double H = mPDFRange[i_bin + 1] - x0;
        const double B1 = mPDFValues[i_bin];
        const double B2 = mPDFValues[i_bin + 1];
        const double x_within = SampleWithinTrapezoid(H, B1, B2);
        return x0 + x_within;
    }

    void PiecewiseLinearRandomVariable::CalculateTrapezoidProbabilitiesAndNormalize(){
        double total_area = 0.0;
        std::vector<double> areas(mPDFRange.size() - 1);

        // Area under each straight bin
        for (std::size_t i = 0; i < areas.size(); ++i) {
            const double trapezoid_area = 0.5 * (mPDFRange[i + 1] - mPDFRange[i]) * (mPDFValues[i + 1] + mPDFValues[i]);
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
        std::random_device rd;
        std::mt19937 gen(rd());
        return mTrapezoidsDiscreteDistribution(gen);
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
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> uniform_distribution(0, 1);
        const double alpha = uniform_distribution(gen);
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
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> uniform_distribution(0, 1);
        const double x = uniform_distribution(gen);
        return std::sqrt(x);
    }

    double PiecewiseLinearRandomVariable::SampleNegativeSlopingStandardTriangle(){
        return 1.0 - SamplePositiveSlopingStandardTriangle();
    }

}

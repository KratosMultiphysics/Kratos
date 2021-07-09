//
// Authors:
// Guillermo Casas gcasas@cimne-upc.edu
//

#ifndef PIECEWISE_LINEAR_RANDOM_VARIABLE_H
#define PIECEWISE_LINEAR_RANDOM_VARIABLE_H


// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "DEM_application_variables.h"
#include "includes/model_part.h"
#include "random_variable.h"

namespace Kratos {

class KRATOS_API(DEM_APPLICATION) PiecewiseLinearRandomVariable: public RandomVariable {

public:

    KRATOS_CLASS_POINTER_DEFINITION(PiecewiseLinearRandomVariable);

    /// Default constructor
    PiecewiseLinearRandomVariable();
    PiecewiseLinearRandomVariable(const Parameters rParameters);
    PiecewiseLinearRandomVariable(const Parameters rParameters, const int seed);

    double Sample() override;
    double ProbabilityDensity(const double x);
    double GetMean() override;

    /// Turn back information as a stemplate<class T, std::size_t dim> tring.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "PiecewiseLinearRandomVariable" ;
        return buffer.str();

    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PiecewiseLinearRandomVariable";
    }

    /// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {

    }

protected:
    void Check() override;


private:

    void CalculateTrapezoidProbabilitiesAndNormalize();
    std::size_t SampleTrapezoidChoice();
    double SampleWithinTrapezoid(const double H, const double B1, const double B2);
    double SampleWithinStandardTrapezoid(const double b);
    double SamplePositiveSlopingStandardTriangle();
    double SampleNegativeSlopingStandardTriangle();

    /// Assignment operator.
    PiecewiseLinearRandomVariable & operator=(PiecewiseLinearRandomVariable const& rOther);

    double mRelativeClosenessTolerance = 0.0;
    double mMean = 0.0 ;
    bool mMeanHasAlreadyBeenCalculated=false;
    std::vector<double> mPDFValues;
    std::vector<double> mPDFBreakpoints;
    std::mt19937 mRandomNumberGenerator;
    std::discrete_distribution<int> mTrapezoidsDiscreteDistribution;

}; // Class PiecewiseLinearRandomVariable

} // namespace Kratos

#endif // PIECEWISE_LINEAR_RANDOM_VARIABLE_H defined



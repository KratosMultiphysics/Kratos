//
// Authors:
// Guillermo Casas gcasas@cimne-upc.edu
//

#ifndef DISCRETE_RANDOM_VARIABLE_H
#define DISCRETE_RANDOM_VARIABLE_H


// System includes
#include <string>
#include <numeric>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "DEM_application_variables.h"
#include "includes/model_part.h"
#include "random_variable.h"

namespace Kratos {

class KRATOS_API(DEM_APPLICATION) DiscreteRandomVariable: public RandomVariable {

public:

    KRATOS_CLASS_POINTER_DEFINITION(DiscreteRandomVariable);

    /// Default constructor
    DiscreteRandomVariable();
    DiscreteRandomVariable(const Parameters rParameters);
    DiscreteRandomVariable(const Parameters rParameters, const int seed);

    double Sample() override;
    double ProbabilityDensity(const double x);
    double GetMean() override;

    /// Turn back information as a stemplate<class T, std::size_t dim> tring.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "DiscreteRandomVariable" ;
        return buffer.str();

    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DiscreteRandomVariable";
    }

    /// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {

    }

protected:
    void Check() override;


private:

    void Normalize();

    /// Assignment operator.
    DiscreteRandomVariable & operator=(DiscreteRandomVariable const& rOther);

    double mRelativeClosenessTolerance = 0.0;
    std::vector<double> mRelativeFrequencies;
    std::vector<double> mPossibleValues;
    std::mt19937 mRandomNumberGenerator;
    std::discrete_distribution<> mTrapezoidsDiscreteDistribution;

}; // Class DiscreteRandomVariable

} // namespace Kratos

#endif // DISCRETE_RANDOM_VARIABLE_H defined



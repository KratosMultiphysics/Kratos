//
// Authors:
// Guillermo Casas gcasas@cimne-upc.edu
//

#ifndef RANDOM_VARIABLE_H
#define RANDOM_VARIABLE_H


// System includes
#include <string>
#include <iostream>
#include <random>

// Project includes
#include "includes/define.h"
#include "DEM_application_variables.h"
#include "includes/model_part.h"

namespace Kratos {

class KRATOS_API(DEM_APPLICATION) RandomVariable {

public:

    KRATOS_CLASS_POINTER_DEFINITION(RandomVariable);

    /// Default constructor
    RandomVariable();
    RandomVariable(const Parameters rParameters);

    /// Destructor
    virtual ~RandomVariable(){};

    virtual double Sample(){KRATOS_ERROR << "You are calling the 'Sample' function of the abstract class 'RandomVariable'. Please instantiate a specific derived class instead."; return 0.0;};

    virtual double GetMean(){KRATOS_ERROR << "You are calling 'GetMean' function of the abstract class 'RandomVariable'. Please instantiate a specific derived class instead."; return 0.0;};

    const array_1d<double, 2>& GetSupport();

    /// Turn back information as a stemplate<class T, std::size_t dim> tring.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


protected:
    virtual void Check(){};
    double mMean = 0.0;
    bool mMeanHasAlreadyBeenCalculated=false;

    template<typename T>
    void CalculateFirstAndLastIndicesWithNonzeroValue(std::vector<T> values, size_t& low_index, size_t& high_index){
                // finding first and last indices that correspond to nonzero probabilites

        auto it = std::find_if(values.begin(), values.end(), [](const double x) { return x != 0; });
        auto reverse_it = std::find_if(values.rbegin(), values.rend(), [](const double x) { return x != 0; });
        low_index = std::distance(values.begin(), it);
        high_index = std::distance(begin(values), reverse_it.base()) - 1;
    }

    void SetSupport(const double Min, const double Max);

private:
    array_1d<double, 2> mSupport;

    /// Assignment operator.
    RandomVariable & operator=(RandomVariable const& rOther);

}; // Class RandomVariable

} // namespace Kratos

#endif // RANDOM_VARIABLE_H defined



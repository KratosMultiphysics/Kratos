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

    /// Turn back information as a stemplate<class T, std::size_t dim> tring.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


protected:
    virtual void Check(){};

private:

    /// Assignment operator.
    RandomVariable & operator=(RandomVariable const& rOther);

}; // Class RandomVariable

} // namespace Kratos

#endif // RANDOM_VARIABLE_H defined



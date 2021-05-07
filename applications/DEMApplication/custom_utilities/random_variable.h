//
// Authors:
// Guillermo Casas gcasas@cimne-upc.edu
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#ifndef RANDOM_VARIABLE_H
#define RANDOM_VARIABLE_H


// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "../DEM_application_variables.h"
#include "includes/model_part.h"

namespace Kratos {

class KRATOS_API(DEM_APPLICATION) RandomVariable {

public:

    KRATOS_CLASS_POINTER_DEFINITION(RandomVariable);

    /// Default constructor
    RandomVariable();
    RandomVariable(const Parameters rParameters);

    /// Destructor
    ~RandomVariable(){};

    /// Turn back information as a stemplate<class T, std::size_t dim> tring.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


protected:


private:

    /// Assignment operator.
    RandomVariable & operator=(RandomVariable const& rOther);

}; // Class RandomVariable

} // namespace Kratos

#endif // RANDOM_VARIABLE_H defined



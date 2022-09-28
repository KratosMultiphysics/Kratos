// Authors:
// Guillermo Casas gcasas@cimne-upc.edu

#include "random_variable.h"

namespace Kratos {
    RandomVariable::RandomVariable(){
    }
    RandomVariable::RandomVariable(const Parameters rParameters){
    }

    void RandomVariable::SetSupport(const double Min, const double Max){
        mSupport[0] = Min;
        mSupport[1] = Max;
    }

    const array_1d<double, 2>& RandomVariable::GetSupport(){
        return mSupport;
    }

    std::string RandomVariable::Info() const
    {
        std::stringstream buffer;
        buffer << "RandomVariable" ;
        return buffer.str();
    }
    void RandomVariable::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Abstract RandomVariable";
    }
    void RandomVariable::PrintData(std::ostream& rOStream) const {}
}



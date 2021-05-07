#include "random_variable.h"

namespace Kratos {
    RandomVariable::RandomVariable(){
    }
    RandomVariable::RandomVariable(const Parameters rParameters){
    }

    std::string RandomVariable::Info() const { return "Random Variable Object"; }
    void RandomVariable::PrintInfo(std::ostream& rOStream) const {}
    void RandomVariable::PrintData(std::ostream& rOStream) const {}
}

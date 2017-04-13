//   $Main author: Guillermo Casas
//

// Project includes

// System includes
#include <limits>
#include <iostream>
#include <iomanip>

// External includes
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "analytic_model_part_filler.h"

namespace Kratos
{

/// Turn back information as a string.
std::string AnalyticModelPartFiller::Info() const {
        return "AnalyticModelPartFiller";
}

/// Print information about this object.
void AnalyticModelPartFiller::PrintInfo(std::ostream& rOStream) const {}

/// Print object's data.
void AnalyticModelPartFiller::PrintData(std::ostream& rOStream) const {}

} // namespace Kratos

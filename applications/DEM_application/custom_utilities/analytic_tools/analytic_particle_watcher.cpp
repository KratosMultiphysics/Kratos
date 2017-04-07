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
#include "analytic_particle_watcher.h"

namespace Kratos
{

/// Turn back information as a string.
std::string AnalyticParticleWatcher::Info() const {
        return "";
}

/// Print information about this object.
void AnalyticParticleWatcher::PrintInfo(std::ostream& rOStream) const {}

/// Print object's data.
void AnalyticParticleWatcher::PrintData(std::ostream& rOStream) const {}


} // namespace Kratos

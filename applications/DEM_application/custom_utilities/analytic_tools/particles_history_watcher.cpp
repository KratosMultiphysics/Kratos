//   $Author: Guillermo Casas

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
#include "particles_history_watcher.h"
#include "DEM_application.h"

namespace Kratos
{

void ParticlesHistoryWatcher::ClearData()
{
    mIds.clear();
    mX0s.clear();
    mY0s.clear();
    mZ0s.clear();
    mRadii.clear();
    mTimesOfCreation.clear();
}

void ParticlesHistoryWatcher::MakeMeasurements(ModelPart& analytic_model_part)
{

}

void ParticlesHistoryWatcher::ClearList(boost::python::list& my_list)
{
    while (len(my_list)){
        my_list.pop(); // only way I found to remove all entries
    }
}

void ParticlesHistoryWatcher::GetNewParticlesData(boost::python::list ids,
                                                  boost::python::list X0s,
                                                  boost::python::list Y0s,
                                                  boost::python::list Z0s,
                                                  boost::python::list radii,
                                                  boost::python::list times_of_creation)
{

}



/// Turn back information as a string.
std::string ParticlesHistoryWatcher::Info() const {
        return "ParticlesHistoryWatcher";
}

/// Print information about this object.
void ParticlesHistoryWatcher::PrintInfo(std::ostream& rOStream) const {}

/// Print object's data.
void ParticlesHistoryWatcher::PrintData(std::ostream& rOStream) const {}


} // namespace Kratos

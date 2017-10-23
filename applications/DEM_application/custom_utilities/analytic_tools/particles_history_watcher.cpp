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

void ParticlesHistoryWatcher::Record(SphericParticle* p_particle, ModelPart& r_model_part)
{
    Node<3>& node = p_particle->GetGeometry()[0];
    mIds.push_back(node.Id());
    mX0s.push_back(node.X0());
    mY0s.push_back(node.Y0());
    mZ0s.push_back(node.Z0());
    mRadii.push_back(node.FastGetSolutionStepValue(RADIUS));
    mTimesOfCreation.push_back(r_model_part.GetProcessInfo()[TIME]);
}


void ParticlesHistoryWatcher::GetNewParticlesData(boost::python::list ids,
                                                  boost::python::list X0s,
                                                  boost::python::list Y0s,
                                                  boost::python::list Z0s,
                                                  boost::python::list radii,
                                                  boost::python::list times_of_creation)
{
    ParticlesHistoryWatcher::ClearList(ids);
    ParticlesHistoryWatcher::ClearList(X0s);
    ParticlesHistoryWatcher::ClearList(Y0s);
    ParticlesHistoryWatcher::ClearList(Z0s);
    ParticlesHistoryWatcher::ClearList(radii);
    ParticlesHistoryWatcher::ClearList(times_of_creation);


    const int n_particles = mIds.size();
    for (int i = 0; i < n_particles; ++i){
        ids.append(mIds[i]);
        X0s.append(mX0s[i]);
        Y0s.append(mY0s[i]);
        Z0s.append(mZ0s[i]);
        radii.append(mRadii[i]);
        times_of_creation.append(mTimesOfCreation[i]);
    }

    ClearData();
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

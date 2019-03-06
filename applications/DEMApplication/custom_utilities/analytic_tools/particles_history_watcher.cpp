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


void ParticlesHistoryWatcher::GetNewParticlesData(std::list<int> ids,
                                                  std::list<double> X0s,
                                                  std::list<double> Y0s,
                                                  std::list<double> Z0s,
                                                  std::list<double> radii,
                                                  std::list<double> times_of_creation)
{
    ids.clear();
    X0s.clear();
    Y0s.clear();
    Z0s.clear();
    radii.clear();
    times_of_creation.clear();

    const int n_particles = mIds.size();
    for (int i = 0; i < n_particles; ++i){
        ids.push_back(mIds[i]);
        X0s.push_back(mX0s[i]);
        Y0s.push_back(mY0s[i]);
        Z0s.push_back(mZ0s[i]);
        radii.push_back(mRadii[i]);
        times_of_creation.push_back(mTimesOfCreation[i]);
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

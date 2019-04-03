//   $Author: Guillermo Casas
#ifndef PARTICLES_HISTORY_WATCHER
#define PARTICLES_HISTORY_WATCHER

// System includes

#include <limits>
#include <iostream>
#include <iomanip>
#include <list>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "analytic_watcher.h"

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos
{
class KRATOS_API(DEM_APPLICATION) ParticlesHistoryWatcher: public AnalyticWatcher {

public:

KRATOS_CLASS_POINTER_DEFINITION(ParticlesHistoryWatcher);

/// Default constructor

ParticlesHistoryWatcher(){}

/// Destructor

virtual ~ParticlesHistoryWatcher(){}

void ClearData() override;

void GetNewParticlesData(std::list<int> ids,
                         std::list<double> X0s,
                         std::list<double> Y0s,
                         std::list<double> Z0s,
                         std::list<double> radii,
                         std::list<double> times_of_creation);

void MakeMeasurements(ModelPart& analytic_model_part) override;

void Record(SphericParticle* p_particle, ModelPart& r_model_part) override;

/// Turn back information as a string
std::string Info() const override;

/// Print information about this object
void PrintInfo(std::ostream& rOStream) const override;

/// Print object's data
void PrintData(std::ostream& rOStream) const override;


private:

std::vector<int> mIds;
std::vector<double> mX0s;
std::vector<double> mY0s;
std::vector<double> mZ0s;
std::vector<double> mRadii;
std::vector<double> mTimesOfCreation;

/// Assignment operator
ParticlesHistoryWatcher & operator=(ParticlesHistoryWatcher const& rOther);

}; // Class ParticlesHistoryWatcher

} // namespace Kratos.

#endif // PARTICLES_HISTORY_WATCHER

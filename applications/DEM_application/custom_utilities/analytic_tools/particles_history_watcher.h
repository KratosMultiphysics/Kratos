//   $Author: Guillermo Casas
#ifndef PARTICLES_HISTORY_WATCHER
#define PARTICLES_HISTORY_WATCHER

// System includes

#include <limits>
#include <iostream>
#include <iomanip>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "analytic_watcher.h"

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif
#include "boost/python/list.hpp"

namespace Kratos
{
class ParticlesHistoryWatcher: public AnalyticWatcher {

public:

KRATOS_CLASS_POINTER_DEFINITION(ParticlesHistoryWatcher);

/// Default constructor

ParticlesHistoryWatcher(){}

/// Destructor

virtual ~ParticlesHistoryWatcher(){}

static void ClearList(boost::python::list& my_list); // its best to pass empty lists in the first place to avoid this operation

void ClearData();

void GetNewParticlesData(boost::python::list ids,
                         boost::python::list X0s,
                         boost::python::list Y0s,
                         boost::python::list Z0s,
                         boost::python::list radii,
                         boost::python::list times_of_creation);

void MakeMeasurements(ModelPart& analytic_model_part);

void Record(SphericParticle* p_particle, ModelPart& r_model_part);

/// Turn back information as a string
std::string Info() const;

/// Print information about this object
void PrintInfo(std::ostream& rOStream) const;

/// Print object's data
void PrintData(std::ostream& rOStream) const;


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

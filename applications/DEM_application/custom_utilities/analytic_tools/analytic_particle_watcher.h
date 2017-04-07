#ifndef ANALYTIC_PARTICLE_WATCHER_H
#define ANALYTIC_PARTICLE_WATCHER_H

// System includes

#include <limits>
#include <iostream>
#include <iomanip>

// Project includes

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos
{

class AnalyticParticleWatcher {

KRATOS_CLASS_POINTER_DEFINITION(AnalyticParticleWatcher);

public:

/// Default constructor

AnalyticParticleWatcher() {}

/// Destructor

virtual ~AnalyticParticleWatcher();

/// Turn back information as a string
virtual std::string Info() const;

/// Print information about this object
virtual void PrintInfo(std::ostream& rOStream) const;

/// Print object's data
virtual void PrintData(std::ostream& rOStream) const;

protected:

    vector<unsigned int> mElementPartition;

private:

    array_1d<double, 3> mInitialCenterOfMassAndMass;
    double mInitialMass;

/// Assignment operator
AnalyticParticleWatcher & operator=(AnalyticParticleWatcher const& rOther);

}; // Class AnalyticParticleWatcher

} // namespace Kratos.

#endif // ANALYTIC_PARTICLE_WATCHER_H

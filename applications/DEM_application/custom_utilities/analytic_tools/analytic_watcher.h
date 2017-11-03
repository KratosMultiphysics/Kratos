//   $Author: Guillermo Casas
#ifndef ANALYTIC_WATCHER
#define ANALYTIC_WATCHER

// System includes

#include <limits>
#include <iostream>
#include <iomanip>

// Project includes
#include "includes/define.h"
#include "../../custom_elements/spheric_particle.h"

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif
#include "boost/python/list.hpp"

namespace Kratos
{
class AnalyticWatcher {

public:

KRATOS_CLASS_POINTER_DEFINITION(AnalyticWatcher);

/// Default constructor

AnalyticWatcher(){}

/// Destructor

virtual ~AnalyticWatcher(){}

static void ClearList(boost::python::list& my_list) // its best to pass empty lists in the first place to avoid this operation
{
    while (len(my_list)){
        my_list.pop(); // only way I found to remove all entries
    }
}

virtual void ClearData(){}

virtual void MakeMeasurements(ModelPart& r_model_part){}

/// Turn back information as a string
virtual std::string Info() const {return "AnalyticWatcher";}

/// Print information about this object
virtual void PrintInfo(std::ostream& rOStream) const {}

/// Print object's data
virtual void PrintData(std::ostream& rOStream) const {}

virtual void Record(SphericParticle* p_particle, ModelPart& r_model_part){}

/// Assignment operator
AnalyticWatcher & operator=(AnalyticWatcher const& rOther);

}; // Class AnalyticWatcher
} // namespace Kratos.

#endif // ANALYTIC_WATCHER

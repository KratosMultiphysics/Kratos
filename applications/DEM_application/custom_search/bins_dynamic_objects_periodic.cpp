#include "bins_dynamic_objects_periodic.h"
#include "../custom_utilities/discrete_particle_configure.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos
{


typedef DiscreteParticleConfigure <3> DiscreteParticleConfigure3D;
template class BinsObjectDynamicPeriodic< DiscreteParticleConfigure3D >;
}  // namespace Kratos.

//Explicit Instantiation


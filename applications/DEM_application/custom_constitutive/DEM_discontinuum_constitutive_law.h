
#if !defined(DEM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_CONSTITUTIVE_LAW_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


/* Project includes */
#include "includes/define.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "custom_utilities/GeometryFunctions.h"
#include "../custom_elements/discrete_element.h"
#include "DEM_application.h"
#include "../custom_elements/Particle_Contact_Element.h"
#include "containers/vector_component_adaptor.h"
#include "containers/array_1d.h"


namespace Kratos
{

/**
 * Base class of constitutive laws.
 */
class DEMDiscontinuumConstitutiveLaw
{

public:

  virtual void CalculateContactForces(double LocalElasticContactForce[3],double indentation,SphericParticle *neighbour_iterator);

protected:

    vector< array_1d<double, 4> > mHistory;

    vector< array_1d<double, 4> > DEMConstitutiveLaw::GetmHistory() { return mHistory;  }

};

} /* namespace Kratos.*/
#endif /* DEM_CONSTITUTIVE_LAW_H_INCLUDED  defined */


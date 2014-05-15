
#if !defined(DEM_CONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_CONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED

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
class DEMContinuumConstitutiveLaw
{

public:

  virtual void CalculateContactForces (double mRadius,double mSqrtOfRealMass, double other_radius,double otherSqrtMass,double distance,double initial_delta,ProcessInfo& rCurrentProcessInfo,PropertiesProxy *myProperties,PropertiesProxy *neighbourProperties,int mapping_new_ini,int mapping_new_cont,int &mNeighbourFailureId_i,double LocalElasticContactForce[3]);

protected:

    vector< array_1d<double, 4> > mHistory;

    vector< array_1d<double, 4> > DEMContinuumConstitutiveLaw::GetmHistory() { return mHistory;  }


};

} /* namespace Kratos.*/
#endif /* DEM_CONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED  defined */


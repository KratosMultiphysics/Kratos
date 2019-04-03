/*
 * File:   particle_data.h
 * Author: gcasas
 *
 * Created on 28 de junio de 2011, 15:31
 */

#if !defined(KRATOS_PARTICLE_DATA_H)
#define  KRATOS_PARTICLE_DATA_H

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>
// External includes


// Project includes
#include "includes/define.h"
#include "includes/define.h"
#include "includes/model_part.h"

//Database includes
#include "spatial_containers/spatial_containers.h"
#include "containers/vector_map.h"
#include "containers/pointer_vector_set.h"
#include "containers/variables_list_data_value_container.h"
#include "utilities/indexed_object.h"
#include "containers/weak_pointer_vector.h"

namespace Kratos
{

class Particle_Data
{
const double mRadius;
const double mStiffness;
const double mDensity;
const double mPi;
const double mProximity_Tol;

Particle_Data(double radius, double stiffness, double density, double prox_tol)
    : mRadius(radius),
      mStiffness(stiffness),
      mDensity(density),
      mPi(3.141592653589793238462643383279),
      mProximity_Tol(prox_tol){};

virtual ~Particle_Data(){};

};
}  // namespace Kratos.

#endif // KRATOS_PARTICLE_DATA_H  defined 


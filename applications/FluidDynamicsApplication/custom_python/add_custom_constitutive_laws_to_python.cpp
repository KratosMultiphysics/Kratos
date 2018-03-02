//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Ruben Zorrilla
//


// System includes
#include <boost/python.hpp>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"

//Application includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"

// 3D constitutive laws
#include "custom_constitutive/euler_3d_law.h"
#include "custom_constitutive/bingham_3d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_constitutive/herschel_bulkey_3d_law.h"

// 2D constitutive laws
#include "custom_constitutive/euler_2d_law.h"
#include "custom_constitutive/newtonian_2d_law.h"

namespace Kratos
{

namespace Python
{

using namespace boost::python;


void  AddCustomConstitutiveLawsToPython()
{

    class_< Euler2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "Euler2DLaw",  init<>() );

    class_< Euler3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "Euler3DLaw",  init<>() );
    
    class_< Bingham3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "Bingham3DLaw",init<>() );

    class_< Newtonian2DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "Newtonian2DLaw",  init<>() );
     
    class_< Newtonian3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "Newtonian3DLaw",  init<>() );
  
    class_< HerschelBulkey3DLaw, bases< ConstitutiveLaw >, boost::noncopyable >( "HerschelBulkey3DLaw",  init<>() );
    
}

}  // namespace Python.
}  // namespace Kratos.

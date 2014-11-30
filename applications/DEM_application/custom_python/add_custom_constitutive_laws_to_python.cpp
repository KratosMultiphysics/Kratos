//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            MSantasusana $
//   Date:                $Date:                Maig 2014 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

 // Project includes
#include "includes/define.h" 

#include "custom_constitutive/DEM_discontinuum_constitutive_law.h"    
#include "custom_constitutive/DEM_continuum_constitutive_law.h" 
#include "custom_constitutive/DEM_Dempack1_CL.h"



namespace Kratos
{

namespace Python
{

using namespace boost::python;


void  AddCustomConstitutiveLawsToPython()
{
    

    //ConstitutiveLaw

    class_< DEMDiscontinuumConstitutiveLaw, DEMDiscontinuumConstitutiveLaw::Pointer, boost::noncopyable >  //bases< ConstitutiveLawBaseType >
    ( "DEMDiscontinuumConstitutiveLaw",
      init<>() )
    .def("Clone",&DEMDiscontinuumConstitutiveLaw::Clone)
    .def("SetConstitutiveLawInProperties",&DEMDiscontinuumConstitutiveLaw::SetConstitutiveLawInProperties)
    ;
   
    class_<Variable<DEMDiscontinuumConstitutiveLaw::Pointer>, boost::noncopyable >( "DEMDiscontinuumConstitutiveLawPointerVariable", no_init )
    .def( self_ns::str( self ) )
    ;
    
    
    
    class_< DEMContinuumConstitutiveLaw, DEMContinuumConstitutiveLaw::Pointer, boost::noncopyable >  //bases< ConstitutiveLawBaseType >
    ( "DEMContinuumConstitutiveLaw",
      init<>() )
    .def("Clone",&DEMContinuumConstitutiveLaw::Clone)
    .def("SetConstitutiveLawInProperties",&DEMContinuumConstitutiveLaw::SetConstitutiveLawInProperties)
    ;
   
    class_<Variable<DEMContinuumConstitutiveLaw::Pointer>, boost::noncopyable >( "DEMContinuumConstitutiveLawPointerVariable", no_init )
    .def( self_ns::str( self ) )
    ;

    class_<DEM_Dempack1, bases< DEMContinuumConstitutiveLaw >, boost::noncopyable >
    ( "DEM_Dempack1",
      init<>() )
    ;




}

}  // namespace Python.
}  // namespace Kratos.

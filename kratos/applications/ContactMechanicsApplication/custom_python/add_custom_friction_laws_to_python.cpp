//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                August 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes
#include <boost/python.hpp>

// External includes

// Project includes

// Application includes
#include "custom_python/add_custom_friction_laws_to_python.h"

//friction laws
#include "custom_friction/friction_law.hpp"
#include "custom_friction/coulomb_adhesion_friction_law.hpp"
#include "custom_friction/hardening_coulomb_friction_law.hpp"

namespace Kratos
{

  namespace Python
  {

    using namespace boost::python;

    typedef Properties::Pointer                    PropertiesPointer;
    typedef FrictionLaw                          FrictionLawBaseType;
    typedef FrictionLaw::Pointer                  FrictionLawPointer;
    typedef std::vector<FrictionLaw::Pointer>   FrictionLawContainer;

    typedef CoulombAdhesionFrictionLaw           CoulombFrictionLawBaseType;
    typedef CoulombAdhesionFrictionLaw::Pointer   CoulombFrictionLawPointer;

    void Push_Back_Friction_Laws( FrictionLawContainer& ThisFrictionLawContainer,
    				  FrictionLawPointer ThisFrictionLaw )
    {
      ThisFrictionLawContainer.push_back( ThisFrictionLaw );
    }

    void  AddCustomFrictionLawsToPython()
    {

       class_< FrictionLawContainer >( "FrictionLawContainer", init<>() )
      	.def( "PushBack", Push_Back_Friction_Laws )
      	;

       class_<Variable<FrictionLaw::Pointer>, bases<VariableData>, boost::noncopyable >( "FrictionLawVariable", no_init )
	 .def( self_ns::str( self ) )
	 ;
       
       //Friction laws
       class_< FrictionLaw, FrictionLaw::Pointer, boost::noncopyable >
      	( "FrictionLaw",  init<>() )
	 .def("Clone",&FrictionLaw::Clone)
	 ;

       class_< CoulombAdhesionFrictionLaw, bases< FrictionLawBaseType >, boost::noncopyable >
      	( "CoulombAdhesionFrictionLaw",  init<>() )
      	;

       class_< HardeningCoulombFrictionLaw, bases< CoulombFrictionLawBaseType >, boost::noncopyable >
       	( "HardeningCoulombFrictionLaw",  init<>() )
       	;

    }

  }  // namespace Python.

}  // namespace Kratos.

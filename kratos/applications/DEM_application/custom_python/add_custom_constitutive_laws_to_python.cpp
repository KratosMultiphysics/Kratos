//
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
#include "custom_constitutive/DEM_D_Linear_viscous_Coulomb_CL.h"
#include "custom_constitutive/DEM_D_Hertz_viscous_Coulomb_CL.h"
#include "custom_constitutive/DEM_D_Linear_viscous_Coulomb_2D_CL.h"
#include "custom_constitutive/DEM_D_Hertz_viscous_Coulomb_2D_CL.h"
#include "custom_constitutive/DEM_D_JKR_cohesive_law.h"
#include "custom_constitutive/DEM_D_DMT_cohesive_law.h"

#include "custom_constitutive/DEM_continuum_constitutive_law.h" 
#include "custom_constitutive/DEM_Dempack_CL.h"
#include "custom_constitutive/DEM_Dempack_2D_CL.h"
#include "custom_constitutive/DEM_KDEM_CL.h"
#include "custom_constitutive/DEM_ExponentialHC_CL.h"




namespace Kratos {

    namespace Python {

        using namespace boost::python;

        void AddCustomConstitutiveLawsToPython() {


//            DEM Discontinuum Constitutive Laws :  

            class_< DEMDiscontinuumConstitutiveLaw, DEMDiscontinuumConstitutiveLaw::Pointer, boost::noncopyable > //bases< ConstitutiveLawBaseType >
                    ("DEMDiscontinuumConstitutiveLaw",
                    init<>())
                    .def("Clone", &DEMDiscontinuumConstitutiveLaw::Clone)
                    .def("SetConstitutiveLawInProperties", &DEMDiscontinuumConstitutiveLaw::SetConstitutiveLawInProperties)
                    ;

            class_<Variable<DEMDiscontinuumConstitutiveLaw::Pointer>, boost::noncopyable >("DEMDiscontinuumConstitutiveLawPointerVariable", no_init)
                    .def(self_ns::str(self))
                    ;

            class_<DEM_D_Linear_viscous_Coulomb, bases< DEMDiscontinuumConstitutiveLaw >, boost::noncopyable >("DEM_D_Linear_viscous_Coulomb",init<>())
                    ;          
            
            class_<DEM_D_Linear_viscous_Coulomb2D, bases< DEM_D_Linear_viscous_Coulomb >, boost::noncopyable >("DEM_D_Linear_viscous_Coulomb2D",init<>())
                    ;  
          
            class_<DEM_D_Hertz_viscous_Coulomb, bases< DEMDiscontinuumConstitutiveLaw >, boost::noncopyable >("DEM_D_Hertz_viscous_Coulomb",init<>())
                    ;

            class_<DEM_D_Hertz_viscous_Coulomb2D, bases< DEM_D_Hertz_viscous_Coulomb >, boost::noncopyable >("DEM_D_Hertz_viscous_Coulomb2D",init<>())
                    ;  
          
            class_<DEM_compound_constitutive_law<DEM_D_Hertz_viscous_Coulomb, DEM_D_JKR_Cohesive_Law>, bases< DEM_D_Hertz_viscous_Coulomb >, boost::noncopyable >("DEM_D_Hertz_viscous_Coulomb_JKR", init<>())
                    ;

            class_<DEM_compound_constitutive_law<DEM_D_Hertz_viscous_Coulomb, DEM_D_DMT_Cohesive_Law>, bases< DEM_D_Hertz_viscous_Coulomb >, boost::noncopyable >("DEM_D_Hertz_viscous_Coulomb_DMT", init<>())
                    ;
            
            class_<DEM_compound_constitutive_law<DEM_D_Linear_viscous_Coulomb, DEM_D_JKR_Cohesive_Law>, bases< DEM_D_Linear_viscous_Coulomb >, boost::noncopyable >("DEM_D_Linear_viscous_Coulomb_JKR", init<>())
                    ;
            
            class_<DEM_compound_constitutive_law<DEM_D_Linear_viscous_Coulomb, DEM_D_DMT_Cohesive_Law>, bases< DEM_D_Linear_viscous_Coulomb >, boost::noncopyable >("DEM_D_Linear_viscous_Coulomb_DMT", init<>())
                    ;
            
//            DEM Continuum Constitutive Laws :  

            class_< DEMContinuumConstitutiveLaw, DEMContinuumConstitutiveLaw::Pointer, boost::noncopyable > //bases< ConstitutiveLawBaseType >
                    ("DEMContinuumConstitutiveLaw",
                    init<>())
                    .def("Clone", &DEMContinuumConstitutiveLaw::Clone)
                    .def("SetConstitutiveLawInProperties", &DEMContinuumConstitutiveLaw::SetConstitutiveLawInProperties)
                    ;

            class_<Variable<DEMContinuumConstitutiveLaw::Pointer>, boost::noncopyable >("DEMContinuumConstitutiveLawPointerVariable", no_init)
                    .def(self_ns::str(self))
                    ;

            class_<DEM_Dempack, bases< DEMContinuumConstitutiveLaw >, boost::noncopyable >("DEM_Dempack",init<>())
                    ;

            class_<DEM_Dempack2D, bases< DEM_Dempack >, boost::noncopyable >("DEM_Dempack2D",init<>())
                    ;

            class_<DEM_KDEM, bases< DEMContinuumConstitutiveLaw >, boost::noncopyable >("DEM_KDEM",init<>())
                    ;
            
            class_<DEM_ExponentialHC, bases< DEMContinuumConstitutiveLaw >, boost::noncopyable >("DEM_ExponentialHC",init<>())
                    ;
        }

    } // namespace Python.
} // namespace Kratos.

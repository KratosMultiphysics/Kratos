//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// System includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "includes/define.h" 

#include "../custom_constitutive/DEM_discontinuum_constitutive_law.h"    
#include "../custom_constitutive/DEM_continuum_constitutive_law.h" 
#include "../custom_constitutive/DEM_compound_constitutive_law.h" 

#include "../custom_constitutive/DEM_D_Linear_viscous_Coulomb_CL.h"
#include "../custom_constitutive/DEM_D_Hertz_viscous_Coulomb_CL.h"
#include "../custom_constitutive/DEM_D_Bentonite_Colloid_CL.h"
#include "../custom_constitutive/DEM_D_Linear_viscous_Coulomb_2D_CL.h"
#include "../custom_constitutive/DEM_D_Hertz_viscous_Coulomb_2D_CL.h"
#include "../custom_constitutive/DEM_D_JKR_cohesive_law.h"
#include "../custom_constitutive/DEM_D_DMT_cohesive_law.h"

#include "../custom_constitutive/DEM_Dempack_CL.h"
#include "../custom_constitutive/DEM_Dempack_2D_CL.h"
#include "../custom_constitutive/DEM_KDEM_CL.h"
#include "../custom_constitutive/DEM_KDEM_Rankine_CL.h"
#include "../custom_constitutive/DEM_KDEM_Mohr_Coulomb_CL.h"
#include "../custom_constitutive/DEM_sintering_continuum_CL.h"
#include "../custom_constitutive/DEM_KDEM_fabric_CL.h"
#include "../custom_constitutive/DEM_ExponentialHC_CL.h"
#include "../custom_constitutive/DEM_Dempack_torque_CL.h"
#include "../custom_constitutive/DEM_Dempack_dev_CL.h"
#include "../custom_constitutive/DEM_Dempack_2D_dev_CL.h"
#include "../custom_constitutive/dem_d_linear_custom_constants_cl.h"
#include "../custom_constitutive/DEM_D_Conical_damage_CL.h"
#include "../custom_constitutive/dem_kdem_2d_cl.h"
#include "../custom_constitutive/dem_kdem_fabric_2d_cl.h"

namespace Kratos {

    namespace Python {

        using namespace boost::python;

        void AddCustomConstitutiveLawsToPython() {

            // DEM Discontinuum Constitutive Laws:  

            class_< DEMDiscontinuumConstitutiveLaw, DEMDiscontinuumConstitutiveLaw::Pointer, boost::noncopyable > //bases< ConstitutiveLawBaseType >
                    ("DEMDiscontinuumConstitutiveLaw",
                    init<>())
                    .def("Clone", &DEMDiscontinuumConstitutiveLaw::Clone)
                    .def("SetConstitutiveLawInProperties", &DEMDiscontinuumConstitutiveLaw::SetConstitutiveLawInProperties)
                    .def("GetTypeOfLaw", &DEMDiscontinuumConstitutiveLaw::GetTypeOfLaw)
                    ;

            class_<Variable<DEMDiscontinuumConstitutiveLaw::Pointer>, boost::noncopyable >("DEMDiscontinuumConstitutiveLawPointerVariable", no_init)
                    .def(self_ns::str(self))
                    ;

            class_<DEM_D_Linear_viscous_Coulomb, bases< DEMDiscontinuumConstitutiveLaw >, boost::noncopyable >("DEM_D_Linear_viscous_Coulomb",init<>())
                    ;          

            class_<DEM_D_Bentonite_Colloid, bases< DEMDiscontinuumConstitutiveLaw >, boost::noncopyable >("DEM_D_Bentonite_Colloid",init<>())
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
                    
            class_<DEM_D_Linear_Custom_Constants, bases< DEM_D_Linear_viscous_Coulomb >, boost::noncopyable >("DEM_D_Linear_Custom_Constants",init<>())
                    ;
                    
            class_<DEM_D_Conical_damage, bases< DEMDiscontinuumConstitutiveLaw >, boost::noncopyable >("DEM_D_Conical_damage",init<>())
                    ;
            
            // DEM Continuum Constitutive Laws:  

            class_< DEMContinuumConstitutiveLaw, DEMContinuumConstitutiveLaw::Pointer, boost::noncopyable > //bases< ConstitutiveLawBaseType >
                    ("DEMContinuumConstitutiveLaw",
                    init<>())
                    .def("Clone", &DEMContinuumConstitutiveLaw::Clone)
                    .def("SetConstitutiveLawInProperties", &DEMContinuumConstitutiveLaw::SetConstitutiveLawInProperties)
                    .def("GetTypeOfLaw", &DEMContinuumConstitutiveLaw::GetTypeOfLaw)
                    ;

            class_<Variable<DEMContinuumConstitutiveLaw::Pointer>, boost::noncopyable >("DEMContinuumConstitutiveLawPointerVariable", no_init)
                    .def(self_ns::str(self))
                    ;

            class_<DEM_Dempack, bases< DEMContinuumConstitutiveLaw >, boost::noncopyable >("DEM_Dempack",init<>())
                    ;

            class_<DEM_Dempack2D, bases< DEM_Dempack >, boost::noncopyable >("DEM_Dempack2D",init<>())
                    ;

            class_<DEM_Dempack_torque, bases< DEM_Dempack >, boost::noncopyable >("DEM_Dempack_torque",init<>())
                    ;

            class_<DEM_Dempack_dev, bases< DEM_Dempack >, boost::noncopyable >("DEM_Dempack_dev",init<>())
                    ;

            class_<DEM_Dempack2D_dev, bases< DEM_Dempack_dev >, boost::noncopyable >("DEM_Dempack2D_dev",init<>())
                    ;

            class_<DEM_KDEM, bases< DEMContinuumConstitutiveLaw >, boost::noncopyable >("DEM_KDEM",init<>())
                    ;

            class_<DEM_sintering_continuum, bases< DEMContinuumConstitutiveLaw >, boost::noncopyable >("DEM_sintering_continuum", init<>())
                    ;
            
            class_<DEM_KDEMFabric, bases< DEM_KDEM >, boost::noncopyable >("DEM_KDEMFabric",init<>())
                    ;
            
            class_<DEM_KDEM_Rankine, bases< DEM_KDEM >, boost::noncopyable >("DEM_KDEM_Rankine",init<>())
                    ;
            
            class_<DEM_KDEM_Mohr_Coulomb, bases< DEM_KDEM_Rankine >, boost::noncopyable >("DEM_KDEM_Mohr_Coulomb",init<>())
                    ;
            
            class_<DEM_KDEM2D, bases< DEM_KDEM >, boost::noncopyable >("DEM_KDEM2D",init<>())
                    ;
            
            class_<DEM_KDEMFabric2D, bases< DEM_KDEM2D >, boost::noncopyable >("DEM_KDEMFabric2D",init<>())
                    ;
            
            class_<DEM_ExponentialHC, bases< DEMContinuumConstitutiveLaw >, boost::noncopyable >("DEM_ExponentialHC",init<>())
                    ;
        }

    } // namespace Python.
} // namespace Kratos.

/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Nelson $
//   Date:                $Date: 2009-02-02 14:03:23 $
//   Revision:            $Revision: 1.5 $
//

// System includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "python/pointer_vector_set_python_interface.h"
#include "python/variable_indexing_python.h"
#include "custom_python/add_fluency_criteria_to_python.h"
#include "fluency_criteria/fluency_criteria.h"
#include "python/vector_python_interface.h"
#include "fluency_criteria/energy_yield_function.h"
#include "fluency_criteria/isotropic_rankine_yield_function.h"
//#include "fluency_criteria/tresca_yield_function.h"
#include "fluency_criteria/von_mises_yield_function.h"     
#include "fluency_criteria/modified_morh_coulomb_yield_function.h"
#include "fluency_criteria/morh_coulomb_yield_function.h"
#include "fluency_criteria/standard_morh_coulomb_yield_function.h"

//#include "fluency_criteria/rankine_yield_function.h"
//#include "fluency_criteria/drucker_prager_yield_function.h"

#include "soft_hard_behavior/softening_hardening_criteria.h"
#include "soft_hard_behavior/exponencial_softening.h"
#include "soft_hard_behavior/linear_softening.h"
#include "soft_hard_behavior/cohesion_softening.h"
#include "soft_hard_behavior/friction_softening.h"
#include "soft_hard_behavior/dilatancy_softening.h"


#include "spaces/ublas_space.h"





namespace Kratos
{

    namespace Python
    {
	    using namespace boost::python;
            typedef FluencyCriteria                   FluencyCriteriaType; 
	    typedef SofteningHardeningCriteria        SofteningHardeningType;
	    typedef Morh_Coulomb_Yield_Function       MorhCoulombType;
	    typedef Isotropic_Rankine_Yield_Function  RankineType;
	    
	    
	    typedef FluencyCriteria::Pointer            FluencyPointerType;
            typedef SofteningHardeningCriteria::Pointer SofteningHardeningPointerType;
	    typedef MorhCoulombType::Pointer            MorhCoulombPointerType;
	    typedef RankineType::Pointer                RankinePointerType;
	   
	        void  AddFluencyCriteriaToPython()
			      {

			      class_< FluencyCriteriaType, boost::noncopyable >
			      ("FluencyCriteriaType",
				init<>() )
			      ;

			      class_< SofteningHardeningType, boost::noncopyable >
			      ("SofteningHardeningType",
				init<>() )
			      ;

			      enum_<myState>("State")
				 .value("Plane_Stress", Plane_Stress)
				 .value("Plane_Strain", Plane_Strain)
                                 .value("Tri_D", Tri_D)
				  ;

			      enum_<myPotencialPlastic>("PotencialPlastic")
				 .value("Not_Associated", Not_Associated)
				 .value("Associated", Associated)
                                  ;

           
	                     /*			 
			       class_<Rankine_Yield_Function, bases< FluencyCriteriaType >, boost::noncopyable >
			      ("RankineYieldFunction",
			      init<myState> () )
			      ;  
			      */
				 
                              class_<Isotropic_Rankine_Yield_Function, bases< FluencyCriteriaType >, boost::noncopyable >
			      ("IsotropicRankineYieldFunction",
			      init<SofteningHardeningPointerType,
			      myState> () )
			      ;  
 
                              /*
			      class_<Tresca_Yield_Function, bases< FluencyCriteriaType >, boost::noncopyable >
			      ("TrescaYieldFunction",
			      init<myState> () )
			      ;  
                              */
			      
			      class_<Von_Mises_Yield_Function, bases< FluencyCriteriaType >, boost::noncopyable >
			      ("VonMisesYieldFunction",
			       init<const myState&, const SofteningHardeningPointerType&> () )
			       //init<myState, myPotencialPlastic> () )
			      ;  

			      class_<Modified_Morh_Coulomb_Yield_Function, bases< FluencyCriteriaType >, boost::noncopyable >
			      ("ModifiedMorhCoulombYieldFunction",
			      init<myState, MorhCoulombPointerType, RankinePointerType> () )
			      /* init<const SofteningHardeningPointerType&,
			           const SofteningHardeningPointerType&,
			           const SofteningHardeningPointerType&,
			           const SofteningHardeningPointerType&,
			           const myState, 
			           const myPotencialPlastic> () )
			           */
			      ;  
			      

			      class_<Morh_Coulomb_Yield_Function, bases< FluencyCriteriaType >, boost::noncopyable >
			      ("MorhCoulombYieldFunction",
			      init<const SofteningHardeningPointerType&,
			           const SofteningHardeningPointerType&,
			           const SofteningHardeningPointerType&,
			           const myState, 
			           const myPotencialPlastic> () )
			      ;  
			      
			      class_<Standard_Morh_Coulomb_Yield_Function, bases< FluencyCriteriaType >, boost::noncopyable >
			      ("StandardMorhCoulombYieldFunction",
			      init<const SofteningHardeningPointerType&,
			           const myState> () )
			      ;  
			      
			      /*
			      class_<Drucker_Prager_Yield_Function, bases< FluencyCriteriaType >, boost::noncopyable >
			      ("DruckerPragerYieldFunction",
			      init<myState> () )
			      ;  
                               */
			      
			      
			      class_<Energy_Yield_Function, bases< FluencyCriteriaType >, boost::noncopyable >
			      ("EnergyYieldFunction",
			      init<myState> () )
			      ;  
                              
			      
			      class_<Exponential_Softening, bases< SofteningHardeningType>, boost::noncopyable >
			      ("ExponentialSoftening",
			      init<> () )
			      ; 

			      class_<Linear_Softening, bases< SofteningHardeningType>, boost::noncopyable >
			      ("LinearSoftening",
			      init<> () )
			      ; 
			      
			      class_<Cohesion_Softening, bases< SofteningHardeningType>, boost::noncopyable >
			      ("CohesionSoftening",
			      init<> () )
			      ; 
			      
			      class_<Friction_Softening, bases< SofteningHardeningType>, boost::noncopyable >
			      ("FrictionSoftening",
			      init<> () )
			      ; 
			      
			      class_<Dilatancy_Softening, bases< SofteningHardeningType>, boost::noncopyable >
			      ("DilatancySoftening",
			      init<> () )
			      ; 
                              

                                }
                

	}  // namespace Python.
  
}  // namespace Kratos.


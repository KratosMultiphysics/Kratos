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
#include "fluency_criteria/energy_yield_function.h"
#include "fluency_criteria/rankine_yield_function.h"
#include "fluency_criteria/tresca_yield_function.h"


#include "spaces/ublas_space.h"





namespace Kratos
{

    namespace Python
    {
	    using namespace boost::python;
            typedef FluencyCriteria  FluencyCriteriaBaseType; 
	   
	        void  AddFluencyCriteriaToPython()
			      {

			      class_< FluencyCriteriaBaseType, boost::noncopyable >
			      ("FluencyCriteriaBaseType",
				init<>() )
			      ;

			     class_<Energy_Yield_Function, bases< FluencyCriteriaBaseType >, boost::noncopyable >
			      ("EnergyYieldFunction",
			      init<int, double > () ) // dimesion, limite elastico
			      ;  
                             
			       class_<Rankine_Yield_Function, bases< FluencyCriteriaBaseType >, boost::noncopyable >
			      ("RankineYieldFunction",
			      init<int, double > () )
			      ;  

			      class_<Tresca_Yield_Function, bases< FluencyCriteriaBaseType >, boost::noncopyable >
			      ("TrescaYieldFunction",
			      init<int, double > () )
			      ;  


                                }
                

	}  // namespace Python.
  
}  // namespace Kratos.


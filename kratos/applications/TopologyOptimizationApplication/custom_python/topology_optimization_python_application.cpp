/*
==============================================================================
KratosTopologyOptimizationApplication 
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
//   Project Name:        $Project:     Topology_Optimization_Application $
//   Last modified by:    $Author:      Malfavón Farías, Baumgärtner      $
//   Date:                $Date:        May 2016                          $
//   Revision:            $Revision:    0.0                               $
//


// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>


// Project includes 
#include "includes/define.h"
#include "topology_optimization_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"

 
namespace Kratos
{

namespace Python
{

  using namespace boost::python;


  
  BOOST_PYTHON_MODULE(KratosTopologyOptimizationApplication)
  {

	  class_<KratosTopologyOptimizationApplication, 
                          KratosTopologyOptimizationApplication::Pointer, 
                          bases<KratosApplication>, boost::noncopyable >("KratosTopologyOptimizationApplication")
                         ;

    AddCustomStrategiesToPython();
    AddCustomUtilitiesToPython();

    //Registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( E_MIN )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( E_0 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( PENAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( X_PHYS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( X_PHYS_OLD )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DCDX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( DVDX )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( SOLID_VOID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE( LOCAL_STRAIN_ENERGY )


  }
  
  
}  // namespace Python.
  
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined

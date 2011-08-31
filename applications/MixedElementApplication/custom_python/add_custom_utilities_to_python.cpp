/*
==============================================================================
KratosTestApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
//#include "linear_solvers/linear_solver.h"
#include "custom_utilities/plot_damage.h"
#include "custom_utilities/check_constitutive_law.h"

#include "includes/model_part.h"
#include "includes/constitutive_law.h"

namespace Kratos
{
	
namespace Python
{

	
  void  AddCustomUtilitiesToPython()
  {
	using namespace boost::python;


//		typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
//		typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
//		typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

        class_<PlotUtility,  boost::noncopyable > ("PlotUtils", init<>())
                    .def("PlotVariable", &PlotUtility::PlotVariable)
                    ;

        class_<CheckConstitutiveLawUtility,  boost::noncopyable > ("CheckConstitutiveLawUtils", init<Properties::Pointer , ProcessInfo&, double >())
                    .def("CheckConstitutiveLaw", &CheckConstitutiveLawUtility::CheckConstitutiveLaw)
                    .def("GetScalarVar",&CheckConstitutiveLawUtility::GetScalarVar)
                    .def("InitializeSolutionStep",&CheckConstitutiveLawUtility::InitializeSolutionStep)
                    .def("FinalizeSolutionStep",&CheckConstitutiveLawUtility::FinalizeSolutionStep)
                    ;

  }
	




}  // namespace Python.

} // Namespace Kratos


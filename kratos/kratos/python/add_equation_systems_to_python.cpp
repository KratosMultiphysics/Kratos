/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "python/add_equation_systems_to_python.h" 
#include "equation_systems/equation_system.h"
#include "includes/dof.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
	
namespace Python
{
  void  AddEquationSystemsToPython()
  {
typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
typedef EquationSystem<SpaceType,  LocalSpaceType, Dof<double> > EquationSystemType;

	void (EquationSystemType::*pointer_to_initalize)() = &EquationSystemType::Initialize;
	const EquationSystemType::SystemMatrixType& (EquationSystemType::*pointer_to_get_system_matrix)() const = &EquationSystemType::GetSystemMatrix;
	const EquationSystemType::SystemVectorType& (EquationSystemType::*pointer_to_get_results)() const = &EquationSystemType::GetResults;
	const EquationSystemType::SystemVectorType& (EquationSystemType::*pointer_to_get_rhs)() const = &EquationSystemType::GetRightHandSide;
	const EquationSystemType::SystemMatrixType& (EquationSystemType::*pointer_to_get_dirichlet_matrix)() const = &EquationSystemType::GetDirichletMatrix;


  using namespace boost::python;

  class_<EquationSystemType, EquationSystemType::Pointer>("EquationSystem")
		  //.def(init<std::string const&>())
		  .def(init<EquationSystemType::SizeType>())
		  .def("Initialize",pointer_to_initalize)
		  .def("ApplyDirichletConditions",&EquationSystemType::ApplyDirichletConditions)
		  .def("Size",&EquationSystemType::Size)
		  .def("DirichletSize",&EquationSystemType::DirichletSize)
		  .add_property("SystemMatrix", make_function(pointer_to_get_system_matrix, return_internal_reference<>()),&EquationSystemType::SetSystemMatrix)
		  .add_property("Results", make_function(pointer_to_get_results, return_internal_reference<>()),&EquationSystemType::SetResults)
		  .add_property("RightHandSide", make_function(pointer_to_get_rhs, return_internal_reference<>()),&EquationSystemType::SetRightHandSide)
		  .add_property("DirichletMatrix", make_function(pointer_to_get_dirichlet_matrix, return_internal_reference<>()))
		  //.def("",&EquationSystemType::)
		  .def(self_ns::str(self))
		  ;
  }
	
}  // namespace Python.

} // Namespace Kratos


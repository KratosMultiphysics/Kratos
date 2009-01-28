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
#include "includes/tracer.h"
#include "linear_solvers_python.h" 
#include "linear_solvers/kratos_linear_solvers.h"
#include "vectorial_spaces/kratos_vectorial_spaces.h"


namespace Kratos
{

namespace Python
{

typedef Kratos::LinearSolver<Kratos::CSRSpace<double>, Kratos::DenseSpace<double> > KratosLinearSolverType;

typedef Kratos::DirectSolver<Kratos::CSRSpace<double>, Kratos::DenseSpace<double> > KratosDirectSolverType;

typedef Kratos::IterativeSolver<Kratos::CSRSpace<double>, Kratos::DenseSpace<double> > KratosIterativeSolverType;

#define KRATOS_REGISTER_LINEAR_SOLVERS_BEGIN using namespace Kratos;

#define KRATOS_REGISTER_LINEAR_SOLVER typedef 
#define KRATOS_AS_LINEAR_SOLVER_NAMED(name) KRATOS_TYPE_NAME_OF(name); \
    bool (name##Type::*name##_solve_method_pointer)(name##Type::SparseMatrixType&, name##Type::VectorType&, name##Type::VectorType&) = &name##Type::Solve; \
       class_<name##Type, bases<name##Type::BaseType> >(#name) \
	 .def("Solve", name##_solve_method_pointer) \
	 .def(self_ns::str(self)) \

#define KRATOS_REGISTER_PRECONDITIONER typedef 
#define KRATOS_AS_PRECONDITIONER_NAMED(name) KRATOS_TYPE_NAME_OF(name); \
    void (name##Type::*name##_initialize_method_pointer)(name##Type::SparseMatrixType&, name##Type::VectorType&, name##Type::VectorType&) = &name##Type::Initialize; \
       class_<name##Type, bases<name##Type::BaseType> >(#name) \
	 .def("Initialize", name##_initialize_method_pointer) \
	 .def(self_ns::str(self)) \

const KratosIterativeSolverType::PreconditionerType::Pointer  (KratosIterativeSolverType::*IterativeSolver_get_preconditioner_method_pointer)(void) const = &KratosIterativeSolverType::GetPreconditioner;

#define KRATOS_AS_ITERATIVE_LINEAR_SOLVER_NAMED(name) \
KRATOS_AS_LINEAR_SOLVER_NAMED(name) \
KRATOS_LINEAR_SOLVER_ADDITIONAL_METHOD(.def(init<double>())) \
KRATOS_LINEAR_SOLVER_ADDITIONAL_METHOD(.def(init<double,unsigned int>())) \
KRATOS_LINEAR_SOLVER_ADDITIONAL_METHOD(.def(init<double,unsigned int, Preconditioner<CSRSpace<double>, DenseSpace<double> >::Pointer >()))

#define KRATOS_LINEAR_SOLVER_ADDITIONAL_METHOD(method) method

//   void Reorder(CSRMatrix<double>& rA, Vector<double>& rX, Vector<double>& rB)
using namespace boost::python;

  void  AddLinearSolversToPython()
  {
       class_<Kratos::Preconditioner<Kratos::CSRSpace<double>, Kratos::DenseSpace<double> > >("Preconditioner", no_init)
	 .def(self_ns::str(self))
        ;

       class_<Kratos::CuthillMcKeeReorderer<Kratos::CSRSpace<double>, Kratos::DenseSpace<double> > >("CuthillMcKeeReorderer")
	 .def("Initialize" , &Kratos::CuthillMcKeeReorderer<Kratos::CSRSpace<double>, Kratos::DenseSpace<double> >::Initialize)
	 .def("Reorder" , &Kratos::CuthillMcKeeReorderer<Kratos::CSRSpace<double>, Kratos::DenseSpace<double> >::Reorder)
	 .def("InverseReorder" , &Kratos::CuthillMcKeeReorderer<Kratos::CSRSpace<double>, Kratos::DenseSpace<double> >::InverseReorder)
        ;

       class_<KratosLinearSolverType>("LinearSolver", no_init)
	 .def(self_ns::str(self))
        ;

       class_<KratosDirectSolverType, bases<KratosLinearSolverType> >("DirectSolver",no_init)
	 .def(self_ns::str(self))
        ;
    
       class_<KratosIterativeSolverType, bases<KratosLinearSolverType> >("IterativeSolver",no_init)
 	 .add_property("Preconditioner", IterativeSolver_get_preconditioner_method_pointer, &KratosIterativeSolverType::SetPreconditioner)
 	 .add_property("MaxIterationsNumber", &KratosIterativeSolverType::GetMaxIterationsNumber, &KratosIterativeSolverType::SetMaxIterationsNumber)
 	 .add_property("IterationsNumber", &KratosIterativeSolverType::GetIterationsNumber, &KratosIterativeSolverType::SetIterationsNumber)
 	 .add_property("Tolerance", &KratosIterativeSolverType::GetTolerance, &KratosIterativeSolverType::SetTolerance)
 	 .add_property("ResidualNorm", &KratosIterativeSolverType::GetResidualNorm, &KratosIterativeSolverType::SetResidualNorm)
 	 .def("IterationNeeded", &KratosIterativeSolverType::IterationNeeded)
 	 .def("IsConverged", &KratosIterativeSolverType::IsConverged)
	 .def(self_ns::str(self))
        ;

#include "components/linear_solvers.h"

  }


}  // namespace Python.
  
} // Namespace Kratos


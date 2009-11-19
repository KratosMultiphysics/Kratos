/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Farshid Mossaiby
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
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


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "python/add_equation_systems_to_python.h" 
#include "spaces/ublas_space.h"

#include "linear_solvers/direct_solver.h"
#include "linear_solvers/iterative_solver.h"
#include "gpu_bicgstab_solver.h"
#include "gpu_cg_solver.h"
#include "amg_solver.h"

#include "GPUPreconditioner.h"
#include "AMGpreconditioner.h"
#include "Kratos_AMGpreconditioner.h"
#include "Diagonalpreconditioner.h"

namespace Kratos
{
    
namespace Python
{
    void  AddLinearSolversToPython()
    {
        typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
        typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
        typedef DirectSolver<SpaceType,  LocalSpaceType> DirectSolverType;
        typedef IterativeSolver<SpaceType, LocalSpaceType> IterativeSolverType;
        typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;
        typedef GPUBICGSTABSolver GPUBICGSTABSolverType;
	typedef GPUCGSolver<SpaceType,  LocalSpaceType> GPUCGSolverType;
	typedef AMGSolver<SpaceType,  LocalSpaceType> AMGSolverType;

	typedef GPUPreconditioner GPUPreconditionerType;
	typedef AMGpreconditioner AMGPreconditionerType;
	typedef Kratos_AMGpreconditioner Kratos_AMGPreconditionerType;
	typedef Diagonalpreconditioner GpuDiagonalPreconditionerType;

        
        using namespace boost::python; 
        
        //***************************************************************************
        //linear solvers
        //***************************************************************************
        
        class_<GPUBICGSTABSolverType, GPUBICGSTABSolverType::Pointer,
        bases<LinearSolverType> >( "GPUBICGSTABSolver" )
                .def(init<double, unsigned int>() )
                .def(init<double, unsigned int, Kratos_AMGPreconditionerType& >() )
		.def(init<double, unsigned int, GpuDiagonalPreconditionerType& >() )
                .def(self_ns::str(self))
                ;
	
	class_<GPUPreconditionerType >( "GPUPreconditioner" )
                .def(init<>() )
                ;

	class_<AMGPreconditionerType, bases<GPUPreconditionerType> >( "AMGPreconditioner" )
		//.def(init<>() )
                .def(init<double, size_t, bool, size_t, size_t, size_t*, size_t*, bool>() )
                ;
	class_<GpuDiagonalPreconditionerType, bases<GPUPreconditionerType> >( "GPUDiagonalPreconditioner" )
		.def(init<>() )
	        ;

	class_<Kratos_AMGPreconditionerType, bases<AMGPreconditionerType> >( "KratosAMGPreconditioner" )
                .def(init<double, size_t, bool, size_t, size_t, const Vector, const Vector, bool>() )
                ;


	class_<GPUCGSolverType, GPUCGSolverType::Pointer,
        bases<LinearSolverType> >( "GPUCGSolver" )
                .def(init<double, unsigned int>() )
                .def(init<double, unsigned int, Kratos_AMGPreconditionerType& >() )
		.def(init<double, unsigned int, GpuDiagonalPreconditionerType& >() )
                .def(self_ns::str(self))
                ;
	
	class_<AMGSolverType, AMGSolverType::Pointer,
        bases<LinearSolverType> >( "AMGSolver" )
                .def(init<double, unsigned int, double, unsigned int, bool, unsigned int, unsigned int, const Vector, const Vector>() )
                .def(self_ns::str(self))
                ;	
  }
	
}  // namespace Python.

} // Namespace Kratos


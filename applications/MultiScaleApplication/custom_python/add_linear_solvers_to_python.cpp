/*
==============================================================================
KratosMultiScaleApplication
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
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:37:00 $
//   Revision:            $Revision: 1.00 $
//
//

// System includes

// External includes
#include <boost/python.hpp>
#include "spaces/ublas_space.h"

// Project includes
#include "add_linear_solvers_to_python.h"

#ifdef MULTISCALE_APPLICATION_USE_MUMPS
#include "custom_linear_solvers/mumps_linear_solver.h"  
#endif // MULTISCALE_APPLICATION_USE_MUMPS

#ifdef MULTISCALE_APPLICATION_USE_EIGEN
#include "custom_linear_solvers/eigenlib_lu_linear_solver.h"  
#endif // MULTISCALE_APPLICATION_USE_EIGEN

#include "custom_linear_solvers/skyline_lu_linear_solver_v2.h"

namespace Kratos
{
namespace Python
{

using namespace boost::python;

void AddLinearSolversToPython()
{
	typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
	typedef Reorderer<SpaceType,  LocalSpaceType > ReordererType;
    typedef LinearSolver<SpaceType,  LocalSpaceType> LinearSolverType;
	typedef DirectSolver<SpaceType,  LocalSpaceType, ReordererType > DirectSolverType;


#ifdef MULTISCALE_APPLICATION_USE_MUMPS
	typedef MUMPSLinearSolver<SpaceType, LocalSpaceType, ReordererType > MUMPSLinearSolverType;
	class_<MUMPSLinearSolverType, MUMPSLinearSolverType::Pointer, bases<DirectSolverType> >("MUMPSSolver")
		.def(init< >())
		.def(init< int > ())
		.def(self_ns::str(self))
		;  
#endif // MULTISCALE_APPLICATION_USE_MUMPS

#ifdef MULTISCALE_APPLICATION_USE_EIGEN
	typedef EigenlibLULinearSolver<SpaceType, LocalSpaceType, ReordererType > EigenlibLULinearSolverType;
	class_<EigenlibLULinearSolverType, EigenlibLULinearSolverType::Pointer, bases<DirectSolverType> >("EigenlibLULinearSolver")
		.def(init< >())
		.def(self_ns::str(self))
		;  
#endif // MULTISCALE_APPLICATION_USE_EIGEN

	typedef SkylineLUFactorizationLinearSolverV2<SpaceType, LocalSpaceType, ReordererType > SkylineLUFactorizationLinearSolverV2Type;
	class_<SkylineLUFactorizationLinearSolverV2Type, SkylineLUFactorizationLinearSolverV2Type::Pointer, bases<DirectSolverType> >
		("SkylineLUFactorizationSolverV2")
		.def(init< >())
		.def(self_ns::str(self))
		;

}


}

}
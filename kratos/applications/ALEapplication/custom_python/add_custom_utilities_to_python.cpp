/*
==============================================================================
KratosPFEMApplication
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

/* ****************************************************************************
 *  Projectname:         $KratosALEApplication
 *  Last Modified by:    $Author: A.Winterstein@tum.de $
 *  Date:                $Date: June 2016 $
 *  Revision:            $Revision: 1.5 $
 * ***************************************************************************/


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/ball_vertex_meshmoving.h"
#include "custom_utilities/ball_vertex_meshmoving3D.h"
#include "custom_utilities/move_mesh_utilities.h"


namespace Kratos
{

namespace Python
{


void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;


    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    class_< BallVertexMeshMoving< 2, SparseSpaceType, LinearSolverType >,  boost::noncopyable >	("BallVertexMeshMoving2D", init<	>() )
    .def("ConstructSystem",&BallVertexMeshMoving< 2, SparseSpaceType, LinearSolverType >::ConstructSystem)
    .def("BuildAndSolveSystem",&BallVertexMeshMoving< 2, SparseSpaceType, LinearSolverType >::BuildAndSolveSystem)
    .def("ClearSystem",&BallVertexMeshMoving< 2, SparseSpaceType, LinearSolverType >::ClearSystem)
    ;


    class_< BallVertexMeshMoving3D< 3, SparseSpaceType, LinearSolverType >,  boost::noncopyable >	("BallVertexMeshMoving3D", init<	>() )
    .def("ConstructSystem",&BallVertexMeshMoving3D< 3, SparseSpaceType, LinearSolverType >::ConstructSystem)
    .def("BuildAndSolveSystem",&BallVertexMeshMoving3D< 3, SparseSpaceType, LinearSolverType >::BuildAndSolveSystem)
    .def("ClearSystem",&BallVertexMeshMoving3D< 3, SparseSpaceType, LinearSolverType >::ClearSystem)
    ;


    class_< MoveMeshUtilities,  boost::noncopyable >	("MoveMeshUtilities", init<	>() )
    .def("BDF_MoveMesh",&MoveMeshUtilities::BDF_MoveMesh)
    ;

}





}  // namespace Python.

} // Namespace Kratos


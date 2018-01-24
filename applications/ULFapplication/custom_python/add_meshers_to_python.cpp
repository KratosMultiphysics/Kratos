/*
==============================================================================
KratosULFApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
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
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2008-10-23 12:50:01 $
//   Revision:            $Revision: 1.10 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_meshers_to_python.h"

#include "external_includes/tetgen_glass.h"
//#include "external_includes/tetgen_mesh_suite_optimized.h"
//#include "external_includes/trigen_mesh_suite.h"
//#include "external_includes/trigen_refine.h"


namespace Kratos
{

namespace Python
{

///////////////////////////////////////////////////////////////////////////////////////////
//											//
//				ADAPTIVE 3D MESHER					//
//											//
//////////////////////////////////////////////////////////////////////////////////////////

//tetgen pfem refine
void TetRegenerateMesh(TetGenGlassModeler& Mesher, char* ElementName, char* ConditionName, ModelPart& model_part, NodeEraseProcess& node_erase,bool rem_nodes, bool add_nodes, double alpha_shape, double h_factor)
{
    Mesher.ReGenerateMesh(model_part,
                          KratosComponents<Element>::Get(ElementName),
                          KratosComponents<Condition>::Get(ConditionName),node_erase,rem_nodes, add_nodes,  alpha_shape, h_factor	);
}


void  AddMeshersToPython()
{

    using namespace boost::python;
    //class that allows 3D adaptive remeshing (inserting and erasing nodes)
    class_<TetGenGlassModeler >("TetGenGlassModeler",
                               init< >())
    .def("ReGenerateMesh",TetRegenerateMesh)
    .def("ReGenerateMesh",&TetGenGlassModeler::ReGenerateMesh)    ;

  

}

}  // namespace Python.

} // Namespace Kratos


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
//   Last modified by:    $Author: jcotela $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "add_modeler_to_python.h"
#include "modeler/modeler.h"
#include "modeler/edge_swapping_2d_modeler.h"
#include "modeler/mpi_connectivity_preserve_modeler.h"
//#include "sources/mpi_connectivity_preserve_modeler.cpp"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

void GenerateModelPart(Modeler& GM, ModelPart& origin_model_part, ModelPart& destination_model_part, const char* ElementName, const char* ConditionName)
{
    GM.GenerateModelPart(origin_model_part, destination_model_part,
                         KratosComponents<Element>::Get(ElementName),
                         KratosComponents<Condition>::Get(ConditionName));

}

void GenerateMesh(Modeler& GM, ModelPart& model_part, const char* ElementName, const char* ConditionName)
{
    GM.GenerateMesh(model_part,
                    KratosComponents<Element>::Get(ElementName),
                    KratosComponents<Condition>::Get(ConditionName));

}


void  AddModelerToPython()
{
    class_<Modeler, Modeler::Pointer, boost::noncopyable>("Modeler")
            .def(init<>())
            .def("GenerateModelPart",&GenerateModelPart)
            .def("GenerateMesh",&GenerateMesh)
            .def("GenerateNodes",&Modeler::GenerateNodes)
    .def(self_ns::str(self))
    ;

    class_<MPIConnectivityPreserveModeler,MPIConnectivityPreserveModeler::Pointer,bases<Modeler>,boost::noncopyable>("MPIConnectivityPreserveModeler")
            ;


    class_< EdgeSwapping2DModeler, EdgeSwapping2DModeler::Pointer, bases<Modeler>, boost::noncopyable  >("EdgeSwapping2DModeler",init< >())
            .def("ReGenerateMesh",&EdgeSwapping2DModeler::Remesh)
    ;
}

}  // namespace Python.

} // Namespace Kratos

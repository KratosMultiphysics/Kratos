// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



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

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


// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
//         -        Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//         -        Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
//                  in the documentation and/or other materials provided with the distribution.
//         -        All advertising materials mentioning features or use of this software must display the following acknowledgement:
//                         This product includes Kratos Multi-Physics technology.
//         -        Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
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
// #include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

// #include "spaces/ublas_space.h"
// #include "linear_solvers/linear_solver.h"

#include "custom_utilities/mapping_utilities.h"

#include "custom_utilities/mapper.h"
#include "custom_utilities/nearest_neighbor_mapper.h"
#include "custom_utilities/iterative_mortar_mapper.h"


namespace Kratos
{

namespace Python
{
    // Wrapper functions for taking the boolean default argument "add_value"
    void Map(Mapper& dummy, const Variable<double>& origin_variable,
                     const Variable<double>& destination_variable) {
        dummy.Map(origin_variable, destination_variable, false);
    }

    void Map(Mapper& dummy, const Variable< array_1d<double,3> >& origin_variable,
                     const Variable< array_1d<double,3> >& destination_variable) {
        dummy.Map(origin_variable, destination_variable, false);
    }

    void InverseMap(Mapper& dummy, const Variable<double>& origin_variable,
                            const Variable<double>& destination_variable) {
        dummy.InverseMap(origin_variable, destination_variable, false);
    }

    void InverseMap(Mapper& dummy, const Variable< array_1d<double,3> >& origin_variable,
                            const Variable< array_1d<double,3> >& destination_variable) {
               dummy.InverseMap(origin_variable, destination_variable, false);
    }

  void  AddCustomUtilitiesToPython()
  {
	using namespace boost::python;

		// typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
		// typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
		// typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
        // class Foo
        // {
        //   ...
        //    void m1(A& a);
        //    void m1(B& b);
        //
        // void (Foo::*m1_a)(A&) = &Foo::m1;
        // void (Foo::*m1_b)(B&) = &Foo::m1;
        //
        // boost::python::class_<Foo>("Foo")
        //     .def("m1", m1_a)
        //     .def("m1", m1_b)

    //     BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Map_overloads, Map, 2, 3)
    //     BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(InverseMap_overloads, InverseMap, 2, 3)
       //
    //     class_<Mapper, boost::noncopyable>("Mapper", no_init)
    //    .def("UpdateInterface",  &Mapper::UpdateInterface)
    //    .def("Map",              &Mapper::Map, Map_overloads())
    //    .def("InverseMap",       &Mapper::InverseMap, InverseMap_overloads())
    //    ;

        void (*MapS1)(Mapper&, const Variable<double> &, const Variable<double> &)                                    = &Map;
        void (*MapV1)(Mapper&, const Variable< array_1d<double,3> > &, const Variable< array_1d<double,3> > &)        = &Map;
        void (*InverseMapS1)(Mapper&, const Variable<double> &, const Variable<double> &)                             = &InverseMap;
        void (*InverseMapV1)(Mapper&, const Variable< array_1d<double,3> > &, const Variable< array_1d<double,3> > &) = &InverseMap;

        void (Mapper::*MapS2)(const Variable<double> &, const Variable<double> &, const bool)                                    = &Mapper::Map;
        void (Mapper::*MapV2)(const Variable< array_1d<double,3> > &, const Variable< array_1d<double,3> > &, const bool)        = &Mapper::Map;
        void (Mapper::*InverseMapS2)(const Variable<double> &, const Variable<double> &, const bool)                             = &Mapper::InverseMap;
        void (Mapper::*InverseMapV2)(const Variable< array_1d<double,3> > &, const Variable< array_1d<double,3> > &, const bool) = &Mapper::InverseMap;

        class_<Mapper, boost::noncopyable>("Mapper", no_init)
       .def("UpdateInterface",  & Mapper::UpdateInterface)
       .def("Map",              MapS1)
       .def("Map",              MapV1)
       .def("InverseMap",       InverseMapS1)
       .def("InverseMap",       InverseMapV1)
       .def("Map",              MapS2)
       .def("Map",              MapV2)
       .def("InverseMap",       InverseMapS2)
       .def("InverseMap",       InverseMapV2)
       ;

    //    class_<Mapper, boost::noncopyable>("Mapper", no_init)
    //   .def("UpdateInterface",  & Mapper::UpdateInterface)
    //   .def("Map",              & Mapper::Map)
    //   .def("InverseMap",       & Mapper::InverseMap)
    //   ;

       class_<  NearestNeighborMapper, bases<Mapper> >("NearestNeighborMapper", init<ModelPart&,ModelPart&>())
       .def(init<ModelPart&,ModelPart&,double,int>()) // second constructor
       // .def("UpdateInterface",     &  CustomNearestNeighborMapperProcess::UpdateInterface)
       // .def("Map",        &  CustomNearestNeighborMapperProcess::Map)
       // .def("InverseMap", &  CustomNearestNeighborMapperProcess::InverseMap)
       ;

       class_<  IterativeMortarMapper, bases<Mapper> >("IterativeMortarMapper", init<ModelPart&,ModelPart&>())
      //  .def(init<ModelPart&,ModelPart&,double,int>()) // second constructor
       // .def("UpdateInterface",     &  CustomNearestNeighborMapperProcess::UpdateInterface)
       // .def("Map",        &  CustomNearestNeighborMapperProcess::Map)
       // .def("InverseMap", &  CustomNearestNeighborMapperProcess::InverseMap)
       ;

       def("ComputeSearchRadius", ComputeSearchRadius);



  }





}  // namespace Python.

} // Namespace Kratos

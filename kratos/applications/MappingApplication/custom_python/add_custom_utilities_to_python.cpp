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

#include "custom_utilities/mapper_flags.h"

#include "custom_utilities/mapper_factory.h"

#include "custom_utilities/mapper.h"
#include "custom_utilities/nearest_neighbor_mapper.h"
// #include "custom_utilities/nearest_element_mapper.h"
// #include "custom_utilities/approximate_mortar_mapper.h"


namespace Kratos
{

namespace Python
{

  // Wrapper functions for taking a default argument for the flags
  void UpdateInterface(MapperFactory& dummy) {
      Kratos::Flags dummy_flags = Kratos::Flags();
      double dummy_search_radius = -1.0f;
      dummy.UpdateInterface(dummy_flags, dummy_search_radius);
  }

  void UpdateInterface(MapperFactory& dummy, Kratos::Flags& options) {
      double dummy_search_radius = -1.0f;
      dummy.UpdateInterface(options, dummy_search_radius);
  }

  void UpdateInterface(MapperFactory& dummy, double search_radius) {
      Kratos::Flags dummy_flags = Kratos::Flags();
      dummy.UpdateInterface(dummy_flags, search_radius);
  }


  void Map(MapperFactory& dummy,
           const Variable<double>& origin_variable,
           const Variable<double>& destination_variable) {
      Kratos::Flags dummy_flags = Kratos::Flags();
      dummy.Map(origin_variable, destination_variable, dummy_flags);
  }

  void Map(MapperFactory& dummy,
           const Variable< array_1d<double,3> >& origin_variable,
           const Variable< array_1d<double,3> >& destination_variable) {
      Kratos::Flags dummy_flags = Kratos::Flags();
      dummy.Map(origin_variable, destination_variable, dummy_flags);
  }

  void InverseMap(MapperFactory& dummy,
                  const Variable<double>& origin_variable,
                  const Variable<double>& destination_variable) {
      Kratos::Flags dummy_flags = Kratos::Flags();
      dummy.InverseMap(origin_variable, destination_variable, dummy_flags);
  }

  void InverseMap(MapperFactory& dummy,
                  const Variable< array_1d<double,3> >& origin_variable,
                  const Variable< array_1d<double,3> >& destination_variable) {
      Kratos::Flags dummy_flags = Kratos::Flags();
      dummy.InverseMap(origin_variable, destination_variable, dummy_flags);
  }

  void  AddCustomUtilitiesToPython()
  {
	using namespace boost::python;

		// typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
		// typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
		// typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

      void (*pUpdateInterface)(MapperFactory &)
                               = &UpdateInterface;

      void (*pUpdateInterfaceOptions)(MapperFactory &, Kratos::Flags &)
                               = &UpdateInterface;

      void (*pUpdateInterfaceSearchRadius)(MapperFactory &, double)
                               = &UpdateInterface;

      void (*pMapScalar)(MapperFactory &,
                         const Variable<double> &,
                         const Variable<double> &)
                         = &Map;

      void (*pMapVector)(MapperFactory &,
                         const Variable< array_1d<double,3> > &,
                         const Variable< array_1d<double,3> > &)
                         = &Map;

      void (*pInverseMapScalar)(MapperFactory &,
                                const Variable<double> &,
                                const Variable<double> &)
                                = &InverseMap;

      void (*pInverseMapVector)(MapperFactory &,
                                const Variable< array_1d<double,3> > &,
                                const Variable< array_1d<double,3> > &)
                                = &InverseMap;


      void (MapperFactory::*pUpdateInterfaceFull)(Kratos::Flags &, double)
                                                  = &MapperFactory::UpdateInterface;

      void (MapperFactory::*pMapScalarOptions)(const Variable<double> &,
                                               const Variable<double> &,
                                               Kratos::Flags &)
                                               = &MapperFactory::Map;

      void (MapperFactory::*pMapVectorOptions)(const Variable< array_1d<double,3> > &,
                                               const Variable< array_1d<double,3> > &,
                                               Kratos::Flags &)
                                               = &MapperFactory::Map;

      void (MapperFactory::*pInverseMapScalarOptions)(const Variable<double> &,
                                                      const Variable<double> &,
                                                      Kratos::Flags &)
                                                      = &MapperFactory::InverseMap;

      void (MapperFactory::*pInverseMapVectorOptions)(const Variable< array_1d<double,3> > &,
                                                      const Variable< array_1d<double,3> > &,
                                                      Kratos::Flags &)
                                                      = &MapperFactory::InverseMap;


      class_< MapperFactory > mapper_factory = class_<MapperFactory>("MapperFactory", init<ModelPart&,ModelPart&,Parameters&>())
      .def("UpdateInterface",  pUpdateInterface)
      .def("UpdateInterface",  pUpdateInterfaceOptions)
      .def("UpdateInterface",  pUpdateInterfaceSearchRadius)
      .def("Map",              pMapScalar)
      .def("Map",              pMapVector)
      .def("InverseMap",       pInverseMapScalar)
      .def("InverseMap",       pInverseMapVector)

      .def("UpdateInterface",  pUpdateInterfaceFull)
      .def("Map",              pMapScalarOptions)
      .def("Map",              pMapVectorOptions)
      .def("InverseMap",       pInverseMapScalarOptions)
      .def("InverseMap",       pInverseMapVectorOptions)
      ;

      mapper_factory.attr("SWAP_SIGN") = MapperFlags::SWAP_SIGN;
      mapper_factory.attr("ADD_VALUES") = MapperFlags::ADD_VALUES;
      mapper_factory.attr("CONSERVATIVE") = MapperFlags::CONSERVATIVE;
    	mapper_factory.attr("REMESHED") = MapperFlags::REMESHED;


      // void (Mapper::*MapScalar)(const Variable<double> &, const Variable<double> &, const bool, const bool, const bool)                             = &Mapper::Map;
      // void (Mapper::*MapVector)(const Variable< array_1d<double,3> > &, const Variable< array_1d<double,3> > &, const bool, const bool, const bool) = &Mapper::Map;
      //
      // class_<Mapper, boost::noncopyable> mapper_interface = class_<Mapper, boost::noncopyable>("Mapper", no_init)
      // .def("UpdateInterface",  & Mapper::UpdateInterface)
      // .def("Map",              MapScalar)
      // .def("Map",              MapVector)
      // ;
      //
      // mapper_interface.attr("SWAP_SIGN") = MapperFlags::SWAP_SIGN;
      // mapper_interface.attr("ADD_VALUES") = MapperFlags::ADD_VALUES;
      // mapper_interface.attr("POINT_WISE_VALUES") = MapperFlags::POINT_WISE_VALUES;
      // mapper_interface.attr("CONSERVATIVE") = MapperFlags::CONSERVATIVE;
    	// mapper_interface.attr("REMESHED") = MapperFlags::REMESHED;
      //
      // // class_< NearestNeighborMapper, bases<Mapper> >("NearestNeighborMapper", init<ModelPart&,ModelPart&,double,int>());
      //
      // // class_< NearestElementMapper, bases<Mapper> >("NearestElementMapper", init<ModelPart&,ModelPart&,double,int>());
      // //
      // // class_< ApproximateMortarMapper, bases<Mapper> >("ApproximateMortarMapper", init<ModelPart&,ModelPart&,double,int,double,int>());
      //
      //
      // class_< MapperUtilities >("MapperUtilities", no_init);
      //
      // def("ComputeNumberOfNodes", ComputeNumberOfNodes);
      // def("ComputeNumberOfConditions", ComputeNumberOfConditions);
      // def("ComputeSearchRadius", ComputeSearchRadius);
  }

}  // namespace Python.

} // Namespace Kratos

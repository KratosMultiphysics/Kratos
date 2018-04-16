//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "spaces/ublas_space.h"

#include "custom_utilities/mapper_flags.h"
#include "custom_python/add_custom_mappers_to_python.h"
#include "custom_utilities/mapper_factory.h"

// Mapper base class
#include "custom_mappers/mapper.h"

// Matrix-free Mappers
#include "custom_mappers/nearest_neighbor_mapper.h"
#include "custom_mappers/nearest_element_mapper.h"


namespace Kratos
{

namespace Python
{

// Wrapper functions for taking a default argument for the flags // TODO inline? Jordi
template<class TSparseSpace, class TDenseSpace>
inline void UpdateInterfaceWithoutArgs(Mapper<TSparseSpace, TDenseSpace>& dummy)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    double dummy_search_radius = -1.0f;
    dummy.UpdateInterface(dummy_flags, dummy_search_radius);
}

template<class TSparseSpace, class TDenseSpace>
inline void UpdateInterfaceWithOptions(Mapper<TSparseSpace, TDenseSpace>& dummy, Kratos::Flags options)
{
    double dummy_search_radius = -1.0f;
    dummy.UpdateInterface(options, dummy_search_radius);
}

template<class TSparseSpace, class TDenseSpace>
inline void UpdateInterfaceWithSearchRadius(Mapper<TSparseSpace, TDenseSpace>& dummy, double search_radius)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.UpdateInterface(dummy_flags, search_radius);
}

template<class TSparseSpace, class TDenseSpace>
inline void MapWithoutOptionsScalar(Mapper<TSparseSpace, TDenseSpace>& dummy,
         const Variable<double>& origin_variable,
         const Variable<double>& destination_variable)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.Map(origin_variable, destination_variable, dummy_flags);
}

template<class TSparseSpace, class TDenseSpace>
inline void MapWithoutOptionsVector(Mapper<TSparseSpace, TDenseSpace>& dummy,
         const Variable< array_1d<double, 3> >& origin_variable,
         const Variable< array_1d<double, 3> >& destination_variable)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.Map(origin_variable, destination_variable, dummy_flags);
}

template<class TSparseSpace, class TDenseSpace>
inline void InverseMapWithoutOptionsScalar(Mapper<TSparseSpace, TDenseSpace>& dummy,
                const Variable<double>& origin_variable,
                const Variable<double>& destination_variable)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.InverseMap(origin_variable, destination_variable, dummy_flags);
}

template<class TSparseSpace, class TDenseSpace>
inline void InverseMapWithoutOptionsVector(Mapper<TSparseSpace, TDenseSpace>& dummy,
                const Variable< array_1d<double, 3> >& origin_variable,
                const Variable< array_1d<double, 3> >& destination_variable)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.InverseMap(origin_variable, destination_variable, dummy_flags);
}

void  AddCustomMappersToPython(pybind11::module& m)
{
    using namespace pybind11;

    // typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    // typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;

    /* void (Mapper::*pMapScalarOptions)(const Variable<double> &,
            const Variable<double> &,
            Kratos::Flags)
        = &Mapper::Map;

    void (Mapper::*pMapVectorOptions)(const Variable< array_1d<double, 3> > &,
            const Variable< array_1d<double, 3> > &,
            Kratos::Flags)
        = &Mapper::Map;

    void (Mapper::*pInverseMapScalarOptions)(const Variable<double> &,
            const Variable<double> &,
            Kratos::Flags)
        = &Mapper::InverseMap;

    void (Mapper::*pInverseMapVectorOptions)(const Variable< array_1d<double, 3> > &,
            const Variable< array_1d<double, 3> > &,
            Kratos::Flags)
        = &Mapper::InverseMap;

    // Exposing the base class of the Mappers to Python, but without constructor
    class_< Mapper, Mapper::Pointer > mapper
        = class_< Mapper, Mapper::Pointer >(m, "Mapper")
            .def("UpdateInterface",  UpdateInterfaceWithoutArgs)
            .def("UpdateInterface",  UpdateInterfaceWithOptions)
            .def("UpdateInterface",  UpdateInterfaceWithSearchRadius)
            .def("Map",              MapWithoutOptionsScalar)
            .def("Map",              MapWithoutOptionsVector)
            .def("InverseMap",       InverseMapWithoutOptionsScalar)
            .def("InverseMap",       InverseMapWithoutOptionsVector)

            .def("UpdateInterface",  &Mapper::UpdateInterface)
            .def("Map",              pMapScalarOptions)
            .def("Map",              pMapVectorOptions)
            .def("InverseMap",       pInverseMapScalarOptions)
            .def("InverseMap",       pInverseMapVectorOptions)

            .def("__repr__",         &Mapper::Info)
            ;

    // Adding the flags that can be used while mapping
    mapper.attr("SWAP_SIGN")        = MapperFlags::SWAP_SIGN;
    mapper.attr("ADD_VALUES")       = MapperFlags::ADD_VALUES;
    mapper.attr("CONSERVATIVE")     = MapperFlags::CONSERVATIVE;
    mapper.attr("REMESHED")         = MapperFlags::REMESHED;

    // Jordi is it possible to expose the mappers without a constructor and use them only through the factory?
    // This would circumvent problems with the wrong space being selected

    // Exposing the Mappers
    // class_< NearestNeighborMapper, bases<Mapper>, boost::noncopyable>
    // ("NearestNeighborMapper", init<ModelPart&, ModelPart&, Parameters, bool>());
    // class_< NearestElementMapper, bases<Mapper>, boost::noncopyable>
    // ("NearestElementMapper", init<ModelPart&, ModelPart&, Parameters, bool>());

    // class_< NearestElementMapper, NearestElementMapper::Pointer, Mapper>
    // (m, "NearestElementMapper")
    //     .def( init<ModelPart&, ModelPart&, Parameters>() );

    // Exposing the MapperFactory
    class_< MapperFactory, MapperFactory::Pointer>(m, "MapperFactory")
        .def_static("CreateMapper", &MapperFactory::CreateMapper); */
}

}  // namespace Python.

} // Namespace Kratos

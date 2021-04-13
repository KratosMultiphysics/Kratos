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
#include "includes/define_python.h"
#include "custom_python/add_custom_mappers_to_python.h"
#include "custom_mappers/mapper.h"
#include "custom_utilities/mapper_factory.h"
#include "custom_utilities/mapper_flags.h"

namespace Kratos {
namespace Python {

// Wrapper functions for taking a default argument for the flags
template<bool TIsDistributed>
inline void UpdateInterfaceWithoutArgs(Mapper<TIsDistributed>& dummy)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    double dummy_search_radius = -1.0f;
    dummy.UpdateInterface(dummy_flags, dummy_search_radius);
}

template<bool TIsDistributed>
inline void UpdateInterfaceWithOptions(Mapper<TIsDistributed>& dummy, Kratos::Flags options)
{
    double dummy_search_radius = -1.0f;
    dummy.UpdateInterface(options, dummy_search_radius);
}

template<bool TIsDistributed>
inline void UpdateInterfaceWithSearchRadius(Mapper<TIsDistributed>& dummy, double search_radius)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.UpdateInterface(dummy_flags, search_radius);
}


template<bool TIsDistributed>
inline void MapWithoutOptionsScalar(Mapper<TIsDistributed>& dummy,
         const Variable<double>& origin_variable,
         const Variable<double>& destination_variable)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.Map(origin_variable, destination_variable, dummy_flags);
}

template<bool TIsDistributed>
inline void MapWithoutOptionsVector(Mapper<TIsDistributed>& dummy,
         const Variable< array_1d<double, 3> >& origin_variable,
         const Variable< array_1d<double, 3> >& destination_variable)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.Map(origin_variable, destination_variable, dummy_flags);
}

template<bool TIsDistributed>
inline void InverseMapWithoutOptionsScalar(Mapper<TIsDistributed>& dummy,
                const Variable<double>& origin_variable,
                const Variable<double>& destination_variable)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.InverseMap(origin_variable, destination_variable, dummy_flags);
}

template<bool TIsDistributed>
inline void InverseMapWithoutOptionsVector(Mapper<TIsDistributed>& dummy,
                const Variable< array_1d<double, 3> >& origin_variable,
                const Variable< array_1d<double, 3> >& destination_variable)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.InverseMap(origin_variable, destination_variable, dummy_flags);
}


template<bool TIsDistributed>
inline void MapWithOptionsScalar(Mapper<TIsDistributed>& dummy,
         const Variable<double>& origin_variable,
         const Variable<double>& destination_variable,
         Kratos::Flags MappingOptions)
{
    dummy.Map(origin_variable, destination_variable, MappingOptions);
}

template<bool TIsDistributed>
inline void MapWithOptionsVector(Mapper<TIsDistributed>& dummy,
         const Variable< array_1d<double, 3> >& origin_variable,
         const Variable< array_1d<double, 3> >& destination_variable,
         Kratos::Flags MappingOptions)
{
    dummy.Map(origin_variable, destination_variable, MappingOptions);
}

template<bool TIsDistributed>
inline void InverseMapWithOptionsScalar(Mapper<TIsDistributed>& dummy,
                const Variable<double>& origin_variable,
                const Variable<double>& destination_variable,
                Kratos::Flags MappingOptions)
{
    dummy.InverseMap(origin_variable, destination_variable, MappingOptions);
}

template<bool TIsDistributed>
inline void InverseMapWithOptionsVector(Mapper<TIsDistributed>& dummy,
                const Variable< array_1d<double, 3> >& origin_variable,
                const Variable< array_1d<double, 3> >& destination_variable,
                Kratos::Flags MappingOptions)
{
    dummy.InverseMap(origin_variable, destination_variable, MappingOptions);
}


template<bool TIsDistributed>
void ExposeMapperToPython(pybind11::module& m, const std::string& rName)
{
    typedef Mapper<TIsDistributed> MapperType;
    namespace py = pybind11;
    // Exposing the base class of the Mappers to Python, but without constructor
    const auto mapper
        = py::class_< MapperType, typename MapperType::Pointer >(m, rName.c_str())
            .def("UpdateInterface",     UpdateInterfaceWithoutArgs<TIsDistributed>)
            .def("UpdateInterface",     UpdateInterfaceWithOptions<TIsDistributed>)
            .def("UpdateInterface",     UpdateInterfaceWithSearchRadius<TIsDistributed>)
            .def("UpdateInterface",     &MapperType::UpdateInterface) // with options & search-radius

            .def("Map",                 MapWithoutOptionsScalar<TIsDistributed>)
            .def("Map",                 MapWithoutOptionsVector<TIsDistributed>)
            .def("Map",                 MapWithOptionsScalar<TIsDistributed>)
            .def("Map",                 MapWithOptionsVector<TIsDistributed>)

            .def("InverseMap",          InverseMapWithoutOptionsScalar<TIsDistributed>)
            .def("InverseMap",          InverseMapWithoutOptionsVector<TIsDistributed>)
            .def("InverseMap",          InverseMapWithOptionsScalar<TIsDistributed>)
            .def("InverseMap",          InverseMapWithOptionsVector<TIsDistributed>)

            .def("GetMappingMatrix",    &MapperType::GetMappingMatrix, py::return_value_policy::reference_internal)
            .def("GetInterfaceModelPartOrigin", &MapperType::GetInterfaceModelPartOrigin, py::return_value_policy::reference_internal)
            .def("GetInterfaceModelPartDestination", &MapperType::GetInterfaceModelPartDestination, py::return_value_policy::reference_internal)

            .def("AreMeshesConforming", &MapperType::AreMeshesConforming)

            .def("__str__",             PrintObject<MapperType>)
            ;

    // Adding the flags that can be used for mapping
    mapper.attr("SWAP_SIGN")           = MapperFlags::SWAP_SIGN;
    mapper.attr("ADD_VALUES")          = MapperFlags::ADD_VALUES;
    mapper.attr("REMESHED")            = MapperFlags::REMESHED;
    mapper.attr("USE_TRANSPOSE")       = MapperFlags::USE_TRANSPOSE;
    mapper.attr("TO_NON_HISTORICAL")   = MapperFlags::TO_NON_HISTORICAL;
    mapper.attr("FROM_NON_HISTORICAL") = MapperFlags::FROM_NON_HISTORICAL;
}

void  AddCustomMappersToPython(pybind11::module& m)
{
    namespace py = pybind11;

    ExposeMapperToPython<false>(m, "Mapper");
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    ExposeMapperToPython<true>(m, "MPIMapper");
#endif

    // Exposing the MapperFactory
    py::class_< MapperFactory, MapperFactory::Pointer>(m, "MapperFactory")
        .def_static("CreateMapper", &MapperFactory::CreateMapper<false>)
        .def_static("HasMapper", &MapperFactory::HasMapper<false>)
        .def_static("GetRegisteredMapperNames", &MapperFactory::GetRegisteredMapperNames<false>)
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
        .def_static("CreateMPIMapper", &MapperFactory::CreateMapper<true>)
        .def_static("HasMPIMapper", &MapperFactory::HasMapper<true>)
        .def_static("GetRegisteredMPIMapperNames", &MapperFactory::GetRegisteredMapperNames<true>)
#endif
    ;
}

}  // namespace Python.
} // Namespace Kratos

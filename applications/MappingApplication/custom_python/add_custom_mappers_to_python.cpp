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
#include "custom_utilities/mapper_typedefs.h"


// Matrix-free Mappers
#include "custom_mappers/nearest_neighbor_mapper.h"
#include "custom_mappers/nearest_element_mapper.h"


namespace Kratos {
namespace Python {

// Wrapper functions for taking a default argument for the flags
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

template<class TSparseSpace, class TDenseSpace>
void (Mapper<TSparseSpace, TDenseSpace>::*pMapScalarOptions)(const Variable<double> &,
        const Variable<double> &,
        Kratos::Flags)
    = &Mapper<TSparseSpace, TDenseSpace>::Map;

template<class TSparseSpace, class TDenseSpace>
void (Mapper<TSparseSpace, TDenseSpace>::*pMapVectorOptions)(const Variable< array_1d<double, 3> > &,
        const Variable< array_1d<double, 3> > &,
        Kratos::Flags)
    = &Mapper<TSparseSpace, TDenseSpace>::Map;

template<class TSparseSpace, class TDenseSpace>
void (Mapper<TSparseSpace, TDenseSpace>::*pInverseMapScalarOptions)(const Variable<double> &,
        const Variable<double> &,
        Kratos::Flags)
    = &Mapper<TSparseSpace, TDenseSpace>::InverseMap;

template<class TSparseSpace, class TDenseSpace>
void (Mapper<TSparseSpace, TDenseSpace>::*pInverseMapVectorOptions)(const Variable< array_1d<double, 3> > &,
        const Variable< array_1d<double, 3> > &,
        Kratos::Flags)
    = &Mapper<TSparseSpace, TDenseSpace>::InverseMap;

template<class TSparseSpace, class TDenseSpace>
void ExposeMapperToPython(pybind11::module& m, const std::string& rName)
{
    typedef Mapper<TSparseSpace, TDenseSpace> MapperType;
    namespace py = pybind11;
    // Exposing the base class of the Mappers to Python, but without constructor
    const auto mapper
        = py::class_< MapperType, typename MapperType::Pointer >(m, rName.c_str())
            .def("UpdateInterface",  UpdateInterfaceWithoutArgs<TSparseSpace, TDenseSpace>)
            .def("UpdateInterface",  UpdateInterfaceWithOptions<TSparseSpace, TDenseSpace>)
            .def("UpdateInterface",  UpdateInterfaceWithSearchRadius<TSparseSpace, TDenseSpace>)
            .def("Map",              MapWithoutOptionsScalar<TSparseSpace, TDenseSpace>)
            .def("Map",              MapWithoutOptionsVector<TSparseSpace, TDenseSpace>)
            .def("InverseMap",       InverseMapWithoutOptionsScalar<TSparseSpace, TDenseSpace>)
            .def("InverseMap",       InverseMapWithoutOptionsVector<TSparseSpace, TDenseSpace>)

            .def("UpdateInterface",  &MapperType::UpdateInterface)
            .def("Map",              pMapScalarOptions<TSparseSpace, TDenseSpace>)
            .def("Map",              pMapVectorOptions<TSparseSpace, TDenseSpace>)
            .def("InverseMap",       pInverseMapScalarOptions<TSparseSpace, TDenseSpace>)
            .def("InverseMap",       pInverseMapVectorOptions<TSparseSpace, TDenseSpace>)

            .def("__str__",          PrintObject<MapperType>)
            ;

    // Adding the flags that can be used for mapping
    mapper.attr("SWAP_SIGN")        = MapperFlags::SWAP_SIGN;
    mapper.attr("ADD_VALUES")       = MapperFlags::ADD_VALUES;
    mapper.attr("REMESHED")         = MapperFlags::REMESHED;
    mapper.attr("USE_TRANSPOSE")    = MapperFlags::REMESHED;
}

void  AddCustomMappersToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef MapperDefinitions::DenseSpaceType DenseSpaceType;
    typedef MapperDefinitions::SparseSpaceType SparseSpaceType;
    ExposeMapperToPython<SparseSpaceType, DenseSpaceType>(m, "Mapper");
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    typedef MapperDefinitions::MPISparseSpaceType MPISparseSpaceType;
    ExposeMapperToPython<MPISparseSpaceType, DenseSpaceType>(m, "MPIMapper");
#endif

    // Exposing the MapperFactory
    py::class_< MapperFactory, MapperFactory::Pointer>(m, "MapperFactory")
        .def_static("CreateMapper", &MapperFactory::CreateMapper<SparseSpaceType, DenseSpaceType>)
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
        .def_static("CreateMPIMapper", &MapperFactory::CreateMapper<MPISparseSpaceType, DenseSpaceType>)
#endif
    ;
}

}  // namespace Python.
} // Namespace Kratos

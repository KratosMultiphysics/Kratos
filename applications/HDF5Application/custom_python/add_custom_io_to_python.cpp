//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/io.h"

// Application includes
#include "custom_io/hdf5_io.h"
#include "custom_io/hdf5_file.h"
#include "custom_io/hdf5_model_part_io.h"
#ifdef KRATOS_USING_MPI
#include "custom_io/hdf5_file_mpi.h"
#include "custom_io/hdf5_partitioned_model_part_io.h"
#endif

namespace Kratos
{

namespace Python
{

void AddCustomIOToPython()
{
    using namespace boost::python;

    class_<HDF5IO, HDF5IO::Pointer, bases<IO>, boost::noncopyable >("HDF5IO", init<std::string, Flags>())
        .def("WriteModelPart",&HDF5IO::WriteModelPart)
    ;

    class_<HDF5File, HDF5File::Pointer, boost::noncopyable >("HDF5File", init<Parameters&>())
    .def("HasPath",&HDF5File::HasPath)
    .def("IsGroup",&HDF5File::IsGroup)
    .def("IsDataSet",&HDF5File::IsDataSet)
    .def("CreateGroup",&HDF5File::CreateGroup)
    .def("AddPath",&HDF5File::AddPath)
    .def("GetDataDimensions",&HDF5File::GetDataDimensions)
    .def("HasIntDataType",&HDF5File::HasIntDataType)
    .def("HasFloatDataType",&HDF5File::HasFloatDataType)
    .def("Flush",&HDF5File::Flush)
    .def("GetFileSize",&HDF5File::GetFileSize)
    .def("GetFileName",&HDF5File::GetFileName)
    ;
    
    class_<HDF5ModelPartIO, HDF5ModelPartIO::Pointer, bases<IO>, boost::noncopyable>(
        "HDF5ModelPartIO", init<Parameters&, HDF5File::Pointer>())
    ;

#ifdef KRATOS_USING_MPI
    class_<HDF5FileMPI, HDF5FileMPI::Pointer, bases<HDF5File>, boost::noncopyable>(
        "HDF5FileMPI", init<Parameters&>());
    class_<HDF5PartitionedModelPartIO, HDF5PartitionedModelPartIO::Pointer, bases<IO>, boost::noncopyable>(
        "HDF5PartitionedModelPartIO", init<Parameters&, HDF5File::Pointer>());
#endif
    
}

} // namespace Python.

} // Namespace Kratos

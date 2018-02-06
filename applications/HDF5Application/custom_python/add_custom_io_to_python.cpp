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
#include "custom_io/hdf5_file.h"
#include "custom_io/hdf5_file_serial.h"
#include "custom_io/hdf5_model_part_io.h"
#include "custom_io/hdf5_nodal_solution_step_data_io.h"
#ifdef KRATOS_USING_MPI
#include "custom_io/hdf5_file_parallel.h"
#include "custom_io/hdf5_partitioned_model_part_io.h"
#endif

namespace Kratos
{

namespace Python
{

void AddCustomIOToPython()
{
    using namespace boost::python;

    class_<HDF5::File, HDF5::File::Pointer, boost::noncopyable >("HDF5File", no_init)
    .def("HasPath",&HDF5::File::HasPath)
    .def("IsGroup",&HDF5::File::IsGroup)
    .def("IsDataSet",&HDF5::File::IsDataSet)
    .def("CreateGroup",&HDF5::File::CreateGroup)
    .def("AddPath",&HDF5::File::AddPath)
    .def("GetDataDimensions",&HDF5::File::GetDataDimensions)
    .def("HasIntDataType",&HDF5::File::HasIntDataType)
    .def("HasFloatDataType",&HDF5::File::HasFloatDataType)
    .def("Flush",&HDF5::File::Flush)
    .def("GetFileSize",&HDF5::File::GetFileSize)
    .def("GetFileName",&HDF5::File::GetFileName)
    ;

    class_<HDF5::FileSerial, HDF5::FileSerial::Pointer, bases<HDF5::File>, boost::noncopyable>(
        "HDF5FileSerial", init<Parameters&>())
    ;
    
    class_<HDF5::ModelPartIO, HDF5::ModelPartIO::Pointer, bases<IO>, boost::noncopyable>(
        "HDF5ModelPartIO", init<Parameters, HDF5::File::Pointer>())
    ;

    class_<HDF5::NodalSolutionStepDataIO, HDF5::NodalSolutionStepDataIO::Pointer, boost::noncopyable>(
        "HDF5NodalSolutionStepDataIO", init<Parameters&, HDF5::File::Pointer>())
        .def("GetPrefix", &HDF5::NodalSolutionStepDataIO::GetPrefix)
        .def("SetPrefix", &HDF5::NodalSolutionStepDataIO::SetPrefix)
        .def("WriteNodalResults", &HDF5::NodalSolutionStepDataIO::WriteNodalResults)
        .def("ReadNodalResults", &HDF5::NodalSolutionStepDataIO::ReadNodalResults)
    ;

#ifdef KRATOS_USING_MPI
    class_<HDF5::FileParallel, HDF5::FileParallel::Pointer, bases<HDF5::File>, boost::noncopyable>(
        "HDF5FileParallel", init<Parameters&>())
    ;

    class_<HDF5::PartitionedModelPartIO, HDF5::PartitionedModelPartIO::Pointer, bases<IO>, boost::noncopyable>(
        "HDF5PartitionedModelPartIO", init<Parameters, HDF5::File::Pointer>())
    ;
#endif
    
}

} // namespace Python.

} // Namespace Kratos

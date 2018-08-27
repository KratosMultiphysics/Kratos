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
#include "pybind11/pybind11.h"


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/io.h"
#include "includes/kratos_parameters.h"
#include "includes/communicator.h"

// Application includes
#include "custom_io/hdf5_file.h"
#include "custom_io/hdf5_file_serial.h"
#include "custom_io/hdf5_model_part_io.h"
#include "custom_io/hdf5_nodal_solution_step_data_io.h"
#include "custom_io/hdf5_element_data_value_io.h"
#include "custom_io/hdf5_nodal_solution_step_bossak_io.h"
#include "custom_io/hdf5_nodal_data_value_io.h"
#include "custom_io/hdf5_data_value_container_io.h"
#ifdef KRATOS_USING_MPI
#include "custom_io/hdf5_file_parallel.h"
#include "custom_io/hdf5_partitioned_model_part_io.h"
#endif

namespace Kratos
{

namespace Python
{

void AddCustomIOToPython(pybind11::module& m)
{
    using namespace pybind11;

    m.def("WriteDataValueContainer", &HDF5::Internals::WriteDataValueContainer, "");
    m.def("ReadDataValueContainer", &HDF5::Internals::ReadDataValueContainer, "");

    class_<HDF5::File, HDF5::File::Pointer >(m,"HDF5File")
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

    class_<HDF5::FileSerial, HDF5::FileSerial::Pointer, HDF5::File>(m,"HDF5FileSerial")
    .def(init<Parameters&>())
    ;

    class_<HDF5::ModelPartIO, HDF5::ModelPartIO::Pointer, IO>(m,"HDF5ModelPartIO")
    .def(init<HDF5::File::Pointer, std::string const&>())
    ;

    class_<HDF5::NodalSolutionStepDataIO, HDF5::NodalSolutionStepDataIO::Pointer>(
        m,"HDF5NodalSolutionStepDataIO")
        .def(init<Parameters, HDF5::File::Pointer>())
        .def("WriteNodalResults", &HDF5::NodalSolutionStepDataIO::WriteNodalResults)
        .def("ReadNodalResults", &HDF5::NodalSolutionStepDataIO::ReadNodalResults)
    ;

    class_<HDF5::NodalSolutionStepBossakIO, HDF5::NodalSolutionStepBossakIO::Pointer>(
        m,"HDF5NodalSolutionStepBossakIO")
        .def(init<Parameters, HDF5::File::Pointer>())
        .def("WriteNodalResults", &HDF5::NodalSolutionStepBossakIO::WriteNodalResults)
        .def("ReadNodalResults", &HDF5::NodalSolutionStepBossakIO::ReadNodalResults)
        .def("SetAlphaBossak", &HDF5::NodalSolutionStepBossakIO::SetAlphaBossak)
    ;

    class_<HDF5::ElementDataValueIO, HDF5::ElementDataValueIO::Pointer>(
        m,"HDF5ElementDataValueIO")
        .def(init<Parameters, HDF5::File::Pointer>())
        .def("WriteElementResults", &HDF5::ElementDataValueIO::WriteElementResults)
        .def("ReadElementResults", &HDF5::ElementDataValueIO::ReadElementResults)
    ;

    class_<HDF5::NodalDataValueIO, HDF5::NodalDataValueIO::Pointer>(
        m,"HDF5NodalDataValueIO")
        .def(init<Parameters, HDF5::File::Pointer>())
        .def("WriteNodalResults", &HDF5::NodalDataValueIO::WriteNodalResults)
        .def("ReadNodalResults", &HDF5::NodalDataValueIO::ReadNodalResults)
    ;

#ifdef KRATOS_USING_MPI
    class_<HDF5::FileParallel, HDF5::FileParallel::Pointer, HDF5::File>(m,"HDF5FileParallel")
    .def(init<Parameters&>())
    ;

    class_<HDF5::PartitionedModelPartIO, HDF5::PartitionedModelPartIO::Pointer, HDF5::ModelPartIO>
        (m,"HDF5PartitionedModelPartIO")
        .def(init<HDF5::File::Pointer, std::string const&>())
    ;
#endif

}

} // namespace Python.

} // Namespace Kratos

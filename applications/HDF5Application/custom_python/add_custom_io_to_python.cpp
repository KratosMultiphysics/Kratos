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
#include "custom_io/hdf5_element_flag_value_io.h"
#include "custom_io/hdf5_condition_data_value_io.h"
#include "custom_io/hdf5_condition_flag_value_io.h"
#include "custom_io/hdf5_nodal_solution_step_bossak_io.h"
#include "custom_io/hdf5_nodal_data_value_io.h"
#include "custom_io/hdf5_nodal_flag_value_io.h"
#include "custom_io/hdf5_data_value_container_io.h"
#ifdef KRATOS_USING_MPI
#include "custom_io/hdf5_file_parallel.h"
#include "custom_io/hdf5_partitioned_model_part_io.h"
#endif

namespace Kratos {
namespace Python {

void AddCustomIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    m.def("WriteDataValueContainer", &HDF5::Internals::WriteDataValueContainer, "");
    m.def("ReadDataValueContainer", &HDF5::Internals::ReadDataValueContainer, "");

    py::class_<HDF5::File, HDF5::File::Pointer >(m,"HDF5File")
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

    py::class_<HDF5::FileSerial, HDF5::FileSerial::Pointer, HDF5::File>(m,"HDF5FileSerial")
        .def(py::init<Parameters&>())
        ;

    py::class_<HDF5::ModelPartIO, HDF5::ModelPartIO::Pointer, IO>(m,"HDF5ModelPartIO")
        .def(py::init<HDF5::File::Pointer, std::string const&>())
        ;

    py::class_<HDF5::NodalSolutionStepDataIO, HDF5::NodalSolutionStepDataIO::Pointer>(
        m,"HDF5NodalSolutionStepDataIO")
        .def(py::init<Parameters, HDF5::File::Pointer>())
        .def("WriteNodalResults", &HDF5::NodalSolutionStepDataIO::WriteNodalResults)
        .def("ReadNodalResults", &HDF5::NodalSolutionStepDataIO::ReadNodalResults)
        ;

    py::class_<HDF5::NodalSolutionStepBossakIO, HDF5::NodalSolutionStepBossakIO::Pointer>(
        m,"HDF5NodalSolutionStepBossakIO")
        .def(py::init<Parameters, HDF5::File::Pointer>())
        .def("WriteNodalResults", &HDF5::NodalSolutionStepBossakIO::WriteNodalResults)
        .def("ReadNodalResults", &HDF5::NodalSolutionStepBossakIO::ReadNodalResults)
        .def("SetAlphaBossak", &HDF5::NodalSolutionStepBossakIO::SetAlphaBossak)
        ;

    py::class_<HDF5::ElementFlagValueIO, HDF5::ElementFlagValueIO::Pointer>(
        m,"HDF5ElementFlagValueIO")
        .def(py::init<Parameters, HDF5::File::Pointer>())
        .def("WriteElementFlags", &HDF5::ElementFlagValueIO::WriteElementFlags)
        .def("ReadElementFlags", &HDF5::ElementFlagValueIO::ReadElementFlags)
        ; 

    py::class_<HDF5::ElementDataValueIO, HDF5::ElementDataValueIO::Pointer>(
        m,"HDF5ElementDataValueIO")
        .def(py::init<Parameters, HDF5::File::Pointer>())
        .def("WriteElementResults", &HDF5::ElementDataValueIO::WriteElementResults)
        .def("ReadElementResults", &HDF5::ElementDataValueIO::ReadElementResults)
        ;

    py::class_<HDF5::ConditionFlagValueIO, HDF5::ConditionFlagValueIO::Pointer>(
        m,"HDF5ConditionFlagValueIO")
        .def(py::init<Parameters, HDF5::File::Pointer>())
        .def("WriteConditionFlags", &HDF5::ConditionFlagValueIO::WriteConditionFlags)
        .def("ReadConditionFlags", &HDF5::ConditionFlagValueIO::ReadConditionFlags)
        ; 

    py::class_<HDF5::ConditionDataValueIO, HDF5::ConditionDataValueIO::Pointer>(
        m,"HDF5ConditionDataValueIO")
        .def(py::init<Parameters, HDF5::File::Pointer>())
        .def("WriteConditionResults", &HDF5::ConditionDataValueIO::WriteConditionResults)
        .def("ReadConditionResults", &HDF5::ConditionDataValueIO::ReadConditionResults)
        ;

    py::class_<HDF5::NodalFlagValueIO, HDF5::NodalFlagValueIO::Pointer>(
        m,"HDF5NodalFlagValueIO")
        .def(py::init<Parameters, HDF5::File::Pointer>())
        .def("WriteNodalFlags", &HDF5::NodalFlagValueIO::WriteNodalFlags)
        .def("ReadNodalFlags", &HDF5::NodalFlagValueIO::ReadNodalFlags)
        ; 

    py::class_<HDF5::NodalDataValueIO, HDF5::NodalDataValueIO::Pointer>(
        m,"HDF5NodalDataValueIO")
        .def(py::init<Parameters, HDF5::File::Pointer>())
        .def("WriteNodalResults", &HDF5::NodalDataValueIO::WriteNodalResults)
        .def("ReadNodalResults", &HDF5::NodalDataValueIO::ReadNodalResults)
        ;

#ifdef KRATOS_USING_MPI
    py::class_<HDF5::FileParallel, HDF5::FileParallel::Pointer, HDF5::File>(m,"HDF5FileParallel")
        .def(py::init<Parameters&>())
        ;

    py::class_<HDF5::PartitionedModelPartIO, HDF5::PartitionedModelPartIO::Pointer, HDF5::ModelPartIO>
        (m,"HDF5PartitionedModelPartIO")
        .def(py::init<HDF5::File::Pointer, std::string const&>())
        ;
#endif

}

} // namespace Python.

} // Namespace Kratos

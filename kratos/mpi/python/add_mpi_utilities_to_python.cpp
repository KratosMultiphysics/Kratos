//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/parallel_environment.h"
#include "add_mpi_utilities_to_python.h"
#include "mpi/utilities/model_part_communicator_utilities.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "mpi/utilities/data_communicator_factory.h"
#include "mpi/utilities/gather_modelpart_utility.h"
#include "mpi/utilities/mpi_normal_calculation_utilities.h"
#include "mpi/utilities/distributed_model_part_initializer.h"

namespace Kratos::Python {

const DataCommunicator& CreateFromListOfRanks(
    const DataCommunicator& rReferenceComm,
    pybind11::list& rRanks,
    const std::string& rNewName)
{
    namespace py = pybind11;

    std::size_t list_length = py::len(rRanks);
    std::vector<int> rank_vector(list_length);
    for (std::size_t i = 0; i < list_length; i++)
    {
        rank_vector[i] = py::int_(rRanks[i]);
    }
    return DataCommunicatorFactory::CreateFromRanksAndRegister(rReferenceComm, rank_vector, rNewName);
}

void AddMPIUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<ModelPartCommunicatorUtilities>(m,"ModelPartCommunicatorUtilities")
    .def_static("SetMPICommunicator", [](ModelPart& rModelPart, const DataCommunicator& rDataCommunicator){
        ModelPartCommunicatorUtilities::SetMPICommunicator(rModelPart, rDataCommunicator);
    })
    .def("SetMPICommunicatorRecursively", &ModelPartCommunicatorUtilities::SetMPICommunicatorRecursively)
    .def_static("SetMPICommunicator", [](ModelPart& rModelPart){
        KRATOS_WARNING("ModelPartCommunicatorUtilities") << "This function is deprecated, please use the one that accepts a DataCommunicator" << std::endl;
        // passing the default data comm as a temp solution until this function is removed to avoid the deprecation warning in the interface
        ModelPartCommunicatorUtilities::SetMPICommunicator(rModelPart, ParallelEnvironment::GetDefaultDataCommunicator());
    })
    ;

    py::class_<ParallelFillCommunicator, ParallelFillCommunicator::Pointer, FillCommunicator>(m,"ParallelFillCommunicator")
    .def(py::init([](ModelPart& rModelPart){
        KRATOS_WARNING("ParallelFillCommunicator") << "Using deprecated constructor. Please use constructor with data communicator!";
        return Kratos::make_shared<ParallelFillCommunicator>(rModelPart, ParallelEnvironment::GetDefaultDataCommunicator());
    }) )
    .def(py::init<ModelPart&, const DataCommunicator& >() )
    ;

    m.def_submodule("DataCommunicatorFactory")
    .def("DuplicateAndRegister", &DataCommunicatorFactory::DuplicateAndRegister, py::return_value_policy::reference)
    .def("SplitAndRegister", &DataCommunicatorFactory::SplitAndRegister, py::return_value_policy::reference)
    .def("CreateFromRanksAndRegister", &CreateFromListOfRanks, py::return_value_policy::reference)
    .def("CreateUnionAndRegister", &DataCommunicatorFactory::CreateUnionAndRegister, py::return_value_policy::reference)
    .def("CreateIntersectionAndRegister", &DataCommunicatorFactory::CreateIntersectionAndRegister, py::return_value_policy::reference)
    ;

    py::class_<GatherModelPartUtility>(m, "GatherModelPartUtility")
        .def(py::init<int, ModelPart&, ModelPart&>())
        .def(py::init<int, ModelPart&, int, ModelPart&>())
        .def("GatherOnMaster", &GatherModelPartUtility::GatherOnMaster<double>)
        .def("GatherOnMaster", &GatherModelPartUtility::GatherOnMaster<array_1d<double, 3>>)
        .def("ScatterFromMaster", &GatherModelPartUtility::ScatterFromMaster<double>)
        .def("ScatterFromMaster",
             &GatherModelPartUtility::ScatterFromMaster<array_1d<double, 3>>);

    py::class_<MPINormalCalculationUtils, MPINormalCalculationUtils::Pointer>(m,"MPINormalCalculationUtils")
        .def(py::init<>())
        .def("Check",&MPINormalCalculationUtils::Check)
        .def("OrientFaces",&MPINormalCalculationUtils::OrientFaces)
        .def("CalculateOnSimplex",&MPINormalCalculationUtils::CalculateOnSimplex)
        ;

    py::class_<DistributedModelPartInitializer>(m, "DistributedModelPartInitializer")
        .def(py::init<ModelPart&, const DataCommunicator&, int>())
        .def("CopySubModelPartStructure", &DistributedModelPartInitializer::CopySubModelPartStructure)
        .def("Execute", &DistributedModelPartInitializer::Execute)
        ;

}

} // namespace Kratos::Python


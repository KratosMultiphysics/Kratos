//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

// System includes
#include <pybind11/stl.h>

// External includes

// Project includes
#include "custom_utilities/mpi_utilities.h"
#include "custom_utilities/mpm_mpi_search.h"
#include "includes/model_part.h"

namespace Kratos::Python {

    void SearchElementMPIAccordingToDimension(
            ModelPart& rBackgroundGridModelPart,
            ModelPart& rMPMModelPart,
            const std::size_t MaxNumberOfResults,
            const double Tolerance)
    {
        const auto dimension = rBackgroundGridModelPart.GetProcessInfo()[DOMAIN_SIZE];
        if(dimension == 2)
            MPM_MPI_SEARCH<2>::SearchElementMPI(rBackgroundGridModelPart, rMPMModelPart, MaxNumberOfResults, Tolerance);
        else if(dimension == 3)
            MPM_MPI_SEARCH<3>::SearchElementMPI(rBackgroundGridModelPart, rMPMModelPart, MaxNumberOfResults, Tolerance);
    }

    void  AddCustomMPIUtilitiesToPython(pybind11::module& m)
    {
        m.def("SearchElementMPI",SearchElementMPIAccordingToDimension);

        namespace py = pybind11;
        py::class_<MPM_MPI_Utilities, MPM_MPI_Utilities::Pointer>(m, "MPM_MPI_Utilities")
            .def_static("TransferElements", &MPM_MPI_Utilities::TransferElements)
            .def_static("TransferConditions", &MPM_MPI_Utilities::TransferConditions)
            .def_static("SetMPICommunicator",&MPM_MPI_Utilities::SetMPICommunicator)
            .def_static("ClearLocalElementsFromCommunicator", &MPM_MPI_Utilities::ClearLocalElementsFromCommunicator)
            .def_static("WriteGlobalParticlesToFile", &MPM_MPI_Utilities::WriteGlobalParticlesToFile );
    }

} // Namespace Kratos::Python

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:          BSD License
//  Original author:  Riccardo Rossi, Ruben Zorrilla
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "spaces/ublas_space.h"
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/FSI_utils.h"
#include "custom_utilities/aitken_utils.h"
#include "custom_utilities/partitioned_fsi_utilities.hpp"

namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;
    typedef UblasSpace<double, Matrix, Vector > TSpace;

    class_<FSIUtils>("FSIUtils", init<>())
//		.def("FSIUtils",&FSIUtils::GenerateCouplingElements)
        .def("CheckPressureConvergence",&FSIUtils::CheckPressureConvergence)
        .def("StructuralPressurePrediction",&FSIUtils::StructuralPressurePrediction)
        ;

    class_<AitkenUtils>("AitkenUtils", init<>())
        .def("ComputeAitkenFactor",&AitkenUtils::ComputeAitkenFactor)
        .def("ComputeRelaxedDisplacement",&AitkenUtils::ComputeRelaxedDisplacement)
        ;

    class_<PartitionedFSIUtilities <TSpace,2>, boost::noncopyable >("PartitionedFSIUtilities2D", init< >())
        .def("GetInterfaceArea",&PartitionedFSIUtilities<TSpace,2>::GetInterfaceArea)
        .def("GetInterfaceResidualSize",&PartitionedFSIUtilities<TSpace,2>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues",&PartitionedFSIUtilities<TSpace,2>::UpdateInterfaceValues)
        .def("ComputeInterfaceVectorResidual",&PartitionedFSIUtilities<TSpace,2>::ComputeInterfaceVectorResidual)
        .def("ComputeFluidInterfaceMeshVelocityResidualNorm",&PartitionedFSIUtilities<TSpace,2>::ComputeFluidInterfaceMeshVelocityResidualNorm)
        .def("ComputeCorrectedInterfaceDisplacementDerivatives",&PartitionedFSIUtilities<TSpace,2>::ComputeCorrectedInterfaceDisplacementDerivatives)
        ;

    class_<PartitionedFSIUtilities <TSpace,3>, boost::noncopyable >("PartitionedFSIUtilities3D", init< >())
        .def("GetInterfaceArea",&PartitionedFSIUtilities<TSpace,3>::GetInterfaceArea)
        .def("GetInterfaceResidualSize",&PartitionedFSIUtilities<TSpace,3>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues",&PartitionedFSIUtilities<TSpace,3>::UpdateInterfaceValues)
        .def("ComputeInterfaceVectorResidual",&PartitionedFSIUtilities<TSpace,3>::ComputeInterfaceVectorResidual)
        .def("ComputeFluidInterfaceMeshVelocityResidualNorm",&PartitionedFSIUtilities<TSpace,3>::ComputeFluidInterfaceMeshVelocityResidualNorm)
        .def("ComputeCorrectedInterfaceDisplacementDerivatives",&PartitionedFSIUtilities<TSpace,3>::ComputeCorrectedInterfaceDisplacementDerivatives)
        ;

}

}  // namespace Python.

} // Namespace Kratos

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi, Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "spaces/ublas_space.h"
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/FSI_utils.h"
#include "custom_utilities/aitken_utils.h"
#include "custom_utilities/partitioned_fsi_utilities.hpp"
#include "custom_utilities/nodal_update_utilities.h"
#include "custom_utilities/variable_redistribution_utility.h"

namespace Kratos
{

namespace Python
{

void AddCustomUtilitiesToPython(pybind11::module &m)
{
    using namespace pybind11;

    typedef UblasSpace<double, Matrix, Vector > TSpace;
    typedef NodalUpdateBaseClass< 2 > NodalUpdateBaseClass2DType;
    typedef NodalUpdateBaseClass< 3 > NodalUpdateBaseClass3DType;

    class_<FSIUtils>(m,"FSIUtils")
        .def(init<>())
        .def("CheckPressureConvergence",&FSIUtils::CheckPressureConvergence)
        .def("StructuralPressurePrediction",&FSIUtils::StructuralPressurePrediction)
        ;

    class_<AitkenUtils>(m,"AitkenUtils")
        .def(init<>())
        .def("ComputeAitkenFactor",&AitkenUtils::ComputeAitkenFactor)
        .def("ComputeRelaxedDisplacement",&AitkenUtils::ComputeRelaxedDisplacement)
        ;

    class_<PartitionedFSIUtilities<TSpace,2>>(m,"PartitionedFSIUtilities2D")
        .def(init<>())
        .def("GetInterfaceArea",&PartitionedFSIUtilities<TSpace,2>::GetInterfaceArea)
        .def("GetInterfaceResidualSize",&PartitionedFSIUtilities<TSpace,2>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues",&PartitionedFSIUtilities<TSpace,2>::UpdateInterfaceValues)
        .def("ComputeInterfaceVectorResidual",&PartitionedFSIUtilities<TSpace,2>::ComputeInterfaceVectorResidual)
        .def("ComputeFluidInterfaceMeshVelocityResidualNorm",&PartitionedFSIUtilities<TSpace,2>::ComputeFluidInterfaceMeshVelocityResidualNorm)
        .def("ComputeAndPrintFluidInterfaceNorms",&PartitionedFSIUtilities<TSpace,2>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms",&PartitionedFSIUtilities<TSpace,2>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid",&PartitionedFSIUtilities<TSpace,2>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure",&PartitionedFSIUtilities<TSpace,2>::CheckCurrentCoordinatesStructure)
        ;

    class_<PartitionedFSIUtilities<TSpace, 3>>(m,"PartitionedFSIUtilities3D")
        .def(init<>())
        .def("GetInterfaceArea", &PartitionedFSIUtilities<TSpace, 3>::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &PartitionedFSIUtilities<TSpace, 3>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues", &PartitionedFSIUtilities<TSpace, 3>::UpdateInterfaceValues)
        .def("ComputeInterfaceVectorResidual", &PartitionedFSIUtilities<TSpace, 3>::ComputeInterfaceVectorResidual)
        .def("ComputeFluidInterfaceMeshVelocityResidualNorm", &PartitionedFSIUtilities<TSpace, 3>::ComputeFluidInterfaceMeshVelocityResidualNorm)
        .def("ComputeAndPrintFluidInterfaceNorms", &PartitionedFSIUtilities<TSpace, 3>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &PartitionedFSIUtilities<TSpace, 3>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &PartitionedFSIUtilities<TSpace, 3>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &PartitionedFSIUtilities<TSpace, 3>::CheckCurrentCoordinatesStructure);

    class_<NodalUpdateBaseClass<2>>(m,"BaseNodalUpdate2D")
        .def(init<>())
        .def("UpdateMeshTimeDerivatives", &NodalUpdateBaseClass<2>::UpdateMeshTimeDerivatives)
        .def("SetMeshTimeDerivativesOnInterface", &NodalUpdateBaseClass<2>::SetMeshTimeDerivativesOnInterface);

    class_<NodalUpdateBaseClass<3>>(m,"BaseNodalUpdate3D")
        .def(init<>())
        .def("UpdateMeshTimeDerivatives", &NodalUpdateBaseClass<3>::UpdateMeshTimeDerivatives)
        .def("SetMeshTimeDerivativesOnInterface", &NodalUpdateBaseClass<3>::SetMeshTimeDerivativesOnInterface);

    class_<NodalUpdateNewmark<2>, NodalUpdateBaseClass2DType>(m,"NodalUpdateNewmark2D")
        .def(init<const double>())
        .def("UpdateMeshTimeDerivatives", &NodalUpdateNewmark<2>::UpdateMeshTimeDerivatives)
        .def("SetMeshTimeDerivativesOnInterface", &NodalUpdateNewmark<2>::SetMeshTimeDerivativesOnInterface);

    class_<NodalUpdateNewmark<3>, NodalUpdateBaseClass3DType>(m,"NodalUpdateNewmark3D")
        .def(init<const double>())
        .def("UpdateMeshTimeDerivatives", &NodalUpdateNewmark<3>::UpdateMeshTimeDerivatives)
        .def("SetMeshTimeDerivativesOnInterface", &NodalUpdateNewmark<3>::SetMeshTimeDerivativesOnInterface);

    typedef void (*DistributePointDoubleType)(ModelPart&, const Variable< double >&, const Variable< double >&, double, unsigned int);
    typedef void (*DistributePointArrayType)(ModelPart&, const Variable< array_1d<double,3> >&, const Variable< array_1d<double,3> >&,double, unsigned int);

    DistributePointDoubleType DistributePointDouble = &VariableRedistributionUtility::DistributePointValues;
    DistributePointArrayType  DistributePointArray  = &VariableRedistributionUtility::DistributePointValues;

    typedef void (*ConvertDistributedDoubleType)(ModelPart&, const Variable< double >&, const Variable< double >&);
    typedef void (*ConvertDistributedArrayType)(ModelPart&, const Variable< array_1d<double,3> >&, const Variable< array_1d<double,3> >&);

    ConvertDistributedDoubleType ConvertDistributedDouble = &VariableRedistributionUtility::ConvertDistributedValuesToPoint;
    ConvertDistributedArrayType  ConvertDistributedArray  = &VariableRedistributionUtility::ConvertDistributedValuesToPoint;

    // Note: The StaticMethod thing should be done only once for each set of overloads
    class_< VariableRedistributionUtility >(m,"VariableRedistributionUtility")
    .def_static("DistributePointValues",DistributePointDouble)
    .def_static("DistributePointValues",DistributePointArray)
    .def_static("ConvertDistributedValuesToPoint",ConvertDistributedDouble)
    .def_static("ConvertDistributedValuesToPoint",ConvertDistributedArray)
    ;
}

}  // namespace Python.

} // Namespace Kratos

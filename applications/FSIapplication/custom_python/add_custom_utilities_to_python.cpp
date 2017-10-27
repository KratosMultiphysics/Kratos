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
#include "custom_utilities/nodal_update_utilities.h"
#include "custom_utilities/variable_redistribution_utility.h"

namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;
    typedef UblasSpace<double, Matrix, Vector > TSpace;
    typedef NodalUpdateBaseClass< 2 > NodalUpdateBaseClass2DType;
    typedef NodalUpdateBaseClass< 3 > NodalUpdateBaseClass3DType;

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
        .def("ComputeAndPrintFluidInterfaceNorms",&PartitionedFSIUtilities<TSpace,2>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms",&PartitionedFSIUtilities<TSpace,2>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid",&PartitionedFSIUtilities<TSpace,2>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure",&PartitionedFSIUtilities<TSpace,2>::CheckCurrentCoordinatesStructure)
        ;

    class_<PartitionedFSIUtilities <TSpace,3>, boost::noncopyable >("PartitionedFSIUtilities3D", init< >())
        .def("GetInterfaceArea",&PartitionedFSIUtilities<TSpace,3>::GetInterfaceArea)
        .def("GetInterfaceResidualSize",&PartitionedFSIUtilities<TSpace,3>::GetInterfaceResidualSize)
        .def("UpdateInterfaceValues",&PartitionedFSIUtilities<TSpace,3>::UpdateInterfaceValues)
        .def("ComputeInterfaceVectorResidual",&PartitionedFSIUtilities<TSpace,3>::ComputeInterfaceVectorResidual)
        .def("ComputeFluidInterfaceMeshVelocityResidualNorm",&PartitionedFSIUtilities<TSpace,3>::ComputeFluidInterfaceMeshVelocityResidualNorm)
        .def("ComputeAndPrintFluidInterfaceNorms",&PartitionedFSIUtilities<TSpace,3>::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms",&PartitionedFSIUtilities<TSpace,3>::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid",&PartitionedFSIUtilities<TSpace,3>::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure",&PartitionedFSIUtilities<TSpace,3>::CheckCurrentCoordinatesStructure)
        ;

    class_<NodalUpdateBaseClass<2>, boost::noncopyable>("BaseNodalUpdate2D", init<>())
        .def("UpdateMeshTimeDerivatives", &NodalUpdateBaseClass<2>::UpdateMeshTimeDerivatives)
        .def("SetMeshTimeDerivativesOnInterface", &NodalUpdateBaseClass<2>::SetMeshTimeDerivativesOnInterface);

    class_<NodalUpdateBaseClass<3>, boost::noncopyable>("BaseNodalUpdate3D", init<>())
        .def("UpdateMeshTimeDerivatives", &NodalUpdateBaseClass<3>::UpdateMeshTimeDerivatives)
        .def("SetMeshTimeDerivativesOnInterface", &NodalUpdateBaseClass<3>::SetMeshTimeDerivativesOnInterface);

    class_<NodalUpdateNewmark<2>, bases<NodalUpdateBaseClass2DType>, boost::noncopyable>("NodalUpdateNewmark2D", init<const double>())
        .def("UpdateMeshTimeDerivatives", &NodalUpdateNewmark<2>::UpdateMeshTimeDerivatives)
        .def("SetMeshTimeDerivativesOnInterface", &NodalUpdateNewmark<2>::SetMeshTimeDerivativesOnInterface);

    class_<NodalUpdateNewmark<3>, bases<NodalUpdateBaseClass3DType>, boost::noncopyable>("NodalUpdateNewmark3D", init<const double>())
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
    class_< VariableRedistributionUtility, boost::noncopyable >("VariableRedistributionUtility", no_init)
    .def("DistributePointValues",DistributePointDouble)
    .def("DistributePointValues",DistributePointArray).staticmethod("DistributePointValues")
    .def("ConvertDistributedValuesToPoint",ConvertDistributedDouble)
    .def("ConvertDistributedValuesToPoint",ConvertDistributedArray).staticmethod("ConvertDistributedValuesToPoint")
    ;
}

}  // namespace Python.

} // Namespace Kratos

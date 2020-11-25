// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 license: CoSimulationApplication/license.txt
//
//  Main authors:
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "spaces/ublas_space.h"
#include "custom_utilities/feti_dynamic_coupling_utilities.h"

namespace Kratos{

namespace Python{

    void  AddCustomUtilitiesToPython(pybind11::module& m)
    {
        typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

        typedef FetiDynamicCouplingUtilities<SparseSpaceType, LocalSpaceType> FetiDynamicCouplingUtilitiesType;
        typedef FetiDynamicCouplingUtilitiesType::SolverIndex FetiSolverIndexType;

        pybind11::class_< FetiDynamicCouplingUtilitiesType>(m, "FetiDynamicCouplingUtilities")
            .def(pybind11::init<ModelPart&, ModelPart&, Parameters>())
            .def("SetOriginAndDestinationDomainsWithInterfaceModelParts",
                &FetiDynamicCouplingUtilitiesType::SetOriginAndDestinationDomainsWithInterfaceModelParts)
            .def("SetEffectiveStiffnessMatrixImplicit",
                &FetiDynamicCouplingUtilitiesType::SetEffectiveStiffnessMatrixImplicit)
            .def("SetEffectiveStiffnessMatrixExplicit",
                &FetiDynamicCouplingUtilitiesType::SetEffectiveStiffnessMatrixExplicit)
            .def("EquilibrateDomains",
                &FetiDynamicCouplingUtilitiesType::EquilibrateDomains)
            .def("SetOriginInitialKinematics",
                &FetiDynamicCouplingUtilitiesType::SetOriginInitialKinematics)
            .def("SetMappingMatrix",
                &FetiDynamicCouplingUtilitiesType::SetMappingMatrix)
            .def("SetLinearSolver",
                &FetiDynamicCouplingUtilitiesType::SetLinearSolver)
            ;

        pybind11::enum_< FetiSolverIndexType>(m, "FetiSolverIndexType")
            .value("Origin", FetiSolverIndexType::Origin)
            .value("Destination", FetiSolverIndexType::Destination)
            ;
    }

}  // namespace Python.
} // Namespace Kratos

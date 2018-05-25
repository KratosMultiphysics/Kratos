//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/NurbsBrepApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "custom_utilities/NurbsBrepModeler.h"
#include "custom_utilities/BrepModelGeometryReader.h"

namespace Kratos
{

namespace Python
{

	void  AddCustomUtilitiesToPython(pybind11::module& m)
	{
		using namespace pybind11;


		typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
		typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
		typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;


		class_<NurbsBrepModeler, typename NurbsBrepModeler::Pointer>(m, "NurbsBrepModeler")
			.def(init<ModelPart::Pointer>())
			.def("LoadGeometry", &NurbsBrepModeler::LoadGeometry)
			.def("CreateIntegrationDomain", &NurbsBrepModeler::CreateIntegrationDomain)
			.def("ApplyGeometryRefinement", &NurbsBrepModeler::ApplyGeometryRefinement)
			.def("ComputeArea", &NurbsBrepModeler::ComputeArea)
			.def("MapNode", &NurbsBrepModeler::MapNode)
			.def("GetInterfaceConditions", &NurbsBrepModeler::GetInterfaceConditions)
			;

		class_<BrepModelGeometryReader, typename BrepModelGeometryReader::Pointer>(m, "BrepModelGeometryReader")
			.def(init<Parameters&>())
			.def("ReadGeometry", &BrepModelGeometryReader::ReadGeometry)
			.def("WriteGaussPoints", &BrepModelGeometryReader::WriteGaussPoints)
			.def("WriteGaussPointsJson", &BrepModelGeometryReader::WriteGaussPointsJson)
			;
	}

}  // namespace Python.

} // Namespace Kratos

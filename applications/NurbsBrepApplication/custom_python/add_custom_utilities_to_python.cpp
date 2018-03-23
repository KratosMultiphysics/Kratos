//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/NurbBrepApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Riccardo Rossi
//

// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "custom_utilities/brep/BrepModel.h"
#include "custom_utilities/NurbsBrepModeler.h"
#include "custom_utilities/BrepModelGeometryReader.h"

namespace Kratos
{

namespace Python
{


	void  AddCustomUtilitiesToPython()
	{
		using namespace boost::python;


		typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
		typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
		typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

		typedef std::vector<BrepModel> BrepModelVector;
		typedef std::vector<BrepFace>  BrepFacesVector;
		typedef std::vector<BrepEdge>  BrepEdgesVector;

		class_<NurbsBrepModeler, boost::noncopyable>("NurbsBrepModeler", init<ModelPart&>())
			.def("LoadGeometry", &NurbsBrepModeler::LoadGeometry)
			.def("CreateIntegrationDomain", &NurbsBrepModeler::CreateIntegrationDomain)
			.def("ApplyGeometryRefinement", &NurbsBrepModeler::ApplyGeometryRefinement)
			.def("ComputeArea", &NurbsBrepModeler::ComputeArea)
			.def("MapNode", &NurbsBrepModeler::MapNode)
			;

		class_<BrepModelGeometryReader, boost::noncopyable>("BrepModelGeometryReader", init<Parameters&>())
			.def("ReadGeometry", &BrepModelGeometryReader::ReadGeometry)
			.def("WriteGaussPoints", &BrepModelGeometryReader::WriteGaussPoints)
			;
	}

}  // namespace Python.

} // Namespace Kratos

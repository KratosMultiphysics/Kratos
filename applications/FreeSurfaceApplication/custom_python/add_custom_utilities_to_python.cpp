//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Antonia Larese
//

// System includes

// External includes

// Project includes
#include "add_custom_utilities_to_python.h"

namespace Kratos
{

	namespace Python
	{

		void AddCustomUtilitiesToPython(pybind11::module &pymodule)
		{
			namespace py = pybind11;

			typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;

			py::class_<MatrixContainer<2, SparseSpaceType>>(pymodule, "MatrixContainer2D")
				.def(py::init<>())
				.def("ConstructCSRVector", &MatrixContainer<2, SparseSpaceType>::ConstructCSRVector)
				.def("BuildCSRData", &MatrixContainer<2, SparseSpaceType>::BuildCSRData)
				.def("Clear", &MatrixContainer<2, SparseSpaceType>::Clear);

			py::class_<MatrixContainer<3, SparseSpaceType>>(pymodule, "MatrixContainer3D")
				.def(py::init<>())
				.def("ConstructCSRVector", &MatrixContainer<3, SparseSpaceType>::ConstructCSRVector)
				.def("BuildCSRData", &MatrixContainer<3, SparseSpaceType>::BuildCSRData)
				.def("Clear", &MatrixContainer<3, SparseSpaceType>::Clear);

            py::class_<MatrixContainerC2C<2, SparseSpaceType>>(pymodule, "MatrixContainerC2C2D")
                .def(py::init<>())
                .def("ConstructCSRVector", &MatrixContainerC2C<2, SparseSpaceType>::ConstructCSRVector)
                .def("BuildCSRData", &MatrixContainerC2C<2, SparseSpaceType>::BuildCSRData)
                .def("Clear", &MatrixContainerC2C<2, SparseSpaceType>::Clear);

            py::class_<MatrixContainerC2C<3, SparseSpaceType>>(pymodule, "MatrixContainerC2C3D")
                .def(py::init<>())
                .def("ConstructCSRVector", &MatrixContainerC2C<3, SparseSpaceType>::ConstructCSRVector)
                .def("BuildCSRData", &MatrixContainerC2C<3, SparseSpaceType>::BuildCSRData)
                .def("Clear", &MatrixContainerC2C<3, SparseSpaceType>::Clear);

			py::class_<EdgebasedLevelsetAuxiliaryUtils <2>>(pymodule, "EdgebasedLevelsetAuxiliaryUtils2D")
				.def(py::init<>())
				.def("CalculateDistances", &EdgebasedLevelsetAuxiliaryUtils <2> ::CalculateDistances)
				.def("FindMaximumEdgeSize", &EdgebasedLevelsetAuxiliaryUtils <2> ::FindMaximumEdgeSize);

			py::class_<EdgebasedLevelsetAuxiliaryUtils <3>>(pymodule,"EdgebasedLevelsetAuxiliaryUtils3D")
				.def(py::init<>())
				.def("CalculateDistances", &EdgebasedLevelsetAuxiliaryUtils <3> ::CalculateDistances)
				.def("FindMaximumEdgeSize", &EdgebasedLevelsetAuxiliaryUtils <3> ::FindMaximumEdgeSize);
		}

	} // namespace Python

} // Namespace Kratos

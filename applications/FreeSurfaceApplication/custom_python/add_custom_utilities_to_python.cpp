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

        void  AddCustomUtilitiesToPython(pybind11::module& pymodule)
        {
	        namespace py = pybind11;

		    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;

    	    py::class_< MatrixContainer < 2, SparseSpaceType> > (pymodule,"MatrixContainer2D")
            .def(py::init< >())
    	    .def("ConstructCSRVector", &MatrixContainer < 2, SparseSpaceType >::ConstructCSRVector)
    	    .def("BuildCSRData", &MatrixContainer < 2, SparseSpaceType >::BuildCSRData)
    	    .def("Clear", &MatrixContainer < 2, SparseSpaceType >::Clear)
    	    ;

    	    py::class_< MatrixContainer < 3, SparseSpaceType> > (pymodule,"MatrixContainer3D")
            .def(py::init< >())
    	    .def("ConstructCSRVector", &MatrixContainer < 3, SparseSpaceType >::ConstructCSRVector)
    	    .def("BuildCSRData", &MatrixContainer < 3, SparseSpaceType >::BuildCSRData)
    	    .def("Clear", &MatrixContainer < 3, SparseSpaceType >::Clear)
    	    ;
        }

    }  // namespace Python

} // Namespace Kratos

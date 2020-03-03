//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Raul Bravo
//
//


// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/rom_residuals_utility.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos {
namespace Python {


void convert_to_numpy(RomResidualsUtility& RomResidualsUtilityObject, const Matrix & KratosMatrix, pybind11::object NumpyMatrix)
{
    PyObject* pobj = NumpyMatrix.ptr();
    Py_buffer pybuf;
    PyObject_GetBuffer(pobj, &pybuf, PyBUF_SIMPLE);
    void *buf = pybuf.buf;
    double *p = (double*)buf;
    Py_XDECREF(pobj);

    unsigned int n_rows = KratosMatrix.size1();
    unsigned int n_cols = KratosMatrix.size2();
    for (unsigned int i = 0; i < n_rows; i++)
    {
        for (unsigned int j = 0; j < n_cols; j++)
        {
            p[i*n_cols+j] = KratosMatrix(i,j);
        }
    }
}


using namespace pybind11;

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    class_<RomResidualsUtility, typename RomResidualsUtility::Pointer>(m, "RomResidualsUtility")
    .def(init<ModelPart&, Parameters, BaseSchemeType::Pointer>()) // 
    .def("GetResiduals",&RomResidualsUtility::Calculate) //
    .def("Kratos2Numpy",convert_to_numpy)
    ;     

}

} // namespace Python.
} // Namespace Kratos

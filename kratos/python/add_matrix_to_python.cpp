//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//


// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "add_matrix_to_python.h"
// #include "python/matrix_python_interface.h"
// #include "python/matrix_scalar_operator_python.h"
// #include "python/matrix_scalar_assignment_operator_python.h"
// #include "python/matrix_matrix_operator_python.h"

namespace Kratos
{
namespace Python
{

using namespace pybind11;

template< typename TMatrixType > class_< TMatrixType > CreateMatrixInterface(pybind11::module& m, std::string Name )
        {
        class_< TMatrixType > binder(m,Name.c_str());

        //binder.def(init<std::TMatrixType& >())
        binder.def("Size1", [](const TMatrixType& self){return self.size1();} );
        binder.def("Size2", [](const TMatrixType& self){return self.size2();} );
        binder.def("Resize", [](TMatrixType& self, typename TMatrixType::size_type new_size1, typename TMatrixType::size_type new_size2)
                            {if(self.size1() != new_size1 || self.size2() != new_size2) 
                                self.resize(new_size1, new_size2, false);} );
        binder.def("__setitem__", [](TMatrixType& self, const tuple& index, const double value)
                                    {
                                        auto index_i = cast<typename TMatrixType::size_type>(index[0]);
                                        auto index_j = cast<typename TMatrixType::size_type>(index[1]);
                                        self(index_i,index_j) = value;
                                    } );
        binder.def("__getitem__", [](TMatrixType& self, const tuple& index)
                                    {
                                        auto index_i = cast<typename TMatrixType::size_type>(index[0]);
                                        auto index_j = cast<typename TMatrixType::size_type>(index[1]);
                                        return self(index_i, index_j);} );
        binder.def("__repr__", [](const TMatrixType& self)
                { std::stringstream out;
                  out << ( self );
                  return out.str();
                });         
        return binder;
        }



void  AddMatrixToPython(pybind11::module& m)
{
    auto matrix_binder = CreateMatrixInterface< Matrix >(m,"Matrix");
    matrix_binder.def(init<matrix<double>::size_type, matrix<double>::size_type>());
    matrix_binder.def(init<matrix<double>::size_type, matrix<double>::size_type, matrix<double>::value_type >());
    /*		  .def(MatrixScalarOperatorPython<matrix<double>, double>())
      .def(MatrixScalarAssignmentOperatorPython<matrix<double>, double>())
          .def(MatrixMatrixOperatorPython<matrix<double>, zero_matrix<double>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, identity_matrix<double>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, scalar_matrix<double>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, banded_matrix<double>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, triangular_matrix<double, upper>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, triangular_matrix<double, lower>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, symmetric_matrix<double, upper>, matrix<double> >())
    #if defined KRATOS_ADD_HERMITIAN_MATRIX_INTERFACE
          .def(MatrixMatrixOperatorPython<matrix<double>, hermitian_matrix<double, upper>, matrix<double> >())
    #endif
          .def(MatrixMatrixOperatorPython<matrix<double>, mapped_matrix<double>, matrix<double> >())
          .def(MatrixMatrixOperatorPython<matrix<double>, compressed_matrix<double>, matrix<double> >())
    #if defined KRATOS_ADD_COORDINATE_MATRIX_INTERFACE
          .def(MatrixMatrixOperatorPython<matrix<double>, coordinate_matrix<double>, matrix<double> >())
    #endif*/
    ;
}

}  // namespace Python.

} // Namespace Kratos


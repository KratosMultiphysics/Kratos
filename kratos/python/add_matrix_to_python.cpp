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
#include <pybind11/complex.h>

// Project includes
#include "includes/define_python.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"
#include "add_matrix_to_python.h"
#include "containers/array_1d.h"
// #include "python/matrix_python_interface.h"
// #include "python/matrix_scalar_operator_python.h"
// #include "python/matrix_scalar_assignment_operator_python.h"
// #include "python/matrix_matrix_operator_python.h"

namespace Kratos
{
namespace Python
{
    namespace py = pybind11;

    template< typename TMatrixType > py::class_< TMatrixType > CreateMatrixInterface(pybind11::module& m, std::string Name )
    {
        py::class_< TMatrixType, Kratos::shared_ptr<TMatrixType> > binder(m,Name.c_str(),py::buffer_protocol());
        binder.def(py::init<>());

        //binder.def(py::init<std::TMatrixType& >())
        binder.def("Size1", [](const TMatrixType& self){return self.size1();} );
        binder.def("Size2", [](const TMatrixType& self){return self.size2();} );
        binder.def("Resize", [](TMatrixType& self, const typename TMatrixType::size_type new_size1, const typename TMatrixType::size_type new_size2)
                            {if(self.size1() != new_size1 || self.size2() != new_size2)
                                self.resize(new_size1, new_size2, false);} );
        binder.def("__setitem__", [](TMatrixType& self, const std::pair<int,int> index, const  typename TMatrixType::value_type value)
                                    {
                                        const int index_i = index.first;
                                        const int index_j = index.second;
                                        self(index_i,index_j) = value;
                                    } );
        binder.def("__getitem__", [](const TMatrixType& self, const std::pair<int,int> index)
                                    {
                                        const int index_i = index.first;
                                        const int index_j = index.second;
                                        return self(index_i, index_j);} );

        //out of place versions
        binder.def("__add__", [](const TMatrixType& m1, const TMatrixType& m2){ return TMatrixType(m1+m2);}, py::is_operator());
        binder.def("__sub__", [](const TMatrixType& m1, const TMatrixType& m2){ return TMatrixType(m1-m2);}, py::is_operator());
        binder.def("__mul__", [](const TMatrixType& m1, const TMatrixType& m2){ return TMatrixType(prod(m1,m2));}, py::is_operator());
        binder.def("__mul__", [](const TMatrixType& m1, const typename TMatrixType::value_type& value){ return TMatrixType(m1*value);}, py::is_operator());
        binder.def("__truediv__", [](const TMatrixType& m1, const typename TMatrixType::value_type& value){ return TMatrixType(m1/value);}, py::is_operator());

        //inplace versions
        binder.def("__imul__", [](TMatrixType& m1, const typename TMatrixType::value_type& value){ m1*=value; return m1;}, py::is_operator());
        binder.def("__itruediv__", [](TMatrixType& m1, const typename TMatrixType::value_type& value){ m1/=value; return m1;}, py::is_operator());

        binder.def("__str__", PrintObject<TMatrixType>);
        return std::move(binder);
    }

    void  AddMatrixToPython(pybind11::module& m)
    {
        //here we add the dense matrix
        auto matrix_binder = CreateMatrixInterface< DenseMatrix<double> >(m,"Matrix");
        matrix_binder.def(py::init<const DenseMatrix<double>::size_type, const DenseMatrix<double>::size_type>());
      #ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        // This constructor is not supported by AMatrix
        //matrix_binder.def(py::init<const DenseMatrix<double>::size_type, const DenseMatrix<double>::size_type, const DenseMatrix<double>::value_type >());
        matrix_binder.def("fill", [](DenseMatrix<double>& self, const typename DenseMatrix<double>::value_type value) { self.fill(value); });
        matrix_binder.def("fill_identity", [](DenseMatrix<double>& self) { self.fill_identity(); });
      #else
        matrix_binder.def(py::init<const DenseMatrix<double>::size_type, const DenseMatrix<double>::size_type, const DenseMatrix<double>::value_type >());
        matrix_binder.def("fill", [](DenseMatrix<double>& self, const typename DenseMatrix<double>::value_type value) { noalias(self) = DenseMatrix<double>(self.size1(),self.size2(),value); });
        matrix_binder.def("fill_identity", [](DenseMatrix<double>& self) { noalias(self) = IdentityMatrix(self.size1()); });
      #endif // KRATOS_USE_AMATRIX
        matrix_binder.def(py::init<const DenseMatrix<double>& >());
        matrix_binder.def("__mul__", [](const DenseMatrix<double>& m1, const Vector& v){ return Vector(prod(m1,v));}, py::is_operator());
        matrix_binder.def("__mul__", [](const DenseMatrix<double>& m1, const array_1d<double,3>& v){ if(m1.size2() != 3) KRATOS_ERROR << "matrix size2 is not 3!" << std::endl; return Vector(prod(m1,v));}, py::is_operator());
        matrix_binder.def_buffer( [](DenseMatrix<double>& self)-> py::buffer_info{
                                                                    return py::buffer_info(
                                                                    self.data().begin(),
                                                                    sizeof(double),
                                                                    py::format_descriptor<double>::format(),
                                                                    2,
                                                                    {self.size1(),self.size2()},
                                                                    {sizeof(double)*self.size2(),
                                                                     sizeof(double) }
                                                                    );
                                                                    });

        ;

        //here we add the sparse matrix
        auto compressed_matrix_binder = CreateMatrixInterface< CompressedMatrix >(m,"CompressedMatrix");
        compressed_matrix_binder.def(py::init<const CompressedMatrix::size_type, const CompressedMatrix::size_type>());
        compressed_matrix_binder.def(py::init<const CompressedMatrix& >());
        compressed_matrix_binder.def("value_data", [](const CompressedMatrix& rA) ->  std::vector<double> 
                                                    {return std::vector<double>(
                                                        rA.value_data().begin(),
                                                        rA.value_data().end()
                                                        ) ;});
        compressed_matrix_binder.def("index1_data", [](const CompressedMatrix& rA) -> std::vector<std::size_t> 
                                                    {return std::vector<std::size_t>(
                                                        rA.index1_data().begin(),
                                                        rA.index1_data().end()
                                                        ) ;});
        compressed_matrix_binder.def("index2_data", [](const CompressedMatrix& rA) -> std::vector<std::size_t> 
                                                    {return std::vector<std::size_t>(
                                                        rA.index2_data().begin(),
                                                        rA.index2_data().end()
                                                        ) ;});

        //here we add the complex dense matrix
        auto cplx_matrix_binder = CreateMatrixInterface< ComplexMatrix >(m,"ComplexMatrix");
        cplx_matrix_binder.def(py::init<const ComplexMatrix::size_type, const ComplexMatrix::size_type>());
      #ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        // This constructor is not supported by AMatrix
        //cplx_matrix_binder.def(py::init<const ComplexMatrix::size_type, const ComplexMatrix::size_type, const ComplexMatrix::value_type >());
        cplx_matrix_binder.def("fill", [](ComplexMatrix& self, const typename ComplexMatrix::value_type value) { self.fill(value); });
        cplx_matrix_binder.def("fill_identity", [](ComplexMatrix& self) { self.fill_identity(); });
      #else
        cplx_matrix_binder.def(py::init<const ComplexMatrix::size_type, const ComplexMatrix::size_type, const ComplexMatrix::value_type >());
        cplx_matrix_binder.def("fill", [](ComplexMatrix& self, const typename ComplexMatrix::value_type value) { noalias(self) = ComplexMatrix(self.size1(),self.size2(),value); });
        cplx_matrix_binder.def("fill_identity", [](ComplexMatrix& self) { noalias(self) = IdentityMatrix(self.size1()); });
      #endif // KRATOS_USE_AMATRIX
        cplx_matrix_binder.def(py::init<const ComplexMatrix& >());
        cplx_matrix_binder.def("__mul__", [](const ComplexMatrix& m1, const ComplexVector& v){ return ComplexVector(prod(m1,v));}, py::is_operator());

        ;
        //here we add the complex sparse matrix
        auto cplx_compressed_matrix_binder = CreateMatrixInterface< ComplexCompressedMatrix >(m,"ComplexCompressedMatrix");
        cplx_compressed_matrix_binder.def(py::init< const ComplexCompressedMatrix::size_type, const ComplexCompressedMatrix::size_type >());
        cplx_compressed_matrix_binder.def(py::init< const CompressedMatrix& >());
        cplx_compressed_matrix_binder.def(py::init< const ComplexCompressedMatrix& >());
        cplx_compressed_matrix_binder.def("value_data", [](const ComplexCompressedMatrix& rA) ->  std::vector<std::complex<double>>
                                                    {return std::vector<std::complex<double>>(
                                                        rA.value_data().begin(),
                                                        rA.value_data().end()
                                                        ) ;});
        cplx_compressed_matrix_binder.def("index1_data", [](const ComplexCompressedMatrix& rA) -> std::vector<std::size_t> 
                                                    {return std::vector<std::size_t>(
                                                        rA.index1_data().begin(),
                                                        rA.index1_data().end()
                                                        ) ;});
        cplx_compressed_matrix_binder.def("index2_data", [](const ComplexCompressedMatrix& rA) -> std::vector<std::size_t> 
                                                    {return std::vector<std::size_t>(
                                                        rA.index2_data().begin(),
                                                        rA.index2_data().end()
                                                        ) ;});

    }

}  // namespace Python.

} // Namespace Kratos


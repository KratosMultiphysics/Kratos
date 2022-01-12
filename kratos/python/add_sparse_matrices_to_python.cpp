
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
#include "includes/define_python.h"
#include "includes/ublas_interface.h"
#include "python/add_sparse_matrices_to_python.h"
#include "containers/sparse_graph.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "containers/csr_matrix.h"
#include "containers/system_vector.h"

namespace Kratos
{

namespace Python
{
namespace py = pybind11;

void AddSparseMatricesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef std::size_t IndexType;

    py::class_<SparseGraph<IndexType>, SparseGraph<IndexType>::Pointer >(m, "SparseGraph")
    .def(py::init<IndexType>())
    .def("GetComm", &SparseGraph<IndexType>::GetComm)
    .def("pGetComm", &SparseGraph<IndexType>::pGetComm)
    .def("Size", &SparseGraph<IndexType>::Size)
    .def("IsEmpty", &SparseGraph<IndexType>::IsEmpty)
    .def("AddEnty", &SparseGraph<IndexType>::AddEntry)
    .def("Finalize", &SparseGraph<IndexType>::Finalize)
    .def("GetGraph", &SparseGraph<IndexType>::GetGraph)
    .def("ExportSingleVectorRepresentation", &SparseGraph<IndexType>::ExportSingleVectorRepresentation)
    .def("__str__", PrintObject<SparseGraph<IndexType>>);

    py::class_<SparseContiguousRowGraph<IndexType>, SparseContiguousRowGraph<IndexType>::Pointer >(m, "SparseGraph")
    .def(py::init<IndexType>())
    .def("GetComm", &SparseContiguousRowGraph<IndexType>::GetComm)
    .def("pGetComm", &SparseContiguousRowGraph<IndexType>::pGetComm)
    .def("Size", &SparseContiguousRowGraph<IndexType>::Size)
    .def("AddEnty", &SparseContiguousRowGraph<IndexType>::AddEntry)
    .def("Finalize", &SparseContiguousRowGraph<IndexType>::Finalize)
    .def("GetGraph", &SparseContiguousRowGraph<IndexType>::GetGraph)
    .def("ExportSingleVectorRepresentation", &SparseContiguousRowGraph<IndexType>::ExportSingleVectorRepresentation)
    .def("__str__", PrintObject<SparseContiguousRowGraph<IndexType>>)
    ;

    py::class_<CsrMatrix<double,IndexType>, CsrMatrix<double,IndexType>::Pointer >(m, "CsrMatrix")
    .def(py::init<SparseGraph<IndexType>>())
    .def(py::init<SparseContiguousRowGraph<IndexType>>())
    .def("GetComm", &CsrMatrix<double,IndexType>::GetComm)
    .def("pGetComm", &CsrMatrix<double,IndexType>::pGetComm)
    .def("SetValue", &CsrMatrix<double,IndexType>::SetValue)
    .def("size1", &CsrMatrix<double,IndexType>::size1)
    .def("size2", &CsrMatrix<double,IndexType>::size2)
    .def("nnz", &CsrMatrix<double,IndexType>::nnz)
    .def("index1_data", [](CsrMatrix<double,IndexType>& rA){
        const auto& data = rA.index1_data();
        std::vector<IndexType> v(data.begin(), data.end());
        return v;
    })
    .def("index2_data", [](CsrMatrix<double,IndexType>& rA){
        const auto& data = rA.index2_data();
        std::vector<IndexType> v(data.begin(), data.end());
        return v;
    })
    .def("value_data", [](CsrMatrix<double,IndexType>& rA){
        const auto& data = rA.value_data();
        std::vector<double> v(data.begin(), data.end());
        return v;
    })
    .def("SpMV", [](CsrMatrix<double,IndexType>& rA,SystemVector<double,IndexType>& x, SystemVector<double,IndexType>& y){
        rA.SpMV(x,y);
    })
    .def("TransposeSpMV", [](CsrMatrix<double,IndexType>& rA,SystemVector<double,IndexType>& x, SystemVector<double,IndexType>& y){
        rA.TransposeSpMV(x,y);
    })
    .def("BeginAssemble", &CsrMatrix<double,IndexType>::BeginAssemble)
    .def("FinalizeAssemble", &CsrMatrix<double,IndexType>::FinalizeAssemble)
    .def("Assemble", [](CsrMatrix<double,IndexType>& rA, 
                        Matrix& values, 
                        std::vector<IndexType>& indices){
       rA.Assemble(values,indices);
    })
    .def("Assemble", [](CsrMatrix<double,IndexType>& rA, 
                        Matrix& values, 
                        std::vector<IndexType>& row_indices,
                        std::vector<IndexType>& col_indices){
       rA.Assemble(values,row_indices, col_indices);
    })
    .def("AssembleEntry", [](CsrMatrix<double,IndexType>& rA, double value, IndexType I, IndexType J){
        rA.AssembleEntry(value,I,J);
    })
    .def("__str__", PrintObject<CsrMatrix<double,IndexType>>);


}

} // namespace Python.
} // Namespace Kratos

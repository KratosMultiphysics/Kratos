
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
#include "utilities/amgcl_csr_spmm_utilities.h"

namespace Kratos
{

namespace Python
{
namespace py = pybind11;

void AddSparseMatricesToPython(pybind11::module& m)
{
    namespace py = pybind11;
    typedef std::size_t IndexType;


    //************************************************************************************************************
    //************************************************************************************************************
    py::class_<SparseGraph<IndexType>, SparseGraph<IndexType>::Pointer >(m, "SparseGraph")
    .def(py::init<>())
    .def(py::init<IndexType>())
    .def("GetComm", &SparseGraph<IndexType>::GetComm)
    .def("pGetComm", &SparseGraph<IndexType>::pGetComm)
    .def("Size", &SparseGraph<IndexType>::Size)
    .def("IsEmpty", &SparseGraph<IndexType>::IsEmpty)
    .def("AddEntry", &SparseGraph<IndexType>::AddEntry)
    .def("Finalize", &SparseGraph<IndexType>::Finalize)
    .def("GetGraph", &SparseGraph<IndexType>::GetGraph)
    .def("AddEntries", [](SparseGraph<IndexType>& self, std::vector<IndexType>& indices){
        self.AddEntries(indices);
    })
    .def("AddEntries", [](SparseGraph<IndexType>& self, std::vector<IndexType>& row_indices, std::vector<IndexType>& col_indices){
        self.AddEntries(row_indices, col_indices);
    })
    .def("ExportSingleVectorRepresentation", &SparseGraph<IndexType>::ExportSingleVectorRepresentation)
    .def("__str__", PrintObject<SparseGraph<IndexType>>);


    //************************************************************************************************************
    //************************************************************************************************************
    py::class_<SparseContiguousRowGraph<IndexType>, SparseContiguousRowGraph<IndexType>::Pointer >(m, "SparseContiguousRowGraph")
    .def(py::init<IndexType>())
    .def("GetComm", &SparseContiguousRowGraph<IndexType>::GetComm)
    .def("pGetComm", &SparseContiguousRowGraph<IndexType>::pGetComm)
    .def("Size", &SparseContiguousRowGraph<IndexType>::Size)
    .def("AddEntry", &SparseContiguousRowGraph<IndexType>::AddEntry)
    .def("AddEntries", [](SparseContiguousRowGraph<IndexType>& self, std::vector<IndexType>& indices){
        self.AddEntries(indices);
    })
    .def("AddEntries", [](SparseContiguousRowGraph<IndexType>& self, std::vector<IndexType>& row_indices, std::vector<IndexType>& col_indices){
        self.AddEntries(row_indices, col_indices);
    })
    .def("Finalize", &SparseContiguousRowGraph<IndexType>::Finalize)
    .def("GetGraph", &SparseContiguousRowGraph<IndexType>::GetGraph)
    .def("ExportSingleVectorRepresentation", &SparseContiguousRowGraph<IndexType>::ExportSingleVectorRepresentation)
    .def("__str__", PrintObject<SparseContiguousRowGraph<IndexType>>)
    ;


    //************************************************************************************************************
    //************************************************************************************************************
    py::class_<CsrMatrix<double,IndexType>, CsrMatrix<double,IndexType>::Pointer >(m, "CsrMatrix")
    .def(py::init<SparseGraph<IndexType>&>())
    .def(py::init<SparseContiguousRowGraph<IndexType>&>())
    .def("GetComm", &CsrMatrix<double,IndexType>::GetComm)
    .def("pGetComm", &CsrMatrix<double,IndexType>::pGetComm)
    .def("SetValue", &CsrMatrix<double,IndexType>::SetValue)
    .def("Size1", &CsrMatrix<double,IndexType>::size1)
    .def("Size2", &CsrMatrix<double,IndexType>::size2)
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
    .def("__setitem__", [](CsrMatrix<double,IndexType>& self, const std::pair<int,int> index, const  double value)
        {
            const int index_i = index.first;
            const int index_j = index.second;
            self(index_i,index_j) = value;
        })
    .def("__getitem__", [](const CsrMatrix<double,IndexType>& self, const std::pair<int,int> index)
        {
            const int index_i = index.first;
            const int index_j = index.second;
            return self(index_i, index_j);
        })        
    .def("SpMV", [](CsrMatrix<double,IndexType>& rA,SystemVector<double,IndexType>& x, SystemVector<double,IndexType>& y){
        rA.SpMV(x,y);
    })
    // .def("__matmul__", [](CsrMatrix<double,IndexType>& rA,SystemVector<double,IndexType>& x, SystemVector<double,IndexType>& y){
    //     rA.SpMV(x,y);
    // })
    .def("TransposeSpMV", [](CsrMatrix<double,IndexType>& rA,SystemVector<double,IndexType>& x, SystemVector<double,IndexType>& y){
        rA.TransposeSpMV(x,y);
    })
    .def("SpMM", [](CsrMatrix<double,IndexType>& rA,CsrMatrix<double,IndexType>& rB){
        return AmgclCSRSpMMUtilities::SparseMultiply(rA,rB);
    })
    .def("__matmul__", [](CsrMatrix<double,IndexType>& rA,CsrMatrix<double,IndexType>& rB){
        return AmgclCSRSpMMUtilities::SparseMultiply(rA,rB);
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



    //************************************************************************************************************
    //************************************************************************************************************
    py::class_<SystemVector<double,IndexType>, SystemVector<double,IndexType>::Pointer >(m, "SystemVector", py::buffer_protocol())
    .def(py::init<IndexType>())
    .def(py::init<IndexType, DataCommunicator&>())
    .def(py::init<SparseGraph<IndexType>&>())
    .def(py::init<SparseContiguousRowGraph<IndexType>&>())
    .def("GetComm", &SystemVector<double,IndexType>::GetComm)
    .def("pGetComm", &SystemVector<double,IndexType>::pGetComm)
    .def("Size", &SystemVector<double,IndexType>::size)
    .def("size", &SystemVector<double,IndexType>::size)
    .def("Clear", &SystemVector<double,IndexType>::Clear)
    .def("SetValue", &SystemVector<double,IndexType>::SetValue)
    .def("__setitem__", [](SystemVector<double,IndexType>& self, const IndexType i, const double value){self[i] = value;} )
    .def("__getitem__", [](const SystemVector<double,IndexType>& self, const IndexType i){return self[i];} )
    .def("Add", &SystemVector<double,IndexType>::Add)
    .def("BeginAssemble", &SystemVector<double,IndexType>::BeginAssemble)
    .def("FinalizeAssemble", &SystemVector<double,IndexType>::FinalizeAssemble)
    .def("BeginAssemble", &SystemVector<double,IndexType>::BeginAssemble)
    .def("data", [](SystemVector<double,IndexType>& self){
       return self.data();
    })
    .def_buffer([](SystemVector<double,IndexType>& self) -> py::buffer_info {
        return py::buffer_info(
            &(self.data())[0],                               /* Pointer to buffer */
            sizeof(double),                          /* Size of one scalar */
            py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
            1,                                      /* Number of dimensions */
            { self.size() },                 /* Buffer dimensions */
            { sizeof(double)  }             /* Strides (in bytes) for each index */
        );
    })
    .def("Assemble", [](SystemVector<double,IndexType>& self, 
                        Vector& values, 
                        std::vector<IndexType>& indices){
       self.Assemble(values,indices);
    })
    //inplace
    .def("__iadd__", [](SystemVector<double,IndexType>& self, const SystemVector<double,IndexType>& other_vec){self += other_vec; }, py::is_operator())
    .def("__isub__", [](SystemVector<double,IndexType>& self, const SystemVector<double,IndexType>& other_vec){self -= other_vec;  }, py::is_operator())
    .def("__imul__", [](SystemVector<double,IndexType>& self, const double& value){ self*=value; }, py::is_operator())
    .def("__itruediv__", [](SystemVector<double,IndexType>& self, const double& value){ self/=value; }, py::is_operator())
    //out of place
    //.def("__add__", [](const SystemVector<double,IndexType>& vec1, const SystemVector<double,IndexType>& vec2){SystemVector<double,IndexType> aux(vec1); aux += vec2; return std::move(aux);}, py::is_operator())
    // .def("__sub__", [](const SystemVector<double,IndexType>& vec1, const SystemVector<double,IndexType>& vec2){SystemVector<double,IndexType> aux(vec1); aux -= vec2; return aux;}, py::is_operator())
    //access operators
    .def("__setitem__", [](SystemVector<double,IndexType>& self, const unsigned int i, const double value){self[i] = value;} )
    .def("__getitem__", [](const SystemVector<double,IndexType>& self, const unsigned int i){return self[i];} )
    .def("__str__", PrintObject<SystemVector<double,IndexType>>)
    ;
}

} // namespace Python.
} // Namespace Kratos

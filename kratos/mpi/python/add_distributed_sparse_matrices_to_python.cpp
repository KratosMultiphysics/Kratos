//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "containers/distributed_system_vector.h"
#include "containers/distributed_sparse_graph.h"
#include "containers/distributed_numbering.h"
#include "containers/distributed_vector_importer.h"
#include "containers/distributed_vector_exporter.h"
#include "containers/distributed_csr_matrix.h"
#include "add_distributed_sparse_matrices_to_python.h"
#include "mpi/utilities/amgcl_distributed_csr_spmm_utilities.h"


namespace Kratos {
namespace Python {

void AddDistributedSparseMatricesToPython(pybind11::module& m)
{
    namespace py = pybind11;
    typedef std::size_t IndexType;
    typedef std::size_t MPIIndexType;

    py::class_<DistributedNumbering<IndexType>,DistributedNumbering<IndexType>::Pointer>(m,"DistributedNumbering")
        .def(py::init<const DataCommunicator&, IndexType>())
        .def(py::init<const DataCommunicator&, IndexType, MPIIndexType>())
        .def("GetComm", &DistributedNumbering<IndexType>::GetComm)
        .def("LocalSize", &DistributedNumbering<IndexType>::LocalSize)
        .def("Size", &DistributedNumbering<IndexType>::Size)
        .def("IsLocal", &DistributedNumbering<IndexType>::IsLocal)
        .def("LocalId", &DistributedNumbering<IndexType>::LocalId)
        .def("GlobalId", &DistributedNumbering<IndexType>::GlobalId)
        .def("RemoteLocalId", &DistributedNumbering<IndexType>::RemoteLocalId)
        .def("OwnerRank", &DistributedNumbering<IndexType>::OwnerRank)
        .def("GetCpuBounds", &DistributedNumbering<IndexType>::GetCpuBounds)
        ;

    py::class_<DistributedVectorExporter<IndexType>,DistributedVectorExporter<IndexType>::Pointer>(m,"DistributedVectorExporter")
        .def(py::init<const DataCommunicator&, const std::vector<IndexType>&, const DistributedNumbering<IndexType>&>())
        .def("GetComm", &DistributedVectorExporter<IndexType>::GetComm)
        .def("ApplyAddition", &DistributedVectorExporter<IndexType>::Apply<DistributedSystemVector<IndexType>, DenseVector<double>, std::plus<double>>)
        .def("ApplySubtraction", &DistributedVectorExporter<IndexType>::Apply<DistributedSystemVector<IndexType>, DenseVector<double>, std::minus<double>>)
        ;

    py::class_<DistributedVectorImporter<double,IndexType>,DistributedVectorImporter<double,IndexType>::Pointer>(m,"DistributedVectorImporter")
        .def(py::init<const DataCommunicator&, const std::vector<IndexType>&, const DistributedNumbering<IndexType>&>())
        .def("GetComm", &DistributedVectorImporter<double,IndexType>::GetComm)
        .def("ImportData", &DistributedVectorImporter<double,IndexType>::ImportData)
        ;

    py::class_<DistributedSparseGraph<IndexType>, DistributedSparseGraph<IndexType>::Pointer>(m,"DistributedSparseGraph")
        .def(py::init<const IndexType, DataCommunicator&>())
        .def("GetComm", &DistributedSparseGraph<IndexType>::GetComm)
        .def("Size", &DistributedSparseGraph<IndexType>::Size)
        .def("LocalSize", &DistributedSparseGraph<IndexType>::LocalSize)
        .def("Has", &DistributedSparseGraph<IndexType>::Has)
        .def("GetRowNumbering", &DistributedSparseGraph<IndexType>::Has)
        .def("ComputeLocalMinMaxColumnIndex", &DistributedSparseGraph<IndexType>::ComputeLocalMinMaxColumnIndex)
        .def("ComputeMaxGlobalColumnIndex", &DistributedSparseGraph<IndexType>::ComputeMaxGlobalColumnIndex)
        .def("AddEntry", &DistributedSparseGraph<IndexType>::AddEntry)
        .def("Finalize", &DistributedSparseGraph<IndexType>::Finalize)
        .def("GetLocalGraph", &DistributedSparseGraph<IndexType>::GetLocalGraph)
        .def("AddEntries", [](DistributedSparseGraph<IndexType>& self, std::vector<IndexType>& indices){
            self.AddEntries(indices);
        })
        .def("AddEntries", [](DistributedSparseGraph<IndexType>& self, std::vector<IndexType>& row_indices, std::vector<IndexType>& col_indices){
            self.AddEntries(row_indices, col_indices);
        })
        .def("__str__", PrintObject<DistributedSparseGraph<IndexType>>)
        ;

    py::class_<DistributedSystemVector<double,IndexType>,  DistributedSystemVector<double,IndexType>::Pointer>(m,"DistributedSystemVector")
        .def(py::init<const DistributedSparseGraph<IndexType>&>())
        .def(py::init<DistributedSystemVector<double,IndexType>&>())
        .def(py::init<const DistributedNumbering<IndexType>&>())
        .def("GetComm", &DistributedSystemVector<double,IndexType>::GetComm)
        .def("LocalSize", [](const DistributedSystemVector<double,IndexType>& self)
            {return self.LocalSize();})
        .def("GetLocalData", [](DistributedSystemVector<double,IndexType>& self)
             {return self.GetLocalData();})
        .def("GetNumbering", &DistributedSystemVector<double,IndexType>::GetNumbering)
        .def("Size", &DistributedSystemVector<double,IndexType>::Size)
        .def("size", &DistributedSystemVector<double,IndexType>::Size)
        .def("Clear", &DistributedSystemVector<double,IndexType>::Clear)
        .def("SetValue", &DistributedSystemVector<double,IndexType>::SetValue)
        .def("__setitem__", [](DistributedSystemVector<double,IndexType>& self, const IndexType i, const double value){self[i] = value;} )
        .def("__getitem__", [](const DistributedSystemVector<double,IndexType>& self, const IndexType i){return self[i];} )
        .def("Add", &DistributedSystemVector<double,IndexType>::Add)
        .def("BeginAssemble", &DistributedSystemVector<double,IndexType>::BeginAssemble)
        .def("FinalizeAssemble", &DistributedSystemVector<double,IndexType>::FinalizeAssemble)
        .def("Assemble", [](DistributedSystemVector<double,IndexType>& self,
                            Vector& values,
                            std::vector<IndexType>& indices){
        self.Assemble(values,indices);
        })
        .def("Dot", &DistributedSystemVector<double,IndexType>::Dot)
        //inplace
        .def("__iadd__", [](DistributedSystemVector<double,IndexType>& self, const DistributedSystemVector<double,IndexType>& other_vec){self += other_vec; }, py::is_operator())
        .def("__isub__", [](DistributedSystemVector<double,IndexType>& self, const DistributedSystemVector<double,IndexType>& other_vec){self -= other_vec;  }, py::is_operator())
        .def("__imul__", [](DistributedSystemVector<double,IndexType>& self, const double& value){ self*=value; }, py::is_operator())
        .def("__itruediv__", [](DistributedSystemVector<double,IndexType>& self, const double& value){ self/=value; }, py::is_operator())
        //access operators
        .def("__setitem__", [](DistributedSystemVector<double,IndexType>& self, const unsigned int i, const double value){self[i] = value;} )
        .def("__getitem__", [](const DistributedSystemVector<double,IndexType>& self, const unsigned int i){return self[i];} )
        .def("__str__", PrintObject<DistributedSystemVector<double,IndexType>>)
        ;


    py::class_<DistributedCsrMatrix<double,IndexType>, DistributedCsrMatrix<double,IndexType>::Pointer>(m,"DistributedCsrMatrix")
        .def(py::init<const DistributedSparseGraph<IndexType>&>())
        .def("GetComm", &DistributedCsrMatrix<double,IndexType>::GetComm)
        .def("GetRowNumbering", &DistributedCsrMatrix<double,IndexType>::GetRowNumbering)
        .def("GetColNumbering", &DistributedCsrMatrix<double,IndexType>::GetColNumbering)
        .def("SetValue", &DistributedCsrMatrix<double,IndexType>::SetValue)
        .def("local_size1", &DistributedCsrMatrix<double,IndexType>::local_size1)
        .def("size2", &DistributedCsrMatrix<double,IndexType>::size2)
        .def("Clear", &DistributedCsrMatrix<double,IndexType>::Clear)
        .def("GetDiagonalBlock", [](DistributedCsrMatrix<double,IndexType>& self) -> CsrMatrix<double,IndexType>&
             {return self.GetDiagonalBlock();})
        .def("GetOffDiagonalBlock", [](DistributedCsrMatrix<double,IndexType>& self) -> CsrMatrix<double,IndexType>&
             {return self.GetOffDiagonalBlock();})
        .def("GetOffDiagonalLocalIds", [](DistributedCsrMatrix<double,IndexType>& self)
             {return self.GetOffDiagonalLocalIds();})
        .def("GetOffDiagonalGlobalIds", [](DistributedCsrMatrix<double,IndexType>& self)
            {
                return std::vector<IndexType>(self.GetOffDiagonalGlobalIds().begin(), self.GetOffDiagonalGlobalIds().end());
            })
        .def("GetOffDiagonalBlockLocalId", &DistributedCsrMatrix<double,IndexType>::GetOffDiagonalBlockLocalId)
        .def("GetOffDiagonalBlockGlobalId", &DistributedCsrMatrix<double,IndexType>::GetOffDiaGlobalId)
        .def("__setitem__", [](DistributedCsrMatrix<double,IndexType>& self,
            const IndexType i,
            const IndexType j,
            const double value){
            self.GetLocalDataByGlobalId(i,j) = value;} )
        .def("__getitem__", [](DistributedCsrMatrix<double,IndexType>& self,
            const IndexType i,
            const IndexType j){
            return self.GetLocalDataByGlobalId(i,j);} )
        .def("GetDiagonalIndex2DataInGlobalNumbering", &DistributedCsrMatrix<double,IndexType>::GetDiagonalIndex2DataInGlobalNumbering)
        .def("GetOffDiagonalIndex2DataInGlobalNumbering", &DistributedCsrMatrix<double,IndexType>::GetOffDiagonalIndex2DataInGlobalNumbering)
        .def("SpMV", [](DistributedCsrMatrix<double,IndexType>& rA,
                        DistributedSystemVector<double,IndexType>& x,
                        DistributedSystemVector<double,IndexType>& y){
            rA.SpMV(x,y);
        })
        .def("TransposeSpMV", [](DistributedCsrMatrix<double,IndexType>& rA,
                        DistributedSystemVector<double,IndexType>& x,
                        DistributedSystemVector<double,IndexType>& y){
            rA.TransposeSpMV(x,y);
        })
        .def("SpMM", [](const DistributedCsrMatrix<double,IndexType>& self,
                        const DistributedCsrMatrix<double,IndexType>& rB){
            return DistributedAmgclCSRSpMMUtilities::SparseMultiply(self,rB);
        })
        .def("NormFrobenius", &DistributedCsrMatrix<double,IndexType>::NormFrobenius)
        .def("BeginAssemble", &DistributedCsrMatrix<double,IndexType>::BeginAssemble)
        .def("FinalizeAssemble", &DistributedCsrMatrix<double,IndexType>::FinalizeAssemble)
        .def("Assemble", [](DistributedCsrMatrix<double,IndexType>& rA,
                            Matrix& values,
                            std::vector<IndexType>& indices){
            rA.Assemble(values,indices);
            })
        .def("Assemble", [](DistributedCsrMatrix<double,IndexType>& rA,
                            Matrix& values,
                            std::vector<IndexType>& row_indices,
                            std::vector<IndexType>& col_indices){
            rA.Assemble(values,row_indices, col_indices);
            })
        .def("AssembleEntry", [](DistributedCsrMatrix<double,IndexType>& rA, double value, IndexType I, IndexType J){
            rA.AssembleEntry(value,I,J);
            })
        .def("__str__", PrintObject<DistributedCsrMatrix<double,IndexType>>);
}

} // namespace Python
} // namespace Kratos


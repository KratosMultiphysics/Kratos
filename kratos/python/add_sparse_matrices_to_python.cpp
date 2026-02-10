//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes
#include <pybind11/numpy.h>

// Project includes
#include "includes/define_python.h"
#include "includes/ublas_interface.h"
#include "python/add_sparse_matrices_to_python.h"
#include "containers/nd_data.h"
#include "containers/sparse_graph.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "containers/csr_matrix.h"
#include "containers/system_vector.h"
#include "utilities/amgcl_csr_spmm_utilities.h"

namespace Kratos::Python
{

namespace py = pybind11;

/**
 * @brief Anonymous namespace to store the assembly helper functions.
 * @note This namespace is not intended to be used outside of this file as these functions are tailored to enable the efficient Python assembly.
 */
namespace
{

    using IndexType = std::size_t;

    using CsrMatrixType = CsrMatrix<double, IndexType>;

    using EquationIdVectorType = std::vector<std::size_t>;

    /**
     * @brief Get the Equation Id Csr Indices object
     * @details This function computes the CSR matrix value vector indices corresponding to the equation ids of the entities in the provided container.
     * @tparam TContainerType Type of the container of the entities.
     * @param rCsrMatrix The CSR matrix to get the indices from.
     * @param rContainer The container of the entities.
     * @param rProcessInfo The process info.
     * @param pNDData Pointer to the NDData object to store the indices.
     */
    template<class TContainerType>
    NDData<int> GetEquationIdCsrIndices(
        const CsrMatrixType& rCsrMatrix,
        const TContainerType& rContainer,
        const ProcessInfo& rProcessInfo)
    {
        // Get equation ids size from the first entity (assuming all entities have the same size)
        EquationIdVectorType equation_ids;
        rContainer.begin()->EquationIdVector(equation_ids, rProcessInfo);
        const std::size_t local_size = equation_ids.size();

        // Assign the input NDData to have shape: number of entities * local_size * local_size
        const std::size_t n_entities = std::distance(rContainer.begin(), rContainer.end());
        DenseVector<unsigned int> nd_data_shape(3);
        nd_data_shape[0] = n_entities;
        nd_data_shape[1] = local_size;
        nd_data_shape[2] = local_size;
        auto eq_ids_data = NDData<int>(nd_data_shape);

        // Loop over the container
        auto eq_ids_data_view = eq_ids_data.ViewData();
        IndexPartition<std::size_t>(n_entities).for_each([&](std::size_t i) {
            // Get current entity
            auto it = rContainer.begin() + i;
            const std::size_t it_pos = i * (local_size * local_size);

            // Get current entity equation ids
            EquationIdVectorType equation_ids;
            it->EquationIdVector(equation_ids, rProcessInfo);

            // Loop over the DOFs
            for (unsigned int i_local = 0; i_local < local_size; ++i_local) {
                const unsigned int i_global = equation_ids[i_local]; // Row global equation id
                for(unsigned int j_local = 0; j_local < local_size; ++j_local) {
                    const unsigned int j_global = equation_ids[j_local]; // Column global equation id
                    const unsigned int csr_index = rCsrMatrix.FindValueIndex(i_global,j_global); // Index in the CSR matrix values vector
                    eq_ids_data_view[it_pos + i_local * local_size + j_local] = csr_index;
                }
            }
        });

        // Return the data container
        return eq_ids_data;
    }

    /**
     * @brief Assemble the global stiffness matrix from the entities contributions and the equation ids indices.
     * @param rLeftHandSideElementContributions The left hand side ontributions to the global stiffness matrix for each entity.
     * @param rElementEqIdsCsrIndices The CSR values vector indices for each left hand side entry.
     * @param rLeftHandSide The global stiffness matrix in which the contributions will be added.
     */
    void AssembleWithCsrIndices(
        const NDData<double>& rLeftHandSideElementContributions,
        const NDData<int>& rElementEqIdsCsrIndices,
        CsrMatrixType& rLeftHandSide)
    {
        // Get and check the provided data shapes
        const auto& shape = rElementEqIdsCsrIndices.Shape();
        const std::size_t n_entities = shape[0];
        const std::size_t local_size_1 = shape[1];
        const std::size_t local_size_2 = shape[2];
        KRATOS_ERROR_IF(shape.size() != rLeftHandSideElementContributions.Shape().size()) << "The shapes of the equation ids indices and left hand side contributions do not match." << std::endl;
        for (std::size_t i = 0; i < shape.size(); ++i) {
            KRATOS_ERROR_IF(shape[i] != rLeftHandSideElementContributions.Shape()[i]) << "The shapes of the equation ids indices and left hand side contributions do not match in component " << i << "." << std::endl;
        }

        // Get the left hand side data
        auto& r_lhs_data = rLeftHandSide.value_data();
        const auto& r_idx_data = rElementEqIdsCsrIndices.ViewData();
        const auto& r_lhs_contribution_data = rLeftHandSideElementContributions.ViewData();

        // Loop over the entities
        IndexPartition<std::size_t>(n_entities).for_each([&](std::size_t i) {
            const std::size_t entity_pos = i * (local_size_1 * local_size_2);
            for (std::size_t i_local = 0; i_local < local_size_1; ++i_local) {
                for (std::size_t j_local = 0; j_local < local_size_2; ++j_local) {
                    const std::size_t aux_idx = entity_pos + i_local * local_size_1 + j_local; // Position in the contributions and equation ids data
                    const int csr_index = r_idx_data[aux_idx]; // Index in the CSR matrix values vector
                    const double lhs_contribution = r_lhs_contribution_data[aux_idx]; // Scalar contribution to the left hand side
                    r_lhs_data[csr_index] += lhs_contribution;
                }
            }
        });
    }

    /**
     * @brief Assemble the global stiffness matrix from the entities contributions and the equation ids.
     * @param rLeftHandSideElementContributions The left hand side contributions to the global stiffness matrix for each entity.
     * @param rElementEqIds The equation ids of each entity.
     * @param rLeftHandSide The global stiffness matrix in which the contributions will be added.
     */
    void AssembleWithEquationIds(
        const NDData<double>& rLeftHandSideElementContributions,
        const NDData<int>& rElementEqIds,
        CsrMatrixType& rLeftHandSide)
    {
        // Get and check the provided data shapes
        const auto& eq_ids_shape = rElementEqIds.Shape();
        const std::size_t n_entities = eq_ids_shape[0];
        const std::size_t local_size = eq_ids_shape[1];
        KRATOS_ERROR_IF(n_entities != rLeftHandSideElementContributions.Shape()[0]) << "The shapes of the equation ids and left hand side contributions refer to different number of entites." << std::endl;
        KRATOS_ERROR_IF(local_size != rLeftHandSideElementContributions.Shape()[1]) << "The shapes of the equation ids and left hand side contributions do not match in component 1." << std::endl;
        KRATOS_ERROR_IF(local_size != rLeftHandSideElementContributions.Shape()[2]) << "The shapes of the equation ids and left hand side contributions do not match in component 2." << std::endl;

        // Get the left hand side data
        const auto& r_eq_ids_data = rElementEqIds.ViewData();
        const auto& r_lhs_contribution_data = rLeftHandSideElementContributions.ViewData();

        // Loop over the entities
        IndexPartition<std::size_t>(n_entities).for_each([&](std::size_t i) {
            const std::size_t eq_id_pos = i * local_size;
            const std::size_t lhs_pos = i * (local_size * local_size);
            for (std::size_t i_local = 0; i_local < local_size; ++i_local) {
                for (std::size_t j_local = 0; j_local < local_size; ++j_local) {
                    const std::size_t aux_idx = lhs_pos + i_local * local_size + j_local; // Position in the LHS contributions data
                    const double lhs_contribution = r_lhs_contribution_data[aux_idx]; // Scalar contribution to the left hand side
                    rLeftHandSide.AssembleEntry(lhs_contribution, r_eq_ids_data[eq_id_pos + i_local], r_eq_ids_data[eq_id_pos + j_local]);
                }
            }
        });
    }

}

void AddSparseMatricesToPython(pybind11::module& m)
{
    namespace py = pybind11;
    typedef std::size_t IndexType;

    using ElementsContainerType = typename ModelPart::ElementsContainerType;

    using ConditionsContainerType = typename ModelPart::ConditionsContainerType;


    //************************************************************************************************************
    //************************************************************************************************************
    py::class_<SparseGraph<IndexType>, SparseGraph<IndexType>::Pointer >(m, "SparseGraph")
    .def(py::init<>())
    .def(py::init<IndexType>())
    .def("GetComm", &SparseGraph<IndexType>::GetComm, py::return_value_policy::reference_internal)
    .def("Size", &SparseGraph<IndexType>::Size)
    .def("IsEmpty", &SparseGraph<IndexType>::IsEmpty)
    .def("AddEntry", &SparseGraph<IndexType>::AddEntry)
    .def("Finalize", &SparseGraph<IndexType>::Finalize)
    .def("GetGraph", &SparseGraph<IndexType>::GetGraph, py::return_value_policy::reference_internal)
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
    .def("GetComm", &SparseContiguousRowGraph<IndexType>::GetComm, py::return_value_policy::reference_internal)
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
    .def(py::init<>())
    .def(py::init<const DataCommunicator&>())
    .def(py::init<const CsrMatrix<double,IndexType>&>())
    .def(py::init<SparseGraph<IndexType>&>())
    .def(py::init<SparseContiguousRowGraph<IndexType>&>())
    .def("copy", [](const CsrMatrix<double,IndexType>& self){
        return Kratos::make_shared<CsrMatrix<double,IndexType>>(self);
    })
    .def("GetComm", &CsrMatrix<double,IndexType>::GetComm, py::return_value_policy::reference_internal)
    .def("SetValue", &CsrMatrix<double,IndexType>::SetValue)
    .def("Size1", &CsrMatrix<double,IndexType>::size1)
    .def("Size2", &CsrMatrix<double,IndexType>::size2)
    .def("size1", &CsrMatrix<double,IndexType>::size1)
    .def("size2", &CsrMatrix<double,IndexType>::size2)
    .def("nnz", &CsrMatrix<double,IndexType>::nnz)
    .def("index1_data", [](CsrMatrix<double,IndexType>& self){
            return py::array_t<IndexType>(
                self.index1_data().size(),
                self.index1_data().begin(),
                py::cast(self) //this is fundamental to avoid copying
            );
    }, py::return_value_policy::reference_internal)
    .def("index2_data", [](CsrMatrix<double,IndexType>& self){
            return py::array_t<IndexType>(
                self.index2_data().size(),
                self.index2_data().begin(),
                py::cast(self) //this is fundamental to avoid copying
            );
    }, py::return_value_policy::reference_internal)
    .def("value_data", [](CsrMatrix<double,IndexType>& self){
            return py::array_t<double>(
                self.value_data().size(),
                self.value_data().begin(),
                py::cast(self) //this is fundamental to avoid copying
            );
    }, py::return_value_policy::reference_internal)
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
    .def("SpMV", [](const CsrMatrix<double,IndexType>& rA,const SystemVector<double,IndexType>& x, SystemVector<double,IndexType>& y){
        rA.SpMV(x,y);
    })
    .def("SpMV", [](const CsrMatrix<double,IndexType>& rA,const Vector& x, Vector& y){
        rA.SpMV(x,y);
    })
    .def("__matmul__", [](const CsrMatrix<double,IndexType>& rA,const SystemVector<double,IndexType>& x){
        auto py  = std::make_shared<SystemVector<double,IndexType>>(rA.size1());
        py->SetValue(0.0);
        rA.SpMV(x,*py);
        return py;
    }, py::is_operator())
    .def("__matmul__", [](const CsrMatrix<double,IndexType>& rA,const Vector& x){
        auto py  = std::make_shared<SystemVector<double,IndexType>>(rA.size1());
        py->SetValue(0.0);
        rA.SpMV(x,*py);
        return py;
    }, py::is_operator())
    .def("TransposeSpMV", [](CsrMatrix<double,IndexType>& rA,SystemVector<double,IndexType>& x, SystemVector<double,IndexType>& y){
        rA.TransposeSpMV(x,y);
    })
    .def("SpMM", [](CsrMatrix<double,IndexType>& rA,CsrMatrix<double,IndexType>& rB){
        return AmgclCSRSpMMUtilities::SparseMultiply(rA,rB);
    })
    .def("__matmul__", [](CsrMatrix<double,IndexType>& rA,CsrMatrix<double,IndexType>& rB){
        return AmgclCSRSpMMUtilities::SparseMultiply(rA,rB);
    }, py::is_operator())
    .def("GetEquationIdCsrIndices", [](CsrMatrix<double,IndexType>& rA, const ElementsContainerType& rElements, const ProcessInfo& rProcessInfo) {
        return GetEquationIdCsrIndices(rA, rElements, rProcessInfo);
    }, py::return_value_policy::move)
    .def("GetEquationIdCsrIndices", [](CsrMatrix<double,IndexType>& rA, const ConditionsContainerType& rConditions, const ProcessInfo& rProcessInfo) {
        return GetEquationIdCsrIndices(rA, rConditions, rProcessInfo);
    }, py::return_value_policy::move)
    .def("Transpose", [](CsrMatrix<double,IndexType>& rA){
        return AmgclCSRConversionUtilities::Transpose<double,IndexType>(rA);
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
    .def("AssembleWithCsrIndices", [](CsrMatrix<double,IndexType>& rA, const NDData<double>& rLeftHandSideContributions, const NDData<int>& rEqIdsCsrIndices) {
        AssembleWithCsrIndices(rLeftHandSideContributions, rEqIdsCsrIndices, rA);
    })
    .def("AssembleWithEquationIds", [](CsrMatrix<double,IndexType>& rA, const NDData<double>& rLeftHandSideContributions, const NDData<int>& rEqIds) {
        AssembleWithEquationIds(rLeftHandSideContributions, rEqIds, rA);
    })
    .def("ApplyHomogeneousDirichlet", &CsrMatrix<double,IndexType>::ApplyHomogeneousDirichlet<Vector>)
    .def("ApplyHomogeneousDirichlet", &CsrMatrix<double,IndexType>::ApplyHomogeneousDirichlet<SystemVector<double,IndexType>>)
    .def("__str__", PrintObject<CsrMatrix<double,IndexType>>);



    //************************************************************************************************************
    //************************************************************************************************************
    py::class_<SystemVector<double,IndexType>, SystemVector<double,IndexType>::Pointer >(m, "SystemVector", py::buffer_protocol())
    .def(py::init<IndexType>())
    .def(py::init<IndexType, DataCommunicator&>())
    .def(py::init<SparseGraph<IndexType>&>())
    .def(py::init<SparseContiguousRowGraph<IndexType>&>())
    .def(py::init<Vector&>())
    .def("copy", [](const SystemVector<double,IndexType>& self){
        return Kratos::make_shared<SystemVector<double,IndexType>>(self);
    })
    .def("GetComm", &SystemVector<double,IndexType>::GetComm, py::return_value_policy::reference_internal)
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
    .def("Data", [](SystemVector<double,IndexType>& self) -> DenseVector<double>&{
       return self.data();
    }, py::return_value_policy::reference_internal)
    .def("SetData", [](SystemVector<double,IndexType>& self, const DenseVector<double>& other_data) { noalias(self.data()) = other_data; })
    .def("Assign", [](SystemVector<double,IndexType>& self, const SystemVector<double,IndexType>& other_vec){self = other_vec; } )
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
    .def("__mul__", [](const SystemVector<double,IndexType>& self, const double factor)
        {auto paux = std::make_shared<SystemVector<double,IndexType>>(self);
         (*paux) *= factor;
         return paux;}, py::is_operator())
    .def("__rmul__", [](const SystemVector<double,IndexType>& self, const double factor)
        {auto paux = std::make_shared<SystemVector<double,IndexType>>(self);
         (*paux) *= factor;
         return paux;}, py::is_operator())
    .def("__iadd__", [](SystemVector<double,IndexType>& self, const SystemVector<double,IndexType>& other_vec) -> SystemVector<double,IndexType>& {self += other_vec; return self;}, py::is_operator())
    .def("__isub__", [](SystemVector<double,IndexType>& self, const SystemVector<double,IndexType>& other_vec) -> SystemVector<double,IndexType>& {self -= other_vec;  return self;}, py::is_operator())
    .def("__imul__", [](SystemVector<double,IndexType>& self, const double& value) -> SystemVector<double,IndexType>& { self*=value; return self;}, py::is_operator())
    .def("__itruediv__", [](SystemVector<double,IndexType>& self, const double& value) -> SystemVector<double,IndexType>& { self/=value; return self;}, py::is_operator())
    //out of place
    .def("__add__", [](const SystemVector<double,IndexType>& vec1, const SystemVector<double,IndexType>& vec2){ //implies an internal copy
            auto paux = std::make_shared<SystemVector<double,IndexType>>(vec1);
            *paux += vec2;
            return paux;}
            , py::is_operator())
    .def("__sub__", [](const SystemVector<double,IndexType>& vec1, const SystemVector<double,IndexType>& vec2){ //implies an internal copy
            auto paux = std::make_shared<SystemVector<double,IndexType>>(vec1);
            *paux -= vec2;
            return paux;}
            , py::is_operator())
    //access operators
    .def("__setitem__", [](SystemVector<double,IndexType>& self, const unsigned int i, const double value){self[i] = value;} )
    .def("__getitem__", [](const SystemVector<double,IndexType>& self, const unsigned int i){return self[i];} )
    .def("__str__", PrintObject<SystemVector<double,IndexType>>)
    ;
}

} // namespace Kratos::Python.

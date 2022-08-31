//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_ADDITIVE_SCHWARZ_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_ADDITIVE_SCHWARZ_PRECONDITIONER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/preconditioner.h"
#include "utilities/variable_utils.h"
#include "utilities/math_utils.h"
#include "utilities/atomic_utilities.h"
#include "utilities/builtin_timer.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class AdditiveSchwarzPreconditioner
 * @todo Use dense solver for computation of inverse/ pseudo inverse.
 *       Check values for relaxation parameter.
 *       Only apply additive schwarz preconditioner for trimmed elements (badly conditioned elements).
 *       Apply Jacobi preconditioning for non-trimmed elements.
 *       Migrate preconditioner into Kratos Core.
 **/

template<class TSparseSpaceType, class TDenseSpaceType>
class AdditiveSchwarzPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of AdditiveSchwarzPreconditioner
    KRATOS_CLASS_POINTER_DEFINITION(AdditiveSchwarzPreconditioner);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename ModelPart::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AdditiveSchwarzPreconditioner()
    {
        mMatrixIsInitializedFlag = false;
    }

    /// Copy constructor.
    AdditiveSchwarzPreconditioner(const AdditiveSchwarzPreconditioner& Other) = delete;

    /// Destructor.
    ~AdditiveSchwarzPreconditioner() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AdditiveSchwarzPreconditioner& operator=(const AdditiveSchwarzPreconditioner& Other) = delete;

    ///@}
    ///@name Get/Set functions
    ///@{

    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    void ProvideAdditionalData(SparseMatrixType& rA, VectorType& rX, VectorType& rB,
                            DofsArrayType& rdof_set, ModelPart& r_model_part ) override
    {
        BuiltinTimer timer;

        if( ! mMatrixIsInitializedFlag ){
            InitializeMatrix(rA);
            mMatrixIsInitializedFlag = true;
        } else {
            TSparseSpaceType::SetToZero(*mpS);
        }

        std::vector<ModelPart::ElementIterator> elements_to_consider;
        elements_to_consider.reserve(r_model_part.NumberOfElements());

        // Find all unique elements. In case of NURBS patches, multiple elements (Quadrature Point Geometries)
        // might be located within the same knot span.
        VariableUtils().SetFlag(VISITED, false, r_model_part.Nodes());
        const auto element_it_begin = r_model_part.ElementsBegin();
        const auto geometry_type = GeometryData::KratosGeometryType::Kratos_Quadrature_Point_Geometry;
        for( IndexType i = 0; i < r_model_part.NumberOfElements(); ++i){
            auto element_it = element_it_begin + i;
            auto& r_geometry = element_it->GetGeometry();
            if( element_it->GetGeometry().GetGeometryType() == geometry_type ){
                if( r_geometry[0].IsNot(VISITED) ){
                    r_geometry[0].Set(VISITED,true);
                    elements_to_consider.push_back(element_it);
                }
            } else {
                elements_to_consider.push_back(element_it);
            }
        }

        // Assemble additive schwarz blocks
        const auto element_to_consider_it_begin = elements_to_consider.begin();
        IndexPartition<IndexType>(elements_to_consider.size()).for_each([&](IndexType i){
            auto element_it = element_to_consider_it_begin + i;

            auto& r_geometry = (*element_it)->GetGeometry();
            const SizeType number_nodes = r_geometry.size();
            const SizeType number_dofs_per_node = r_geometry[0].GetDofs().size();

            std::vector<std::size_t> equation_ids;
            equation_ids.reserve(number_nodes*number_dofs_per_node);
            for( IndexType j = 0; j < r_geometry.size(); ++j){
                for(auto& dof : r_geometry[j].GetDofs()){
                    equation_ids.push_back(dof->EquationId());
                }
            }
            // Relaxation parameter
            // Values from: Jomo et al., Hierarchical multigrid approaches for the finite cell
            // method on uniform and multi-level hp-refined grids, CMAME, 2010.
            // @TODO: Values must be adopted for overlapping "IGA additive schwarz blocks".
            double omega = 1;
            if( r_geometry.LocalSpaceDimension() == 2 ){
                omega = 1.0/6.0;
            } else if (r_geometry.LocalSpaceDimension() == 3){
                omega = 1.0/18.0;
            }
            // Assemble preconditioning matrix S.
            AssembleContributions(rA, *mpS, equation_ids, omega);
        });

        KRATOS_INFO("AdditiveSchwarzPreconditioner:") << "Build Matrix Time: " << timer.ElapsedSeconds() << std::endl;
    }

    ///@}
    ///@name Operations
    ///@{

    VectorType& ApplyInverseRight(VectorType& rX) override
    {
        return rX;
    }

    void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY) override
    {

        VectorType tmp = rX;
        TSparseSpaceType::Mult(rA, tmp, rY);

        ApplyLeft(rY);

    }

    VectorType& ApplyLeft(VectorType& rX) override
    {
        VectorType tmp = rX;
        TSparseSpaceType::Mult(*mpS, tmp, rX);

        return rX;
    }

    VectorType& Finalize(VectorType& rX) override{

        return rX;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Return information about this object.
    std::string Info() const override
    {
        return "AdditiveSchwarzPreconditioner";
    }

    /// Print information about this object.
    void  PrintInfo(std::ostream& OStream) const override
    {
        OStream << "AdditiveSchwarzPreconditioner";
    }

    void PrintData(std::ostream& OStream) const override
    {
    }

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{
    SparseMatrixPointerType mpS;
    bool mMatrixIsInitializedFlag = false;

    ///@}

private:

    //@name Private operations

    void InitializeMatrix(SparseMatrixType& rA){

        mpS = TSparseSpaceType::CreateEmptyMatrixPointer();
        SparseMatrixType& S = *mpS;

        S.resize(rA.size1(), rA.size1(), false);
        S.reserve(rA.nnz_capacity(), false);

        double* Avalues = S.value_data().begin();
        std::size_t* Arow_indices = S.index1_data().begin();
        std::size_t* Acol_indices = S.index2_data().begin();

        std::size_t* Arow_indices_old = rA.index1_data().begin();
        std::size_t* Acol_indices_old = rA.index2_data().begin();

        // Copy row indices
        IndexPartition<std::size_t>(rA.index1_data().size()).for_each([&](IndexType i){
            Arow_indices[i] = Arow_indices_old[i];
        });
        // Copy column indices
        IndexPartition<std::size_t>(rA.index2_data().size()).for_each([&](IndexType i){
            Acol_indices[i] = Acol_indices_old[i];
            Avalues[i] =  0.0;
        });

        S.set_filled(rA.filled1(),rA.filled2());
    }

    void AssembleContributions(
        SparseMatrixType& rA,
        SparseMatrixType& rS,
        Element::EquationIdVectorType& rEquationId,
        double omega
    )
    {
        const SizeType local_size = rEquationId.size();
        std::sort(rEquationId.begin(), rEquationId.end() );

        DenseMatrixType M = ZeroMatrix(local_size, local_size);
        for (IndexType i_local = 0; i_local < local_size; ++i_local) {
            const IndexType i_global = rEquationId[i_local];
            GetRowEntries(rA, M, i_global, i_local, rEquationId);
        }

        DenseMatrixType M2;
        M2.resize(local_size, local_size);
        for( IndexType i = 0; i < local_size; ++i){
            for( IndexType j = 0; j < local_size; ++j){
                M2(i,j) = rA(rEquationId[i],rEquationId[j]);
            }
        }

        double M_det = 0.0;
        DenseMatrixType M_invert = ZeroMatrix(local_size, local_size);
        M_invert = omega * M_invert;
        //@TODO Replace this by dense solver.
        MathUtils<double>::InvertMatrix(M, M_invert, M_det, -1e-6);

        for (IndexType i_local = 0; i_local < local_size; ++i_local) {
            const IndexType i_global = rEquationId[i_local];
            AssembleRowEntries(rS, M_invert, i_global, i_local, rEquationId);
        }
    }

    void AssembleRowEntries(SparseMatrixType& rA, const DenseMatrixType& rAlocal, const unsigned int i, const unsigned int i_local, std::vector<std::size_t>& rEquationId){
        double* values_vector = rA.value_data().begin();
        std::size_t* index1_vector = rA.index1_data().begin();
        std::size_t* index2_vector = rA.index2_data().begin();

        size_t left_limit = index1_vector[i];

        //find the first entry
        size_t last_pos = ForwardFind(rEquationId[0],left_limit,index2_vector);
        size_t last_found = rEquationId[0];

        double& r_a = values_vector[last_pos];
        const double& v_a = rAlocal(i_local,0);
        AtomicAdd(r_a,  v_a);

        //now find all of the other entries
        size_t pos = 0;
        for (unsigned int j=1; j<rEquationId.size(); j++) {
            unsigned int id_to_find = rEquationId[j];
            if(id_to_find > last_found) {
                pos = ForwardFind(id_to_find,last_pos+1,index2_vector);
            } else if(id_to_find < last_found) {
                pos = BackwardFind(id_to_find,last_pos-1,index2_vector);
            } else {
                pos = last_pos;
            }

            double& r = values_vector[pos];
            const double& v = rAlocal(i_local,j);
            AtomicAdd(r,  v);

            last_found = id_to_find;
            last_pos = pos;
        }
    }

    void GetRowEntries(SparseMatrixType& rA, DenseMatrixType& rAlocal, const unsigned int i, const unsigned int i_local, std::vector<std::size_t>& rEquationId){
        double* values_vector = rA.value_data().begin();
        std::size_t* index1_vector = rA.index1_data().begin();
        std::size_t* index2_vector = rA.index2_data().begin();

        size_t left_limit = index1_vector[i];

        // Find the first entry
        size_t last_pos = ForwardFind(rEquationId[0],left_limit,index2_vector);
        size_t last_found = rEquationId[0];

        const double& r_a = values_vector[last_pos];
        double& v_a = rAlocal(i_local,0);
        AtomicAdd(v_a, r_a);

        // Now find all of the other entries
        size_t pos = 0;
        for (unsigned int j=1; j<rEquationId.size(); j++) {
            unsigned int id_to_find = rEquationId[j];
            if(id_to_find > last_found) {
                pos = ForwardFind(id_to_find,last_pos+1,index2_vector);
            } else if(id_to_find < last_found) {
                pos = BackwardFind(id_to_find,last_pos-1,index2_vector);
            } else {
                pos = last_pos;
            }

            const double& r = values_vector[pos];
            double& v = rAlocal(i_local,j);
            AtomicAdd(v,  r);

            last_found = id_to_find;
            last_pos = pos;
        }
    }

    inline unsigned int ForwardFind(const unsigned int id_to_find,
                                    const unsigned int start,
                                    const std::size_t* index_vector)
    {
        unsigned int pos = start;
        while(id_to_find != index_vector[pos]) pos++;
        return pos;
    }

    inline unsigned int BackwardFind(const unsigned int id_to_find,
                                     const unsigned int start,
                                     const std::size_t* index_vector)
    {
        unsigned int pos = start;
        while(id_to_find != index_vector[pos]) pos--;
        return pos;
    }

}; // Class AdditiveSchwarzPreconditioner

///@}
///@}

///@name Input and output
///@{

/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::istream& operator >> (std::istream& IStream,
                                  AdditiveSchwarzPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const AdditiveSchwarzPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_ADDITIVE_SCHWARZ_PRECONDITIONER_H_INCLUDED  defined


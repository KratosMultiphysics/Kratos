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

namespace Kratos
{

///@name  Preconditioners
///@{

/// AdditiveSchwarzPreconditioner class.
/**   */
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

    typedef typename TSparseSpace::MatrixPointerType SparseMatrixPointerType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename ModelPart::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AdditiveSchwarzPreconditioner()
    {

    }


    /// Copy constructor.
    AdditiveSchwarzPreconditioner(const AdditiveSchwarzPreconditioner& Other) {}


    /// Destructor.
    ~AdditiveSchwarzPreconditioner() override
    {
        //@TODO
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AdditiveSchwarzPreconditioner& operator=(const AdditiveSchwarzPreconditioner& Other)
    {

    }

    void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {



        //ApplyLeft(rY);
    }
    ///@}
    ///@name Get/Set functions
    ///@{
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    void GetEntries(
        SparseMatrixType& rA,
        SparseMatrixType& rS,
        Element::EquationIdVectorType& rEquationId
    )
    {
        const SizeType local_size = rEquationId.size();
        std::sort(rEquationId.begin(), rEquationId.end() );

        DenseMatrixType M = ZeroMatrix(local_size, local_size);
        for (IndexType i_local = 0; i_local < local_size; i_local++) {
            const IndexType i_global = rEquationId[i_local];
            GetRowEntries(rA, M, i_global, i_local, rEquationId);
        }

        DenseMatrixType M2;
        M2.resize(local_size, local_size);
        for( int i = 0; i < local_size; ++i){
            for( int j = 0; j < local_size; ++j){
                M2(i,j) = rA(rEquationId[i],rEquationId[j]);
            }
        }

        double M_det = 0.0; //MathUtils<double>::Det(M);
        DenseMatrixType M_invert = ZeroMatrix(local_size, local_size);
        MathUtils<double>::InvertMatrix(M, M_invert, M_det, -1e-6);

        for (IndexType i_local = 0; i_local < local_size; i_local++) {
            const IndexType i_global = rEquationId[i_local];
            AssembleRowEntries(rS, M_invert, i_global, i_local, rEquationId);
        }
    }

    void GetRowEntries(SparseMatrixType& rA, DenseMatrixType& Alocal, const unsigned int i, const unsigned int i_local, std::vector<std::size_t>& EquationId){
        double* values_vector = rA.value_data().begin();
        std::size_t* index1_vector = rA.index1_data().begin();
        std::size_t* index2_vector = rA.index2_data().begin();

        size_t left_limit = index1_vector[i];
//    size_t right_limit = index1_vector[i+1];

        //find the first entry
        size_t last_pos = ForwardFind(EquationId[0],left_limit,index2_vector);
        size_t last_found = EquationId[0];

        const double& r_a = values_vector[last_pos];
        double& v_a = Alocal(i_local,0);
        AtomicAdd(v_a, r_a);

        //now find all of the other entries
        size_t pos = 0;
        for (unsigned int j=1; j<EquationId.size(); j++) {
            unsigned int id_to_find = EquationId[j];
            if(id_to_find > last_found) {
                pos = ForwardFind(id_to_find,last_pos+1,index2_vector);
            } else if(id_to_find < last_found) {
                pos = BackwardFind(id_to_find,last_pos-1,index2_vector);
            } else {
                pos = last_pos;
            }

            const double& r = values_vector[pos];
            double& v = Alocal(i_local,j);
            AtomicAdd(v,  r);

            last_found = id_to_find;
            last_pos = pos;
        }
    }

    void AssembleRowEntries(SparseMatrixType& rA, const DenseMatrixType& Alocal, const unsigned int i, const unsigned int i_local, std::vector<std::size_t>& EquationId){
        double* values_vector = rA.value_data().begin();
        std::size_t* index1_vector = rA.index1_data().begin();
        std::size_t* index2_vector = rA.index2_data().begin();

        size_t left_limit = index1_vector[i];
//    size_t right_limit = index1_vector[i+1];

        //find the first entry
        size_t last_pos = ForwardFind(EquationId[0],left_limit,index2_vector);
        size_t last_found = EquationId[0];

        double& r_a = values_vector[last_pos];
        const double& v_a = Alocal(i_local,0);
        AtomicAdd(r_a,  v_a);

        //now find all of the other entries
        size_t pos = 0;
        for (unsigned int j=1; j<EquationId.size(); j++) {
            unsigned int id_to_find = EquationId[j];
            if(id_to_find > last_found) {
                pos = ForwardFind(id_to_find,last_pos+1,index2_vector);
            } else if(id_to_find < last_found) {
                pos = BackwardFind(id_to_find,last_pos-1,index2_vector);
            } else {
                pos = last_pos;
            }

            double& r = values_vector[pos];
            const double& v = Alocal(i_local,j);
            AtomicAdd(r,  v);

            last_found = id_to_find;
            last_pos = pos;
        }
    }

    void InitializeMatrix(SparseMatrixType& rA){

        mpS = TSparseSpaceType::CreateEmptyMatrixPointer();
        SparseMatrixType& S = *mpS;
        S.resize(rA.size1(), rA.size2(), rA.nnz());

        double* Avalues = S.value_data().begin();
        std::size_t* Arow_indices = S.index1_data().begin();
        std::size_t* Acol_indices = S.index2_data().begin();

        std::size_t* Arow_indices_old = rA.index1_data().begin();
        std::size_t* Acol_indices_old = rA.index2_data().begin();
        //filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
        Arow_indices[0] = 0;
        for (int i = 0; i < static_cast<int>(rA.size1()); i++) {
            Arow_indices[i+1] = Arow_indices_old[i+1];
        }

        IndexPartition<std::size_t>(rA.size1()).for_each([&](std::size_t i){
            const unsigned int row_begin = Arow_indices_old[i];
            //const unsigned int row_end = Arow_indices[i+1];
            unsigned int k = row_begin;
            for (int i = 0; i < rA.index1_data().size() ; ++i) {
                Acol_indices[k] = Acol_indices_old[k];
                Avalues[k] = 0.0;
                k++;
            }

            // indices[i].clear(); //deallocating the memory

            // //std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);

        });

        S.set_filled(rA.size1()+1, rA.nnz());
    }

    void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        DofsArrayType& rdof_set,
        ModelPart& r_model_part
    ) override
    {

        InitializeMatrix(rA);
        //mpS = rA; //SparseMatrixType(rA.size1(), rA.size2(), rA.nnz_capacity() );
        TSparseSpaceType::SetToZero(mpS);

        VariableUtils().SetFlag(VISITED, false, r_model_part.Nodes());
        auto element_it_begin = r_model_part.ElementsBegin();
        int count = 0;
        for( IndexType i = 0; i < r_model_part.NumberOfElements(); ++i){
            auto element_it = element_it_begin + i;
            //if( element_it->GetGeometry().IsQuadraturePoint() ){
                auto& r_geometry = element_it->GetGeometry();
                std::vector<std::size_t> new_equation_ids;
                if( r_geometry[0].IsNot(VISITED) ){
                    r_geometry[0].Set(VISITED,true);
                    count++;
                    const SizeType dimension = r_geometry.WorkingSpaceDimension();
                    bool rotation_dofs =r_geometry[0].HasDofFor(ROTATION_X);
                    for( IndexType j = 0; j < r_geometry.size(); ++j){
                        new_equation_ids.push_back( r_geometry[j].GetDof(DISPLACEMENT_X).EquationId() );
                        new_equation_ids.push_back( r_geometry[j].GetDof(DISPLACEMENT_Y).EquationId() );
                        if( dimension == 3 ){
                            new_equation_ids.push_back( r_geometry[j].GetDof(DISPLACEMENT_Z).EquationId() );
                        }
                        if( rotation_dofs ){
                            new_equation_ids.push_back( r_geometry[j].GetDof(ROTATION_X).EquationId() );
                            new_equation_ids.push_back( r_geometry[j].GetDof(ROTATION_Y).EquationId() );
                            if( dimension == 3 ){
                                new_equation_ids.push_back( r_geometry[j].GetDof(ROTATION_Z).EquationId() );
                            }
                        }
                    }
                    GetEntries(rA, mpS, new_equation_ids);
                    //dof_blocks.push_back(new_equation_ids);
                }
        }
    }

    ///@}
    ///@name Operations
    ///@{

    VectorType& ApplyInverseRight(VectorType& rX) override
    {
        // VectorType z = rX;
        // TSparseSpaceType::Mult(mpS, z, rX);
        return rX;
    }


    void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY) override
    {

        VectorType zz = rX;
        TSparseSpaceType::Mult(rA, zz, rY);

        ApplyLeft(rY);

    }

    void TransposeMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY) override
    {
        // VectorType z = rX;
        // ApplyTransposeLeft(z);
        // TSparseSpaceType::TransposeMult(rA,z, rY);
    }

    /** multiply first rX by L^-1 and store result in temp
        then multiply temp by U^-1 and store result in rX
        @param rX  Unknows of preconditioner suystem
    */
    VectorType& ApplyLeft(VectorType& rX) override
    {
        // for( int i = 0; i < mpS.size1(); ++i){
        //     std::cout << mpS(i,i) << std::endl;
        // }
        VectorType z = rX;
        TSparseSpaceType::Mult(mpS, z, rX);

        return rX;
    }

    VectorType& Finalize(VectorType& rX) override{
        TSparseSpaceType::Clear(mpS);
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
    std::vector<std::vector<size_t>> dof_blocks;
    SparseMatrixPointerType mpS;

    ///@}

private:

    //@name Private operations

    inline unsigned int ForwardFind(const unsigned int id_to_find,
                                    const unsigned int start,
                                    const size_t* index_vector)
    {
        unsigned int pos = start;
        while(id_to_find != index_vector[pos]) pos++;
        return pos;
    }

    inline unsigned int BackwardFind(const unsigned int id_to_find,
                                     const unsigned int start,
                                     const size_t* index_vector)
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


//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes
#include <cmath>
#include <complex>

// External includes

// Project includes
#include "factories/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class ScalingSolver
 * @ingroup KratosCore
 * @brief This solvers rescales in order to improve the conditioning of the system
 * @details Rescales the matrix, and uses a given linear solver
 * @author Riccardo Rossi
 * @tparam TSparseSpaceType The sparse space definition
 * @tparam TDenseSpaceType The dense space definition
 * @tparam TReordererType The reorder considered
 */
template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class ScalingSolver
    : public LinearSolver<TSparseSpaceType, TDenseSpaceType,  TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ScalingSolver
    KRATOS_CLASS_POINTER_DEFINITION(ScalingSolver);

    /// Definition of the base type
    using BaseType = LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>;

    /// The definition of the spaces (sparse matrix)
    using SparseMatrixType = typename TSparseSpaceType::MatrixType;

    /// The definition of the spaces (vector)
    using VectorType = typename TSparseSpaceType::VectorType;

    /// The definition of the spaces (dense matrix)
    using DenseMatrixType = typename TDenseSpaceType::MatrixType;

    /// The definition of the linear solver factory type
    using LinearSolverFactoryType = LinearSolverFactory<TSparseSpaceType, TDenseSpaceType>;

    /// The index type definition to be consistent
    using IndexType = typename TSparseSpaceType::IndexType;

    /// Definition of the index iterator type
    using IndexIteratorType = typename boost::numeric::ublas::compressed_matrix<typename TSparseSpaceType::DataType>::index_array_type::iterator;

    /// Definition of the const index iterator type
    using ConstIndexIteratorType = typename boost::numeric::ublas::compressed_matrix<typename TSparseSpaceType::DataType>::index_array_type::const_iterator;

    /// Definition of the value iterator type
    using ValueIteratorType = typename boost::numeric::ublas::compressed_matrix<typename TSparseSpaceType::DataType>::value_array_type::iterator;

    /// Definition of the const value iterator type
    using ConstValueIteratorType = typename boost::numeric::ublas::compressed_matrix<typename TSparseSpaceType::DataType>::value_array_type::const_iterator;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ScalingSolver()
    {
    }

    /**
     * @brief Constructor without parameters
     * @param pLinearSolver The linear solver to be scaled
     * @param SymmetricScaling If the scaling is symmetric (true by default)
     */
    ScalingSolver(
        typename BaseType::Pointer pLinearSolver,
        const bool SymmetricScaling = true
        ) : BaseType (),
            mpLinearSolver(pLinearSolver),
            mSymmetricScaling(SymmetricScaling)
    {
    }

    /**
     * @brief Constructor with parameters
     * @param ThisParameters The configuration parameters of the linear solver
     */
    ScalingSolver(Parameters ThisParameters)
        : BaseType ()
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(ThisParameters.Has("solver_type")) << "Solver_type must be specified to construct the ScalingSolver" << std::endl;

        mpLinearSolver = LinearSolverFactoryType().Create(ThisParameters);

        mSymmetricScaling = ThisParameters.Has("symmetric_scaling") ? ThisParameters["symmetric_scaling"].GetBool() : true;

        KRATOS_CATCH("")
    }

    /// Copy constructor.
    ScalingSolver(const ScalingSolver& Other) : BaseType(Other) {}


    /// Destructor.
    ~ScalingSolver() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ScalingSolver& operator=(const ScalingSolver& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
    * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
    * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
    * which require knowledge on the spatial position of the nodes associated to a given dof.
    * This function tells if the solver requires such data
    */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return mpLinearSolver->AdditionalPhysicalDataIsNeeded();
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
    void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    ) override
    {
        mpLinearSolver->ProvideAdditionalData(rA,rX,rB,rdof_set,r_model_part);
    }

    void InitializeSolutionStep (SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        mpLinearSolver->InitializeSolutionStep(rA,rX,rB);
    }

    /** This function is designed to be called at the end of the solve step.
     * for example this is the place to remove any data that we do not want to save for later
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    void FinalizeSolutionStep (SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        mpLinearSolver->FinalizeSolutionStep(rA,rX,rB);
    }

    /** This function is designed to clean up all internal data in the solver.
     * Clear is designed to leave the solver object as if newly created.
     * After a clear a new Initialize is needed
     */
    void Clear() override
    {
        mpLinearSolver->Clear();
    }

    /** Normal solve method.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

        VectorType scaling_vector(rX.size());

        //obtain the scaling matrix
        GetScalingWeights(rA,scaling_vector);

        //scale system
        if(mSymmetricScaling == false)
        {
            KRATOS_THROW_ERROR(std::logic_error,"not yet implemented","")
        }
        else
        {
            IndexPartition<std::size_t>(scaling_vector.size()).for_each([&](std::size_t Index){
                scaling_vector[Index] = sqrt(std::abs(scaling_vector[Index]));
            });

            SymmetricScaling(rA,scaling_vector);
        }

        //scale RHS
        IndexPartition<std::size_t>(scaling_vector.size()).for_each([&](std::size_t Index){
            rB[Index] /= scaling_vector[Index];
        });

        //solve the problem
        bool is_solved = mpLinearSolver->Solve(rA,rX,rB);

        //backscale the solution
        if(mSymmetricScaling == true)
        {
            IndexPartition<std::size_t>(scaling_vector.size()).for_each([&](std::size_t Index){
                rX[Index] /= scaling_vector[Index];
            });
        }

        return is_solved;
    }

    ///@}
    ///@name Access
    ///@{

    IndexType GetIterationsNumber() override
    {
        return mpLinearSolver->GetIterationsNumber();
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Composite Linear Solver. Uses internally the following linear solver " << mpLinearSolver->Info();
        return  buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    typename LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>::Pointer mpLinearSolver; /// Pointer to the internal linear solver
    bool mSymmetricScaling; /// Flag to indicate if the scaling is symmetric or not

    ///@}
    ///@name Private Operators
    ///@{
    static void SymmetricScaling( SparseMatrixType& A, const VectorType& aux)
    {
        IndexPartition<std::size_t>(A.size1()).for_each([&](std::size_t Index){
            // Get the row iterator, index2 iterator and value iterator for the current row (Index)
            auto row_begin_it = A.index1_data().begin() + Index;
            auto index2_begin_it = A.index2_data().begin() + (*row_begin_it);
            auto value_begin_it = A.value_data().begin() + (*row_begin_it);
            int number_of_entries_in_row = *(row_begin_it+1) - *row_begin_it;

            // Call perform_matrix_scaling_row for the current row
            perform_matrix_scaling_row(
                Index, // Current row index
                number_of_entries_in_row,
                row_begin_it,
                index2_begin_it, // Iterator to the beginning of the current row's data in index2_data
                value_begin_it, // Iterator to the beginning of the current row's data in value_data
                aux
            );
        });
    }

    static void perform_matrix_scaling_row(
        std::size_t CurrentRowIndex,
        int number_of_entries_in_row,
        IndexIteratorType row_it, // Should be index1_data().begin() + CurrentRowIndex
        IndexIteratorType index2_begin,
        ValueIteratorType value_begin,
        const VectorType& weights
    )
    {
        const typename TSparseSpaceType::DataType row_weight = weights[CurrentRowIndex];

        for(int i = 0; i < number_of_entries_in_row; i++) {
            const typename TSparseSpaceType::DataType col_weight = weights[*index2_begin];
            typename TSparseSpaceType::DataType t = (*value_begin);
            t /= (row_weight*col_weight);
            (*value_begin) = t;
            value_begin++;
            index2_begin++;
        }
    }

    static void GetScalingWeights( const SparseMatrixType& A, VectorType& aux)
    {
        IndexPartition<std::size_t>(A.size1()).for_each([&A, &aux](std::size_t Index){
            ConstIndexIteratorType row_begin_it = A.index1_data().begin() + Index;
            ConstIndexIteratorType index2_begin_it = A.index2_data().begin() + (*row_begin_it); // Not strictly needed for GS2weights_row, but kept for consistency
            ConstValueIteratorType value_begin_it = A.value_data().begin() + (*row_begin_it);
            int number_of_entries_in_row = *(row_begin_it+1) - *row_begin_it;

            GS2weights_row(
                Index, // Current row index
                number_of_entries_in_row,
                row_begin_it,
                index2_begin_it, // Not used in current GS2weights_row implementation
                value_begin_it,
                aux
            );
        });
    }

    static void GS2weights_row(
        std::size_t CurrentRowIndex,
        int number_of_entries_in_row,
        ConstIndexIteratorType row_it,
        ConstIndexIteratorType index2_begin, // Not used in original GS2weights, but kept for consistency if needed later
        ConstValueIteratorType value_begin,
        VectorType& weights
    )
    {
        double t = 0.0;

        for(int i = 0; i < number_of_entries_in_row; i++)
        {
            double tmp = std::abs(*value_begin);
            t += tmp*tmp;
            value_begin++;
        }
        t = std::sqrt(t);
        weights[CurrentRowIndex] = t;
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; // Class ScalingSolver

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::istream& operator >> (std::istream& IStream,
                                  ScalingSolver<TSparseSpaceType, TDenseSpaceType,
                                  TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const ScalingSolver<TSparseSpaceType, TDenseSpaceType,
                                  TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}

}  // namespace Kratos.
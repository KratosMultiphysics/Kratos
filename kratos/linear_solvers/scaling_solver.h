//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_SCALING_SOLVER_H_INCLUDED )
#define  KRATOS_SCALING_SOLVER_H_INCLUDED

// System includes
#include <cmath>
#include <complex>

// External includes

// Project includes
#include "includes/define.h"
#include "factories/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"

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
    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    /// The definition of the spaces (sparse matrix)
    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    /// The definition of the spaces (vector)
    typedef typename TSparseSpaceType::VectorType VectorType;

    /// The definition of the spaces (dense matrix)
    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /// The definition of the linear solver factory type
    typedef LinearSolverFactory<TSparseSpaceType,TDenseSpaceType> LinearSolverFactoryType;

    /// The index type definition to be consistent
    typedef typename TSparseSpaceType::IndexType IndexType;

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
            #pragma omp parallel for
            for(int i=0; i< static_cast<int>(scaling_vector.size()); i++)
                scaling_vector[i] = sqrt(std::abs(scaling_vector[i]));

            SymmetricScaling(rA,scaling_vector);


        }

        //scale RHS
        #pragma omp parallel for
        for(int i=0; i< static_cast<int>(scaling_vector.size()); i++)
            rB[i] /= scaling_vector[i];


        //solve the problem
        bool is_solved = mpLinearSolver->Solve(rA,rX,rB);

        //backscale the solution
        if(mSymmetricScaling == true)
        {
            #pragma omp parallel for
            for(int i=0; i< static_cast<int>(scaling_vector.size()); i++)
                rX[i] /= scaling_vector[i];
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
    typename LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>::Pointer mpLinearSolver;
    bool mSymmetricScaling;

    ///@}
    ///@name Private Operators
    ///@{
    static void SymmetricScaling( SparseMatrixType& A, const VectorType& aux)
    {
      //typedef  unsigned int size_type;
        //typedef  double value_type;

        //create partition
        OpenMPUtils::PartitionVector partition;
        int number_of_threads = ParallelUtilities::GetNumThreads();
        OpenMPUtils::DivideInPartitions(A.size1(),number_of_threads,  partition);
        //parallel loop

        #pragma omp parallel
        {
            int thread_id = OpenMPUtils::ThisThread();
            int number_of_rows = partition[thread_id+1] - partition[thread_id];
            typename boost::numeric::ublas::compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::iterator row_iter_begin = A.index1_data().begin()+partition[thread_id];
            typename boost::numeric::ublas::compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::iterator index_2_begin = A.index2_data().begin()+*row_iter_begin;
            typename boost::numeric::ublas::compressed_matrix<typename TDenseSpaceType::DataType>::value_array_type::iterator value_begin = A.value_data().begin()+*row_iter_begin;

            perform_matrix_scaling(    number_of_rows,
                                       row_iter_begin,
                                       index_2_begin,
                                       value_begin,
                                       partition[thread_id],
                                       aux
                                  );
        }
    }

    /**
     * calculates partial product resetting to Zero the output before
     */
    static void perform_matrix_scaling(
        int number_of_rows,
        typename boost::numeric::ublas::compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::iterator row_begin,
        typename boost::numeric::ublas::compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::iterator index2_begin,
        typename boost::numeric::ublas::compressed_matrix<typename TDenseSpaceType::DataType>::value_array_type::iterator value_begin,
        unsigned int output_begin_index,
        const VectorType& weights
    )
    {
        int row_size;
        typename SparseMatrixType::index_array_type::const_iterator row_it = row_begin;
        int kkk = output_begin_index;
        for(int k = 0; k < number_of_rows; k++)
        {
            row_size= *(row_it+1)-*row_it;
            row_it++;
            const typename TDenseSpaceType::DataType row_weight = weights[kkk++];

            for(int i = 0; i<row_size; i++)
            {
                const typename TDenseSpaceType::DataType col_weight = weights[*index2_begin];
                typename TDenseSpaceType::DataType t = (*value_begin);
                t /= (row_weight*col_weight);
                (*value_begin) = t; //check if this is correcct!!
                value_begin++;
                index2_begin++;
            }

        }
    }


    static void GetScalingWeights( const SparseMatrixType& A, VectorType& aux)
    {
      //typedef  unsigned int size_type;
      //typedef  double value_type;

        //create partition
        OpenMPUtils::PartitionVector partition;
        int number_of_threads = ParallelUtilities::GetNumThreads();
        OpenMPUtils::DivideInPartitions(A.size1(),number_of_threads,  partition);
        //parallel loop

        #pragma omp parallel
        {
            int thread_id = OpenMPUtils::ThisThread();
            int number_of_rows = partition[thread_id+1] - partition[thread_id];
            typename boost::numeric::ublas::compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::const_iterator row_iter_begin = A.index1_data().begin()+partition[thread_id];
            typename boost::numeric::ublas::compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::const_iterator index_2_begin = A.index2_data().begin()+*row_iter_begin;
            typename boost::numeric::ublas::compressed_matrix<typename TDenseSpaceType::DataType>::value_array_type::const_iterator value_begin = A.value_data().begin()+*row_iter_begin;

            GS2weights(    number_of_rows,
                           row_iter_begin,
                           index_2_begin,
                           value_begin,
                           partition[thread_id],
                           aux
                      );
        }
    }

    /**
     * calculates partial product resetting to Zero the output before
     */
    static void GS2weights(
        int number_of_rows,
        typename boost::numeric::ublas::compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::const_iterator row_begin,
        typename boost::numeric::ublas::compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::const_iterator index2_begin,
        typename boost::numeric::ublas::compressed_matrix<typename TDenseSpaceType::DataType>::value_array_type::const_iterator value_begin,
        unsigned int output_begin_index,
        VectorType& weights
    )
    {
        int row_size;
        typename SparseMatrixType::index_array_type::const_iterator row_it = row_begin;
        int kkk = output_begin_index;
        for(int k = 0; k < number_of_rows; k++)
        {
            row_size= *(row_it+1)-*row_it;
            row_it++;
            double t = 0.0;

            for(int i = 0; i<row_size; i++)
            {
                double tmp = std::abs(*value_begin);
                t += tmp*tmp;
                value_begin++;
            }
            t = sqrt(t);
            weights[kkk++] = t;
        }
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

#endif // KRATOS_SCALING_SOLVER_H_INCLUDED  defined



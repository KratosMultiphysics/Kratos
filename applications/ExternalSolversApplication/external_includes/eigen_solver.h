#if !defined(KRATOS_EIGEN_SOLVER_H_INCLUDED)
#define  KRATOS_EIGEN_SOLVER_H_INCLUDED

// External includes
#include "boost/smart_ptr.hpp"

#include "includes/ublas_interface.h"
#include "Eigen/Core"
#include "Eigen/Sparse"

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"

#include <chrono>
using namespace std::chrono;

namespace ublas = boost::numeric::ublas;

namespace Kratos
{

struct EigenSolverTypes
{
    typedef typename Eigen::SparseMatrix<double, Eigen::RowMajor, int> SparseMatrix;
        
    typedef typename Eigen::SparseLU<SparseMatrix> SparseLU;

    typedef typename Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper> ConjugateGradient;
};

template<
    class TSolverType,
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>
>
class EigenSolver
: public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of EigenSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(EigenSolver);

    typedef DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
    
    
    EigenSolver() {}

    /**
     * Default constructor
     */
    EigenSolver(Parameters settings): BaseType(settings) {}

    /**
     * Destructor
     */
    ~EigenSolver() override {}

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {        
        // --- begin time measure code
        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        std::vector<int> index1_vector(rA.index1_data().size());
        std::vector<int> index2_vector(rA.index2_data().size());
        // --- end time measure code

        for (int unsigned i = 0; i < rA.index1_data().size(); i++) {
            index1_vector[i] = (int)rA.index1_data()[i];
        }

        for (unsigned int i = 0; i < rA.index2_data().size(); i++) {
            index2_vector[i] = (int)rA.index2_data()[i];
        }
    
        Eigen::Map<EigenSolverTypes::SparseMatrix> a(rA.size1(), rA.size2(), rA.nnz(), index1_vector.data(), index2_vector.data(), rA.value_data().begin());
        Eigen::Map<Eigen::VectorXd> x(rX.data().begin(), rX.size());
        Eigen::Map<Eigen::VectorXd> b(rB.data().begin(), rB.size());

        TSolverType solver;
        solver.compute(a);
        x = solver.solve(b);

        // --- begin time measure code
        high_resolution_clock::time_point t2 = high_resolution_clock::now();

        double solver_time = duration_cast<microseconds>( t2 - t1 ).count() / 1000000.0;
        
        KRATOS_WATCH(solver_time);
        // --- end time measure code
               
        bool success = (solver.info() == Eigen::Success);
 
        return success;
    }

    /**
     * Multi solve method for solving a set of linear systems with same coefficient matrix.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        bool is_solved = true;

        return is_solved;
    }

    /**
     * Print information about this object.
     */
    void  PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Eigen solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const override
    {
    }

private:

    /**
     * Assignment operator.
     */
    EigenSolver& operator=(const EigenSolver& Other);

    /**
     * Copy constructor.
     */
    EigenSolver(const EigenSolver& Other);

}; // Class SkylineLUFactorizationSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, EigenSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const EigenSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_EIGEN_SOLVER_H_INCLUDED  defined 



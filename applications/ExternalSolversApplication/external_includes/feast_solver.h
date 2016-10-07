//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:		 BSD License
//   Kratos default license: kratos/license.txt
//
//   Project Name:        $ExternalSolversApplication   $
//   Last modified by:    $Author: michael.andre@tum.de $
//   Date:                $Date:         September 2016 $
//   Revision:            $Revision:                0.0 $
//
//

// System includes
#include <iostream>
#include <complex>
#include <vector>
#include <unordered_set>
#include <algorithm>

// External includes
#include <boost/smart_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/value_type/complex.hpp>
#include <amgcl/solver/skyline_lu.hpp>
extern "C" {
#include "feast.h"
}

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"

#if !defined(KRATOS_FEAST_SOLVER )
#define  KRATOS_FEAST_SOLVER

namespace Kratos {

///@name Kratos Classes
///@{

/// Sparse system matrix required for FEAST's inner iterations.
/**
 *  Matrix arrays are accessed the same as for ublas compressed_matrix.
 */
struct FEASTSystemMatrix
{
    typedef std::complex<double> ScalarType;
    typedef boost::numeric::ublas::compressed_matrix<double> RealMatrixType;

    /**
     * The crs matrix structure remains constant during the inner iterations.
     * It is initialized here.
     */
    void Initialize(const RealMatrixType& M, const RealMatrixType& K)
    {
        const std::size_t dimension = M.size1();

        row_ptr.resize(dimension + 1);

        std::vector<std::unordered_set<std::size_t> > indices(dimension);

        row_ptr[0] = 0;
        for (std::size_t i = 0; i < dimension; i++)
        {
            std::size_t row_begin, row_end;
            indices[i].reserve(40);

            row_begin = M.index1_data()[i];
            row_end = M.index1_data()[i + 1];

            indices[i].insert(M.index2_data().begin() + row_begin,
                    M.index2_data().begin() + row_end);

            row_begin = K.index1_data()[i];
            row_end = K.index1_data()[i + 1];

            indices[i].insert(K.index2_data().begin() + row_begin,
                    K.index2_data().begin() + row_end);

            row_ptr[i + 1] = row_ptr[i] + indices[i].size();
        }

        col_idx.resize(row_ptr[dimension]);

        size_t k = 0;
        for (std::size_t i = 0; i < dimension; i++)
        {
            std::for_each(std::begin(indices[i]), std::end(indices[i]),
                    [&](std::size_t j) {col_idx[k++] = j;});

            indices[i].clear();

            std::sort(col_idx.begin() + row_ptr[i],
                    col_idx.begin() + row_ptr[i + 1]);
        }

        values.resize(row_ptr[dimension]);
    }

    /**
     * Similar to FEAST's zdaddcsr subroutine. Assumes crs matrix structure is
     * already initialized.
     */
    void Calculate(const RealMatrixType& M, const RealMatrixType& K, ScalarType z)
    {
        std::size_t jm, jk;
        double mij, kij;
        const std::size_t dimension = M.size1();

        std::size_t ptr = 0;
        for (std::size_t i = 0; i < dimension; i++)
        {
            std::size_t m_ptr = M.index1_data()[i];
            std::size_t k_ptr = K.index1_data()[i];
            while (m_ptr < M.index1_data()[i + 1] || k_ptr < K.index1_data()[i + 1])
            {
                jm = (m_ptr < M.index1_data()[i + 1]) ?
                        M.index2_data()[m_ptr] : dimension;
                jk = (k_ptr < K.index1_data()[i + 1]) ?
                        K.index2_data()[k_ptr] : dimension;

                if (jm < jk)
                {
                    mij = M(i, jm);
                    values[ptr] = z * mij;
                    m_ptr++;
                }
                else if (jm > jk)
                {
                    kij = K(i, jk);
                    values[ptr] = -kij;
                    k_ptr++;
                }
                else
                { // jm == jk
                    mij = M(i, jm);
                    kij = K(i, jk);
                    values[ptr] = z * mij - kij;
                    m_ptr++;
                    k_ptr++;
                }
                ptr++;
            }
        }
    }

    boost::numeric::ublas::unbounded_array<std::ptrdiff_t>&
    index1_data()
    {
        return row_ptr;
    }

    boost::numeric::ublas::unbounded_array<std::ptrdiff_t>&
    index2_data()
    {
        return col_idx;
    }

    boost::numeric::ublas::unbounded_array<ScalarType>&
    value_data()
    {
        return values;
    }

    std::size_t size1() const
    {
        return row_ptr.size() - 1;
    }

    std::size_t size2() const
    {
        return size1();
    }

    boost::numeric::ublas::unbounded_array<std::ptrdiff_t> row_ptr;
    boost::numeric::ublas::unbounded_array<std::ptrdiff_t> col_idx;
    boost::numeric::ublas::unbounded_array<ScalarType> values;
};

template<typename TScalarType>
struct SkylineLUSolver
{
    typedef typename amgcl::backend::builtin<TScalarType>::matrix MatrixType;
    typedef amgcl::solver::skyline_lu<TScalarType> SolverType;

    boost::shared_ptr<MatrixType> pBuiltinMatrix;

    boost::shared_ptr<SolverType> pSolver;

    ~SkylineLUSolver()
    {
        Clear();
    }

    template<class Matrix>
    void Factorize(Matrix &A)
    {
        Clear();

        pBuiltinMatrix = amgcl::adapter::zero_copy(
                A.size1(),
                A.index1_data().begin(),
                A.index2_data().begin(),
                A.value_data().begin());

        pSolver = boost::make_shared<SolverType>(*pBuiltinMatrix);
    }

    void Solve(std::vector<TScalarType> &b, std::vector<TScalarType> &x)
    {
        (*pSolver)(b, x);
    }

    void Clear()
    {
        pSolver.reset();
        pBuiltinMatrix.reset();
    }
};

/// Adapter to FEAST eigenvalue problem solver.
template<class TSparseSpaceType, class TDenseSpaceType,
        class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class FEASTSolver: public LinearSolver<TSparseSpaceType, TDenseSpaceType,
        TReordererType> {

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( FEASTSolver );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType SparseVectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /**
     * Parameters let the user control the settings of the FEAST library.
     */
    FEASTSolver(Parameters::Pointer pParam)
    : mpParam(pParam)
    {
        Parameters default_params(R"(
        {
            "solver_type": "FEAST",
            "print_feast_output": false,
            "perform_stochastic_estimate": true,
            "solve_eigenvalue_problem": true,
            "lambda_min": 0.0,
            "lambda_max": 1.0,
            "number_of_eigenvalues": 0,
            "search_dimension": 10,
            "linear_solver_settings": {
                "solver.type": "skyline_lu"
            }
        })");

        mpParam->RecursivelyValidateAndAssignDefaults(default_params);
    }

    /// Deleted copy constructor.
    FEASTSolver(const FEASTSolver& Other) = delete;

    /// Destructor.
    virtual ~FEASTSolver() {}

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    FEASTSolver& operator=(const FEASTSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Solve the generalized eigenvalue problem.
    /**
     * K is a symmetric matrix. M is a symmetric positive-definite matrix.
     */
    virtual void Solve(
            SparseMatrixType& K,
            SparseMatrixType& M,
            DenseVectorType& Eigenvalues,
            DenseMatrixType& Eigenvectors)
    {
        const auto SystemSize = K.size1();

        Parameters& FEAST_Settings = *mpParam;
        const double EigenvalueRangeMin = FEAST_Settings["lambda_min"].GetDouble();
        const double EigenvalueRangeMax = FEAST_Settings["lambda_max"].GetDouble();

        int SearchDimension = FEAST_Settings["search_dimension"].GetInt();
        int NumEigenvalues = FEAST_Settings["number_of_eigenvalues"].GetInt();

        Eigenvalues.resize(SearchDimension,false);
        Eigenvectors.resize(SearchDimension,SystemSize,false);

        if (FEAST_Settings["perform_stochastic_estimate"].GetBool())
        {
            // this estimates the number of eigenvalues in the interval [lambda_min, lambda_max]
            Calculate(M,K,EigenvalueRangeMin,EigenvalueRangeMax,SearchDimension,
                    NumEigenvalues,Eigenvalues,Eigenvectors,true);

            std::cout << "Estimated number of eigenvalues = " << NumEigenvalues << std::endl;

            // recommended estimate of search dimension from FEAST documentation
            SearchDimension = NumEigenvalues + NumEigenvalues/2 + 1;
            FEAST_Settings["search_dimension"].SetInt(SearchDimension);
        }
        if (FEAST_Settings["solve_eigenvalue_problem"].GetBool())
        {
            // this attempts to solve the generalized eigenvalue problem
            Calculate(M,K,EigenvalueRangeMin,EigenvalueRangeMax,SearchDimension,
                    NumEigenvalues,Eigenvalues,Eigenvectors,false);

            Eigenvalues.resize(NumEigenvalues,true);
            Eigenvectors.resize(NumEigenvalues,SystemSize,true);

        }
        FEAST_Settings["number_of_eigenvalues"].SetInt(NumEigenvalues);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FEAST solver.";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Parameters::Pointer mpParam;

    ///@}
    ///@name Private Operations
    ///@{

    /// Wrapper for FEAST library.
    void Calculate(
            SparseMatrixType& rMassMatrix,
            SparseMatrixType& rStiffnessMatrix,
            double EigenvalueRangeMin,
            double EigenvalueRangeMax,
            int SearchDimension,
            int& rNumEigenvalues,
            DenseVectorType& rEigenvalues,
            DenseMatrixType& rEigenvectors,
            bool PerformStochasticEstimate)
    {
        KRATOS_TRY

        int FEAST_Params[64] = {};
        int NumIter, Info, SystemSize;
        double Epsout;
        DenseVectorType Residual(SearchDimension);
        std::vector<std::complex<double> > IntegrationNodes, IntegrationWeights;
        SystemSize = static_cast<int>(rMassMatrix.size1());
        matrix<double,column_major> work(SystemSize,SearchDimension);
        matrix<std::complex<double>,column_major> zwork(SystemSize,SearchDimension);
        matrix<double,column_major> Aq(SearchDimension,SearchDimension);
        matrix<double,column_major> Bq(SearchDimension,SearchDimension);
        std::complex<double> Ze;
        FEASTSystemMatrix Az;
        std::vector<std::complex<double> > b(SystemSize);
        std::vector<std::complex<double> > x(SystemSize);

        Az.Initialize(rMassMatrix,rStiffnessMatrix);

        Parameters& FEAST_Settings = *mpParam;

        // warning: if an iterative solver is used, very small tolerances (~1e-15)
        // may be needed for FEAST to work properly.
        SkylineLUSolver<std::complex<double> > Solver;

        // initialize FEAST eigenvalue solver (see FEAST documentation for details)
        feastinit(FEAST_Params);
        if (FEAST_Settings["print_feast_output"].GetBool())
            FEAST_Params[0] = 1;
        FEAST_Params[2] = 8; // stopping convergence criteria 10^-FEAST_Params[2]
        FEAST_Params[28] = 1;// not sure if this is needed
        if (PerformStochasticEstimate)
        {
            FEAST_Params[1] = 4; // number of quadrature points (default: 8)
            FEAST_Params[13] = 2;
        }

        IntegrationNodes.resize(FEAST_Params[1]);
        IntegrationWeights.resize(FEAST_Params[1]);

        // get quadrature nodes and weights
        zfeast_contour(&EigenvalueRangeMin,
                &EigenvalueRangeMax,
                &FEAST_Params[1],
                &FEAST_Params[15],
                &FEAST_Params[17],
                (double *)IntegrationNodes.data(),
                (double *)IntegrationWeights.data());

        int ijob = -1;
        // solve the eigenvalue problem
        while (ijob != 0)
        {
            // FEAST's reverse communication interface
            dfeast_srcix(&ijob,&SystemSize,(double *)&Ze,(double *)work.data().begin(),
                    (double *)zwork.data().begin(),(double *)Aq.data().begin(),
                    (double *)Bq.data().begin(),FEAST_Params,&Epsout,&NumIter,
                    &EigenvalueRangeMin,&EigenvalueRangeMax,&SearchDimension,
                    (double *)rEigenvalues.data().begin(),
                    (double *)rEigenvectors.data().begin(),
                    &rNumEigenvalues,(double *)Residual.data().begin(),&Info,
                    (double *)IntegrationNodes.data(),
                    (double *)IntegrationWeights.data());

            switch (ijob)
            {
                case 10:
                {
                    // set up quadrature matrix (ZeM-K) and solver
                    Az.Calculate(rMassMatrix,rStiffnessMatrix,Ze);
                    Solver.Factorize(Az);
                } break;
                case 11:
                {
                    // solve the linear system for one quadrature point
                    for (int j=0; j < FEAST_Params[22]; j++)
                    {
                        for (int i=0; i < SystemSize; i++)
                            b[i] = zwork(i,j);
                        Solver.Solve(b,x);
                        for (int i=0; i < SystemSize; i++)
                            zwork(i,j) = x[i];
                    }
                } break;
                case 30:
                {
                    // multiply Kx
                    for (int i=0; i < FEAST_Params[24]; i++)
                    {
                        int k = FEAST_Params[23]-1+i;
                        noalias(column(work,k)) = prod(rStiffnessMatrix,row(rEigenvectors,k));
                    }
                } break;
                case 40:
                {
                    // multiply Mx
                    for (int i=0; i < FEAST_Params[24]; i++)
                    {
                        int k = FEAST_Params[23]-1+i;
                        noalias(column(work,k)) = prod(rMassMatrix,row(rEigenvectors,k));
                    }
                }
            } // switch
        } // while

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class FEASTSolver

///@}

///@name Input and output
///@{

/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(std::istream& rIStream,
        FEASTSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    return rIStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(std::ostream& rOStream,
        const FEASTSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}// namespace Kratos.

#endif // KRATOS_FEAST_SOLVER  defined

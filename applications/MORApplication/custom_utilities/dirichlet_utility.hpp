#if !defined(KRATOS_DIRICHLET_UTILITY_H_INCLUDED)
#define KRATOS_DIRICHLET_UTILITY_H_INCLUDED

// System includes

// External includes
// #include <Eigen/QR>

// Project includes
#include "includes/define.h"

namespace Kratos
{
namespace DirichletUtility
{

    /**
     * @brief Finds all fixed dofs
     * @details vector<bool> is not thread safe, so we use unsigned int here
     * @param rDofSet the problem's dof set
     * @param rFixedDofs array where 1 specifies a fixed dof´
     */
    template<typename DofsArrayType>
    void GetDirichletConstraints(DofsArrayType& rDofSet, std::vector<unsigned int>& rFixedDofs)
    {
        const size_t n_dofs = rDofSet.size();

        if( rFixedDofs.size() != n_dofs )
            rFixedDofs.resize( n_dofs, false);

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(n_dofs)
        for( int k = 0; k < static_cast<int>(n_dofs); k++ )
        {
            auto dof_iterator = std::begin(rDofSet) + k;
            rFixedDofs[k] = dof_iterator->IsFixed() ? 1 : 0;
        }
    }

    /**
     * @brief Applies the dirichlet conditions. This operation may be very heavy or completely
     * unexpensive depending on the implementation choosen and on how the System Matrix is built.
     * This should be part of the builder and solver.
     * @details For explanation of how it works for a particular implementation the user
     * should refer to the particular Builder And Solver choosen
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    template<typename SparseSpaceType>
    void ApplyDirichletConditions(
        typename SparseSpaceType::MatrixType& A,
        typename SparseSpaceType::VectorType& b,
        const std::vector<unsigned int>& FixedDofSet,
        double Factor)
    {
        std::size_t system_size = A.size1();
        double* Avalues = A.value_data().begin();
        std::size_t* Arow_indices = A.index1_data().begin();
        std::size_t* Acol_indices = A.index2_data().begin();

        //detect if there is a line of all zeros and set the diagonal to a 1 if this happens
        if( Factor != 0.0 )
        {
            #pragma omp parallel for firstprivate(system_size)
            for (int k = 0; k < static_cast<int>(system_size); ++k){
                std::size_t col_begin = Arow_indices[k];
                std::size_t col_end = Arow_indices[k+1];
                bool empty = true;
                for (std::size_t j = col_begin; j < col_end; ++j)
                {
                    if(Avalues[j] != 0.0)
                    {
                        empty = false;
                        break;
                    }
                }

                if(empty == true)
                {
                    A(k,k) = 1.0;
                    b[k] = 0.0;
                }
            }
        }

        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(system_size); ++k)
        {
            std::size_t col_begin = Arow_indices[k];
            std::size_t col_end = Arow_indices[k+1];
            const unsigned int is_fixed = FixedDofSet[k];
            if( is_fixed == 1 )
            {
                // zero out the whole row, except the diagonal
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if (static_cast<int>(Acol_indices[j]) != k )
                        Avalues[j] = 0.0;
                    else
                        Avalues[j] *= Factor;

                // zero out the RHS
                b[k] = 0.0;
            }
            else
            {
                // zero out the column which is associated with the zero'ed row
                for( std::size_t j = col_begin; j < col_end; ++j )
                    if( FixedDofSet[ Acol_indices[j] ] )
                        Avalues[j] = 0.0;
            }
        }
    }

    // template <typename DenseSpaceType>
    // void OrthogonalizeQR(typename DenseSpaceType::MatrixType& rA)
    // {
    //     typedef typename DenseSpaceType::DataType ScalarType;

    //     Eigen::Map<Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>> A(rA.data().begin(), rA.size1(), rA.size2());
    //     Eigen::HouseholderQR<Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>> qr(A);

    //     Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> thinQ;
    //     Eigen::MatrixXd I(Eigen::MatrixXd::Identity(rA.size1(),rA.size2()));

    //     thinQ = qr.householderQ() * I;
    //     A = thinQ;
    // }

} // namespace DirichletUtility

} // namespace Kratos

#endif
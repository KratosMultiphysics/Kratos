/*          
* =======================================================================*
* kkkk   kkkk  kkkkkkkkkk   kkkkk    kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkk  kkkk   kkkk   kkkk  kkkkkk   kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkkkkkkk    kkkk   kkkk  kkkkkkk     kkkk    kkk    kkk  kkkk         *
* kkkkkkkkk    kkkkkkkkkkk  kkkk kkk	kkkk    kkk    kkk    kkkk       *
* kkkk  kkkk   kkkk  kkkk   kkkk kkkk   kkkk    kkk    kkk      kkkk     *
* kkkk   kkkk  kkkk   kkkk  kkkk  kkkk  kkkk    kkkkkkkkkk  kkkkkkkkkk   *
* kkkk    kkkk kkkk    kkkk kkkk   kkkk kkkk    kkkkkkkkkk  kkkkkkkkkk 	 *
*                                                                        *
* krATos: a fREe opEN sOURce CoDE for mULti-pHysIC aDaptIVe SoLVErS,     *
* aN extEnsIBLe OBjeCt oRiEnTEd SOlutION fOR fInITe ELemEnt fORmULatIONs *
* Copyleft by 2003 ciMNe                                                 *
* Copyleft by 2003 originary authors Copyleft by 2003 your name          *
* This library is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License as         *
* published by the Free Software Foundation; either version 2.1 of       *
* the License, or any later version.                                     *
*                                                                        *
* This library is distributed in the hope that it will be useful, but    *
* WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
* See the GNU Lesser General Public License for more details.            *
*                                                                        *
* You should have received a copy of the GNU Lesser General Public       *
* License along with this library; if not, write to International Centre *
* for Numerical Methods in Engineering (CIMNE),                          *
* Edifici C1 - Campus Nord UPC, Gran Capità s/n, 08034 Barcelona.        *
*                                                                        *
* You can also contact us to the following email address:                *
* kratos@cimne.upc.es                                                    *
* or fax number: +34 93 401 65 17                                        *
*                                                                        *
* Created at Institute for Structural Mechanics                          *
* Ruhr-University Bochum, Germany                                        *
* Last modified by:    $Author: janosch $  				 *
* Date:                $Date: 2008-07-23 14:46:50 $			 *
* Revision:            $Revision: 1.1 $ 				 *
*========================================================================*
* International Center of Numerical Methods in Engineering - CIMNE	 *
* Barcelona - Spain 							 *
*========================================================================*
*/

#if !defined(KRATOS_PARALLEL_SUPERLU_SOLVER_H_INCLUDED )
#define  KRATOS_PARALLEL_SUPERLU_SOLVER_H_INCLUDED

// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

// External includes 
#include <omp.h>
#include <boost/timer.hpp>

#include "boost/smart_ptr.hpp"
// #include "utilities/superlu_interface.h"
#include "includes/ublas_interface.h"
// #include "external_includes/bindings/superlu/superlu.hpp"
// #include "external_includes/bindings/traits/traits.hpp"
// #include <boost/numeric/bindings/traits/ublas_sparse.hpp>
// #include "boost/numeric/bindings/traits/ublas_matrix.hpp"
// #include "structural_application/custom_utilities/ublas_matrix.hpp"
// #include "boost/numeric/bindings/traits/ublas_vector2.hpp"

// #include <boost/numeric/bindings/traits/sparse_traits.hpp>
// #include <boost/numeric/bindings/traits/ublas_matrix.hpp>
// #include <boost/numeric/bindings/traits/ublas_sparse.hpp>
// #include <boost/numeric/bindings/traits/ublas_vector.hpp>

#include <boost/numeric/bindings/traits/sparse_traits.hpp>

#include "external_includes/superlu/superlu.hpp"
// #include "utilities/superlu.hpp"
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"

namespace slu = boost::numeric::bindings::superlu;
namespace ublas = boost::numeric::ublas;

namespace Kratos
{
    template< class TSparseSpaceType, class TDenseSpaceType, 
              class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
        class ParallelSuperLUSolver : public DirectSolver< TSparseSpaceType, 
              TDenseSpaceType, TReordererType>
    {
        public:
            /**
             * Counted pointer of SuperLUSolver
             */
            typedef boost::shared_ptr<ParallelSuperLUSolver> Pointer;
            
            typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 
            
            typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
            
            typedef typename TSparseSpaceType::VectorType VectorType;
            
            typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
            
            /**
             * Default constructor
             */
            ParallelSuperLUSolver(){}
            
            /**
             * Destructor
             */
            virtual ~ParallelSuperLUSolver(){}
            
            /** 
             * Normal solve method.
             * Solves the linear system Ax=b and puts the result on SystemVector& rX.
             * rX is also th initial guess for iterative methods.
             * @param rA. System matrix
             * @param rX. Solution vector.
             * @param rB. Right hand side vector.
             */
            bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
            {
                std::cout << "this is parallel SuperLU solver" << std::endl;
                std::cout << "matrix size in solver: " << rA.size1() << std::endl;
                std::cout << "RHS size in solver: " << rB.size() << std::endl;
                typedef boost::numeric::bindings::traits::sparse_matrix_traits<SparseMatrixType> matraits;
                typedef typename matraits::value_type val_t;
                
                typedef ublas::compressed_matrix<double, ublas::column_major, 0,
                ublas::unbounded_array<std::size_t>, ublas::unbounded_array<double> > cm_t;
                typedef ublas::matrix<double, ublas::column_major> m_t;
                
                boost::timer copy_time;
                
                //TODO: reorder A matrix from row_major format to column_major format
                 std::cout << "number_of_non_zeros: " << rA.nnz() << std::endl;
                 std::cout << "size of index1_data: " << rA.index1_data().size() << std::endl;
                 std::cout << "size of index2_data: " << rA.index2_data().size() << std::endl;
//                 for( int i=0; i<rA.index1_data().size(); i++ )
//                     std::cout << rA.index1_data()[i] << "  " ;
//                 std::cout << std::endl;
                //cm_t A(rA);
		std::cout << "generating matrix" << std::endl;
		cm_t A(rA.size1(),rA.size2(),rA.nnz());
		std::cout << "copy index1_data" << std::endl;
		A.index1_data() = rA.index1_data();
		std::cout << "copy index2_data" << std::endl;
		A.index2_data() = rA.index2_data();
		typedef typename SparseMatrixType::iterator1 i1_t;
		typedef typename SparseMatrixType::iterator2 i2_t;
		std::cout << "copy values" << std::endl;

		A.set_filled( rA.filled1(), rA.filled2() );
//		A.value_data() = rA.value_data();
		for (i1_t i1 = rA.begin1(); i1 != rA.end1(); ++i1) 
		{
			for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
				A(i2.index1(),i2.index2()) = *i2;		
		}
                std::cout << "Copy Time : " << copy_time.elapsed() << std::endl;
                 std::cout << "dest number_of_non_zeros: " << A.nnz() << std::endl;
                 std::cout << "dest size of index1_data: " << A.index1_data().size() << std::endl;
                 std::cout << "dest size of index2_data: " << A.index2_data().size() << std::endl;
//                 for( int i=0; i<rA.index1_data().size(); i++ )
//                     std::cout << A.index1_data()[i] << "  " ;
//                 std::cout << std::endl;
                
                
                
                if(IsNotConsistent(rA, rX, rB))
                {
                    std::cout << "NOT CONSISTENT!" << std::endl;
                    return false;
                }
                
                boost::timer rhs_time;
                
                //manually create RHS matrix
                m_t b( rB.size(), 1 );
                for( unsigned int i=0; i<rB.size(); i++ )
                {
//                     if( ! (abs(rB[i])< 1.0e-15) )
					b(i,0) = rB[i];
//                     else b[i][0] = 0.0;
                }
                std::cout << "RHS Time : " << rhs_time.elapsed() << std::endl;
//                 KRATOS_WATCH(b);
                
//                 for( int i=0; i< rA.index1_data().size(); i++ )
//                 {
//                     std::cout << rA.index1_data()[i]  << " ";
//                     row_index_vector[i] = rA.index1_data()[i];
//                 }
//                 std::cout << std::endl;
//                 for( int i=0; i< rA.index2_data().size(); i++ )
//                 {
//                     std::cout << rA.index2_data()[i]  << " ";
// //                     row_index_vector[i] = rA.index1_data()[i];
//                 }
//                 std::cout << std::endl;
                
                /**
                 * specify number of threads
                 */
                int number_of_threads = omp_get_max_threads();
                KRATOS_WATCH(number_of_threads);
                //call solver routine
                
                boost::timer solve_time;
                
//                 for( int i=0; i<4; i++ )
//                 {
//                     std::cout << "running solver for iteration no. " << i << std::endl;
                    slu::pgssv (number_of_threads, A, b, slu::atpla_min_degree);
//                 }
                
                std::cout << "solve Time : " << solve_time.elapsed() << std::endl;
                
                boost::timer result_resubstitution_time;
                //resubstitution of results
                for( unsigned int i=0; i<rB.size(); i++ )
                {
                    rX[i] = b(i,0);
                }
                std::cout << "Vector Resubstitution Time : " << result_resubstitution_time.elapsed() << std::endl;
//                 boost::timer matrix_resubstitution_time1;
//                 //resubstituting matrix
//                 std::cout << "copy index1_data" << std::endl;
// 		rA.index1_data() = A.index1_data();
//                 std::cout << "Matrix Resubstitution Time 1: " << matrix_resubstitution_time1.elapsed() << std::endl;
//                 boost::timer matrix_resubstitution_time2;
// 		std::cout << "copy index2_data" << std::endl;
// 		rA.index2_data() = A.index2_data();
//                 std::cout << "Matrix Resubstitution Time 2: " << matrix_resubstitution_time2.elapsed() << std::endl;
//                 boost::timer matrix_resubstitution_time3;
// 		typedef typename cm_t::iterator1 ii1_t;
// 		typedef typename cm_t::iterator2 ii2_t;
// 		std::cout << "copy values" << std::endl;
// 
// 		rA.set_filled( A.filled1(), A.filled2() );
// //		A.value_data() = rA.value_data();
// 		for (ii1_t i1 = A.begin1(); i1 != A.end1(); ++i1) 
// 		{
// 			for (ii2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
// 				rA(i2.index1(),i2.index2()) = *i2;		
// 		}
//                 
//                 std::cout << "Matrix Resubstitution Time 3: " << matrix_resubstitution_time3.elapsed() << std::endl;
                return true;
            }
            
            /** 
             * Multi solve method for solving a set of linear systems with same coefficient matrix.
             * Solves the linear system Ax=b and puts the result on SystemVector& rX. 
             * rX is also th initial guess for iterative methods.
             * @param rA. System matrix
             * @param rX. Solution vector.
             * @param rB. Right hand side vector.
             */
            bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
            {
            /**
             * TODO: 
                 * translate SparseMatrixType into SuperMatrix
                 * call solving routine from SuperLU
             */
//                 slu::gssv ( rA, rB, slu::atpla_min_degree);
                
//                 std::cout<<"Matrix Test:"<<std::endl;
//                 std::cout<<"boost matrix:"<<std::endl;
//                 KRATOS_WATCH( rA );
            //             const int size1 = TDenseSpaceType::Size1(rX);
//             const int size2 = TDenseSpaceType::Size2(rX);

                bool is_solved = true;

//             VectorType x(size1);
//             VectorType b(size1);

            // define an object to store skyline matrix and factorization
//             LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType> myFactorization;
            // copy myMatrix into skyline format
//             myFactorization.copyFromCSRMatrix(rA);
            // factorize it
//             myFactorization.factorize();

//             for(int i = 0 ; i < size2 ; i++)
//             {
//                 TDenseSpaceType::GetColumn(i,rX, x);
//                 TDenseSpaceType::GetColumn(i,rB, b);

                // and back solve
//                 myFactorization.backForwardSolve(size1, b, x);

//                 TDenseSpaceType::SetColumn(i,rX, x);
//                 TDenseSpaceType::SetColumn(i,rB, b);
//             }

                return is_solved;
            }
            
            /**
             * Print information about this object.
             */
            void  PrintInfo(std::ostream& rOStream) const
            {
                rOStream << "Parallel SuperLU solver finished.";
            }
            
            /**
             * Print object's data.
             */
            void  PrintData(std::ostream& rOStream) const 
            {
            }
        
        private:
            
            /**
             * Assignment operator.
             */
            ParallelSuperLUSolver& operator=(const ParallelSuperLUSolver& Other);
            
            /**
             * Copy constructor.
             */
//             ParallelSuperLUSolver(const ParallelSuperLUSolver& Other);
    
    }; // Class ParallelSuperLUSolver

    
    /**
     * input stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
    inline std::istream& operator >> (std::istream& rIStream, ParallelSuperLUSolver< TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
    }
    
    /**
     * output stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    inline std::ostream& operator << (std::ostream& rOStream, 
                                      const ParallelSuperLUSolver<TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
 
  
}  // namespace Kratos.

#endif // KRATOS_PARALLEL_SUPERLU_SOLVER_H_INCLUDED  defined 



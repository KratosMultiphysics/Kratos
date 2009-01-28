/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
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
* Last modified by:    $Author: rrossi $  				 *
* Date:                $Date: 2007-03-06 10:30:33 $			 *
* Revision:            $Revision: 1.2 $ 				 *
*========================================================================*
* International Center of Numerical Methods in Engineering - CIMNE	 *
* Barcelona - Spain 							 *
*========================================================================*
*/

#if !defined(KRATOS_SUPERLU_SOLVER_H_INCLUDED )
#define  KRATOS_SUPERLU_SOLVER_H_INCLUDED

// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

// External includes 
#include "boost/smart_ptr.hpp"
// #include "utilities/superlu_interface.h"
#include "includes/ublas_interface.h"
//#include "boost/numeric/bindings/superlu/superlu.hpp"
// #include "boost/numeric/bindings/traits/traits.hpp"
// #include <boost/numeric/bindings/traits/ublas_sparse.hpp>
// #include "boost/numeric/bindings/traits/ublas_matrix.hpp"
// #include "structural_application/custom_utilities/ublas_matrix.hpp"
// #include "boost/numeric/bindings/traits/ublas_vector2.hpp"

// #include <boost/numeric/bindings/traits/sparse_traits.hpp>
// #include <boost/numeric/bindings/traits/ublas_matrix.hpp>
// #include <boost/numeric/bindings/traits/ublas_sparse.hpp>
// #include <boost/numeric/bindings/traits/ublas_vector.hpp>

#include <numeric/bindings/traits/sparse_traits.hpp>

#include "utilities/superlu.hpp"
// #include <boost/numeric/bindings/superlu/superlu.hpp>
#include <numeric/bindings/traits/ublas_matrix.hpp>
#include <numeric/bindings/traits/ublas_sparse.hpp>
#include <numeric/bindings/traits/ublas_vector.hpp>

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"
namespace slu = boost::numeric::bindings::superlu;
namespace ublas = boost::numeric::ublas;

namespace Kratos
{
    template< class TSparseSpaceType, class TDenseSpaceType, 
              class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
        class SuperLUSolver : public DirectSolver< TSparseSpaceType, 
              TDenseSpaceType, TReordererType>
    {
        public:
            /**
             * Counted pointer of SuperLUSolver
             */
            typedef boost::shared_ptr<SuperLUSolver> Pointer;
            
            typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 
            
            typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
            
            typedef typename TSparseSpaceType::VectorType VectorType;
            
            typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
            
            /**
             * Default constructor
             */
            SuperLUSolver(){}
            
            /**
             * Destructor
             */
            virtual ~SuperLUSolver(){}
            
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
				std::cout << "matrix size in solver: " << rA.size1() << std::endl;
				std::cout << "RHS size in solver: " << rB.size() << std::endl;
                typedef boost::numeric::bindings::traits::sparse_matrix_traits<SparseMatrixType> matraits;
                typedef typename matraits::value_type val_t;
                
                typedef ublas::compressed_matrix<double, ublas::row_major, 0,
                ublas::unbounded_array<int>, ublas::unbounded_array<double> > cm_t;
                typedef ublas::matrix<double, ublas::row_major> m_t;
                
                if(IsNotConsistent(rA, rX, rB))
                    return false;
                
                //manually create RHS matrix
                m_t b( rB.size(), 1 );
                for( int i=0; i<rB.size(); i++ )
                {
//                     if( ! (abs(rB[i])< 1.0e-15) )
					b(i,0) = rB[i];
//                     else b[i][0] = 0.0;
                }
                
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
                
                //call solver routine
                slu::gssv (rA, b, slu::atpla_min_degree);
                
                //resubstitution of results
                for( int i=0; i<rB.size(); i++ )
                {
					rX[i] = b(i,0);
                }
                
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
                rOStream << "SuperLU solver finished.";
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
            SuperLUSolver& operator=(const SuperLUSolver& Other);
            
            /**
             * Copy constructor.
             */
            SuperLUSolver(const SuperLUSolver& Other);
    
    }; // Class SkylineLUFactorizationSolver 

    
    /**
     * input stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
    inline std::istream& operator >> (std::istream& rIStream, SuperLUSolver< TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
    }
    
    /**
     * output stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    inline std::ostream& operator << (std::ostream& rOStream, 
                                      const SuperLUSolver<TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
 
  
}  // namespace Kratos.

#endif // KRATOS_SUPERLU_SOLVER_H_INCLUDED  defined 



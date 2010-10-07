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

#if !defined(KRATOS_GMRES_SOLVER_H_INCLUDED )
#define  KRATOS_GMRES_SOLVER_H_INCLUDED
#include "boost/smart_ptr.hpp"
// #include "utilities/superlu_interface.h"
#include "includes/ublas_interface.h"
// #include "boost/numeric/bindings/superlu/superlu.hpp"
// #include "boost/numeric/bindings/traits/traits.hpp"
// #include <boost/numeric/bindings/traits/ublas_sparse.hpp>
// #include "boost/numeric/bindings/traits/ublas_matrix.hpp"
// #include "structural_application/custom_utilities/ublas_matrix.hpp"
// #include "boost/numeric/bindings/traits/ublas_vector2.hpp"
// #include <numeric/bindings/traits/sparse_traits.hpp>
// // #include <boost/numeric/bindings/superlu/superlu.hpp>
// #include <numeric/bindings/traits/ublas_matrix.hpp>
// #include <numeric/bindings/traits/ublas_sparse.hpp>
// #include <numeric/bindings/traits/ublas_vector.hpp>
// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"

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
  
  /// Short class definition.
  /** Detail class definition.
  */
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>, 
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
    class GMRESSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of GMRESSolver
      typedef boost::shared_ptr<GMRESSolver> Pointer;

      typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      GMRESSolver(){}

      GMRESSolver(double NewTolerance) : BaseType(NewTolerance){}

      GMRESSolver(double NewTolerance, unsigned int NewMaxIterationsNumber) : BaseType(NewTolerance, NewMaxIterationsNumber){}

      GMRESSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) : 
	BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){
		}

      /// Copy constructor.
   	GMRESSolver(const GMRESSolver& Other) : BaseType(Other){}

      /// Destructor.
      virtual ~GMRESSolver(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      GMRESSolver& operator=(const GMRESSolver& Other)
      {
        BaseType::operator=(Other);
	return *this;
      }
      
      ///@}
      ///@name Operations
      ///@{
      
      /** Normal solve method.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  rX is also th initial guess for iterative methods.
	  @param rA. System matrix
	  @param rX. Solution vector. it's also the initial 
	  guess for iterative linear solvers.
 	  @param rB. Right hand side vector.
      */
      bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
	{
	  	if(this->IsNotConsistent(rA, rX, rB))
	    		return false;

		BaseType::mBNorm = TSparseSpaceType::TwoNorm(rB);

	 	bool is_solved= false;
// 
	  //GetTimeTable()->Start(Info());

// 	  BaseType::GetPreconditioner()->Initialize(rA,rX,rB);
//  	  BaseType::GetPreconditioner()->ApplyInverseRight(rX);
// 	  BaseType::GetPreconditioner()->ApplyLeft(rB);

  	  	long double *x;
  	  	long double *rhs;
	  	unsigned int nz_sum=0;

		x= new long double[rB.size()];
		rhs= new long double[rB.size()];

		for(unsigned int i=0; i<rB.size(); i++)
		{
			x[i]= rX(i);
			rhs[i]= rB(i);
		}
		for(unsigned int i=0; i<rA.size1(); i++)
			for(unsigned int j=0; j<rA.size2(); j++)
				if(rA(i,j) != 0.0)
					nz_sum++;
// 
  	  	long double *a;
  	  	int *ia;
  	  	int *ja;
		
		a= new long double[nz_sum];
		ia= new int[nz_sum];
		ja= new int[nz_sum];

	 	unsigned int counter= 0;

		for(unsigned int i=0; i<rA.size1(); i++)
			for(unsigned int j=0; j<rA.size2(); j++)
				if(rA(i,j) != 0.0)
				{
					a[counter]=rA(i,j);
 					ia[counter]= i;
 					ja[counter]= j;
					counter++;
				}
		unsigned int inner_iterations=100;
		if(inner_iterations> rB.size()-1)
			inner_iterations= rB.size()-1;

		mNumberOfRestarts= 0;
// 		for(unsigned int i=0; i<2; i++)
// 		{
// 	 		if(
			is_solved = mgmres (a, ia, ja, x, rhs, 
  				rB.size(), nz_sum, BaseType::GetMaxIterationsNumber(), inner_iterations, 
				BaseType::GetTolerance() );/*)*/
// 			{
// 				break;
// 			}
// 			mNumberOfRestarts++;
// 		}

	  	for(unsigned int i=0; i<rB.size(); i++)
		{
			rX(i)=x[i];
		}

//  	  BaseType::GetPreconditioner()->Finalize(rX);

// 	  //GetTimeTable()->Stop(Info());
// 
		delete[] x;
 		delete[] rhs;
		delete[] a;
		delete[] ia;
		delete[] ja;

		if(!is_solved)
		{
			SkylineLUFactorizationSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>().Solve(rA, rX, rB);
			std::cout<<"####GMRES not converged -> System solved with SkylineLUFactorizationSolver####"<<std::endl;
		}

	  	return is_solved;
	}
      
      /** Multi solve method for solving a set of linear systems with same coefficient matrix.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  rX is also th initial guess for iterative methods.
	  @param rA. System matrix
	  @param rX. Solution vector. it's also the initial 
	  guess for iterative linear solvers.
 	  @param rB. Right hand side vector.
      */
      bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
	{
	  //DOES NOTHING AT THE MOMENT
          bool is_solved = true;

	  //GetTimeTable()->Stop(Info());

	  return is_solved;
	}

      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Return information about this object.
      virtual std::string Info() const
	{
	  std::stringstream buffer;
	  buffer << "GMRES iterative solver [at the moment unfortunately without] with " << BaseType::GetPreconditioner()->Info();
	  return  buffer.str();
	}
      
      /// Print information about this object.
      void  PrintInfo(std::ostream& OStream) const
	{
	  OStream << "GMRES iterative solver [at the moment unfortunately without] with ";
	  BaseType::GetPreconditioner()->PrintInfo(OStream);
	}

      /// Print object's data.
      void  PrintData(std::ostream& OStream) const 
	{
	  OStream << "GMRES "<<mNumberOfRestarts<<" times restarted"<<std::endl;
	  BaseType::PrintData(OStream);
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
        
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
// //******************************************************************************
// //
// //  Purpose:
// //
// //    AX computes A * X for a sparse matrix.
// //
// //  Discussion:
// //
// //    The matrix A is assumed to be sparse.  To save on storage, only
// //    the nonzero entries of A are stored.  For instance, the K-th nonzero
// //    entry in the matrix is stored by:
// //
// //      A(K) = value of entry,
// //      IA(K) = row of entry,
// //      JA(K) = column of entry.
// //
// //  Modified:
// //
// //    08 August 2006
// //
// //  Author:
// //
// //    Lili Ju
// //
// //    C++ translation by John Burkardt
// //
// //  Reference:
// //
// //    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
// //      June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
// //      Charles Romine, Henk van der Vorst,
// //    Templates for the Solution of Linear Systems:
// //      Building Blocks for Iterative Methods,
// //    SIAM, 1994.
// //
// //    Yousef Saad,
// //    Iterative Methods for Sparse Linear Systems,
// //    Second Edition,
// //    SIAM, 20003,
// //    ISBN: 0898715342.
// //
// //  Parameters:
// //
// //    Input, double A[NZ_NUM], the matrix values.
// //
// //    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
// //    of the matrix values.
// //
// //    Input, double X[N], the vector to be multiplied by A.
// //
// //    Output, double W[N], the value of A*X.
// //
// //    Input, int N, the order of the system.
// //
// //    Input, int NZ_NUM, the number of nonzeros.
// //
      void ax ( long double *a, int *ia, int *ja, long double *x, long double *w, int n, int nz_num )
	{
  		int i;
  		int j;
  		int k;

  		for ( i = 0; i < n; i++ )
  		{
    			w[i] = 0.0;
  		}

  		for ( k = 0; k < nz_num; k++ )
  		{
    			i = ia[k];
    			j = ja[k];
    			w[i] = w[i] + a[k] * x[j];
  		}
  		return;
	}
// //******************************************************************************
// //******************************************************************************
// //
// //  Purpose:
// //
// //    DOT_PRODUCT computes the dot product of two vectors.
// //
// //  Modified:
// //
// //    08 August 2006
// //
// //  Author:
// //
// //    Lili Ju
// //
// //    C++ translation by John Burkardt
// //
// //  Parameters:
// //
// //    Input, double V0[N], V1[N], two vectors whose product is desired.
// //
// //    Input, int N, the dimension of the vectors.
// //
// //    Output, double dot_product, the dot product of the two vectors.
// //
	long double dot_product ( long double v0[], long double v1[], int n )
	{
  		int i;
  		long double s;

  		s = 0.0;
  		for ( i = 0; i < n; i++ )
  		{
    			s = s + v0[i] * v1[i];
  		}

  		return s;
	}
// //******************************************************************************
// //******************************************************************************
// //
// //  Purpose:
// //
// //    MGMRES applies the restarted GMRES iteration to a linear system.
// //
// //  Discussion:
// //
// //    The linear system A*X=B is solved iteratively.
// //
// //    The matrix A is assumed to be sparse.  To save on storage, only
// //    the nonzero entries of A are stored.  For instance, the K-th nonzero
// //    entry in the matrix is stored by:
// //
// //      A(K) = value of entry,
// //      IA(K) = row of entry,
// //      JA(K) = column of entry.
// //
// //    The "matrices" H and V are treated as one-dimensional vectors
// //    which store the matrix data in row major form.
// //
// //    This requires that references to H[I][J] be replaced by references
// //    to H[I*MR+J] and references to V[I][J] by V[I*N+J].
// //
// //  Modified:
// //
// //    30 October 2006
// //
// //  Author:
// //
// //    Lili Ju
// //
// //    C++ translation by John Burkardt
// //
// //  Reference:
// //
// //    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
// //      June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
// //      Charles Romine, Henk van der Vorst,
// //    Templates for the Solution of Linear Systems:
// //      Building Blocks for Iterative Methods,
// //    SIAM, 1994.
// //
// //    Yousef Saad,
// //    Iterative Methods for Sparse Linear Systems,
// //    Second Edition,
// //    SIAM, 20003,
// //    ISBN: 0898715342.
// //
// //  Parameters:
// //
// //    Input, double A[NZ_NUM], the matrix values.
// //
// //    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
// //    of the matrix values.
// //
// //    Input/output, double X[N]; on input, an approximation to
// //    the solution.  On output, an improved approximation.
// //
// //    Input, double RHS[N], the right hand side of the linear system.
// //
// //    Input, int N, the order of the linear system.
// //
// //    Input, int NZ_NUM, the number of nonzero matrix values.
// //
// //    Input, int ITRMAX, the maximum number of (outer) iterations to take.
// //
// //    Input, int MR, the maximum number of (inner) iterations to take.
// //    MR must be less than N.
// //
// //    Input, double TOL, a tolerance.
// //
	bool mgmres ( long double a[], int ia[], int ja[],long double x[],long double rhs[], 
  		int n, int nz_num, int itrmax, int mr, long double tol )
	{
  		long double av;
  		long double *c;
  		long double delta = 1.0e-3;
  		long double *g;
  		long double *h;
  		long double htmp;
		int itr;
  		int i;
   		int j;
  		int k;
  		long double mu;
  		long double *r;
  		long double rho = 0.0;
  		long double rho_tol = 0.0;
  		long double *s;
  		long double *v;
  		long double *y;
  		bool isConverged= false;

  		c = new long double[mr+1];
  		g = new long double[mr+1];
  		h = new long double[(mr+1)*mr];
  		r = new long double[n];
  		s = new long double[mr+1];
  		v = new long double[(mr+1)*n];
  		y = new long double[mr+1];

  		for (itr = 0; itr < itrmax; itr++ )
  		{
    			ax ( a, ia, ja, x, r, n, nz_num );

    			for ( i = 0; i < n; i++ )
    			{
      				r[i] = rhs[i] - r[i];
    			}

    			rho = sqrt ( dot_product ( r, r, n ) );

    			if ( itr == 0 ) 
    			{
      				rho_tol = rho * tol;
    			}

    			for (i = 0; i < n; i++)
    			{
      				v[0*n+i] = r[i] / rho;
    			}

    			for ( i = 0; i <= mr; i++ )
    			{
      				g[i] = 0.0;
      				for ( j = 0; j < mr; j++ ) 
      				{
        				h[i*mr+j] = 0.0;
      				}
    			}

    			g[0] = rho;
    			k = 0;

    			while ( ( rho_tol < rho ) && ( k < mr ) ) 
    			{
      				ax ( a, ia, ja, v+k*n, v+(k+1)*n, n, nz_num );
      				av = sqrt ( dot_product ( v+(k+1)*n, v+(k+1)*n, n ) );

      				for ( j = 0; j <= k; j++ )
      				{
        				h[j*mr+k] = dot_product ( v+(k+1)*n, v+j*n, n );
        				for ( i = 0; i < n; i++ ) 
        				{
          					v[(k+1)*n+i] = v[(k+1)*n+i] - h[j*mr+k] * v[j*n+i];
        				}
      				}

      				h[(k+1)*mr+k] = sqrt ( dot_product ( v+(k+1)*n, v+(k+1)*n, n ) );

      				if ( ( av + delta * h[(k+1)*mr+k] ) == av )
      				{
         				for ( j = 0; j <= k; j++ )
         				{
           					htmp = dot_product ( v+(k+1)*n, v+j*n, n );
           					h[j*mr+k] = h[j*mr+k] + htmp;
           					for ( i = 0; i < n; i++ )
           					{
             						v[(k+1)*n+i] = v[(k+1)*n+i] - htmp * v[j*n+i];
          	 				}
         				}
         				h[(k+1)*mr+k] = sqrt ( dot_product ( v+(k+1)*n, v+(k+1)*n, n ) );
      				}

      				for ( i = 0; i < n; i++ ) 
      				{
        				v[(k+1)*n+i] = v[(k+1)*n+i] / h[(k+1)*mr+k];
      				}

      				if ( 0 < k )
      				{
        				for ( i = 0; i <= k+1; i++ )
        				{
          					y[i] = h[i*mr+k];
        				}
        				for ( j = 0; j < k; j++ ) 
        				{
          					mult_givens ( c[j], s[j], j, y );
        				}
        				for ( i = 0; i <= k+1; i++ ) 
        				{
          					h[i*mr+k] = y[i];
        				}
      				}
      				mu = sqrt ( pow ( h[k*mr+k], 2 ) + pow ( h[(k+1)*mr+k], 2 ) );
      				c[k] = h[k*mr+k] / mu;
      				s[k] = -h[(k+1)*mr+k] / mu;
      				h[k*mr+k] = c[k] * h[k*mr+k] - s[k] * h[(k+1)*mr+k];
      				h[(k+1)*mr+k] = 0;
      				mult_givens ( c[k], s[k], k, g );
      				rho = fabs ( g[k] );
      				k = k + 1;
    			}
    			k = k-1;
    			y[k] = g[k] / h[k*mr+k];

    			for ( i = k-1; 0 <= i; i-- )
    			{
      				y[i] = g[i];
      				for ( j = k; i < j; j-- ) 
      				{
        				y[i] = y[i] - h[i*mr+j] * y[j];
      				}
      				y[i] = y[i] / h[i*mr+i];
    			}

    			for ( i = 0; i < n; i++ )
    			{
      				for ( j = 0; j < k; j++ )
      				{
        				x[i] = x[i] + v[j*n+i] * y[j];
      				}
   			}

    			if ( rho <= rho_tol ) 
    			{
				isConverged= true;
      				break;
    			}
  		}

		BaseType::SetIterationsNumber(( itr + 1 ) * mr);
		BaseType::SetResidualNorm(rho);

  		delete [] c;
  		delete [] g;
 		delete [] h;
  		delete [] r;
  		delete [] s;
  		delete [] v;
  		delete [] y;

  		return isConverged;
	}
// //******************************************************************************
// //******************************************************************************
// //
// //  Purpose:
// //
// //    MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
// //
// //  Modified:
// //
// //    08 August 2006
// //
// //  Author:
// //
// //    Lili Ju
// //
// //    C++ translation by John Burkardt
// //
// //  Reference:
// //
// //    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
// //      June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
// //      Charles Romine, Henk van der Vorst,
// //    Templates for the Solution of Linear Systems:
// //      Building Blocks for Iterative Methods,
// //    SIAM, 1994.
// //
// //    Yousef Saad,
// //    Iterative Methods for Sparse Linear Systems,
// //    Second Edition,
// //    SIAM, 20003,
// //    ISBN: 0898715342.
// //
// //  Parameters:
// //
// //    Input, double C, S, the cosine and sine of a Givens
// //    rotation.
// //
// //    Input, int K, indicates the location of the first vector entry.
// //
// //    Input/output, double G[K+2], the vector to be modified.  On output,
// //    the Givens rotation has been applied to entries G(K) and G(K+1).
// //
	void mult_givens ( long double c,long double s, int k, long double g[] )
	{
  		long double g1;
  		long double g2;

 	 	g1 = c * g[k] - s * g[k+1];
  		g2 = s * g[k] + c * g[k+1];

  		g[k]   = g1;
  		g[k+1] = g2;

  		return;
	}
// //******************************************************************************
// //******************************************************************************
// //
// //  Purpose:
// //
// //    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
// //
// //  Discussion:
// //
// //    This routine implements the recursion
// //
// //      seed = 16807 * seed mod ( 2**31 - 1 )
// //      unif = seed / ( 2**31 - 1 )
// //
// //    The integer arithmetic never requires more than 32 bits,
// //    including a sign bit.
// //
// //  Modified:
// //
// //    19 August 2004
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Reference:
// //
// //    Paul Bratley, Bennett Fox, Linus Schrage,
// //    A Guide to Simulation,
// //    Springer Verlag, pages 201-202, 1983.
// //
// //    Bennett Fox,
// //    Algorithm 647:
// //    Implementation and Relative Efficiency of Quasirandom
// //    Sequence Generators,
// //    ACM Transactions on Mathematical Software,
// //    Volume 12, Number 4, pages 362-376, 1986.
// //
// //    Peter Lewis, Allen Goodman, James Miller,
// //    A Pseudo-Random Number Generator for the System/360,
// //    IBM Systems Journal,
// //    Volume 8, pages 136-143, 1969.
// //
// //  Parameters:
// //
// //    Input, int N, the number of entries in the vector.
// //
// //    Input/output, int *SEED, a seed for the random number generator.
// //
// //    Output, double R8VEC_UNIFORM_01[N], the vector of pseudorandom values.
// //
	long double *r8vec_uniform_01 ( int n, int *seed )
	{
  		int i;
  		int k;
  		long double *r;

  		if ( *seed == 0 )
  		{
    			std::cout<< std::endl;
    			std::cout << "R8VEC_UNIFORM_01 - Fatal error!"<< std::endl;
    			std::cout << "  Input value of SEED = 0."<< std::endl;
  		}

  		r = new long double[n];

  		for ( i = 0; i < n; i++ )
  		{
    			k = *seed / 127773;

    			*seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    			if ( *seed < 0 )
    			{
      				*seed = *seed + 2147483647;
    			}

    			r[i] = ( long double ) ( *seed ) * 4.656612875E-10;
  		}
  		return r;
	}

        unsigned int mNumberOfRestarts;

    }; // Class GMRESSolver 
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
  inline std::istream& operator >> (std::istream& rIStream, 
				      GMRESSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
        return rIStream;
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const GMRESSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
    }
//   ///@} 
//   
//   
}  // namespace Kratos.

#endif // KRATOS_GMRES_SOLVER_H_INCLUDED  defined 



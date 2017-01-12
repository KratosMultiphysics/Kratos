/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2015-14-08 12:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

/*
Krylov accelerator utility, using LAPACK subroutine DGELS
to solve the least square problem required by the algorithm.
Reference:
Scott, Michael H., and Gregory L. Fenves. 
"A Krylov subspace accelerated Newton algorithm." 
Proc., 2003 ASCE Structures Congress. 2003.
*/

#if !defined(KRYLOV_SUBSPACE_ACCELERATOR_LAPACK_UTILITIES_H_INCLUDED)
#define KRYLOV_SUBSPACE_ACCELERATOR_LAPACK_UTILITIES_H_INCLUDED

#ifdef MULTISCALE_APPLICATION_LAPACK_HAS_DGELS


extern "C" int  MULTISCALE_APPLICATION_LAPACK_DGELS(
					char *T, int *M, int *N, int *NRHS,
					double *A, int *LDA, double *B, int *LDB,
					double *WORK, int *LWORK, int *INFO);
#define CALL_DGELS MULTISCALE_APPLICATION_LAPACK_DGELS

namespace Kratos
{

namespace KrylovSubspaceAcceleratorUtilties
{

	class KrylovAcceleratorLapack
	{
	
	private:

		// Storage for update vectors
		std::vector<Vector> mv;
		// Storage for subspace vectors
		std::vector<Vector> mAv;
		// Array data sent to LAPACK subroutine
		double * mAvData;
		double * mrData;
		double * mwork;
		// Length of work array
		int mlwork;
		// Size information
		int mnumEqns;
		int mmaxDimension;
		// current subspace size
		int m_krylov_dim;

	public:
		
		KrylovAcceleratorLapack(int max_krylov_dimension = 10)
		{
			this->InitializeData();
			mmaxDimension = std::max(1,max_krylov_dimension);
		}
		
		~KrylovAcceleratorLapack()
		{
			this->DestroyData();
		}
		
	public:

		// to be called before using it
		void InitializeData()
		{
			// Storage for update vectors
			mv.clear();
			// Storage for subspace vectors
			mAv.clear();
			// Array data sent to LAPACK subroutine
			mAvData = NULL;
			mrData  = NULL;
			mwork   = NULL;
			// Length of work array
			mlwork = 0;
			// Size information
			mnumEqns = 0;
			mmaxDimension = 10;
		}

		// to be called to clean data
		void DestroyData()
		{
			// Storage for update vectors
			mv.clear();
			// Storage for subspace vectors
			mAv.clear();
			// Array data sent to LAPACK subroutine
			if (mAvData != NULL)
				delete [] mAvData;
			mAvData = NULL;
			if (mrData != NULL)
				delete [] mrData;
			mrData = NULL;
			if (mwork != NULL)
				delete [] mwork;
			mwork = NULL;
			// Length of work array
			mlwork = 0;
			// Size information
			mnumEqns = 0;
			mmaxDimension = 10;
		}

		void BeginSolutionStep(int new_numEqns)
		{
			if(new_numEqns != mnumEqns)
				this->DestroyData();
			mnumEqns = new_numEqns;
			if (mmaxDimension > mnumEqns) mmaxDimension = mnumEqns;
			if(mv.size() == 0) {
				mv.resize(mmaxDimension+1);
				for (int i = 0; i < mmaxDimension+1; i++)
					mv[i] = Vector(mnumEqns);
			}
			if(mAv.size() == 0) {
				mAv.resize(mmaxDimension+1);
				for (int i = 0; i < mmaxDimension+1; i++)
					mAv[i] = Vector(mnumEqns);
			}
			if (mAvData == NULL)
				mAvData = new double [mmaxDimension*mnumEqns];
			if (mrData == NULL)
				mrData = new double [(mnumEqns > mmaxDimension) ? mnumEqns : mmaxDimension];
			mlwork = 2 * ((mnumEqns < mmaxDimension) ? mnumEqns : mmaxDimension);
			if (mwork == NULL)
				mwork = new double [mlwork];
			
			// current dimension of Krylov subspace
			m_krylov_dim = 0;
		}
		
		bool BuildSystemMatrix()
		{
			if (m_krylov_dim > mmaxDimension) 
			{
				m_krylov_dim = 0;
				return true;
			}
			else 
			{
				return false;
			}
		}
		
		int LeastSquares(const Vector& r)
		{
			int k = m_krylov_dim;
			
			// v_{k+1} = w_{k+1} + q_{k+1}
			noalias(mv[k]) = r;
			noalias(mAv[k]) = r;

			// Subspace is empty
			if (k == 0)
				return 0;

			// Compute Av_k = f(y_{k-1}) - f(y_k) = r_{k-1} - r_k
			noalias(mAv[k-1]) -= r;

			int i,j;

			// Put subspace vectors into AvData
			for (i = 0; i < k; i++) {
				Vector &Ai = mAv[i];
				for (j = 0; j < mnumEqns; j++)
					mAvData[i*mnumEqns + j] = Ai(j);
			}

			// Put residual vector into rData (need to save r for later!)
			for(i = 0; i < mnumEqns; i++)
				mrData[i] = r(i);
  
			// No transpose
			char *trans = "N";
		
			// The number of right hand side vectors
			int nrhs = 1;
		
			// Leading dimension of the right hand side vector
			int ldb = (mnumEqns > k) ? mnumEqns : k;
		
			// Subroutine error flag
			int info = 0;

			// Call the LAPACK least squares subroutine
			CALL_DGELS(trans, &mnumEqns, &k, &nrhs, mAvData, &mnumEqns, mrData, &ldb, mwork, &mlwork, &info);
			
			// Check for error returned by subroutine
			if (info < 0) {
				std::stringstream ss;
				ss << "WARNING Krylov Newton::leastSquares() - \n";
				ss << "error code " << info << " returned by LAPACK dgels\n";
				std::cout << ss.str();
				return info;
			}
  
			// Compute the correction vector
			double cj;
			for (j = 0; j < k; j++) {
			
				// Solution to least squares is written to rData
				cj = mrData[j];
			
				// Compute w_{k+1} = c_1 v_1 + ... + c_k v_k
				noalias(mv[k]) += cj*mv[j];
			
				// Compute least squares residual q_{k+1} = r_k - (c_1 Av_1 + ... + c_k Av_k)
				noalias(mv[k]) -= cj*mAv[j];
			}
		
			return 0;
		}
		
		void AccelerateSolution(Vector& r)
		{
			noalias(r) = mv[m_krylov_dim]; // <-- use the accelerated krylov vector
			// Increase current dimension of Krylov subspace
			m_krylov_dim++;
		}
		
	};

}

} // namespace Kratos



#else

namespace Kratos
{

namespace KrylovSubspaceAcceleratorUtilties
{

	class KrylovAcceleratorLapack
	{

	public:
		
		KrylovAcceleratorLapack(int max_krylov_dimension = 10)
		{
			this->InitializeData();
		}
		
		~KrylovAcceleratorLapack()
		{
			this->DestroyData();
		}

		void InitializeData() {}

		void DestroyData() {}

		void BeginSolutionStep(int new_numEqns) {}
		
		bool BuildSystemMatrix(){ return true; }
		
		int LeastSquares(const Vector& r){ return 0; }
		
		void AccelerateSolution(Vector& r){}
		
	};

}

} // namespace Kratos

#endif // MULTISCALE_APPLICATION_LAPACK_HAS_DGELS

#endif // KRYLOV_SUBSPACE_ACCELERATOR_LAPACK_UTILITIES_H_INCLUDED

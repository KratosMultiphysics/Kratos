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

#if !defined(KRYLOV_SUBSPACE_ACCELERATOR_EIGENLIB_UTILITIES_H_INCLUDED)
#define KRYLOV_SUBSPACE_ACCELERATOR_EIGENLIB_UTILITIES_H_INCLUDED

#include <iostream>
#include "Eigen/Dense"

namespace Kratos
{

namespace KrylovSubspaceAcceleratorUtilties
{

	class KrylovAcceleratorEigenlib
	{
	
	private:

		// Storage for update vectors
		std::vector<Vector> mv;
		// Storage for subspace vectors
		std::vector<Vector> mAv;
		// data for eigen
		Eigen::MatrixXd mAvData;
		Eigen::VectorXd mrData;
		Eigen::VectorXd mSolution;
		// Size information
		int mnumEqns;
		int mmaxDimension;
		// current subspace size
		int m_krylov_dim;

	public:
		
		KrylovAcceleratorEigenlib(int max_krylov_dimension = 10)
		{
			this->InitializeData();
			mmaxDimension = std::max(1,max_krylov_dimension);
		}
		
		~KrylovAcceleratorEigenlib()
		{
			this->DestroyData();
		}
		
	public:

		void InitializeData()
		{
			// Storage for update vectors
			mv.clear();
			// Storage for subspace vectors
			mAv.clear();
			// data for eigen
			mAvData.resize(0,0);
			mrData.resize(0);
			mSolution.resize(0);
			// Size information
			mnumEqns = 0;
			mmaxDimension = 10;
		}

		void DestroyData()
		{
			// Storage for update vectors
			mv.clear();
			// Storage for subspace vectors
			mAv.clear();
			// data for eigen
			mAvData.resize(0,0);
			mrData.resize(0);
			mSolution.resize(0);
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
			mAvData.resize(mnumEqns, mmaxDimension);
			mrData.resize(mnumEqns);
			mSolution.resize(mnumEqns);
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
					mAvData(j,i) = Ai(j);
			}

			// Put residual vector into rData (need to save r for later!)
			for(i = 0; i < mnumEqns; i++)
				mrData(i) = r(i);

			// Normal Equations - the fastest but the least accurate, especially with ill-conditioned matrices
			//mSolution = ( mAvData.block(0,0,mnumEqns,k).transpose() * mAvData.block(0,0,mnumEqns,k) ).ldlt().solve(
			//			  mAvData.block(0,0,mnumEqns,k).transpose() * mrData );

			// SVD - most accurate but the slowest
			//mSolution = mAvData.block(0,0,mnumEqns,k).jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(mrData);

			// QR versions, in between the previous two:
			// 1) HouseHolder QR - fastest, but least accurate

			//mSolution = mAvData.block(0,0,mnumEqns,k).householderQr().solve(mrData);
			// 2) QR with column pivoting - a bit slower but more accurate
			mSolution = mAvData.block(0,0,mnumEqns,k).colPivHouseholderQr().solve(mrData);

			// 3) QR with full pivoting - slowest, most accurate
			//mSolution = mAvData.block(0,0,mnumEqns,k).fullPivHouseholderQr().solve(mrData);
			
  
			// Compute the correction vector
			double cj;
			for (j = 0; j < k; j++) {
			
				cj = mSolution(j);
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

#endif // KRYLOV_SUBSPACE_ACCELERATOR_EIGENLIB_UTILITIES_H_INCLUDED

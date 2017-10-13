/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
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
//   Date:                $Date: 2013-06-06 10:37:00 $
//   Revision:            $Revision: 1.00 $
//
//


#if !defined(KRATOS_EIGENLIB_LU_LINEAR_SOLVER_H_INCLUDED )
#define  KRATOS_EIGENLIB_LU_LINEAR_SOLVER_H_INCLUDED

//#define USE_MUMPS_STATIC_LIB

// External includes

// Project includes
#include "utilities/timer.h"
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"
#include <iostream>

#include "Eigen\Dense"
#include "Eigen\Sparse"

namespace Kratos
{

template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class EigenlibLULinearSolver
    : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( EigenlibLULinearSolver );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
    typedef typename TSparseSpaceType::VectorType VectorType;
    typedef typename TDenseSpaceType::MatrixType  DenseMatrixType;

	typedef typename TSparseSpaceType::DataType  DataType;
	typedef typename TSparseSpaceType::IndexType IndexType;
	
	typedef typename SparseMatrixType::const_iterator1 SparseIteratorType1;
	typedef typename SparseMatrixType::const_iterator2 SparseIteratorType2;

	typedef Eigen::Triplet<DataType, int>                                        TripletType;
	typedef std::vector< TripletType >                                           TripletVectorType;
	typedef Eigen::SparseMatrix<DataType, Eigen::ColMajor, int>                  EigenSparseMatrixType;
	typedef Eigen::Matrix< DataType, -1, 1 >                                     EigenVectorType;
	typedef Eigen::Map< EigenVectorType >                                        EigenVectorMapType;
	typedef Eigen::SparseLU< EigenSparseMatrixType, Eigen::COLAMDOrdering<int> > EigenFactorizationType;

public:
	
    EigenlibLULinearSolver()
		: m_pSolver(NULL)
		, m_failed(false)
	{
	}

    virtual ~EigenlibLULinearSolver()
	{
		this->Clear();
	}

private:

public:
	
	void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
    }

    void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
		this->Clear();
		
		int n = rX.size();

		TripletVectorType tr;
		tr.reserve(rA.nnz());
		for( SparseIteratorType1 it1 = rA.begin1(); it1 != rA.end1(); it1++) 
			for( SparseIteratorType2 it2 = it1.begin(); it2 != it1.end(); it2++) 
				tr.push_back(TripletType(it2.index1(),it2.index2(),*it2));

		EigenSparseMatrixType esm(n,n);
		esm.setFromTriplets(tr.begin(), tr.end());
		esm.makeCompressed();

		m_pSolver = new EigenFactorizationType(esm);
		if(m_pSolver->info() != Eigen::Success)
		{
			m_failed = true;
		}
    }

    void PerformSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
		int n = rX.size();
		if(m_failed)
		{
			TSparseSpaceType::SetToZero(rX);
			return;
		}
		EigenVectorMapType map_rX(&(rX.data()[0]), n);
		EigenVectorMapType map_rB(&(rB.data()[0]), n);
		map_rX = m_pSolver->solve(map_rB);
    }

    void FinalizeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
		this->Clear();
    }

    void Clear()
    {
		if(m_pSolver) {
			delete m_pSolver;
			m_pSolver = NULL;
		}
		m_failed = false;
    }

    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        if(IsNotConsistent(rA, rX, rB)) return false;
		int n = rX.size();
		if(n < 1) return true;

		this->InitializeSolutionStep(rA, rX, rB);
		this->PerformSolutionStep(rA, rX, rB);
		this->FinalizeSolutionStep(rA, rX, rB);

		return true;
    }

    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
		if(IsNotConsistent(rA, rX, rB)) return false;
		int n = rX.size1();
		if(n < 1) return true;
		int nrhs = rX.size2();
		if(nrhs < 1) return true;

		VectorType aux_x(n);
		VectorType aux_b(n);

		this->InitializeSolutionStep(rA, aux_x, aux_b);
		for(int i = 0; i < nrhs; i++)
		{
			noalias(aux_b) = column(rB,i);
			noalias(aux_x) = column(rX,i);
			this->PerformSolutionStep(rA, aux_x, aux_b);
			column(rX,i) = aux_x;
		}
		this->FinalizeSolutionStep(rA, aux_x, aux_b);

		return true;
    }

    void  PrintInfo(std::ostream& rOStream) const
    {
		std::cout << "Eigenlib LU Linear\n";
    }

    void  PrintData(std::ostream& rOStream) const
    {
    }

private:

    EigenlibLULinearSolver& operator=(const EigenlibLULinearSolver& Other);
	
private:

	EigenFactorizationType* m_pSolver;
	bool m_failed;

}; // Class EigenlibLULinearSolver


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, EigenlibLULinearSolver<TSparseSpaceType,
                                  TDenseSpaceType,
                                  TReordererType>& rThis)
{
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const EigenlibLULinearSolver<TSparseSpaceType,
                                  TDenseSpaceType,
                                  TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_EIGENLIB_LU_LINEAR_SOLVER_H_INCLUDED  defined 
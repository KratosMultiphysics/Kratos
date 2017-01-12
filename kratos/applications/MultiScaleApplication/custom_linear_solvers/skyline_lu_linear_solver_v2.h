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


#if !defined(KRATOS_SKYLINE_LINEAR_SOLVER_V2_H_INCLUDED )
#define  KRATOS_SKYLINE_LINEAR_SOLVER_V2_H_INCLUDED

//#define USE_MUMPS_STATIC_LIB

// External includes

// Project includes
#include "utilities/timer.h"
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include <iostream>


namespace Kratos
{


template< class TSparseSpaceType, class TDenseSpaceType > class LUSkylineFactorizationV2
{
public:

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef std::size_t IndexType;


    int     size;
    int*    rowIndex;
    int*    perm;
    double* entriesL;
    double* entriesD;
    double* entriesU;

    void clear()
    {
        // if size==0, we have just created it and the arrays have not been allocated
        if (size!=0)
        {
            delete[] rowIndex;
            delete[] entriesL;
            delete[] entriesD;
            delete[] entriesU;
            delete[] perm;
            size= 0;
            rowIndex= perm= NULL;
            entriesL= entriesD= entriesU= NULL;
        }
    }

    //**********************************************************************************
    //**********************************************************************************

    void copyFromCSRMatrix( SparseMatrixType& A)
    {
        int i, j, newi, newj, indexj, ordering, *invperm;
        double entry;

        // First of all, if there was another factorization stored in this object, erase it
        clear();
        // The size of the factorization will be the same of A
        size = A.size1();

        ordering = 1 ; // 0 - no reordering
        // 1 - Cuthill-McKee ordering
        // 2 - reverse Cuthill-McKee ordering

        // Find the permutation for the reordering.
        perm   = new int[size];
        invperm= new int[size];
        if (ordering==0)
        {
            // no reordering
            for (i=0; i<size; i++) perm[i]=i;
        }
        else
        {
            int initialNode, next, currentLevelSet, nMDICLS, node, soughtDegree;
            int maxDegreeInCurrentLevelSet, maxDegree, firstVal, finalVal, increment;
            std::vector<IndexType> neighbors;
            bool empty, found;

            // Cuthill-McKee ordering (or its reverse)
            initialNode= 0; // node where to start search
            maxDegree=0;
            int *degree        = new int[size];
            int *levelSet      = new int[size];
            int *nextSameDegree= new int[size];
            for (i=0; i<size; i++)
            {
                degree[i]      = TSparseSpaceType::GraphDegree(i, A);
                if (degree[i]>maxDegree) maxDegree=degree[i];
                levelSet[i]    = 0;
                nextSameDegree[i]= -1;
            }
            int  *firstWithDegree = new int[maxDegree+1];
            int *nFirstWithDegree = new int[maxDegree+1];
            for (i=0; i<maxDegree+1; i++) firstWithDegree[i]=-1;

            // The data structure used to sort and traverse the level sets is the following:
            // The current level set is currentLevelSet.
            // In this level set, there are nodes with degrees from 0 (not really useful)
            // to maxDegreeInCurrentLevelSet.
            // firstWithDegree[i] points to a node with degree i, or to -1 if it does not
            // exist. nextSameDegree[firstWithDegree[i]] points to the second node with
            // that degree, etc.
            // While the level set is being traversed, the structure for the next level
            // set is generated; nMDICLS will be the next maxDegreeInCurrentLevelSet
            // and nFirstWithDegree will be firstWithDegree.

            // Initialize the first level set, made up by initialNode alone

            perm[0]= initialNode;
            currentLevelSet= 1;
            levelSet[initialNode]= currentLevelSet;
            maxDegreeInCurrentLevelSet= degree[initialNode];
            firstWithDegree[maxDegreeInCurrentLevelSet]= initialNode;
            next= 1;

            // Main loop

            while (next<size)
            {
                nMDICLS= 0;
                for (i=0; i<maxDegree+1; i++) nFirstWithDegree[i]=-1;
                empty= true; // used to detect different connected components

                if (ordering==1)
                {
                    // Cuthill-McKee ordering
                    firstVal = 0;
                    finalVal = maxDegreeInCurrentLevelSet+1;
                    increment= 1;
                }
                else
                {
                    // reverse Cuthill-McKee ordering
                    firstVal = maxDegreeInCurrentLevelSet;
                    finalVal = -1;
                    increment= -1;
                }


                for (soughtDegree=firstVal; soughtDegree!=finalVal; soughtDegree+=increment)
                {
                    node= firstWithDegree[soughtDegree];
                    while (node!=-1)
                    {
                        // visit neighbors
                        TSparseSpaceType::GraphNeighbors(node, A, neighbors);
                        for (indexj=0; indexj<degree[node]; indexj++)
                        {
                            j    = neighbors[indexj];
                            if (levelSet[j]==0)
                            {
                                levelSet[j]= currentLevelSet+1;
                                perm[next]= j;
                                next++;
                                empty= false; // This level set is not empty
                                nextSameDegree[j]= nFirstWithDegree[degree[j]];
                                nFirstWithDegree[degree[j]]= j;
                                if (nMDICLS<degree[j]) nMDICLS=degree[j];
                            }
                        }
                        node= nextSameDegree[node];
                    }
                }



                currentLevelSet++;
                maxDegreeInCurrentLevelSet= nMDICLS;
                for (i=0; i<nMDICLS+1; i++) firstWithDegree[i]= nFirstWithDegree[i];

                if (empty==true)
                {
                    // The graph contains another connected component that we cannot reach.
                    // Search for a node that has not yet been included in a level set,
                    // and start exploring from it.
                    found= false;
                    for (i=0; i<size; i++)
                    {
                        if (levelSet[i]==0)
                        {
                            perm[next]=i;
                            next++;
                            levelSet[i]= currentLevelSet;
                            maxDegreeInCurrentLevelSet= degree[i];
                            firstWithDegree[maxDegreeInCurrentLevelSet]= i;
                            found= true;
                            break;
                        }
                    }
                    if (found==false)
                    {
                        // There must be a weird problem somewhere, because we cannot find a
                        // new connected component.
                        throw std::runtime_error("This cannot happen! Internal consistency error at LUSkylineFactorization::copyFromCSRMatrix(CSRMatrix A)");
                    }
                }

            }
            delete[] degree;
            delete[] levelSet;
            delete[] nextSameDegree;

            delete[] firstWithDegree;
            delete[] nFirstWithDegree;
        }

        // find inverse permutation
        for (i=0; i<size; i++) invperm[perm[i]]=i;

        // Let us find how large the rows of L and the columns of U should be.
        // Provisionally, we will store in rowIndex[i] the minimum required height
        // of column i over the diagonal, and length of row i below the diagonal.
        // The entry (i,j) in the reordered matrix will be the same as the entry
        // (perm[i],perm[j]) in the original matrix; or, the entry (i,j) in the
        // original matrix will be the same as (invperm[i],invperm[j]) in the
        // reordered matrix.
        rowIndex = new int[size+1];
        for (i=0; i<=size; i++) rowIndex[i]=0;
        // Traverse the matrix finding nonzero elements
        for (typename SparseMatrixType::iterator1 a_iterator = A.begin1();
                a_iterator != A.end1(); a_iterator++)
        {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            for (typename SparseMatrixType::iterator2 row_iterator = a_iterator.begin() ;
                    row_iterator != a_iterator.end() ; ++row_iterator)
            {
#else
            for ( SparseMatrixType::iterator2 row_iterator = begin(a_iterator,
                    boost::numeric::ublas::iterator1_tag());
                    row_iterator != end(a_iterator,
                                        boost::numeric::ublas::iterator1_tag()); ++row_iterator )
            {
#endif
                i = row_iterator.index1();
                j = row_iterator.index2();
                entry = *row_iterator;
                newi = invperm[i];
                newj = invperm[j];
                if( entry != 0 )
                {
                    if (newi>newj)
                    {
                        // row newi needs length at least newi-newj
                        if (rowIndex[newi] < newi-newj)  rowIndex[newi]= newi-newj;
                    }
                    else if (newi<newj)
                    {
                        // column newj needs height at least newj-newi
                        if (rowIndex[newj] < newj-newi)  rowIndex[newj]= newj-newi;
                    }
                }
            }
        }

        // Transform rowIndex so that it doesn't contain the required lengths
        // and heights, but the indexes to the entries
        int oldCarryOut=0;
        for (i=1; i<=size; i++)
        {
            int newCarryOut= rowIndex[i];
            rowIndex[i]= rowIndex[i-1]+oldCarryOut;
            oldCarryOut= newCarryOut;
        }

        // Allocate variables for skyline format entries
        entriesL= new double[rowIndex[size]];
        entriesD= new double[         size ];
        entriesU= new double[rowIndex[size]];
        for (i=0; i<rowIndex[size]; i++)
        {
            entriesL[i]=0.0;
            entriesU[i]=0.0;
        }

        for (i=0; i<size; i++) entriesD[i]=0.0;

        // And finally traverse again the CSR matrix, copying its entries into the
        // correct places in the skyline format
        for (typename SparseMatrixType::iterator1 a_iterator = A.begin1();
                a_iterator != A.end1(); a_iterator++)
        {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            for (typename SparseMatrixType::iterator2 row_iterator = a_iterator.begin() ;
                    row_iterator != a_iterator.end() ; ++row_iterator)
            {
#else
            for ( SparseMatrixType::iterator2 row_iterator = begin(a_iterator,
                    boost::numeric::ublas::iterator1_tag());
                    row_iterator != end(a_iterator,
                                        boost::numeric::ublas::iterator1_tag()); ++row_iterator )
            {
#endif
                i = row_iterator.index1();
                j = row_iterator.index2();
                entry = *row_iterator;
                newi= invperm[i];
                newj= invperm[j];
                if (entry != 0.00)
                {
                    if (newi<newj)
                    {
                        entriesU[ rowIndex[newj+1] + newi - newj ]= entry;
                    }
                    else if (newi==newj)
                    {
                        entriesD[newi]= entry;
                    }
                    else
                    {
                        entriesL[ rowIndex[newi+1] + newj - newi ]= entry;
                    }
                }
            }
        }

        delete[] invperm;

    }//copyFromCSRMatrix

    /**
     * Perform and in-place LU factorization of a skyline matrix by Crout's
     * algorithm. The diagonal of U contains the 1's.
     * The equivalent MATLAB code for a full matrix would be:
     * for k=1:n-1
     * A(1,k+1)=A(1,k+1)/A(1,1);
     * for i=2:k
     * sum=A(i,k+1); for j=1:i-1; sum=sum-A(i,j)*A(j,k+1); end; A(i,k+1)=sum/A(i,i);
     * end
     * for i=2:k
     * sum=A(k+1,i); for j=1:i-1; sum=sum-A(j,i)*A(k+1,j); end; A(k+1,i)=sum;
     * end
     * sum=A(k+1,k+1); for i=1:k; sum=sum-A(k+1,i)*A(i,k+1); end; A(k+1,k+1)=sum;
     * end
     */
    void factorize()
    {
        int i, j, k;
        int indexEntry, indexL, indexU;
        int jBeginRow, jBeginCol, jBeginMult, iBeginCol;
        double sum;


        if(entriesD[0] == 0.00)
        {
            std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LUSkylineFactorization::factorize: Error zero in diagonal" << std::endl;
        } //Roport Error zero in diagonal!!!!!!!!!!!

        for (k=0; k<size-1; k++)
        {
            // check whether A(1,k+1) lies within the skyline structure
            if (rowIndex[k+1] +k+1 == rowIndex[k+2])
            {
                entriesU[rowIndex[k+1]] = entriesU[rowIndex[k+1]] / entriesD[0] ;
            }
            // Compute column k+1 of U
            indexEntry= rowIndex[k+1];
            iBeginCol= k+1-rowIndex[k+2]+rowIndex[k+1];
            i= iBeginCol;
            while (i<=k)
            {
                // if i==0, there is nothing to do
                if (i>0)
                {
                    sum= entriesU[indexEntry]; // this is element U(i,k+1)
                    // Multiply row i of L and Column k+1 of U
                    jBeginRow= i-rowIndex[i+1]+rowIndex[i];
                    jBeginMult= ( iBeginCol > jBeginRow ) ? iBeginCol : jBeginRow ;
                    indexL= rowIndex[i  ] + jBeginMult - jBeginRow;
                    indexU= rowIndex[k+1] + jBeginMult - iBeginCol;
                    for (j=jBeginMult; j<i; j++)
                    {
                        sum= sum - entriesL[indexL]*entriesU[indexU];
                        indexL++;
                        indexU++;
                    }
                    entriesU[indexEntry]= sum/entriesD[i];
                }
                indexEntry++;
                i++;
            }
            // Compute row k+1 of L
            indexEntry= rowIndex[k+1];
            jBeginRow= k+1-rowIndex[k+2]+rowIndex[k+1];
            i= iBeginCol;
            while (i<=k)
            {
                // if i==0, there is nothing to do
                if (i>0)
                {
                    sum=entriesL[indexEntry]; // this is the element L(k+1,i)
                    // multiply row k+1 of L and column i of U
                    jBeginCol = i-rowIndex[i+1]+rowIndex[i];
                    jBeginMult= ( jBeginCol > jBeginRow ) ? jBeginCol : jBeginRow ;
                    indexL= rowIndex[k+1] + jBeginMult - jBeginRow;
                    indexU= rowIndex[i  ] + jBeginMult - jBeginCol;
                    for (j=jBeginMult; j<i; j++)
                    {
                        sum= sum - entriesL[indexL]*entriesU[indexU];
                        indexL++;
                        indexU++;
                    }
                    entriesL[indexEntry]= sum; // !!!!!!!!!! optimizar, no hace falta usar sum
                }
                indexEntry++;
                i++;
            }
            // Find element in diagonal
            sum= entriesD[k+1];
            for (j=rowIndex[k+1]; j<rowIndex[k+2]; j++)
            {
                sum= sum - entriesL[j] * entriesU[j];
            }
            if(sum == 0.00)
            {
                //std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LUSkylineFactorization::factorize: Error zero sum" << std::endl;
            } // Error reporting !!!!!!!!!!!!!!!!!!!

            entriesD[k+1]= sum;
        }
    }


    //***************************************************************************************************
    //*
    //***************************************************************************************************
    /* bugfix janosch void backForwardSolve(int vector_size, Vector<double>& b, Vector<double>& x) // n, b, x*/
    void backForwardSolve(int vector_size, const VectorType& b, VectorType& x) // n, b, x
    {
        // y = L^-1 * perm[b] ;
        // y = U^-1 * y ;
        // x = invperm[y];
        int i, j, indexL, indexU;
        double sum, *y;
        if (this->size != vector_size)
        {
            throw std::runtime_error("matrix and vector have different sizes at LUSkylineFactorization::backForwardSolve");
        }

        y= new double[size];
        for (i=0; i<size; i++)
        {
            j= i-rowIndex[i+1]+rowIndex[i];
            sum= b[perm[i]];
            for (indexL=rowIndex[i]; indexL<rowIndex[i+1]; indexL++)
            {
                sum= sum - entriesL[indexL]*y[j];
                j++;
            }
            y[i]= sum/entriesD[i];
        }
        for (j=size-1; j>=0; j--)
        {
            i= j-rowIndex[j+1]+rowIndex[j];
            for (indexU=rowIndex[j]; indexU<rowIndex[j+1]; indexU++)
            {
                y[i]= y[i] - entriesU[indexU] * y[j] ;
                i++;
            }
        }
        for (i=0; i<size; i++) x[perm[i]]= y[i];

        delete[] y;
    };






    //**************************************************************************
    //**************************************************************************
    LUSkylineFactorizationV2()
    {
        size=0;
        rowIndex= NULL;
        entriesL= entriesD= entriesU= NULL;
    };

    //**************************************************************************
    //**************************************************************************
    ~LUSkylineFactorizationV2()
    {
        clear();
    };
};//class LUSkylineFactorization



template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class SkylineLUFactorizationLinearSolverV2
    : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( SkylineLUFactorizationLinearSolverV2 );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

	typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
	
	typedef LUSkylineFactorizationV2<TSparseSpaceType, TDenseSpaceType> FactorizationType;
	
public:
	
    SkylineLUFactorizationLinearSolverV2()
		: m_pSolver(NULL)
	{
	}

    virtual ~SkylineLUFactorizationLinearSolverV2()
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
		m_pSolver = new FactorizationType();
		m_pSolver->copyFromCSRMatrix(rA);
        m_pSolver->factorize();
    }

    void PerformSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        const int size = TSparseSpaceType::Size(rX);
        m_pSolver->backForwardSolve(size, rB, rX);
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
    }

    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        if(this->IsNotConsistent(rA, rX, rB)) return false;
		int n = TSparseSpaceType::Size(rX);
		if(n < 1) return true;

		this->InitializeSolutionStep(rA, rX, rB);
		this->PerformSolutionStep(rA, rX, rB);
		this->FinalizeSolutionStep(rA, rX, rB);

		return true;
    }

    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
		if(this->IsNotConsistent(rA, rX, rB)) return false;
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
		std::cout << "Skyline LU Factorization Linear\n";
    }

    void  PrintData(std::ostream& rOStream) const
    {
    }

private:

    SkylineLUFactorizationLinearSolverV2& operator=(const SkylineLUFactorizationLinearSolverV2& Other);
	
private:

	FactorizationType* m_pSolver;

}; // Class SkylineLUFactorizationLinearSolverV2


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, SkylineLUFactorizationLinearSolverV2<TSparseSpaceType,
                                  TDenseSpaceType,
                                  TReordererType>& rThis)
{
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SkylineLUFactorizationLinearSolverV2<TSparseSpaceType,
                                  TDenseSpaceType,
                                  TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_SKYLINE_LINEAR_SOLVER_V2_H_INCLUDED  defined 
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Ruben Zorrilla
//

#pragma once

// External includes

// Project includes
#include "future/linear_solvers/direct_solver.h"
#include "future/linear_solvers/reorderer.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"

namespace Kratos::Future
{

template< class TMatrixType, class TVectorType >
class LUSkylineFactorization
{
public:

    using VectorType = TVectorType;

    using SparseMatrixType = TMatrixType;

    using DataType = typename TMatrixType::DataType;

    using IndexType = typename TMatrixType::IndexType;

    using DenseMatrixType = DenseMatrix<DataType>;

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

    void copyFromCSRMatrix(SparseMatrixType& A)
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
                degree[i] = A.GraphDegree(i);
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
                        A.GraphNeighbours(node, neighbors);
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
        for (IndexType i = 0; i < A.size1(); ++i) {
            const IndexType row_begin = A.index1_data()[i];
            const IndexType row_end = A.index1_data()[i+1];
            for (IndexType j = row_begin; j < row_end; ++j) {
                entry = A.value_data()[j];
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

        // And finally traverse again the CSR matrix, copying its entries into the correct places in the skyline format
        for (IndexType i = 0; i < A.size1(); ++i) {
            const IndexType row_begin = A.index1_data()[i];
            const IndexType row_end = A.index1_data()[i+1];
            for (IndexType j = row_begin; j < row_end; ++j) {
                entry = A.value_data()[j];
                newi = invperm[i];
                newj = invperm[j];
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

    }

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
                std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LUSkylineFactorization::factorize: Error zero sum" << std::endl;
            } // Error reporting !!!!!!!!!!!!!!!!!!!

            entriesD[k+1]= sum;
        }
    }


    //***************************************************************************************************
    //*
    //***************************************************************************************************
    /* bugfix janosch void backForwardSolve(int vector_size, Vector<double>& b, Vector<double>& x) // n, b, x*/
    void backForwardSolve(
        int vector_size,
        const VectorType& b,
        VectorType& x) // n, b, x
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
    LUSkylineFactorization()
    {
        size=0;
        rowIndex= NULL;
        entriesL= entriesD= entriesU= NULL;
    };

    //**************************************************************************
    //**************************************************************************
    ~LUSkylineFactorization()
    {
        clear();
    };
};//class LUSkylineFactorization



template<class TMatrixType = CsrMatrix<>, class TVectorType = SystemVector<>, class TReordererType = Reorderer<TMatrixType, TVectorType> >
class SkylineLUFactorizationSolver : public Future::DirectSolver<TMatrixType, TVectorType, TReordererType>
{
public:

    /// Counted pointer of SkylineLUFactorizationSolver
    KRATOS_CLASS_POINTER_DEFINITION(SkylineLUFactorizationSolver);

    typedef Future::DirectSolver<TMatrixType, TVectorType, TReordererType> BaseType;

    using VectorType = TVectorType;

    using SparseMatrixType = TMatrixType;

    using DenseMatrixType = typename BaseType::DenseMatrixType;

    /// Default constructor
    SkylineLUFactorizationSolver() = default;

    SkylineLUFactorizationSolver(Parameters settings)
    : BaseType(settings)
    {
    }

    /// Copy constructor
    SkylineLUFactorizationSolver(const SkylineLUFactorizationSolver& Other) = delete;

    /// Assignment operator
    SkylineLUFactorizationSolver& operator=(const SkylineLUFactorizationSolver& Other) = delete;

    /// Destructor.
    ~SkylineLUFactorizationSolver() override = default;

    /** Normal solve method.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector.
    @param rB. Right hand side vector.
    */
    bool Solve(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB) override
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

        const int size = rX.size();

        // define an object to store skyline matrix and factorization
        LUSkylineFactorization<TMatrixType, TVectorType> myFactorization;

        // copy myMatrix into skyline format
        myFactorization.copyFromCSRMatrix(rA);

        // factorize it
        myFactorization.factorize();

        // and back solve
        myFactorization.backForwardSolve(size, rB, rX);

        return true;
    }

    /** Multi solve method for solving a set of linear systems with same coefficient matrix.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector.
    @param rB. Right hand side vector.
    */
    bool Solve(
        SparseMatrixType& rA,
        DenseMatrixType& rX,
        DenseMatrixType& rB) override
    {
        const int size1 = rX.size1();
        const int size2 = rX.size2();

        bool is_solved = true;

        VectorType x(size1);
        VectorType b(size1);

        // define an object to store skyline matrix and factorization
        LUSkylineFactorization<TMatrixType, TVectorType> myFactorization;
        // copy myMatrix into skyline format
        myFactorization.copyFromCSRMatrix(rA);
        // factorize it
        myFactorization.factorize();

        for(int i = 0 ; i < size2 ; i++)
        {
            for (IndexType i_row = 0; i_row < size1; ++i_row) {
                x[i_row] = rX(i_row,i);
                b[i_row] = rB(i_row,i);
            }

            // and back solve
            myFactorization.backForwardSolve(size1, b, x);

            for (IndexType i_row = 0; i_row < size1; ++i_row) {
                rX(i_row,i) = x[i_row];
                rB(i_row,i) = b[i_row];
            }
        }

        return is_solved;
    }

    /// Print information about this object.
    void  PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "LU factorization solver finished.";
    }

    /// Print object's data.
    void  PrintData(std::ostream& rOStream) const override
    {
    }

}; // Class SkylineLUFactorizationSolver


/// input stream function
template<class TMatrixType, class TVectorType,class TReordererType>
inline std::istream& operator >> (
    std::istream& rIStream,
    SkylineLUFactorizationSolver<TMatrixType,TVectorType,TReordererType>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TMatrixType, class TVectorType, class TReordererType>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const SkylineLUFactorizationSolver<TMatrixType, TVectorType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos::Future.

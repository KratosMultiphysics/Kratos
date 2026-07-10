//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:         BSD License 
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/global_pointer_variables.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"

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

/** 
* @brief This class defines utility for aggregation of node clusters to be used in deflated solvers.
*/
class DeflationUtils
{
public:
    ///@name Type Definitions
    ///@{

    using SparseMatrixType = boost::numeric::ublas::compressed_matrix<double>;

    using SparseVectorType = boost::numeric::ublas::vector<double>;

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///visualize aggregates. This function assumes that neighbours are calculated and builds the connectivity matrix
    ///then writes it to a nodal variable so that it can be used for visualizing it.
    void VisualizeAggregates(ModelPart::NodesContainerType& rNodes, Variable<double>& rVariable, const int max_reduced_size)
    {
        SparseMatrixType A(rNodes.size(),rNodes.size());
        SparseMatrixType Adeflated;

        //first of all build the connectivty matrix
        std::vector< std::vector<int> > index_list(rNodes.size());

        std::size_t total_size = 0;

        //renumber nodes consecutively
        int new_id = 1;
        for(ModelPart::NodesContainerType::iterator in = rNodes.begin();  in!=rNodes.end(); in++)
            in->SetId(new_id++);


        //constructing the system matrix row by row
        std::size_t index_i;
        for(ModelPart::NodesContainerType::iterator in = rNodes.begin();
                in!=rNodes.end(); in++)
        {
            index_i = (in)->Id()-1;
            auto& neighb_nodes = in->GetValue(NEIGHBOUR_NODES);

            std::vector<int>& indices = index_list[index_i];
            indices.reserve(neighb_nodes.size()+1);

            //filling the first neighbours list
            indices.push_back(index_i);
            for( auto i = neighb_nodes.begin();
                    i != neighb_nodes.end(); i++)
            {

                int index_j = i->Id()-1;
                indices.push_back(index_j);
            }

            //sorting the indices and elminating the duplicates
            std::sort(indices.begin(),indices.end());
            std::vector<int>::iterator new_end = std::unique(indices.begin(),indices.end());

            indices.erase(new_end,indices.end());

            total_size += indices.size();

        }

        A.reserve(total_size,false);

        //setting to zero the matrix (and the diagonal matrix)
        for(unsigned int i=0; i<A.size1(); i++)
        {
            std::vector<int>& indices = index_list[i];
            for(unsigned int j=0; j<indices.size(); j++)
            {
                A.push_back(i,indices[j] , 0.00);
            }
        }

        //now call aggregation function to color the nodes
        std::vector<int> w(rNodes.size());
        ConstructW(max_reduced_size, A, w, Adeflated);

        //finally write the color to the nodes so that it can be visualized
        int counter = 0;
        for(ModelPart::NodesContainerType::iterator in=rNodes.begin(); in!=rNodes.end(); in++)
        {
            in->FastGetSolutionStepValue(rVariable) = w[counter++];
        }
    }

    ///this function constructs the structure of a smaller matrix using a technique taken from MIS aggregation
    ///the algorythm is taken from the pyamg lib
    static void ConstructW(const std::size_t max_reduced_size, SparseMatrixType& rA, std::vector<int>& w, SparseMatrixType&  deflatedA)
    {
        KRATOS_TRY

        std::size_t full_size = rA.size1();
        w.resize(full_size,0);

        //call aggregation function to fill mw with "colors"
        std::size_t reduced_size = standard_aggregation<int>(rA.size1(),rA.index1_data().begin(), rA.index2_data().begin(), &w[0]);

        // Non-zero structure of deflatedA

        std::vector<std::set<std::size_t> > deflatedANZ(reduced_size);

        // Loop over non-zero structure of A and build non-zero structure of deflatedA
        SparseMatrixType::iterator1 a_iterator = rA.begin1();

        for (std::size_t i = 0; i < full_size; i++)
        {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            for (SparseMatrixType::iterator2 row_iterator = a_iterator.begin() ;
                    row_iterator != a_iterator.end() ; ++row_iterator)
            {
#else
            for (typename SparseMatrixType::iterator2 row_iterator = begin(a_iterator,
                    boost::numeric::ublas::iterator1_tag());
                    row_iterator != end(a_iterator,
                                        boost::numeric::ublas::iterator1_tag()); ++row_iterator )
            {
#endif
                deflatedANZ[w[a_iterator.index1()]].insert(w[row_iterator.index2()]);
            }

            a_iterator++;
        }

        // Count the number of non-zeros in deflatedA
        std::size_t NZ = 0;
        for (std::size_t i = 0; i < reduced_size; i++)
            NZ += deflatedANZ[i].size();

        deflatedA = SparseMatrixType(reduced_size, reduced_size,NZ);

        // Insert the non-zero structure into deflatedA
        for(std::size_t i = 0 ; i < reduced_size ; i++)
        {
            for(std::set<std::size_t>::iterator j = deflatedANZ[i].begin() ; j != deflatedANZ[i].end() ; j++)
            {
                deflatedA.push_back(i,*j, 0.00);
            }
        }

        if(reduced_size > max_reduced_size)
        {
            SparseMatrixType Areduced;
            std::vector<int> wsmaller;

            //need to reduce further!! - do it recursively
            ConstructW(max_reduced_size, deflatedA, wsmaller, Areduced);

            //now change deflatedA and w on the coarser size
            for(unsigned int i=0; i<full_size; i++)
            {
                int color = w[i];
                int new_color = wsmaller[color];
                w[i] = new_color;
            }
            deflatedA.clear();
            deflatedA = Areduced;

            reduced_size = wsmaller.size();
        }

        KRATOS_CATCH("")
    }

    ///block version of ConstructW. To be used when multiple DOFS are associated to the same node.
    static void ConstructW(const int max_reduced_size, SparseMatrixType& rA, std::vector<int>& w, SparseMatrixType&  deflatedA, const std::size_t block_size)
    {
        if(block_size == 1)
        {
            ConstructW(max_reduced_size,rA, w, deflatedA);
        }
        else
        {
            //simple checks to verify blocks are effectively respected
            KRATOS_ERROR_IF(rA.size1()%block_size != 0 || rA.size2()%block_size != 0) << "The number of rows is not a multiple of block_size. Can not use the block deflation" <<  std::endl;
            KRATOS_ERROR_IF(rA.nnz()%block_size != 0) << "The number of non zeros is not a multiple of block_size. Can not use the block deflation" << std::endl;

            //construct Ascalar
            SparseMatrixType Ascalar;
            ConstructScalarMatrix(rA.size1(),block_size,rA.index1_data().begin(), rA.index2_data().begin(), Ascalar);

            //deflate it using the standard methods
            SparseMatrixType deflatedAscalar;
            std::vector<int> wscalar;
            ConstructW(max_reduced_size/block_size,Ascalar, wscalar, deflatedAscalar);

            //compute w for the block structured problem

            std::vector<int> w(wscalar.size()*block_size);
            for(std::size_t i=0; i<wscalar.size(); i++)
            {
                for(std::size_t j=0; j<block_size; j++)
                {
                    w[i*block_size + j] = wscalar[i]*block_size+j;
                }
            }

            //compute deflatedA
            SparseMatrixType deflatedA(deflatedAscalar.size1()*block_size,deflatedAscalar.size2()*block_size);
            //do reserve!!
            deflatedA.reserve(deflatedAscalar.nnz()*block_size*block_size);
            ExpandScalarMatrix(rA.size1(),block_size,rA.index1_data().begin(), rA.index2_data().begin(), deflatedA);
        }
    }

    //W is the deflation matrix, stored as a single vector of indices
    //y is a vector of "full" size
    //x is a vector of reduced size
    //y = W*x;
    static void ApplyW(const std::vector<int>& w, const SparseVectorType& x, SparseVectorType& y)
    {
        IndexPartition<std::size_t>(w.size()).for_each([&](std::size_t i) {
            y[i] = x[w[i]];
        });
    }

    //W is the deflation matrix, stored as a single vector of indices
    //y is a vector of "reduced" size
    //s is a vector of "full" size
    //y = Wtranspose*x;
    static void ApplyWtranspose(const std::vector<int>& w, const SparseVectorType& x, SparseVectorType& y)
    {
        //first set to zero the destination vector
        IndexPartition<std::size_t>(y.size()).for_each([&](std::size_t i) {
            y[i] = 0.0;
        });

        //now apply the Wtranspose
        IndexPartition<std::size_t>(w.size()).for_each([&](std::size_t i) {
            AtomicAdd(y[w[i]], x[i]);
        });
    }

    //*******************************************************************************
    //*******************************************************************************
    static void FillDeflatedMatrix( const SparseMatrixType& rA, std::vector<int>& w, SparseMatrixType&  Ah)
    {
        KRATOS_TRY

        double* abegin = Ah.value_data().begin();
        const int size = Ah.value_data().size();
        IndexPartition<int>(size).for_each([&](int i) {
            *(abegin+i) = 0.0;
        });
        // TSparseSpaceType::SetToZero(Ah);

        // Now building Ah
        SparseMatrixType::const_iterator1 a_iterator = rA.begin1();
        int full_size = rA.size1();

        for (int i = 0; i < full_size; i++)
        {
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
            for (SparseMatrixType::const_iterator2 row_iterator = a_iterator.begin() ;
                    row_iterator != a_iterator.end() ; ++row_iterator)
            {
#else
            for (typename SparseMatrixType::iterator2 row_iterator = begin(a_iterator,
                    boost::numeric::ublas::iterator1_tag());
                    row_iterator != end(a_iterator,
                                        boost::numeric::ublas::iterator1_tag()); ++row_iterator )
            {
#endif
                Ah(w[a_iterator.index1()], w[row_iterator.index2()]) += *row_iterator;
            }
            a_iterator++;
        }

        KRATOS_CATCH("");
    }


    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Acces
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
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

    /*
    * Compute aggregates for a matrix A stored in CSR format
    *
    * Parameters:
    *   n_row         - number of rows in A
    *   Ap[n_row + 1] - CSR row pointer
    *   Aj[nnz]       - CSR column indices
    *    x[n_row]     - aggregate numbers for each node
    *
    * Returns:
    *  The number of aggregates (== max(x[:]) + 1 )
    *
    * Notes:
    *    It is assumed that A is structurally symmetric.
    *    A may contain diagonal entries (self loops)
    *    Unaggregated nodes are marked with a -1
    */
    template <class I>
    static I standard_aggregation(const I n_row,
                                  const std::size_t Ap[],
                                  const std::size_t Aj[],
                                  I  x[])
    {
        // Bj[n] == -1 means i-th node has not been aggregated
        std::fill(x, x + n_row, 0);

        I next_aggregate = 1; // number of aggregates + 1

        //Pass #1
        for(I i = 0; i < n_row; i++)
        {
            if(x[i])
            {
                continue;    //already marked
            }

            const I row_start = Ap[i];
            const I row_end   = Ap[i+1];

            //Determine whether all neighbors of this node are free (not already aggregates)
            bool has_aggregated_neighbors = false;
            bool has_neighbors            = false;
            for(I jj = row_start; jj < row_end; jj++)
            {
                const I j = Aj[jj];
                if( i != j )
                {
                    has_neighbors = true;
                    if( x[j] )
                    {
                        has_aggregated_neighbors = true;
                        break;
                    }
                }
            }

            if(!has_neighbors)
            {
                //isolated node, do not aggregate
                x[i] = -n_row;
            }
            else if (!has_aggregated_neighbors)
            {
                //Make an aggregate out of this node and its neighbors
                x[i] = next_aggregate;
                for(I jj = row_start; jj < row_end; jj++)
                {
                    x[Aj[jj]] = next_aggregate;
                }
                next_aggregate++;
            }
        }

        //Pass #2
        // Add unaggregated nodes to any neighboring aggregate
        for(I i = 0; i < n_row; i++)
        {
            if(x[i])
            {
                continue;    //already marked
            }

            for(I jj = static_cast<I>(Ap[i]); jj < static_cast<I>(Ap[i+1]); jj++)
            {
                const I j = Aj[jj];

                const I xj = x[j];
                if(xj > 0)
                {
                    x[i] = -xj;
                    break;
                }
            }
        }

        next_aggregate--;

        //Pass #3
        for(I i = 0; i < n_row; i++)
        {
            const I xi = x[i];

            if(xi != 0)
            {
                // node i has been aggregated
                if(xi > 0)
                    x[i] = xi - 1;
                else if(xi == -n_row)
                    x[i] = -1;
                else
                    x[i] = -xi - 1;
                continue;
            }

            // node i has not been aggregated
            const I row_start = Ap[i];
            const I row_end   = Ap[i+1];

            x[i] = next_aggregate;

            for(I jj = row_start; jj < row_end; jj++)
            {
                const I j = Aj[jj];

                if(x[j] == 0)  //unmarked neighbors
                {
                    x[j] = next_aggregate;
                }
            }
            next_aggregate++;
        }


        return next_aggregate; //number of aggregates
    }




    static void ConstructScalarMatrix(const std::size_t n_row, const std::size_t block_size,
                                      const std::size_t Ap[],
                                      const std::size_t Aj[],
                                      SparseMatrixType& Ascalar
                                     )
    {
        Ascalar.resize(n_row/block_size,n_row/block_size,0.0);
        std::size_t scalar_size = (Ap[n_row]-Ap[0])/(block_size*block_size);
        Ascalar.reserve(scalar_size);

        for(std::size_t i = 0; i < n_row; i++)
        {
            if(i%block_size == 0)
            {
                std::size_t iscalar = i/block_size;
                const std::size_t row_start = Ap[i];
                const std::size_t row_end   = Ap[i+1];

                for(std::size_t jj = row_start; jj < row_end; jj++)
                {
                    const std::size_t j = Aj[jj];
                    if(j%block_size == 0)
                    {
                        std::size_t jscalar = j/block_size;
                        Ascalar.push_back(iscalar,jscalar,0.0);
                    }
                }
            }
        }
    }

    static void ExpandScalarMatrix(const std::size_t n_row, const std::size_t block_size,
                                   const std::size_t Ap[],
                                   const std::size_t Aj[],
                                   SparseMatrixType& Aexpanded
                                  )
    {
        Aexpanded.resize(n_row*block_size,n_row*block_size,0.0);
        std::size_t expanded_size = (Ap[n_row]-Ap[0])*block_size*block_size;
        Aexpanded.reserve(expanded_size);

        for(std::size_t i = 0; i < n_row; i++)
        {
            const std::size_t row_start = Ap[i];
            const std::size_t row_end   = Ap[i+1];

            for(std::size_t isub=0; isub<block_size; isub++)
            {
                std::size_t iexpanded = i*block_size + isub;
                for(std::size_t jj = row_start; jj < row_end; jj++)
                {
                    const std::size_t j = Aj[jj];

                    for(std::size_t jsub=0; jsub<block_size; jsub++)
                    {
                        std::size_t jexpanded = j*block_size+jsub;
                        Aexpanded.push_back(iexpanded,jexpanded,0.0);
                    }
                }
            }
        }
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Acces
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; /// class DeflationUtils

///@}
///@name Type Definitions
///@{

///@}

} /// namespace Kratos

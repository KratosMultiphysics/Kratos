
//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Denis Demidov
//                   Riccardo Rossi
//                   
//

#if !defined(KRATOS_CSR_CONVERSION_UTILITIES_H_INCLUDED)
#define  KRATOS_CSR_CONVERSION_UTILITIES_H_INCLUDED

#include "containers/csr_matrix.h"

#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/zero_copy.hpp>

namespace Kratos
{

/**
Utilities to convert the distributed_csr matrix to other libraries
 */
class AmgclCSRConversionUtilities
{

public:

    /**
     This function returns a shared pointer to an Amgcl distributed_matrix
     */
	template< class TDataType, class TIndexType >
	static Kratos::shared_ptr<typename amgcl::backend::builtin<TDataType>::matrix > ConvertToAmgcl(
			const CsrMatrix<TDataType, TIndexType>& rA)
	{
    	Kratos::shared_ptr<typename amgcl::backend::builtin<TDataType>::matrix > pAmgcl = amgcl::adapter::zero_copy(
                rA.size1(),
                rA.index1_data().data().begin(),
                rA.index2_data().data().begin(),
                rA.value_data().data().begin()
            );
        return pAmgcl;
	}

    template< class TDataType, class TIndexType >
	static void ConvertToCsrMatrix(
			typename amgcl::backend::builtin<TDataType>::matrix& rA,
            CsrMatrix<TDataType, TIndexType>& rAconverted
			)
	{
        //TODO: find a way to avoid this copying!! the prolem as of now is that ublas does not allow initializing by moving data
        rAconverted.index1_data().resize(rA.nrows+1,false);
        rAconverted.index2_data().resize(rA.nnz,false);
        rAconverted.value_data().resize(rA.nnz,false);

        rAconverted.SetRowSize(rA.nrows);
        rAconverted.SetColSize(rA.ncols);

        IndexPartition<IndexType>(rAconverted.index1_data().size()).for_each( [&](IndexType i){
            rAconverted.index1_data()[i] = rA.ptr[i];
        });

        IndexPartition<IndexType>(rAconverted.index2_data().size()).for_each( [&](IndexType i){
            rAconverted.index2_data()[i] = rA.col[i];
        });

        IndexPartition<IndexType>(rAconverted.value_data().size()).for_each( [&](IndexType i){
            rAconverted.value_data()[i] = rA.val[i];
        });

        rA.free_data();
        rA.nnz = 0;
        rA.ncols = 0;
        rA.nrows = 0;
	}	



};

}

#endif // KRATOS_CSR_CONVERSION_UTILITIES_H_INCLUDED

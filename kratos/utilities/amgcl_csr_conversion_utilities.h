
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

// System includes

// External includes
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <span/span.hpp>

// Project includes
#include "containers/csr_matrix.h"

namespace Kratos
{

/**
Utilities to convert the distributed_csr matrix to an AMGCL matrix
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
                rA.index1_data().begin(),
                rA.index2_data().begin(),
                rA.value_data().begin()
            );
        pAmgcl->ncols = rA.size2();
        pAmgcl->own_data = false;

        return pAmgcl;
	}

    
    //Note that we deliberately return a unique_ptr as it can be moved to a shared_ptr as needed
    template< class TDataType, class TIndexType >
	static typename CsrMatrix<TDataType, TIndexType>::UniquePointer ConvertToCsrMatrix(
			typename amgcl::backend::builtin<TDataType>::matrix& rA
			)	
    {
        auto pAconverted = Kratos::make_unique<CsrMatrix<TDataType, TIndexType>>();

        if(rA.own_data == false){ //if rA is not the owner, Aconverted cannot be
            pAconverted->SetIsOwnerOfData(false);
        }
        else{ //if rA is the owner, transfer ownership to the csr_matrix
            rA.own_data = false;
            pAconverted->SetIsOwnerOfData(true);
        }

        pAconverted->SetRowSize(rA.nrows);
        pAconverted->SetColSize(rA.ncols);
        pAconverted->AssignIndex1Data((TIndexType*)(rA.ptr), rA.nrows+1);
        pAconverted->AssignIndex2Data((TIndexType*)(rA.col), rA.nnz);
        pAconverted->AssignValueData((TDataType*)(rA.val), rA.nnz);
        return pAconverted;
	}
	
};

}

#endif // KRATOS_CSR_CONVERSION_UTILITIES_H_INCLUDED

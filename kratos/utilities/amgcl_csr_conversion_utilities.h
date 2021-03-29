
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
#include <span/span.hpp>

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

    template< class TDataType, class TIndexType >
	static CsrMatrix<TDataType, TIndexType> ConvertToCsrMatrix(
			typename amgcl::backend::builtin<TDataType>::matrix& rA
			)
	{
        CsrMatrix<TDataType, TIndexType> Aconverted;

        if(rA.own_data == false){ //if rA is not the owner, Aconverted cannot be
            Aconverted.SetIsOwnerOfData(false);
        }
        else{ //if rA is the owner, transfer ownership to the csr_matrix
            rA.own_data = false;
            Aconverted.SetIsOwnerOfData(true);
        }

        Aconverted.SetRowSize(rA.nrows);
        Aconverted.SetColSize(rA.ncols);
        Aconverted.AssignIndex1Data((TIndexType*)(rA.ptr), rA.nrows+1);
        Aconverted.AssignIndex2Data((TIndexType*)(rA.col), rA.nnz);
        Aconverted.AssignValueData((TDataType*)(rA.val), rA.nnz);
        return Aconverted;
	}	



};

}

#endif // KRATOS_CSR_CONVERSION_UTILITIES_H_INCLUDED

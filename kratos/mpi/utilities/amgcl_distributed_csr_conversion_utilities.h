
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

#if !defined(KRATOS_DISTRIBUTED_CSR_CONVERSION_UTILITIES_H_INCLUDED)
#define  KRATOS_DISTRIBUTED_CSR_CONVERSION_UTILITIES_H_INCLUDED

#include "containers/distributed_csr_matrix.h"
#include "amgcl/backend/builtin.hpp"
#include "amgcl/backend/interface.hpp"
#include "amgcl/mpi/distributed_matrix.hpp"
#include "amgcl/adapter/zero_copy.hpp"

namespace Kratos
{

/**
Utilities to convert the distributed_csr matrix to other libraries
 */
class AmgclDistributedCSRConversionUtilities
{

public:

    /**
     This function returns a shared pointer to an Amgcl distributed_matrix
     */
	template< class TDataType, class TIndexType >
	static Kratos::shared_ptr<amgcl::mpi::distributed_matrix<amgcl::backend::builtin<double>>> ConvertToAmgcl(
			const DistributedCsrMatrix<TDataType, TIndexType>& rA,
			DenseVector<TIndexType>& global_index2)
	{

		IndexType chunk = rA.local_size1();

    	auto loc_a = amgcl::adapter::zero_copy(chunk, 
                rA.GetDiagBlock().index1_data().data().begin(),
                rA.GetDiagBlock().index2_data().data().begin(),
                rA.GetDiagBlock().value_data().data().begin()
			);

		auto rem_a = amgcl::adapter::zero_copy(chunk, 
			rA.GetOffDiagBlock().index1_data().data().begin(),
			global_index2.data().begin(),
			rA.GetOffDiagBlock().value_data().data().begin()
		);

		auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator( rA.GetComm());
		amgcl::mpi::communicator comm(raw_mpi_comm);

		auto pAmgcl = Kratos::make_shared<amgcl::mpi::distributed_matrix<amgcl::backend::builtin<double>>>(comm, loc_a, rem_a);
		pAmgcl->move_to_backend();
		return pAmgcl;
	}

};

}

#endif // KRATOS_DISTRIBUTED_CSR_CONVERSION_UTILITIES_H_INCLUDED

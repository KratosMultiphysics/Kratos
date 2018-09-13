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
//
//

// System includes
#include <set>
// External includes


// Project includes
#include "testing/testing.h"
#include "includes/pool_object.h"


namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(FixedSizeMemoryPoolConstruction, KratosCoreFastSuite)
		{
			int max_threads = LockObject::GetNumberOfThreads();
			std::size_t block_size = 5;
			std::size_t default_chunk_size = 1024 * 1024; // 1M
			FixedSizeMemoryPool fixed_size_memory_pool(block_size);
			std::size_t empty_size = sizeof(ThreadFixedSizeMemoryPool) * max_threads + sizeof(FixedSizeMemoryPool);
			KRATOS_CHECK_EQUAL(fixed_size_memory_pool.GetNumberOfAllocatedChunks(), 0);
			KRATOS_CHECK_EQUAL(fixed_size_memory_pool.ChunkSize(), default_chunk_size);
			KRATOS_CHECK_EQUAL(fixed_size_memory_pool.MemoryUsed(), empty_size);
		}

		KRATOS_TEST_CASE_IN_SUITE(FixedSizeMemoryPoolAllocationDeallocation, KratosCoreFastSuite)
		{
			int max_threads = LockObject::GetNumberOfThreads();
			std::size_t block_size = 15;
			std::size_t default_chunk_size = 1024 * 1024; // 1M
			std::size_t number_of_blocks = (default_chunk_size - 2 * max_threads * sizeof(Chunk::SizeType)) / 16;
			FixedSizeMemoryPool fixed_size_memory_pool(block_size);


			std::size_t repeat_number = 2;
			for (std::size_t i_repeat = 1; i_repeat <= repeat_number; i_repeat++)
			{
				std::set<void *> pointer_set;
				for (std::size_t i_chunk = 0; i_chunk < i_repeat; i_chunk++)
					for (std::size_t i_block = 0; i_block < number_of_blocks; i_block++)
					{
						void* p = fixed_size_memory_pool.Allocate();
						KRATOS_ERROR_IF(pointer_set.find(p) != pointer_set.end())
							<< "The allocated position #" << pointer_set.size() + 1 << " was allocated before";
						pointer_set.insert(p);
					}

				for (auto i_pointer = pointer_set.begin(); i_pointer != pointer_set.end(); i_pointer++) {
					fixed_size_memory_pool.Deallocate(*i_pointer);
				}
				std::size_t number_of_chunks = i_repeat;// ((i_repeat / OpenMPUtils::GetNumThreads()) + 1) *  OpenMPUtils::GetNumThreads();
				KRATOS_CHECK_EQUAL(fixed_size_memory_pool.GetNumberOfAllocatedChunks(), number_of_chunks) << " (number_of_chunks = " << number_of_chunks << ")";
				KRATOS_CHECK_EQUAL(fixed_size_memory_pool.ChunkSize(), default_chunk_size);
			}
		}

		//KRATOS_TEST_CASE_IN_SUITE(FixedSizeMemoryPoolStressTest, KratosCoreStressSuite)
		//{
		//	int max_threads = LockObject::GetNumberOfThreads();
		//	std::size_t block_size = 64;
		//	std::size_t default_chunk_size = 1024 * 1024;// 1M
		//	std::size_t number_of_blocks = (default_chunk_size - 2 * max_threads * sizeof(Chunk::SizeType)) / 64;
		//	FixedSizeMemoryPool fixed_size_memory_pool(block_size);

		//	auto repeat_number = 100;
		//	for (auto i_repeat = 1; i_repeat <= repeat_number; i_repeat++)
		//	{
		//		KRATOS_WATCH(i_repeat);
		//		std::vector<void *> pointer_vector;
		//		for (auto i_chunk = 0; i_chunk < i_repeat; i_chunk++)
		//			for (std::size_t i_block = 0; i_block < number_of_blocks; i_block++)
		//			{
		//				void* p = fixed_size_memory_pool.Allocate();
		//				KRATOS_CHECK_NOT_EQUAL(p, nullptr) << " i_repeat : " << i_repeat << " , i_chunk : " << i_chunk << " , i_block : " << i_block;
		//				pointer_vector.push_back(p);

		//			}

		//		for (auto i_pointer = pointer_vector.begin(); i_pointer != pointer_vector.end(); i_pointer++) {
		//			fixed_size_memory_pool.Deallocate(*i_pointer);
		//		}
		//	}

		//}

		//KRATOS_TEST_CASE_IN_SUITE(MemoryPool, KratosCoreFastSuite)
		//{
		//	std::size_t block_size = 61;
		//	auto number_of_blocks = 4096;
		//	MemoryPool::AddPool("TestPool", new FixedSizeMemoryPool(block_size, number_of_blocks));

		//	KRATOS_ERROR << MemoryPool::Info() << std::endl;
		//}

		class Block64byte : public PoolObject{
            public:
                Block64byte() {Block[0] = 0;};
            private:
                double Block[8];
		};

		KRATOS_DISABLED_TEST_CASE_IN_SUITE(PoolObjectParallelNewDelete, KratosCoreFastSuite)
		{
			std::cout << MemoryPool::Info() << std::endl;
			auto repeat_number = 120;
#pragma omp parallel for
			for (auto i_repeat = 0; i_repeat < repeat_number; i_repeat++)
			{
				try {
					std::size_t size = 1024 * 1024;
					std::vector<Block64byte*> the_vector(size);

					for (std::size_t i = 0; i < size; i++)
						the_vector[i] = new Block64byte;

					for (std::size_t i = 0; i < size; i++)
						delete the_vector[i];
				}
				catch(Exception& e){
					std::cout << e << std::endl;
				}
				catch (...) {
					std::cout << "Unknown Error" << std::endl;
				}
			}

			std::cout << MemoryPool::Info() << std::endl;
		}

		//KRATOS_TEST_CASE_IN_SUITE(NormalNewDelete, KratosCoreStressSuite)
		//{
		//	auto repeat_number = 128;
		//	for (auto i_repeat = 1; i_repeat <= repeat_number; i_repeat++)
		//	{
		//		std::vector<Block64byte*> the_vector;

		//		for (auto i_chunk = 0; i_chunk < i_repeat; i_chunk++)
		//			for (auto i_block = 0; i_block < 4096; i_block++)
		//				the_vector.push_back(new Block64byte);

		//		for (std::size_t i = 0; i < the_vector.size(); i++) {
		//			delete the_vector[i];
		//		}
		//	}

		//}


	}
}  // namespace Kratos.

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
			std::size_t block_size = 5;
			unsigned char number_of_blocks = 2;
			FixedSizeMemoryPool fixed_size_memory_pool(block_size, number_of_blocks);
			std::size_t empty_size = 16 + sizeof(Chunk) + sizeof(FixedSizeMemoryPool) + sizeof(std::size_t);

			KRATOS_CHECK_EQUAL(fixed_size_memory_pool.GetNumberOfChunks(), 1);
			KRATOS_CHECK_EQUAL(fixed_size_memory_pool.ChunkSize(), 16);
			KRATOS_CHECK_EQUAL(fixed_size_memory_pool.MemoryUsed(), empty_size);
		}

		KRATOS_TEST_CASE_IN_SUITE(FixedSizeMemoryPoolAllocationDeallocation, KratosCoreFastSuite)
		{
			std::size_t block_size = 5;
			unsigned char number_of_blocks = 2;
			FixedSizeMemoryPool fixed_size_memory_pool(block_size, number_of_blocks);
			std::size_t empty_size = 16 + sizeof(Chunk) + sizeof(FixedSizeMemoryPool) + sizeof(std::size_t);

			auto repeat_number = 13;
			for (std::size_t i_repeat = 1; i_repeat <= repeat_number; i_repeat++)
			{
				std::set<void *> pointer_set;
				for (auto i_chunk = 0; i_chunk < i_repeat; i_chunk++)
					for (auto i_block = 0; i_block < number_of_blocks; i_block++)
					{
						void* p = fixed_size_memory_pool.Allocate();
						KRATOS_ERROR_IF(pointer_set.find(p) != pointer_set.end())
							<< "The allocated position #" << pointer_set.size() + 1 << " was allocated before";
						pointer_set.insert(p);
					}

				for (auto i_pointer = pointer_set.begin(); i_pointer != pointer_set.end(); i_pointer++) {
					fixed_size_memory_pool.Deallocate(*i_pointer);
				}
				KRATOS_CHECK_EQUAL(fixed_size_memory_pool.GetNumberOfChunks(), i_repeat);
				KRATOS_CHECK_EQUAL(fixed_size_memory_pool.ChunkSize(), 16);
				KRATOS_CHECK_EQUAL(fixed_size_memory_pool.MemoryUsed(), empty_size) << fixed_size_memory_pool << std::endl;
				empty_size += sizeof(Chunk) + sizeof(std::size_t);

			}
		}

		KRATOS_TEST_CASE_IN_SUITE(FixedSizeMemoryPoolStressTest, KratosCoreStressSuite)
		{
			std::size_t block_size = 61;
			auto number_of_blocks = 4096;
			FixedSizeMemoryPool fixed_size_memory_pool(block_size, number_of_blocks);

			auto repeat_number = 128;
			for (auto i_repeat = 1; i_repeat <= repeat_number; i_repeat++)
			{
				std::vector<void *> pointer_vector;
				for (auto i_chunk = 0; i_chunk < i_repeat; i_chunk++)
					for (auto i_block = 0; i_block < number_of_blocks; i_block++)
					{
						void* p = fixed_size_memory_pool.Allocate();
						pointer_vector.push_back(p);

					}

				for (auto i_pointer = pointer_vector.begin(); i_pointer != pointer_vector.end(); i_pointer++) {
					fixed_size_memory_pool.Deallocate(*i_pointer);
				}
			}

		}

		//KRATOS_TEST_CASE_IN_SUITE(MemoryPool, KratosCoreFastSuite)
		//{
		//	std::size_t block_size = 61;
		//	auto number_of_blocks = 4096;
		//	MemoryPool::AddPool("TestPool", new FixedSizeMemoryPool(block_size, number_of_blocks));

		//	KRATOS_ERROR << MemoryPool::Info() << std::endl;
		//}

		class Block46byte : public PoolObject{
			double block[8];
		};

		KRATOS_TEST_CASE_IN_SUITE(PoolObjectNew, KratosCoreFastSuite)
		{
			std::cout << MemoryPool::Info() << std::endl;
			Block46byte* p_new = new Block46byte;

			std::cout << MemoryPool::Info() << std::endl;

			delete p_new;
			std::cout << MemoryPool::Info() << std::endl;
		}

		KRATOS_TEST_CASE_IN_SUITE(NormalNewDelete, KratosCoreStressSuite)
		{
			auto repeat_number = 128;
			for (auto i_repeat = 1; i_repeat <= repeat_number; i_repeat++)
			{
				std::vector<Block46byte*> the_vector;

				for (auto i_chunk = 0; i_chunk < i_repeat; i_chunk++)
					for (auto i_block = 0; i_block < 4096; i_block++)
						the_vector.push_back(new Block46byte);

				for (std::size_t i = 0; i < the_vector.size(); i++) {
					delete the_vector[i];
				}
			}

		}


	}
}  // namespace Kratos.

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
#include "utilities/timer.h"


namespace Kratos {
	namespace Testing {


		KRATOS_TEST_CASE_IN_SUITE(ThreadFixedSizeMemoryPoolAllocationDeallocation, KratosCoreFastSuite)
		{
			std::size_t block_size = 15;
			std::size_t default_chunk_size = 1024 * 1024; // 1M
			std::size_t number_of_blocks = (default_chunk_size) / 16;
			std::size_t thread_id = 0;
			ThreadFixedSizeMemoryPool thread_fixed_size_memory_pool(block_size,default_chunk_size, thread_id);


			std::size_t repeat_number = 2;
			for (std::size_t i_repeat = 1; i_repeat <= repeat_number; i_repeat++)
			{
				std::set<void *> pointer_set;
				for (std::size_t i_chunk = 0; i_chunk < i_repeat; i_chunk++)
					for (std::size_t i_block = 0; i_block < number_of_blocks; i_block++)
					{
						void* p = thread_fixed_size_memory_pool.Allocate();
						KRATOS_ERROR_IF(pointer_set.find(p) != pointer_set.end())
							<< "The allocated position #" << pointer_set.size() + 1 << " was allocated before";
						pointer_set.insert(p);
					}

				for (auto i_pointer = pointer_set.begin(); i_pointer != pointer_set.end(); i_pointer++) {
					thread_fixed_size_memory_pool.Deallocate(*i_pointer);
				}
				std::size_t number_of_chunks = i_repeat;// ((i_repeat / ParallelUtilities::GetNumThreads()) + 1) *  ParallelUtilities::GetNumThreads();
				KRATOS_EXPECT_EQ(thread_fixed_size_memory_pool.GetNumberOfChunks(), number_of_chunks) << " (number_of_chunks = " << number_of_chunks << ")";
				KRATOS_EXPECT_EQ(thread_fixed_size_memory_pool.ChunkSize(), default_chunk_size);
			}
		}

		KRATOS_TEST_CASE_IN_SUITE(FixedSizeMemoryPoolConstruction, KratosCoreFastSuite)
		{
			std::size_t block_size = 5;
			std::size_t default_chunk_size = 1024 * 1024; // 1M
			FixedSizeMemoryPool fixed_size_memory_pool(block_size);
			KRATOS_EXPECT_EQ(fixed_size_memory_pool.GetNumberOfAllocatedChunks(), 0);
			KRATOS_EXPECT_EQ(fixed_size_memory_pool.ChunkSize(), default_chunk_size);
		}

		KRATOS_TEST_CASE_IN_SUITE(FixedSizeMemoryPoolAllocationDeallocation, KratosCoreFastSuite)
		{
			int max_threads = OpenMPUtils::GetCurrentNumberOfThreads();
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
				std::size_t number_of_chunks = i_repeat;// ((i_repeat / ParallelUtilities::GetNumThreads()) + 1) *  ParallelUtilities::GetNumThreads();
				KRATOS_EXPECT_EQ(fixed_size_memory_pool.GetNumberOfAllocatedChunks(), number_of_chunks) << " (number_of_chunks = " << number_of_chunks << ")";
				KRATOS_EXPECT_EQ(fixed_size_memory_pool.ChunkSize(), default_chunk_size);
			}
		}

		KRATOS_TEST_CASE_IN_SUITE(FixedSizeMemoryPoolStressTest, KratosCoreStressSuite)
		{
			GTEST_SKIP() << "This test is disabled" << std::endl;
			std::size_t block_size = 64;
			std::size_t default_chunk_size = 128 ;// 1M
			std::size_t number_of_blocks = (default_chunk_size ) / block_size;
			ThreadFixedSizeMemoryPool fixed_size_memory_pool(block_size, default_chunk_size, 0);

			auto repeat_number = 100;
			for (auto i_repeat = 3; i_repeat <= repeat_number; i_repeat++)
			{
				std::vector<void *> pointer_vector;
				for (auto i_chunk = 0; i_chunk < i_repeat; i_chunk++)
					for (std::size_t i_block = 0; i_block < number_of_blocks; i_block++)
					{
						void* p = fixed_size_memory_pool.Allocate();
						KRATOS_EXPECT_NE(p, nullptr) << " i_repeat : " << i_repeat << " , i_chunk : " << i_chunk << " , i_block : " << i_block;
						pointer_vector.push_back(p);

					}
				for (std::size_t i_pointer = 0; i_pointer < pointer_vector.size(); i_pointer++) {
					try{
					fixed_size_memory_pool.Deallocate(pointer_vector[i_pointer]);
					}
					catch(Exception& e){
						std::cout << i_pointer << " : " << pointer_vector[i_pointer] << std::endl;
						std::cout << e.what() << std::endl;
						throw e;
					}
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

namespace{
		class PoolBlock64byte : public PoolObject {
            public:
                PoolBlock64byte() {Block[0] = 0;};
            private:
                double Block[228];
		};

		class Block64byte {
            public:
                Block64byte() {Block[0] = 0;};
            private:
                double Block[228];
		};
}

		KRATOS_TEST_CASE_IN_SUITE(PoolObjectParallelNewDelete, KratosCoreFastSuite)
		{
			std::cout << MemoryPool::Info() << std::endl;
			auto repeat_number = 120;
#pragma omp parallel for
			for (auto i_repeat = 0; i_repeat < repeat_number; i_repeat++)
			{
				try {
					std::size_t size = 1024 * 128;
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

			// std::cout << MemoryPool::Info() << std::endl;
		}

		KRATOS_TEST_CASE_IN_SUITE(NewDeleteComparison, KratosCoreStressSuite)
		{
			GTEST_SKIP() << "This test is disabled" << std::endl;
			auto repeat_number = 256;
			Timer::Start("Pool");
// #pragma omp parallel for
			for (auto i_repeat = 1; i_repeat <= repeat_number; i_repeat++)
			{
				try {
					std::vector<PoolBlock64byte*> the_vector;

					Timer::Start("Pool new");
					for (auto i_chunk = 0; i_chunk < i_repeat; i_chunk++)
						for (auto i_block = 0; i_block < 4096; i_block++)
							the_vector.push_back(new PoolBlock64byte);
					Timer::Stop("Pool new");
					Timer::Start("Pool delete");
					for (std::size_t i = 0; i < the_vector.size(); i++) {
						delete the_vector[i];
					}
					Timer::Stop("Pool delete");
				}
				catch(Exception& e){
					std::cout << e << std::endl;
				}
				catch (...) {
					std::cout << "Unknown Error" << std::endl;
				}

			}
			Timer::Stop("Pool");

			Timer::Start("glib");
// #pragma omp parallel for
			for (auto i_repeat = 1; i_repeat <= repeat_number; i_repeat++)
			{
				std::vector<Block64byte*> the_vector;

				Timer::Start("glib new");
				for (auto i_chunk = 0; i_chunk < i_repeat; i_chunk++)
					for (auto i_block = 0; i_block < 4096; i_block++)
						the_vector.push_back(new Block64byte);

				Timer::Stop("glib new");
				Timer::Start("glib delete");
				for (std::size_t i = 0; i < the_vector.size(); i++) {
					delete the_vector[i];
				}
				Timer::Stop("glib delete");
			}
			Timer::Stop("glib");
			std::cout << Timer() << std::endl;

		}

		KRATOS_TEST_CASE_IN_SUITE(ChunkNewDeleteComparison, KratosCoreStressSuite)
		{
			GTEST_SKIP() << "This test is disabled" << std::endl;
			auto repeat_number = 128;
			std::size_t block_size_in_bytes = 230;
		    std::size_t chunk_size_in_bytes =  1024*1024;
			std::size_t number_of_objects_per_chunk = chunk_size_in_bytes / 8;

			Timer::Start("Chunk");
			Chunk chunk(block_size_in_bytes, chunk_size_in_bytes);
			for (auto i_repeat = 1; i_repeat <= repeat_number; i_repeat++)
			{
				chunk.Initialize();
				std::vector<void*> the_vector(number_of_objects_per_chunk);

				// Timer::Start("Chunk new");
					for (std::size_t i_block = 0; i_block < number_of_objects_per_chunk; i_block++)
						the_vector[i_block] = chunk.Allocate();
				// Timer::Stop("Chunk new");
				// Timer::Start("Chunk delete");
				for (std::size_t i = 0; i < the_vector.size(); i++) {
					chunk.Deallocate(the_vector[i]);
				}
				// Timer::Stop("Chunk delete");

			}
			Timer::Stop("Chunk");

			Timer::Start("glib");
			for (auto i_repeat = 1; i_repeat <= repeat_number; i_repeat++)
			{
				std::vector<Block64byte*> the_vector(number_of_objects_per_chunk);

				// Timer::Start("glib new");
				// for (auto i_chunk = 0; i_chunk < i_repeat; i_chunk++)
					for (std::size_t i_block = 0; i_block < number_of_objects_per_chunk; i_block++)
						the_vector[i_block] = new Block64byte;

				// Timer::Stop("glib new");
				// Timer::Start("glib delete");
				for (std::size_t i = 0; i < the_vector.size(); i++) {
					delete the_vector[i];
				}
				// Timer::Stop("glib delete");
			}
			Timer::Stop("glib");
			std::cout << Timer() << std::endl;

		}
	}
}  // namespace Kratos.

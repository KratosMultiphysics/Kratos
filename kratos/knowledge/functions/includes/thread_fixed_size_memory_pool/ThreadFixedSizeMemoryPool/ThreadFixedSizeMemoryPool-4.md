---
type: C++ Method
title: ThreadFixedSizeMemoryPool
resource: includes/thread_fixed_size_memory_pool.h#L73-L82
generated:
  by: okf-rs/0.7.0
---

# Signature

`ThreadFixedSizeMemoryPool(std::size_t BlockSizeInBytes, SizeType ChunkSize, std::size_t ThreadNumber) : mBlockSizeInBytes(BlockSizeInBytes) , mChunkSize(ChunkSize) , mThreadNumber(ThreadNumber) , mChunks() , mAvailableChunks() , mNumberOfReleasedChunks(0) , mpCurrentChunk(nullptr)`
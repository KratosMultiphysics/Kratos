from math import prod

class BufferPool:
    def __init__(self, sizes, xp, dtype):
        self.xp = xp
        self.dtype = dtype
        self.max_sizes = sizes
        self.buffers = [xp.empty(s, dtype=dtype) for s in sizes]
        self.in_use = [False] * len(sizes)
        # Store base buffer data pointers for address-based lookup
        self._buffer_pointers = []
        for buf in self.buffers:
            ptr = self._get_data_pointer(buf)
            self._buffer_pointers.append(ptr)

    def _get_data_pointer(self, arr):
        """Get the memory data pointer from a numpy or cupy array."""
        if hasattr(arr, '__array_interface__'):
            # NumPy uses __array_interface__['data'][0]
            return arr.__array_interface__['data'][0]
        elif hasattr(arr, '__cuda_array_interface__'):
            # CuPy uses __cuda_array_interface__['data'][0]
            return arr.__cuda_array_interface__['data'][0]
        else:
            raise TypeError(f"Unsupported array type: {type(arr)}")

    def _find_buffer_index_by_pointer(self, arr):
        """Find the buffer index by comparing memory addresses."""
        arr_ptr = self._get_data_pointer(arr)
        for i, buf_ptr in enumerate(self._buffer_pointers):
            # Check if the array pointer falls within this buffer's range
            buf_size = self.max_sizes[i] * arr.itemsize
            if buf_ptr <= arr_ptr < (buf_ptr + buf_size):
                return i
        raise ValueError(
            f"Array pointer {arr_ptr} does not match any buffer in the pool"
        )

    def Get(self, index, shape):
        if self.in_use[index]:
            raise RuntimeError(f"Buffer {index} already in use")

        shape = tuple(int(s) for s in shape)
        required = prod(shape)

        if required > self.max_sizes[index]:
            raise ValueError(
                f"Requested shape {shape} needs {required} elements, "
                f"but buffer {index} has max {self.max_sizes[index]}"
            )

        self.in_use[index] = True
        return self.buffers[index][:required].reshape(shape)

    def Release(self, index_or_array):
        """Release a buffer by index or by passing the array itself.
        
        When an array is passed, the buffer index is automatically detected
        by comparing memory addresses.
        """
        # Check if we received an array (has __array_interface__ or __cuda_array_interface__)
        has_array_interface = hasattr(index_or_array, '__array_interface__')
        has_cuda_array_interface = hasattr(index_or_array, '__cuda_array_interface__')
        
        if has_array_interface or has_cuda_array_interface:
            # Auto-detect index from memory address
            index = self._find_buffer_index_by_pointer(index_or_array)
        else:
            # Use provided index
            index = index_or_array

        if not self.in_use[index]:
            raise RuntimeError(f"Buffer {index} is not in use")
        self.in_use[index] = False

    def ReleaseByIndex(self, index):
        self.in_use[index] = False

    def Status(self):
        return self.in_use
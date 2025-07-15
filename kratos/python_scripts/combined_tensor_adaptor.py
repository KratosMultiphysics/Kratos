import numpy
from typing import Union
import KratosMultiphysics as Kratos

class CombinedTensorAdaptor:
    """A tensor adaptor which may be used to combine given any TensorAdaptor and/or another CombinedTensorAdaptor
       This CombinedTensorAdaptor behaves like a TensorAdaptor even though it is not derived from the TensorAdaptor
       base class.
    """
    def __init__(self, list_of_tensor_adaptors: 'list[Union[Kratos.TensorAdaptors.BoolTensorAdaptor, Kratos.TensorAdaptors.IntTensorAdaptor, Kratos.TensorAdaptors.DoubleTensorAdaptor, CombinedTensorAdaptor]]', axis = 0):
        """Construct a CombinedTensor adaptor which concatenates the given TensorAdaptors and/or CombinedTensorAdaptors
           along the specified axis.

           CombinedTensorAdaptor will have the maximum number of dimensions from all the number of dimensions provided by
           list_of_tensor_adaptors.

           If a provided TensorAdaptor / CombinedTensorAdaptor does not have the same number of dimensions as the previously mentioned maximum
           number of dimensions, new dimensions will be added to the end of the shape having value "1", which is called
           the modified tensor adaptor shape.

        Args:
            list_of_tensor_adaptors (list[Union[Kratos.TensorAdaptors.BoolTensorAdaptor, Kratos.TensorAdaptors.IntTensorAdaptor, Kratos.TensorAdaptors.DoubleTensorAdaptor, &#39;CombinedTensorAdaptor&#39;]]): _description_
            axis (int, optional): _description_. Defaults to 0.

        Raises:
            RuntimeError: If the list_of_tensor_adaptors is empty.
            RuntimeError: If the axis is equal or larger than the maximum number of dimensions found in the list_of_tensor_adaptors.
            RuntimeError: IF list_of_tensor_adaptors contains TensorAdaptors / CombinedTensorAdaptors having different numpy.dtypes.
            RuntimeError: If the modified tensor adaptor shape's number of components in each dimension does not match with other tensor adaptor's same except in the axis dimension.
        """
        self.__list_of_tensor_adaptors = list_of_tensor_adaptors
        self.__list_of_tensor_adaptor_shapes: list[Kratos.DenseVectorUnsignedInt] = []
        self.__axis = axis

        if len(self.__list_of_tensor_adaptors) == 0:
            raise RuntimeError("Cannot construct a combined tensor adaptor with an empty list of tensor adaptors.")

        # get the initial parameters
        shape = self.__list_of_tensor_adaptors[0].Shape()
        dtype = self.__list_of_tensor_adaptors[0].ViewData().dtype

        for current_tensor_adaptor in self.__list_of_tensor_adaptors:
            if shape.Size() < current_tensor_adaptor.Shape().Size():
                shape = current_tensor_adaptor.Shape()

        if self.__axis >= shape.Size():
            raise RuntimeError(f"The axis should be less than the number of dimensions [ axis = {self.__axis}, shape = {shape} ].\n")

        shape[self.__axis] = 0
        for current_tensor_adaptor in self.__list_of_tensor_adaptors:
            if current_tensor_adaptor.ViewData().dtype != dtype:
                raise RuntimeError(f"Only allowed to combine same type of data [ combined dtype =  {dtype.name}, current tensor adaptor dtype = {current_tensor_adaptor.ViewData().dtype.name}, current tensor adaptor = {current_tensor_adaptor} ].\n")

            # we add 1 to the end dimensions of the shape
            # so that they can be properly broadcasted using numpy
            current_ta_shape = Kratos.DenseVectorUnsignedInt(shape.Size(), 1)
            for i in range(current_tensor_adaptor.Shape().Size()):
                current_ta_shape[i] = current_tensor_adaptor.Shape()[i]

            # now check whether the components of the dimensions match
            for i in range(shape.Size()):
                if i != self.__axis and shape[i] != current_ta_shape[i]:
                    raise RuntimeError(f"Number of components in each dimension should match except in the axis dimension [ combined shape = {shape}, axis = {self.__axis}, current_tensor adaptor shape = {current_ta_shape}, tensor adaptor = {current_tensor_adaptor} ].\n")

            shape[self.__axis] += current_ta_shape[self.__axis]
            self.__list_of_tensor_adaptor_shapes.append(current_ta_shape)

        # now allocate the memory
        self.__data: numpy.ndarray = numpy.empty(shape, dtype=dtype)

    def GetTensorAdaptors(self):
        """Returns the list of tensor adaptors.
        """
        return list(self.__list_of_tensor_adaptors)

    def CollectData(self, recursively = True) -> None:
        if recursively:
            for current_tensor_adaptor in self.__list_of_tensor_adaptors:
                if isinstance(current_tensor_adaptor, CombinedTensorAdaptor):
                    current_tensor_adaptor.CollectData(recursively)
                else:
                    current_tensor_adaptor.CollectData()

        list_of_numpy_arrays = [self.__Reshape(ta.ViewData(), self.__list_of_tensor_adaptor_shapes[i]) for i, ta in enumerate(self.__list_of_tensor_adaptors)]

        # now do the numpy concatenation
        self.__data[:] = numpy.concatenate(list_of_numpy_arrays, axis = self.__axis)

    def StoreData(self, recursively = True) -> None:
        slicing = [slice(None)] * len(self.__data.shape)
        slicing[self.__axis] = slice(0, 0)
        for i, ta in enumerate(self.__list_of_tensor_adaptors):
            slicing[self.__axis] = slice(slicing[self.__axis].stop, slicing[self.__axis].stop + self.__list_of_tensor_adaptor_shapes[i][self.__axis])
            ta.SetData(self.__Reshape(self.__data[tuple(slicing)], ta.Shape()))

        if recursively:
            for current_tensor_adaptor in self.__list_of_tensor_adaptors:
                if isinstance(current_tensor_adaptor, CombinedTensorAdaptor):
                    current_tensor_adaptor.StoreData(recursively)
                else:
                    current_tensor_adaptor.StoreData()

    def ViewData(self) -> numpy.ndarray:
        return numpy.array(self.__data, copy=False)

    def SetData(self, values: numpy.ndarray) -> None:
        self.__data[:] = values

    def Shape(self) -> Kratos.DenseVectorUnsignedInt:
        return Kratos.DenseVectorUnsignedInt(list(self.__data.shape))

    def DataShape(self) -> Kratos.DenseVectorUnsignedInt:
        return Kratos.DenseVectorUnsignedInt(list(self.__data.shape[1:]))

    def Size(self) -> int:
        return self.__data.size

    def __Reshape(self, numpy_array: numpy.ndarray, new_shape: 'list[int]'):
        """This is a wrapper method for numpy.reshape

        The copy keywordarg is only introduced in numpy 2.1, hence all the other previous
        versions, we cannot use it. The copy = False make sure that the reshaping of a
        numpy array does not create a new copy of the underlying data. If it cannot reshape,
        then it throws an error. This is a good check to have, eventhough our TensorAdaptor
        will not require any copying for reshaping.
        """
        numpy_version_info = [int(v) for v in numpy.__version__.split(".")]
        if numpy_version_info[0] > 2 or (numpy_version_info[0] == 2 and numpy_version_info[1] >= 1):
            return numpy.reshape(numpy_array, new_shape, copy=False)
        else:
            return numpy.reshape(numpy_array, new_shape)

    @property
    def data(self):
        return self.ViewData()

    @data.setter
    def data(self, values: numpy.ndarray):
        self.__data[:] = values


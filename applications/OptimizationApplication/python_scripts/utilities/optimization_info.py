class OptimizationInfo:
    """Optimization info data container with buffered data ability

    Instances of this can hold (str, any) data pairs in hierachychal data
    structure where each sub item can have their own buffer sizes. Hence,
    Different sub structures can have their own buffers and can be advanced
    seperately or all together.

    The buffered data is stored in a cyclic buffer.
    """
    def __init__(self, buffer_size: int = 1):
        """Creates an instance of OptimizationInfo with specified buffer size

        This creates an instance of Optimization info with given buffer size.

        Args:
            buffer_size (int, optional): Cyclic buffer size. Defaults to 1.
        """
        self.__buffer_index = 0
        self.__buffered_data: 'list[dict[str, any]]' = []
        self.__sub_items: 'dict[str, OptimizationInfo]' = {}

        # sets the buffer
        self.SetBufferSize(buffer_size)


    def AdvanceStep(self, recursive = True) -> None:
        """Advances the buffered data containers

        Advances to the next index of the buffer in the cyclic buffer
        making the current index accesible as a past value indes (step_index = 1)

        Args:
            recursive (bool, optional): Advance step recursively in sub items. Defaults to True.
        """
        # first advance the current instance
        self.__buffer_index = (self.__buffer_index + 1) % self.GetBufferSize()

        if recursive:
            # now advance the sub items
            for sub_item in self.__sub_items.values():
                sub_item.AdvanceStep(recursive)

    def ClearStep(self, step_index: int = 0, recursive = True) -> None:
        """Clear the given step index data holders.

        This makes current instances (and sub items) given step_index data
        holders to be empty.

        Args:
            step_index (int, optional): Step index to be cleared. Defaults to 0.
            recursive (bool, optional): Clear sub_items same step_index. Defaults to True.
        """
        # first clear the current data
        self.__buffered_data[self.__GetBufferIndex(step_index)] = {}

        if recursive:
            # now clear the sub items
            for sub_item in self.__sub_items.values():
                sub_item.ClearStep(step_index, recursive)

    def SetBufferSize(self, buffer_size: int, recursive = False) -> None:
        """Sets the buffer size of the optimization info

        This sets buffer size of the cyclic buffer in this instance (and sub items)
        to specified value.

        This clears the whole buffer of this instance (and sub items).

        Args:
            buffer_size (int): Size of the cyclic buffer
            recursive (bool, optional): Set sub item cyclic buffers. Defaults to False.
        """
        self.__buffer_index = 0
        self.__buffered_data.clear()
        for _ in range(buffer_size):
            self.__buffered_data.append({})

        if recursive:
            for sub_item in self.__sub_items.values():
                sub_item.SetBufferSize(buffer_size, recursive)

    def GetBufferSize(self) -> int:
        """Returns the buffer size of current isntance of optimization info

        Returns:
            int: Buffer size of the cyclic buffer
        """
        return len(self.__buffered_data)

    def HasValue(self, key: str, step_index: int = 0) -> bool:
        """Checks whether key exists in the given step_index

        The key in here should always be relative to the instance it is called upon.
        Hence, this key can have "/" seperations for subitems, then this method
        will navigate through the sub items as well to check whether given key is available.

        Args:
            key (str): Relative path of the key to be checked for.
            step_index (int, optional): Step index to be looked for. Defaults to 0.

        Returns:
            bool: True if the key is found in given step index. False otherwise.
        """
        pos = key.find("/")

        if pos == -1:
            # this is a leaf key, then look for the value
            # first in the value container
            if step_index < self.GetBufferSize() and step_index >= 0:
                if key in self.__buffered_data[self.__GetBufferIndex(step_index)].keys():
                    return True
            # now if it was not in the buffered data, then check it in sub items.
            return key in self.__sub_items.keys()
        else:
            # not the leaf key. hence check for the sub items
            current_key = key[:pos]
            if current_key in self.__sub_items.keys():
                # found the subkey. Call HasValue recursively
                return self.__sub_items[current_key].HasValue(key[pos+1:], step_index)
            else:
                # not found the subkey. Hence the given key is not available.
                return False

    def GetValue(self, key: str, step_index: int = 0) -> any:
        """Get the value given by the key at the specified step_index.

        This method retrieves the value given by the key at the specified step_index.
        The key must be the relative path w.r.t. current instance of the OptimizationInfo.
        It can include "/" seperators to get a value which is in sub items. In this case,
        it will be retrieved by recursive calls.

        Args:
            key (str): Relative path to the value.
            step_index (int, optional): Step index the value should be retrieved. Defaults to 0.

        Raises:
            RuntimeError: If the given key is not found.
            RuntimeError: If parent keys are not found.

        Returns:
            any: Value stored at the key for the specified step_index.
        """
        pos = key.find("/")

        if pos == -1:
            # this is a leaf key, then look for the value
            if not self.HasValue(key, step_index):
                raise RuntimeError(f"The key \"{key}\" not found in the optimization info. OptimizationInfo:\n{self}")

            # first in the value container
            if step_index < self.GetBufferSize() and step_index >= 0:
                buffer_data = self.__buffered_data[self.__GetBufferIndex(step_index)]
                if key in buffer_data.keys():
                    return buffer_data[key]

            # now if it was not in the buffered data, then it is in sub items.
            return self.__sub_items[key]
        else:
            current_key = key[:pos]
            if current_key in self.__sub_items.keys():
                return self.__sub_items[current_key].GetValue(key[pos+1:], step_index)
            else:
                raise RuntimeError(f"The key \"{current_key}\" not found in the optimization info which is a parent key for \"{key}\". OptimizationInfo:\n{self}")

    def SetValue(self, key: str, value: any, step_index: int = 0, overwrite = False) -> None:
        """Sets a value for specified key at specified step_index.

        This method sets the value at the key at the specified step_index.
        The key must be the relative path w.r.t. current instance of the OptimizationInfo.
        It can include "/" seperators to set a value which is in sub items. In this case,
        it will be set by recursive calls.

        Overwrite allows to overwrite existing keys with new values.

        Args:
            key (str): Relative path to the value.
            value (any): value to be set
            step_index (int, optional): Step index the value should be set. Defaults to 0.
            overwrite (bool, optional): Allows overwriting of existing values. Defaults to False.

        Raises:
            RuntimeError: If overwrite is False, and a value is found for the same key.
        """
        pos = key.find("/")
        if pos == -1:
            if not overwrite and self.HasValue(key, step_index):
                raise RuntimeError(f"Trying to overwrite \"{key}\" with value {value} [ Existing value = {self.GetValue(key)}]. OptimizationInfo:\n{self}")

            # this is a leaf, then we check the type of the value and add appropriately
            if isinstance(value, dict):
                # if the given value is a dict, then convert the structure
                # to OptimizationInfo while keeping the sub_item structure.
                sub_item = OptimizationInfo(self.GetBufferSize())
                self.SetValue(key, sub_item)
                for sub_key, sub_value in value.items():
                    sub_item.SetValue(f"{sub_key}", sub_value, step_index, overwrite)
            elif isinstance(value, OptimizationInfo):
                # if the given value is of type OptimizationInfo, then put it to sub_items.
                self.__sub_items[key] = value
            else:
                # if not any of the above, it is a normal value. Then put it to buffer.
                buffer_data = self.__buffered_data[self.__GetBufferIndex(step_index)]
                buffer_data[key] = value
        else:
            # if it is not a leaf key. then look for the sub_item, and recursively
            # call the Set method.
            current_key = key[:pos]
            if not current_key in self.__sub_items.keys():
                # no sub key found. then create it.
                self.__sub_items[current_key] = OptimizationInfo(self.GetBufferSize())

            self.__sub_items[current_key].SetValue(key[pos+1:], value, step_index, overwrite)

    def RemoveValue(self, key: str, step_index: int = 0) -> None:
        """Remove value at the key in the specified step_index

        This method removes value at the key specified at the step_index.
        The key must be the relative path w.r.t. current instance of the OptimizationInfo.
        It can include "/" seperators to remove a value which is in sub items. In this case,
        it will be removed by recursive calls.

        Args:
            key (str): Relative path to the value.
            step_index (int, optional): Step index the value should be removed from. Defaults to 0.

        Raises:
            RuntimeError: If the given key is not found.
            RuntimeError: If the given parents of the key is not found.
        """
        pos = key.find("/")
        if pos == -1:
            # found a leaf key.
            is_reomved = False

            if step_index < self.GetBufferSize() and step_index >= 0:
                buffer_data = self.__buffered_data[self.__GetBufferIndex(step_index)]
                # check whether the key is available in the given step_index
                if key in buffer_data.keys():
                    is_reomved = True
                    del buffer_data[key]

            if not is_reomved:
                # if it is not removed from the buffered data, then check whether it is available
                # in the sub_items
                if key in self.__sub_items.keys():
                    is_reomved = True
                    del self.__sub_items[key]

            if not is_reomved:
                raise RuntimeError(f"\"{key}\" is not found. OptimizationInfo:\n{self}")
        else:
            # it is not a leaf key
            current_key = key[:pos]
            if current_key in self.__sub_items.keys():
                # call recursively the sub_items remove value
                self.__sub_items[current_key].RemoveValue(key[pos+1:], step_index)
            else:
                raise RuntimeError(f"\"{key}\" is not found. OptimizationInfo:\n{self}")

    def PrintData(self, step_index: int = -1, tabbing = "") -> str:
        """Prints containing data in a json like structure.

        This method prints containing data and sub_item data recursively.

        If the step_index == -1 then, it prints all the values on all the steps of the buffer
        in current and all the sub_items recursively.

        If the step_index > 0 then, it will only print values at that specified buffer index for
        current instance as well as all the sub_items.

        Args:
            step_index (int, optional): Step index to be printed. Defaults to -1.
            tabbing (str, optional): Initial tabbing for json output.. Defaults to "".

        Returns:
            str: String will all the containing values formatted like json.
        """
        info = self.__Info(step_index, tabbing)
        info = f"{tabbing}{info}"
        return info

    def __GetBufferIndex(self, step_index: int) -> int:
        if step_index >= self.GetBufferSize() or step_index < 0:
            raise RuntimeError(f"Invalid step index. Allowed step indices are in [0, {self.GetBufferSize()}) [ StepIndex = {step_index} ]. OptimizationInfo:\n{self}")

        return (self.__buffer_index - step_index) % self.GetBufferSize()

    def __Info(self, step_index: int = -1, tabbing = "") -> str:
        info = "{"

        if step_index == -1:
            for i in range(self.GetBufferSize()):
                info += f"\n{tabbing}\t--- Step = {i:2d} ---"
                for key, value in self.__buffered_data[self.__GetBufferIndex(i)].items():
                    info += f"\n{tabbing}\t\"{key}\": {value}"
        else:
            for key, value in self.__buffered_data[self.__GetBufferIndex(step_index)].items():
                info += f"\n{tabbing}\t\"{key}\": {value}"

        if step_index == -1:
            info += f"\n{tabbing}\t--- Sub items ---"

        next_tabbing = f"{tabbing}\t"
        for sub_key, sub_value in self.__sub_items.items():
            info += f"\n{tabbing}\t\"{sub_key}\": {sub_value.__Info(step_index, next_tabbing)}"

        info += f"\n{tabbing}" + "}"
        return info

    def __getitem__(self, key: str) -> any:
        if isinstance(key, str):
            return self.GetValue(key)
        elif isinstance(key, tuple) and len(key) == 2:
            return self.GetValue(key[0], key[1])
        else:
            raise RuntimeError(f"The key should be either a string (representing the key path) or (string, int) tuple (representing key path and step index) [ key = {key}].")

    def __setitem__(self, key: str, value: any) -> None:
        if isinstance(key, str):
            return self.SetValue(key, value)
        elif isinstance(key, tuple) and len(key) == 2:
            return self.SetValue(key[0], value, key[1])
        else:
            raise RuntimeError(f"The key should be either a string (representing the key path) or (string, int) tuple (representing key path and step index) [ key = {key}].")

    def __delitem__(self, key: str) -> None:
        if isinstance(key, str):
            return self.RemoveValue(key)
        elif isinstance(key, tuple) and len(key) == 2:
            return self.RemoveValue(key[0], key[1])
        else:
            raise RuntimeError(f"The key should be either a string (representing the key path) or (string, int) tuple (representing key path and step index) [ key = {key}].")

    def __str__(self) -> str:
        return self.__Info()


from __future__ import annotations
from typing import Any

class BufferedDict:
    """Buffered dict container with buffered data ability

    Instances of this can hold (str, Any) data pairs in hierachychal data
    structure where each sub item can have their own buffer sizes. Hence,
    Different sub structures can have their own buffers and can be advanced
    separately or all together.

    Overwriting of the existing data is not allowed because, this is used
    as a central place to store data from various objects, hence allowing
    overwriting will make difficult to track the changes.

    Buffer size setting can only be done at the construction time
    of the instance to not to allow re sizing because then it is not
    clear how to handle the connected steps.

    The buffered data is stored in a cyclic buffer.
    """
    def __init__(self, buffer_size: int = 1, clear_buffer_when_advancing: bool = True):
        """Creates an instance of BufferedDict with specified buffer size

        This creates an instance of buffered data with given buffer size.

        Args:
            buffer_size (int, optional): Cyclic buffer size. Defaults to 1.
            clear_buffer_when_advancing (bool): Clears the buffer and sub item buffers when advancing.
        """
        self.__parent = None
        self.__buffer_index = 0
        self.__buffered_data: 'list[dict[str, Any]]' = []
        self.__sub_items: 'dict[str, BufferedDict]' = {}
        self.__clear_buffer_when_advancing = clear_buffer_when_advancing

        # sets the buffer
        self.__SetBufferSize(buffer_size)

    def GetBufferSize(self) -> int:
        """Returns the buffer size of current instance of buffered data

        Returns:
            int: Buffer size of the cyclic buffer
        """
        return len(self.__buffered_data)

    def AdvanceStep(self, recursive = True) -> None:
        """Advances the buffered data containers

        Advances to the next index of the buffer in the cyclic buffer
        making the current index accessible as a past value index (step_index = 1)

        Clears the values in the current step after advancement.

        Args:
            recursive (bool, optional): Advance step recursively in sub items. Defaults to True.
        """
        # first advance the current instance
        self.__buffer_index = (self.__buffer_index + 1) % self.GetBufferSize()
        if self.__clear_buffer_when_advancing:
            self.__buffered_data[self.__GetBufferIndex(0)] = {}

        if recursive:
            # now advance the sub items
            for sub_item in self.__sub_items.values():
                sub_item.AdvanceStep(recursive)

    def HasValue(self, key: str, step_index: int = 0) -> bool:
        """Checks whether key exists in the given step_index

        The key in here should always be relative to the instance it is called upon.
        Hence, this key can have "/" separations for subitems, then this method
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

    def GetValue(self, key: str, step_index: int = 0) -> Any:
        """Get the value given by the key at the specified step_index.

        This method retrieves the value given by the key at the specified step_index.
        The key must be the relative path w.r.t. current instance of the BufferedDict.
        It can include "/" separators to get a value which is in sub items. In this case,
        it will be retrieved by recursive calls.

        Args:
            key (str): Relative path to the value.
            step_index (int, optional): Step index the value should be retrieved. Defaults to 0.

        Raises:
            RuntimeError: If the given key is not found.
            RuntimeError: If parent keys are not found.

        Returns:
            Any: Value stored at the key for the specified step_index.
        """
        pos = key.find("/")

        if pos == -1:
            # this is a leaf key, then look for the value
            if not self.HasValue(key, step_index):
                raise RuntimeError(f"The key \"{key}\" not found in the buffered data [ step_index = {step_index} ]. BufferedDict:\n{self}")

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
                raise RuntimeError(f"The key \"{current_key}\" not found in the buffered data which is a parent key for \"{key}\". BufferedDict:\n{self}")

    def SetValue(self, key: str, value: Any, step_index: int = 0, overwrite: bool = False) -> None:
        """Sets a value for specified key at specified step_index.

        This method sets the value at the key at the specified step_index.
        The key must be the relative path w.r.t. current instance of the BufferedDict.
        It can include "/" separators to set a value which is in sub items. In this case,
        it will be set by recursive calls.

        Args:
            key (str): Relative path to the value.
            value (Any): value to be set
            step_index (int, optional): Step index the value should be set. Defaults to 0.
            overwrite (bool): Overwrites if existing key is found with the new value, otherwise creates a key with the value.

        Raises:
            RuntimeError: If a value is overwritten.
        """
        pos = key.find("/")
        if pos == -1:
            # this is a leaf, then we check the type of the value and add appropriately
            if isinstance(value, dict):
                # if the given value is a dict, then convert the structure
                # to BufferedDict while keeping the sub_item structure.
                sub_item = BufferedDict(self.GetBufferSize(), self.__clear_buffer_when_advancing)
                self.__AddSubItem(key, sub_item)

                # now iterate through all the sub_keys and values of the dictionary and
                # add them to sub_item
                for sub_key, sub_value in value.items():
                    sub_item.SetValue(f"{sub_key}", sub_value, step_index)
            elif isinstance(value, BufferedDict):
                # if the given value is of type BufferedDict, then put it to sub_items.
                self.__AddSubItem(key, value)
            else:
                # if not any of the above, it is a normal value. Then put it to buffer.
                self.__AddBufferedValue(key, value, step_index, overwrite)
        else:
            # if it is not a leaf key. then look for the sub_item, and recursively
            # call the Set method.
            current_key = key[:pos]
            if not current_key in self.__sub_items.keys():
                # no existing key found then create it.
                self.__AddSubItem(current_key, BufferedDict(self.GetBufferSize(), self.__clear_buffer_when_advancing))

            self.__sub_items[current_key].SetValue(key[pos+1:], value, step_index, overwrite)

    def RemoveValue(self, key: str, step_index: int = 0) -> None:
        """Remove value at the key in the specified step_index

        This method removes value at the key specified at the step_index.
        The key must be the relative path w.r.t. current instance of the BufferedDict.
        It can include "/" separators to remove a value which is in sub items. In this case,
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
                raise RuntimeError(f"\"{key}\" is not found. BufferedDict:\n{self}")
        else:
            # it is not a leaf key
            current_key = key[:pos]
            if current_key in self.__sub_items.keys():
                # call recursively the sub_items remove value
                self.__sub_items[current_key].RemoveValue(key[pos+1:], step_index)
            else:
                raise RuntimeError(f"\"{key}\" is not found. BufferedDict:\n{self}")

    def GetValueItems(self, step_index: int = 0) -> 'dict[str, Any]':
        return self.__buffered_data[self.__GetBufferIndex(step_index)]

    def GetSubItems(self) -> 'dict[str, BufferedDict]':
        return self.__sub_items

    def GetMap(self, step_index: int = 0) -> 'dict[str, Any]':
        """Generates (str, Any) pair map with recursively collecting all items in the BufferedDict

        This method recursively collects all items in the buffered data and returns a dict
        with (str, Any) pairs. All the sub_items are also added as (str, BufferedDict) pairs.

        Args:
            step_index (int, optional):  Step index the value should be looked at. Defaults to 0.

        Returns:
            dict[str, Any]: (str, Any) map of all the items.
        """
        key_value_pair_map = {}
        self.__AddKeyValuePairsToMap("", key_value_pair_map, step_index)
        return key_value_pair_map

    def GetParent(self) -> BufferedDict:
        """Get the parent of the current buffered data

        Returns:
            BufferedDict: Returns parent of the current instance.
        """
        return self.__parent

    def GetRoot(self) -> BufferedDict:
        """Get the root parent of the current buffered data

        Get the root parent of the current buffered data by recursive calls.

        Returns:
            BufferedDict: Root parent of the current buffered data.
        """
        if self.GetParent() is not None:
            return self.GetParent().GetRoot()
        else:
            return self

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

    def __AddKeyValuePairsToMap(self, current_path: str, key_value_pair_map: 'dict[str, Any]', step_index: int) -> None:
        # first add the buffered data
        if step_index >= 0 and step_index < self.GetBufferSize():
            buffer_data = self.__buffered_data[self.__GetBufferIndex(step_index)]
            for k, v in buffer_data.items():
                key_value_pair_map[f"{current_path}{k}"] = v

        # now add the sub item data recursively
        for k, v in self.__sub_items.items():
            key_value_pair_map[f"{current_path}{k}"] = v
            v.__AddKeyValuePairsToMap(f"{current_path}{k}/", key_value_pair_map, step_index)

    def __AddBufferedValue(self, key: str, value: Any, step_index: int, overwrite: bool) -> None:
        # first check whether this key exists in the specified step_index
        buffer_data = self.__buffered_data[self.__GetBufferIndex(step_index)]
        if not overwrite and key in buffer_data.keys():
            raise RuntimeError(f"Trying to add a buffer value with key = \"{key}\" when already value exists for the key [ Existing value = {buffer_data[key]}]. BufferedDict:\n{self}")

        # now check whether a sub item exists
        if key in self.__sub_items.keys():
            raise RuntimeError(f"Trying to add a buffer value with key = \"{key}\" when a subitem exists with the same key. BufferedDict:\n{self}")

        # now add the buffer value
        buffer_data[key] = value

    def __AddSubItem(self, key: str, value: BufferedDict) -> None:
        # first check if any of the buffered data has the same key.
        # this is because, sub_items are valid for all step_indices.
        if any([key in buffer_data.keys() for buffer_data in self.__buffered_data]):
            raise RuntimeError(f"Trying to add a new a sub item with key = \"{key}\" when a value with the same key exists in buffered data. BufferedDict:\n{self}")

        # now check if the item exists in the sub_items
        if key in self.__sub_items.keys():
            raise RuntimeError(f"Trying to add a new sub_item with key = \"{key}\" when already a sub item exists. BufferedDict:\n{self}")

        # now add the sub_item
        value.__parent = self
        self.__sub_items[key] = value

    def __SetBufferSize(self, buffer_size: int) -> None:
        self.__buffer_index = 0
        self.__buffered_data.clear()
        for _ in range(buffer_size):
            self.__buffered_data.append({})

    def __GetBufferIndex(self, step_index: int) -> int:
        if step_index >= self.GetBufferSize() or step_index < 0:
            raise RuntimeError(f"Invalid step index. Allowed step indices are in [0, {self.GetBufferSize()}) [ StepIndex = {step_index} ]. BufferedDict:\n{self}")

        return (self.__buffer_index - step_index) % self.GetBufferSize()

    def __Info(self, step_index: int, tabbing: str) -> str:
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

    def __getitem__(self, key: str) -> Any:
        if isinstance(key, str):
            return self.GetValue(key)
        elif isinstance(key, tuple) and len(key) == 2:
            return self.GetValue(key[0], key[1])
        else:
            raise RuntimeError(f"The key should be either a string (representing the key path) or (string, int) tuple (representing key path and step index) [ key = {key}].")

    def __setitem__(self, key: str, value: Any) -> None:
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
        return self.__Info(-1, "")

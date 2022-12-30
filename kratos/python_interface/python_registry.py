from enum import Enum

import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import GetListOfAvailableApplications

class RegistryContext(Enum):
    ALL = 0
    CPP = 1
    PYTHON = 2

class PythonRegistryIterator:
    ''' PythonRegistry iterator class
    '''
    def __init__(self, ThisPythonRegistry):
       # Member variables to keep track of iteration
       self.__index = 0
       self.__cpp_root_keys = ThisPythonRegistry._PythonRegistry__cpp_registry.keys()
       self.__py_root_keys = list(ThisPythonRegistry._PythonRegistry__python_registry.keys())

    def __next__(self):
        '''Returns the next object from the PythonRegistry wrapper root level
        '''
        n_py_items = len(self.__py_root_keys)
        n_total_items = len(self.__cpp_root_keys) + n_py_items
        if self.__index < n_total_items:
            # Check if we are at the Python or c++ registry
            if self.__index < n_py_items:
                # Return the item key from the Python registry
                item_key = self.__py_root_keys[self.__index]
            else:
                # Return the item key from the c++ registry
                item_key = self.__cpp_root_keys[self.__index-n_py_items]
            self.__index += 1
            return item_key
        # End of iteration
        raise StopIteration

class PythonRegistry(object):

    #TODO: Change this by a frozenmap() once this is available
    CppGetFunctionNamesMap = {
        "Operations": "GetOperation",
        "Processes": "GetProcess"
    }

    def __init__(self):
        # We get a "private" pointer to the c++ registry
        self.__cpp_registry = KratosMultiphysics.CppRegistry

        # Dictionary to act as Python registry
        self.__python_registry = {}

    def __getitem__(self, Name):
        if self.HasItem(Name, RegistryContext.PYTHON):
            # Get the class from the Python registry
            return self.__GetPythonItem(Name)
        elif self.HasItem(Name, RegistryContext.CPP):
            # Get the class from the c++ registry
            return self.__GetCppItem(Name)
        else:
            err_msg = f"Asking to retrieve '{Name}' non-registered item."
            raise Exception(err_msg)

    def __iter__(self):
       ''' Returns the Iterator object '''
       return PythonRegistryIterator(self)

    def keys(self, Name):
        # Check if Name registry level is registered and if it is a value
        if self.HasItem(Name):
            if not self.HasItems(Name):
                err_mgs = f"Asking for the keys of '{Name}'. '{Name}' item has no subitems."
                raise Exception(err_mgs)
        else:
            err_msg = f"Asking for the keys of '{Name}' non-registered item."
            raise Exception(err_msg)

        # Get the c++ registry item keys
        cpp_keys = None
        if self.HasItems(Name, RegistryContext.CPP):
            cpp_keys = self.__GetCppItem(Name).keys()
        # Get the Python registry item keys
        # Note that this is a view object of the Python dictionary keys
        py_keys = None
        if self.HasItems(Name, RegistryContext.PYTHON):
            py_keys = self.__GetPythonItem(Name).keys()

        # Return the c++ and Python registry keys in a tuple
        # Note that we return None if the item was not present in one of the registries
        return py_keys, cpp_keys

    def values(self, Name):
        # Check if Name registry level is registered and if it is a value
        if self.HasItem(Name):
            if not self.HasItems(Name):
                err_mgs = f"Asking for the values of '{Name}'. '{Name}' item has no subitems."
                raise Exception(err_mgs)
        else:
            err_msg = f"Asking for the values of '{Name}' non-registered item."
            raise Exception(err_msg)

        # Get the c++ registry item values
        cpp_values = None
        if self.HasItems(Name, RegistryContext.CPP):
            cpp_values = self.__GetCppItem(Name).values()
        # Get the Python registry item values
        # Note that this is a view object of the Python dictionary values
        py_values = None
        if self.HasItems(Name, RegistryContext.PYTHON):
            py_values = self.__GetPythonItem(Name).values()

        # Return the c++ and Python registry values in a tuple
        # Note that we return None if the item was not present in one of the registries
        return py_values, cpp_values

    def items(self, Name):
        # Check if Name registry level is registered and if it is a value
        if self.HasItem(Name):
            if not self.HasItems(Name):
                err_mgs = f"Asking for the items of '{Name}'. '{Name}' item has no subitems."
                raise Exception(err_mgs)
        else:
            err_msg = f"Asking for the items of '{Name}' non-registered item."
            raise Exception(err_msg)

        # Get the keys and values
        py_keys, cpp_keys = self.keys(Name)
        py_values, cpp_values = self.values(Name)

        # Loop the keys and values to return the items
        py_items = None
        if py_keys is not None:
            py_items = []
            for py_key, py_value in zip(py_keys, py_values):
                py_items.append((py_key, py_value))

        cpp_items = None
        if cpp_keys is not None:
            cpp_items = []
            for cpp_key, cpp_value in zip(cpp_keys, cpp_values):
                cpp_items.append((cpp_key, cpp_value))

        # Return the c++ and Python registry items in a tuple
        return py_items, cpp_items

    def NumberOfItems(self, Name, Context=RegistryContext.ALL):
        '''Returns the total number of items in the registry
        '''
        py_keys, cpp_keys = self.keys(Name)
        total_items = 0
        if (Context == RegistryContext.ALL or Context == RegistryContext.PYTHON):
            if py_keys is not None:
                total_items += len(py_keys)
        #TODO: Implement this in c++
        if (Context == RegistryContext.ALL or Context == RegistryContext.CPP):
            if cpp_keys is not None:
                total_items += len(cpp_keys)
        return total_items

    def HasItem(self, Name, Context=RegistryContext.ALL):
        if Context == RegistryContext.ALL:
            is_registered_in_cpp = self.__cpp_registry.HasItem(Name)
            is_registered_in_python = self.__HasPythonItem(Name)
            return is_registered_in_cpp or is_registered_in_python
        elif Context == RegistryContext.CPP:
            return self.__cpp_registry.HasItem(Name)
        elif Context == RegistryContext.PYTHON:
            return self.__HasPythonItem(Name)
        else:
            raise Exception("Wrong registry context. Accepted options are 'ALL', 'CPP' and 'PYTHON'.")

    def HasValue(self, Name):
        '''Checks if Name registry item holds a value
        '''
        has_value = False
        if self.HasItem(Name, RegistryContext.PYTHON):
            has_value = not isinstance(self.__GetPythonItem(Name), dict)
        elif self.HasItem(Name, RegistryContext.CPP):
            has_value = self.__cpp_registry.HasValue(Name)
        else:
            err_msg = f"Asking 'HasValue' for the '{Name}' non-registered item."
            raise Exception(err_msg)

        return has_value

    def HasItems(self, Name, Context=RegistryContext.ALL):
        '''Checks if Name registry item saves a value
        '''
        # Check if current item has subitems in Python registry
        has_subitems_py = None
        has_item_py = self.HasItem(Name, RegistryContext.PYTHON)
        if has_item_py:
            has_subitems_py = isinstance(self.__GetPythonItem(Name), dict)

        # Check if current item has subitems in c++ registry
        has_subitems_cpp = None
        has_item_cpp = self.HasItem(Name, RegistryContext.CPP)
        if has_item_cpp:
            has_subitems_cpp = self.__cpp_registry.HasItems(Name)

        # First check that current item is registered
        if not (has_item_py or has_item_cpp):
            err_msg = f"Asking 'HasItems' for the '{Name}' non-registered item."
            raise Exception(err_msg)

        # Return the corresponding result
        if Context == RegistryContext.ALL:
            return has_subitems_py or has_subitems_cpp
        elif Context == RegistryContext.PYTHON:
            return has_subitems_py
        elif Context == RegistryContext.CPP:
            return has_subitems_cpp
        else:
            err_msg = f"Wrong registry context '{Context}'"
            raise Exception(err_msg)

    def AddItem(self, Name, ModuleName):
        # Add current item
        # First check if it is already registered in the c++ and Python registries
        # If not registered, add it to the Python registry (note that we do an "All" registry as well as the module one)
        if self.HasItem(Name, RegistryContext.ALL):
            err_msg = f"Trying to register '{Name}' but it is already registered."
            raise Exception(err_msg)
        else:
            # Add to the corresponding item All block
            self.__AddItemToAll(Name, ModuleName)
            # Add item to the corresponding module
            self.__AddItem(Name, ModuleName)

    def RemoveItem(self, Name):
        if self.HasItem(Name, RegistryContext.PYTHON):
            self.__RemoveItem(Name)
        elif self.HasItem(Name, RegistryContext.CPP):
            #TODO: Discuss about removing stuff from the c++ registry. Do we allow it or not?
            # self.__cpp_registry.RemoveItem(Name)
            err_msg = f"Trying to remove '{Name}' from c++ registry. This operation is forbidden."
            raise Exception(err_msg)
        else:
            err_msg = f"Asking to remove '{Name}' non-registered item."
            raise Exception(err_msg)

    @staticmethod
    def __GetCppRegistryGetFunctionNameFromItemKeyword(ItemKeyword):
        if ItemKeyword in PythonRegistry.CppGetFunctionNamesMap:
            return PythonRegistry.CppGetFunctionNamesMap[ItemKeyword]
        else:
            err_msg = f"'{ItemKeyword}' has not been addded to 'CppGetFunctionNamesMap'."
            raise Exception(err_msg)

    def __HasPythonItem(self, Name):
        split_name = Name.split('.')
        aux_dict = self.__python_registry
        for i_name in split_name[:-1]:
            if not i_name in aux_dict:
                return False
            aux_dict = aux_dict[i_name]
        return split_name[-1] in aux_dict

    def __GetPythonItem(self, Name):
        split_name = Name.split('.')
        aux_dict = self.__python_registry
        for i_name in split_name[:-1]:
            aux_dict = aux_dict[i_name]
        value = aux_dict[split_name[-1]]
        return value

    def __GetCppItem(self, Name):
        # Check if current item has value
        if self.__cpp_registry.HasValue(Name):
            # Return the object from the c++ registry by filtering its type keyword
            item_keyword = Name.split('.')[0]
            get_function_name = self.__GetCppRegistryGetFunctionNameFromItemKeyword(item_keyword)
            if hasattr(self.__cpp_registry, get_function_name):
                return getattr(self.__cpp_registry, get_function_name)(Name)
            else:
                err_msg = f"Retrieving '{item_keyword}' items from c++ registry is not supported."
                raise Exception(err_msg)
        else:
            # Return current subitem
            return self.__cpp_registry.GetItem(Name)

    def __RemoveItem(self, Name):
        split_name = Name.split('.')
        aux_dict = self.__python_registry
        for i_name in split_name[:-1]:
            aux_dict = aux_dict[i_name]
        aux_dict.pop(split_name[-1])

    def __AddItem(self, Name, Value):
        split_name = Name.split('.')
        aux_dict = self.__python_registry
        for i_name in split_name[:-1]:
            if not i_name in aux_dict:
                aux_dict[i_name] = {}
            aux_dict = aux_dict[i_name]
        aux_dict[split_name[-1]] = Value

    def __AddItemToAll(self, Name, Value):
        split_name = Name.split('.')
        item_keyword = split_name[0]
        class_name = split_name[-1]
        all_full_name = f"{item_keyword}.All.{class_name}"
        if self.HasItem(all_full_name, RegistryContext.ALL):
            err_msg = f"Trying to register '{Name}' but there is already an item with the same '{class_name}' name in the '{item_keyword}.All' block."
            raise Exception(err_msg)
        self.__AddItem(all_full_name, Value)

def RegisterPrototype(RegistryPointName: str):
    '''Auxiliary function to be used as a class decorator
    This function is intended to be used as a class decorator in order
    to manually create a registry entry for the class of interest.
    The decorator has only one argument, which is the registry point string.
    The registry point is expected to have the item keyword in front followed
    by the application module or by 'UserDefined' for those modules placed in
    the current work directory. Some examples are:
    - Processes.KratosMultiphysics
    - Processes.KratosMultiphysics.FluidDynamicsApplication
    - Processes.UserDefined
    '''

    def register_wrapper(Class):
        # Get the list of compiled applications
        # This will be used for checking the available submodules
        available_apps = GetListOfAvailableApplications()

        # Check input registry point name
        split_name = RegistryPointName.split('.')
        if len(split_name) < 2:
            err_msg = f"Wrong provided item name '{RegistryPointName}' structure. A structure of the type 'ItemKeyWord.Module.Submodule' is expected."
            raise Exception(err_msg)

        # Check input registry module
        module_keys = split_name[1:]
        available_root_keys = ["KratosMultiphysics","UserDefined"]
        if module_keys[0] not in available_root_keys:
            err_msg = f"Wrong root module '{module_keys[0]}'. This is expected to be either 'KratosMultiphysics' or 'UserDefined'."
            raise Exception(err_msg)

        # Manual application-like registry
        # Check that the corresponding application is compiled
        if module_keys[0] == "KratosMultiphysics":
            if len(module_keys) != 1:
                sub_module_key = module_keys[1]
                if sub_module_key not in available_apps:
                    err_msg = f"Wrong submodule '{sub_module_key}'. Compile the corresponding application."
                    raise Exception(err_msg)

        # Call the Kratos registry to register the current item
        # Note that the item class name is used as
        full_name = f"{RegistryPointName}.{Class.__name__ }"
        prototype_data = {
            "Prototype" : Class
        }
        KratosMultiphysics.Registry.AddItem(full_name, prototype_data)

        return Class
    return register_wrapper

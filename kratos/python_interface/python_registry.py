import KratosMultiphysics

class RegistryContext():
    ALL = 0
    CPP = 1
    PYTHON = 2

class PythonRegistry(object):

    CppGetFunctionNamesMap = {
        "Operations": "GetOperation",
        "Processes": "GetProcess"
    }

    def __init__(self):
        # We get a "private" pointer to the c++ registry
        self.__cpp_registry = KratosMultiphysics.CppRegistry

        # Dictionary to act as Python registry
        self.__python_registry = {}

    def HasItem(self, Name, Context=RegistryContext.ALL):
        if Context == RegistryContext.ALL:
            is_registered_in_cpp = self.__cpp_registry.HasItem(Name)
            is_registered_in_python = self.__InternalHasItem(Name)
            return is_registered_in_cpp or is_registered_in_python
        elif Context == RegistryContext.CPP:
            return self.__cpp_registry.HasItem(Name)
        elif Context == RegistryContext.PYTHON:
            return self.__InternalHasItem(Name)
        else:
            raise Exception("Wrong registry context. Accepted options are 'ALL', 'CPP' and 'PYTHON'.")

    def AddItem(self, Name, Class):
        # Add current item
        # First check if it is already registered in the c++ and Python registries
        # If not registered, add it to the Python registry (note that we do an "All" registry as well as the module one)
        if self.HasItem(Name, RegistryContext.ALL):
            err_msg = f"Trying to register '{Name}' but it is already registered."
            raise Exception(err_msg)
        else:
            # Add to the corresponding item All block
            self.__InternalAddItemToAll(Name, Class)
            # Add item to the corresponding module
            self.__InternalAddItemToModule(Name, Class)

    def RemoveItem(self, Name):
        if self.HasItem(Name, RegistryContext.PYTHON):
            self.__InternalRemoveItem(Name)
        elif self.HasItem(Name, RegistryContext.CPP):
            #TODO: Discuss about removing stuff from the c++ registry. Do we allow it or not?
            # self.__cpp_registry.RemoveItem(Name)
            err_msg = f"Trying to remove '{Name}' from c++ registry. This operation is forbidden."
            raise Exception(err_msg)
        else:
            err_msg = f"Asking to remove '{Name}' non-registered item."
            raise Exception(err_msg)

    def RemoveItems(self, Name, Context=RegistryContext.PYTHON):
        # TODO: Discuss about implementing in here a method to remove a list of items from the Python (not c++) registry
        pass

    def GetItem(self, Name):
        if self.HasItem(Name, RegistryContext.PYTHON):
            # Get the class from the Python registry
            return self.__InternalGetItem(Name)
        elif self.HasItem(Name, RegistryContext.CPP):
            # Split the registry time to get the item keyword
            item_keyword, _, _ = self.__SplitRegistryNameString(Name)

            # Get the object from the c++ registry by filtering its type keyword
            get_function_name = self.__GetCppRegistryGetFunctionNameFromItemKeyword(item_keyword)
            if hasattr(self.__cpp_registry, get_function_name):
                return getattr(self.__cpp_registry, get_function_name)(Name)
            else:
                err_msg = f"Retrieving '{item_keyword}' items from c++ registry is not supported."
                raise Exception(err_msg)
        else:
            err_msg = f"Asking to retrieve '{Name}' non-registered item."
            raise Exception(err_msg)

    def GetItems(self, Name, Context=RegistryContext.ALL):
        # TODO: Implement a GetItems or GetKeys method to return a list with the corresponding keys of that levels, like we do with parameters
        pass

    @staticmethod
    def __SplitRegistryNameString(Name):
        '''Split the registry name
        e.g. for Processes.KratosMultiphysics.ConvectionDiffusionApplication.MyProcess we should get
        - item_keyword: "Processes"
        - module_name: "KratosMultiphysics.ConvectionDiffusionApplication"
        - class_name: "MyProcess"
        '''

        split_name = Name.split('.')
        if len(split_name) < 3:
            err_msg = f"Provided item name '{Name}' is wrong. 'ItemKeyword.Module.ClassName' structure is expected (e.g. 'Processes.KratosMultiphysics.MyProcess')."
            raise Exception(err_msg)
        item_keyword = split_name[0]
        module_name = split_name[1:-1]
        class_name = split_name[-1]

        return item_keyword, module_name, class_name

    @staticmethod
    def __GetCppRegistryGetFunctionNameFromItemKeyword(ItemKeyword):
        if ItemKeyword in PythonRegistry.CppGetFunctionNamesMap:
            return PythonRegistry.CppGetFunctionNamesMap[ItemKeyword]
        else:
            err_msg = f"'{ItemKeyword}' has not been addded to 'CppGetFunctionNamesMap'."
            raise Exception(err_msg)

    def __InternalHasItem(self, Name):
        split_name = Name.split('.')
        aux_dict = self.__python_registry
        for i_name in split_name[:-1]:
            if not i_name in aux_dict:
                return False
            aux_dict = aux_dict[i_name]
        return split_name[-1] in aux_dict

    def __InternalGetItem(self, Name):
        split_name = Name.split('.')
        aux_dict = self.__python_registry
        for i_name in split_name[:-1]:
            aux_dict = aux_dict[i_name]
        return aux_dict[split_name[-1]]

    def __InternalRemoveItem(self, Name):
        split_name = Name.split('.')
        aux_dict = self.__python_registry
        for i_name in split_name[:-1]:
            aux_dict = aux_dict[i_name]
        aux_dict.pop(split_name[-1])

    def __InternalAddItemToAll(self, Name, Class):
        item_keyword, _, class_name = self.__SplitRegistryNameString(Name)
        if self.HasItem(f"{item_keyword}.All.{class_name}", RegistryContext.ALL):
            err_msg = f"Trying to register '{Name}' but there is already an item with the same '{class_name}' name in the '{item_keyword}.All' block."
            raise Exception(err_msg)

        if item_keyword in self.__python_registry:
            item_all_dict = self.__python_registry[item_keyword]["All"]
            item_all_dict[class_name] = Class
        else:
            self.__python_registry[item_keyword] = {}
            self.__python_registry[item_keyword]["All"] = {}
            self.__python_registry[item_keyword]["All"][class_name] = Class

    def __InternalAddItemToModule(self, Name, Class):
        split_name = Name.split('.')
        aux_dict = self.__python_registry
        for i_name in split_name[:-1]:
            if not i_name in aux_dict:
                aux_dict[i_name] = {}
            aux_dict = aux_dict[i_name]
        aux_dict[split_name[-1]] = Class

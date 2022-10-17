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
            is_registered_in_python = Name in self.__python_registry
            return is_registered_in_cpp or is_registered_in_python
        elif Context == RegistryContext.CPP:
            return self.__cpp_registry.HasItem(Name)
        elif Context == RegistryContext.PYTHON:
            return Name in self.__python_registry
        else:
            raise Exception("Wrong registry context. Accepted options are 'ALL', 'CPP' and 'PYTHON'.")

    def AddItem(self, Name, Class):
        # Split the current item registry name
        item_keyword, _, class_name = self.__SplitRegistryNameString(Name)
        all_key = f"{item_keyword}.All.{class_name}"

        # Add current item
        # First check if it is already registered in the c++ and Python registries
        # If not registered, add it to the Python registry (note that we do an "All" registry as well as the module one)
        if self.HasItem(Name, RegistryContext.ALL):
            err_msg = f"Trying to register '{Name}' but it is already registered."
            raise Exception(err_msg)
        else:
            if self.HasItem(all_key, RegistryContext.ALL):
                err_msg = f"Trying to register '{Name}' but there is already an item with the same '{class_name}' name in the '{item_keyword}.All' block."
                raise Exception(err_msg)
            else:
                #FIXME: THIS MUST BE HIERARCHICAL
                self.__python_registry[Name] = Class
                self.__python_registry[all_key] = Class
                #FIXME: THIS MUST BE HIERARCHICAL

    def RemoveItem(self, Name):
        if self.HasItem(Name, RegistryContext.PYTHON):
            self.__python_registry.pop(Name)
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
            return self.__python_registry[Name]
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

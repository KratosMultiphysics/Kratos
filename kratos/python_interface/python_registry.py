import KratosMultiphysics

class PythonRegistry(object):

    def __init__(self):
        # We get a "private" pointer to the c++ registry
        self.__cpp_registry = KratosMultiphysics.Registry

        # Then we point "KratosMultiphysics.Registry" to current Python registry
        # Note that we want to force the user to access the registry through this class
        # Also note that this can be safely done since the c++ registry occurs at compile time
        KratosMultiphysics.Registry = self

        # Dictionary to act as Python registry
        self.__python_registry = {}

    #TODO: Remove this. I don't think we'll manage to throw the exception
    @staticmethod
    def _CppRegistryException(Key):
        err_msg = "c++ registry must be accessed as 'KratosMultiphysics.KratosGlobals.Registry'."
        raise Exception(err_msg)

    def HasItem(self, Name):
        is_registered_in_cpp = self.__cpp_registry.HasItem(Name)
        is_registered_in_python = self.__python_registry.has_key(Name)
        return is_registered_in_cpp or is_registered_in_python

    def AddItem(self, Name, Class):
        # Split the current item registry name
        item_keyword, _, class_name = self.__SplitRegistryNameString(Name)

        # Add current item
        # First check if it is already registered in the c++ and Python registries
        # If not registered, add it to the Python registry (note that we do an "All" registry as well as the module one)
        if self.HasItem(Name):
            err_msg = f"Trying to register '{Name}' but it is already registered."
            raise Exception(err_msg)
        else:
            module_key = Name
            all_key = f"{item_keyword}.All.{class_name}"
            self.__python_registry[all_key] = Class
            self.__python_registry[module_key] = Class

    def RemoveItem(self, Name):
        if self.__cpp_registry.HasItem(Name):
            self.__cpp_registry.RemoveItem(Name)
        if self.__python_registry.has_key(Name):
            self.__python_registry.pop(Name)

    def GetItem(self, Name):
        if self.__python_registry.has_key(Name):
            # Get the class from the Python registry
            return self.__python_registry[Name]
        elif self.__cpp_registry.HasItem(Name):
            # Split the registry time to get the item keyword
            item_keyword, _, _ = self.__SplitRegistryNameString(Name)

            # Get the object from the c++ registry by filtering its type keyword
            if item_keyword == "Operations":
                return self.__cpp_registry.GetOperation(Name)
            elif item_keyword == "Processes":
                return self.__cpp_registry.GetProcess(Name)
            else:
                err_msg = f"Retrieving '{item_keyword}' items from c++ registry is not supported yet."
                raise Exception(err_msg)

    def __SplitRegistryNameString(self, Name):
        '''Split the registry name
        e.g. for Processes.KratosMultiphysics.ConvectionDiffusionApplication.MyProcess we should get
        - item_keyword: "Processes"
        - module_name: "KratosMultiphysics.ConvectionDiffusionApplication"
        - class_name: "MyProcess"
        '''

        split_name = Name.split('.')
        if len(split_name):
            err_msg = f"Provided item name '{Name}' is wrong. 'ItemKeyword.Module.ClassName' structure is expected (e.g. 'Processes.KratosMultiphysics.MyProcess')."
            raise Exception(err_msg)
        item_keyword = split_name[0]
        module_name = split_name[1:-1]
        class_name = split_name[-1]

        return item_keyword, module_name, class_name

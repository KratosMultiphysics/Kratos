import os
import pkgutil
import importlib

# Initialize an empty list for __all__
__all__ = []

# Get the current package name
package_name = __name__

# Iterate over all modules in the current package
for module_info in pkgutil.iter_modules([os.path.dirname(__file__)]):
    module_name = module_info.name

    try:
        # Import the module
        module = importlib.import_module(f".{module_name}", package_name)

        # Ensure the imported object is indeed a Python module
        if not hasattr(module, '__file__'):
            continue

        # Get all public attributes from the module
        public_attributes = [name for name in dir(module) if not name.startswith('_')]

        # Extend __all__ with public attributes
        __all__.extend(public_attributes)

        # Dynamically set the attributes in the current namespace
        globals().update({name: getattr(module, name) for name in public_attributes})

    except ImportError:
        # Skip the module if it cannot be imported
        continue
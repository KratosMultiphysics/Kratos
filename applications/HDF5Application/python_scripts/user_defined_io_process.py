'''Construct a user-defined input/output process.

This factory is for advanced use cases. It allows users to define a custom 
input/output process via json settings which are passed directly to
HDF5Application.core.

license: HDF5Application/license.txt
'''
import KratosMultiphysics
import KratosMultiphysics.HDF5Application.core as _core

def Factory(process_settings, Model):
    """Construct a user-defined input/output process for HDF5."""

    return _core.Factory(process_settings["Parameters"], Model)

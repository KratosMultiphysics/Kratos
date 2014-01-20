from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.OpenCLApplication import *


cl_sources = '../../custom_utilities'

device_group = OpenCLDeviceGroup(cl_device_type.CL_DEVICE_TYPE_GPU, True)
device_group.AddCLSearchPath(cl_sources)

BenchmarkDevice(device_group)

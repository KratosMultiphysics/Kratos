
from KratosMultiphysics import *
from KratosMultiphysics.OpenCLApplication import *


cl_sources = '../../custom_utilities'

device_group = OpenCLDeviceGroup(cl_device_type.CL_DEVICE_TYPE_GPU, True)
device_group.AddCLSearchPath(cl_sources)

BenchmarkDevice(device_group)

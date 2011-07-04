##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = '../../../../libs' ##kratos_root/libs
kratos_applications_path    = '../../../../applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_OpenCLApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosOpenCLApplication import *

cl_sources = '../../custom_utilities'

device_group = OpenCLDeviceGroup(cl_device_type.CL_DEVICE_TYPE_GPU, True)
device_group.AddCLSearchPath(cl_sources)

BenchmarkDevice(device_group)

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Time control starts
import time as timer
print(timer.ctime())
# Measure process time
t0p = timer.clock()
# Measure wall time
t0w = timer.time()

def StartTimeMeasuring():
    # Measure process time
    time_ip = timer.clock()
    return time_ip

def StopTimeMeasuring(time_ip, process, report):
    # Measure process time
    time_fp = timer.clock()
    if( report ):
        used_time = time_fp - time_ip
        print("::[Material Modeling]:: [ %.2f" % round(used_time,2),"s", process," ] ")


# Import kratos core and applications
import KratosMultiphysics

        
#### PARSING THE PARAMETERS ####

# Import input
parameter_file = open("material_parameters.json",'r')
ProjectParameters = KratosMultiphysics.Parameters( parameter_file.read())

# constitutive process
law_test_module  = __import__(ProjectParameters["material_model"]["python_module"].GetString())
material_process = law_test_module.CreateProcess(ProjectParameters["material_model"]["parameters"])

clock_time = StartTimeMeasuring()

material_process.ExecuteInitialize()

clock_time = StartTimeMeasuring()

# time testing (loop for massive calculation):
calls = 0
for i in range(0,calls):
    material_process.Execute()

material_process.ExecuteFinalize()

StopTimeMeasuring(clock_time,"integration time",True)

#### END SOLUTION ####

# Measure process time
tfp = timer.clock()
# Measure wall time
tfw = timer.time()

print("::[Material Modeling]:: [Elapsed Time = %.2f" % (tfp - t0p),"seconds] (%.2f" % (tfw - t0w),"seconds of cpu/s time)")

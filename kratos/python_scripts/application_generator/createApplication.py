from __future__ import print_function, absolute_import, division

import sys

from classes.elementCreator import ElementCreator
from classes.conditionCreator import ConditionCreator
from classes.processCreator import ProcessCreator
from classes.classMemberCreator import ClassMemberCreator
from classes.variableCreator import VariableCreator

from applicationGenerator import ApplicationGenerator

if len(sys.argv) != 2:
    print("Error: Invalid usage.\nPlease call this script with your application name (in camelCase) like this: ")
    print("python createApplication.py MyNewApplication")
    exit()

# Read the application name and generate Camel, Caps and Low
appNameCamel = sys.argv[1]

# Fetch the applications directory
debugApp = ApplicationGenerator(appNameCamel)


# Add KratosVariables
debugApp.AddVariables([
    VariableCreator(name='VELOCITY_MEAN', vtype='double', is3D=True),
    VariableCreator(name='VELOCITY_VARIANCE', vtype='double', is3D=True),
    VariableCreator(name='PRESSURE_MEAN', vtype='double'),
    VariableCreator(name='PRESSURE_VARIANCE', vtype='double')    
])

debugApp.AddProcesses([
    ProcessCreator("StandardDeviationProcess", author="Suneth Warnakulasuriya (https://github.com/sunethwarna)"),
    ProcessCreator("MeanProcess", author="Suneth Warnakulasuriya (https://github.com/sunethwarna)"),
    ProcessCreator("MedianProcess", author="Suneth Warnakulasuriya (https://github.com/sunethwarna)"),
    ProcessCreator("RootMeanSquareProcess", author="Suneth Warnakulasuriya (https://github.com/sunethwarna)"),
    ProcessCreator("MinMaxProcess", author="Suneth Warnakulasuriya (https://github.com/sunethwarna)")
])

debugApp.Generate()

print("Your application has been generated in: applications/{}Application".format(appNameCamel))

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
    VariableCreator(name='PERTURBATION_SIZE', vtype='double'),
    VariableCreator(name='ADAPT_PERTURBATION_SIZE', vtype='bool'),
    VariableCreator(name='HAS_ROTATION_DOFS', vtype='bool'),
    VariableCreator(name='SENSOR_NODE_ID', vtype='int'),
    VariableCreator(name='SENSOR_WEIGHT', vtype='double'),
    VariableCreator(name='SENSOR_DIRECTION', vtype='double', is3D=True),
])

# Add test element
debugApp.Generate()

print(f"Your application has been generated in: applications/{appNameCamel}Application")

import sys

from classes.elementCreator import ElementCreator
from classes.conditionCreator import ConditionCreator
from classes.processCreator import ProcessCreator
from classes.classMemberCreator import ClassMemberCreator
from classes.variableCreator import VariableCreator

from applicationGenerator import ApplicationGenerator

if len(sys.argv) != 2:
    print("Error: Invalid usage.\nPlease call this script with your application name (in camelCase) like this: ")
    print("python laplacian_application_example.py MyNewApplication")
    exit()

# Read the application name and generate Camel, Caps and Low
appNameCamel = sys.argv[1]

# Fetch the applications directory
debugApp = ApplicationGenerator(appNameCamel)

# Add KratosVariables
debugApp.AddVariables([
    VariableCreator(name='MY_SCALAR', vtype='double'),
])

# Add test element
debugApp.AddElements([
    ElementCreator('LaplacianElement')
    .AddDofs(['TEMPERATURE'])
    .AddFlags([])
    .AddClassMemberVariables([])
])

debugApp.AddConditions([
    ConditionCreator('PointSourceCondition')
    .AddDofs(['TEMPERATURE'])
    .AddFlags([])
])

# debugApp.AddProcesses([
#     ProcessCreator('DoSomethingProcess')
# ])

debugApp.Generate()

print("Your application has been generated in: applications/{}Application".format(appNameCamel))

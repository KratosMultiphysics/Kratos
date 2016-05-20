from classes.elementCreator import ElementCreator
from classes.conditionCreator import ConditionCreator
from classes.processCreator import ProcessCreator
from classes.variableCreator import VariableCreator

from applicationGenerator import ApplicationGenerator

# Read the application name and generate Camel, Caps and Low
appCamel = "MyLaplacian"  # DO NOT put "application" in the name, it will be added automatically

# Fetch the applications directory
debugApp = ApplicationGenerator(appCamel)

# Add KratosVariables
debugApp.AddVariables([
    VariableCreator(name='MY_SCALAR', vtype='double'),
    VariableCreator(name='MY_VECTOR', vtype='double', is3D=True),
])

# Add test element
debugApp.AddElements([
    ElementCreator('MyLaplacianElement')
    .AddDofs(['TEMPERATURE'])
])

# debugApp.AddConditions([
#     ConditionCreator('MyFaceCondition')
#     .AddDofs(['TEMPERATURE'])
# ])
#
# debugApp.AddProcesses([
#     ProcessCreator('DoSomethingProcess')
# ])

debugApp.Generate()

print("Your application has been generated in: applications/{}Application".format(appCamel))

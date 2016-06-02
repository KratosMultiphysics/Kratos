import sys

from classes.elementCreator import ElementCreator
from classes.conditionCreator import ConditionCreator
from classes.processCreator import ProcessCreator
from classes.classMemberCreator import ClassMemberCreator
from classes.variableCreator import VariableCreator

from applicationGenerator import ApplicationGenerator

# Read the application name and generate Camel, Caps and Low
appCamel = "MyLaplacian"

# Fetch the applications directory
debugApp = ApplicationGenerator(appCamel)

# Add KratosVariables
debugApp.AddVariables([
    VariableCreator(name='MY_SCALAR', vtype='double'),
    VariableCreator(name='MY_VECTOR', vtype='int', is3D=True),
])

# Add test element
debugApp.AddElements([
    ElementCreator('MyLaplacianElement')
    .AddDofs(['TEMPERATURE'])
    .AddDofs(['MY_SCALAR'])
    .AddFlags([])
    .AddClassMemberVariables([])
])

# debugApp.AddConditions([
#     ConditionCreator('MyFaceCondition')
# 	.AddDofs(['TEMPERATURE'])
#     .AddFlags([])
# ])
#
# debugApp.AddProcesses([
#     ProcessCreator('DoSomethingProcess')
# ])

debugApp.Generate()

print("Your application has been generated in: applications/{}Application".format(appCamel))

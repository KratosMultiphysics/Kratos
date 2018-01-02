import sys

from classes.elementCreator import ElementCreator
from classes.conditionCreator import ConditionCreator
from classes.processCreator import ProcessCreator
from classes.classMemberCreator import ClassMemberCreator
from classes.variableCreator import VariableCreator

from applicationGenerator import ApplicationGenerator

# Set the application name and generate Camel, Caps and Low
appNameCamel = "PureDiffusion"

# Fetch the applications directory
debugApp = ApplicationGenerator(appNameCamel)

# Add KratosVariables
debugApp.AddVariables([
    VariableCreator(name='POINT_HEAT_SOURCE', vtype='double'),
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

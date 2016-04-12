import sys

from classes.elementCreator import ElementCreator
from classes.conditionCreator import ConditionCreator
from classes.processCreator import ProcessCreator
from classes.classMemberCreator import ClassMemberCreator
from classes.variableCreator import VariableCreator

from applicationGenerator import ApplicationGenerator

# Read the application name and generate Camel, Caps and Low
appCamel = sys.argv[1]

# Fetch the applications directory
debugApp = ApplicationGenerator(appCamel)

# Add KratosVariables
debugApp.AddVariables([
    VariableCreator(name='1DVariable', vtype='double'),
    VariableCreator(name='3DVariable', vtype='int', is3D=True),
])

# Add test element
debugApp.AddElements([
    ElementCreator('CustomTestElement')
    .AddDofs(['DOF_1', 'DOF_2'])
    .AddFlags(['FLAG_1', 'FLAG_2'])
    .AddClassMemberVariables([
        ClassMemberCreator(name='Variable1', vtype='double *', default='nullptr'),
        ClassMemberCreator(name='Variable2', vtype='int', default='0'),
        ClassMemberCreator(name='Variable3', vtype='std::string &', default='0'),
        # ClassMemberCreator(name='Warnvar1', vtype='std::string &',)
    ])
])

debugApp.AddConditions([
    ConditionCreator('CustomTestCondition')
])

debugApp.AddProcesses([
    ProcessCreator('CustomTestProcessAlpha'),
    ProcessCreator('CustomTestProcessDelta')
])

debugApp.Generate()

print("Your application has been generated in: applications/{}Application".format(appCamel))

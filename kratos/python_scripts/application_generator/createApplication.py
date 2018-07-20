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
    VariableCreator(name='DOF_1', vtype='double'),
    VariableCreator(name='DOF_2', vtype='double'),
    VariableCreator(name='ScalarVariable', vtype='double'),
    VariableCreator(name='VectorVariable', vtype='double', is3D=True),
])

# Add test element
debugApp.AddElements([
    ElementCreator('CustomTestElement')
    .AddDofs(['DOF_1', 'DOF_2'])
    .AddFlags(['FLAG_1', 'FLAG_2'])
    .AddClassMemberVariables([
        ClassMemberCreator(name='VariableA', vtype='double *', default='nullptr'),
        ClassMemberCreator(name='VariableB', vtype='int', default='0'),
        ClassMemberCreator(name='VariableC', vtype='std::string', default='"Usefull String"'),
    ])
])

debugApp.AddConditions([
    ConditionCreator('CustomTestCondition')
])

debugApp.Generate()

print("Your application has been generated in: applications/{}Application".format(appNameCamel))

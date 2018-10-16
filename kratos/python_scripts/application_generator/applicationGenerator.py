from __future__ import print_function, absolute_import, division

import os
import sys
import shutil
import difflib

from classes.elementCreator import ElementCreator
from classes.conditionCreator import ConditionCreator
from classes.processCreator import ProcessCreator
from classes.variableCreator import VariableCreator

from utils.io import ToUpperFromCamel
from utils.io import ToLowerFromCamel
from utils.io import GetApplicationsDirectory
from utils.io import bcolors, Formatc

from utils.constants import ptab, ctab
from utils.templateRule import TemplateRule


class ApplicationGenerator(TemplateRule):

    def __init__(self, name):

        super(ApplicationGenerator, self).__init__()

        appStrPos = name.lower().find('application')
        maxi = 5
        while appStrPos != -1 and maxi > 0:
            oldname = name
            name = name[0:appStrPos] + name[appStrPos+len('application'):len(name)]

            msg = Formatc([
                {'color': bcolors.WARNING, 'msg': '[WARNING]'},
                {'color': bcolors.CYAN, 'msg': ' {}'.format(oldname)},
                {'color': None, 'msg': ' already contains the substring "'},
                {'color': bcolors.CYAN, 'msg': '{}'.format('Application')},
                {'color': None, 'msg': '". Removing... : '+name},
            ], sys.stderr)

            print(msg, file=sys.stderr)
            appStrPos = name.lower().find('application')
            maxi = maxi -1

        if difflib.SequenceMatcher(None, name.lower(), 'application').ratio() > 0.65:
            msg = Formatc([
                {'color': bcolors.WARNING, 'msg': '[WARNING]'},
                {'color': None, 'msg': ' Your application name contains something wrong after automatic fix, please select another name.'},
            ], sys.stderr)

            print(msg, file=sys.stderr)
            exit()


        self._appDir = GetApplicationsDirectory()

        self._nameCamel = name
        self._nameUpper = ToUpperFromCamel(self._nameCamel)
        self._nameLower = ToLowerFromCamel(self._nameCamel)

        self._elements = []
        self._conditions = []
        self._processes = []
        self._strategies = []
        self._variables = []

        self.isDerived = False

        # Define the application tree rules ( directory names, files, etc...)
        self.rules += [
            {'token': '@{APP_NAME_CAMEL}', 'value': self._nameCamel},
            {'token': '@{APP_NAME_CAPS}', 'value': self._nameUpper},
            {'token': '@{APP_NAME_LOW}', 'value': self._nameLower},
            {'token': '@{APP_DEPEND_LIST}', 'value': ''},
            {'token': '@{APP_SOURCE_LIST}', 'value': ''}
        ]

        # Define the templates files needed for every class
        self._classTemplateFiles = {
            'Elements': ['element_template.cpp', 'element_template.h'],
            'Conditions': ['condition_template.cpp', 'condition_template.h'],
            'Processes': ['process_template.cpp', 'process_template.h'],
            'Variables': [
                self._nameLower+'_application.h',
                self._nameLower+'_application.cpp',
                self._nameLower+'_application_variables.h',
                self._nameLower+'_application_variables.cpp',
                'custom_python/' + self._nameLower+'_python_application.cpp',
            ]
        }

        # Define where the generated classes should be placed
        self._classTemplatePath = {
            'Elements': ['classes', 'custom_elements'],
            'Conditions': ['classes', 'custom_conditions'],
            'Processes': ['classes', 'custom_processes']
        }

    def GenerateFile(self, src, dst, subsMap, removeOriginal=True):
        with open(src) as srcFile, open(dst, 'w+') as dstFile:
            for l in srcFile:
                for pair in subsMap:
                    l = l.replace(pair['token'], pair['value'])
                if removeOriginal:
                    dstFile.write(l)

    def AddElements(self, elements):
        self._addGenericClass(elements, self._elements, ElementCreator)

    def AddConditions(self, conditions):
        self._addGenericClass(conditions, self._conditions, ConditionCreator)

    def AddProcesses(self, processes):
        self._addGenericClass(processes, self._processes, ProcessCreator)

    def AddVariables(self, variables):
        self._addGenericClass(variables, self._variables, VariableCreator)

    def AddCustomStrategies(self, strategies):
        for strategy in strategies:
            self._elements.append(strategy)

    def AddCustomVariables(self, variables):
        for variable in variables:
            self._elements.append(variable)

    def AddCustomConstitutiveLaws(self, constituiveLaws):
        for constitutiveLaw in constituiveLaws:
            self._elements.append(constitutiveLaw)

    def Generate(self):
        ''' Generates the application by copying it from the template app
        '''

        root = os.path.dirname(os.path.realpath(__file__)) + "/"

        tpldir = root + "../../templates/" + "template_application"
        appdir = self._appDir + self._nameCamel + 'Application'

        # TODO: Catch the exception
        # Create the directory
        shutil.copytree(tpldir, appdir, symlinks=False, ignore=None)

        # Copy all files from the template application
        filesToRename = {
            "template_application.cpp.in":
                self._nameLower + "_application.cpp.in",
            "template_application.h.in":
                self._nameLower + "_application.h.in",
            "template_application_variables.cpp.in":
                self._nameLower + "_application_variables.cpp.in",
            "template_application_variables.h.in":
                self._nameLower + "_application_variables.h.in",
            "TemplateApplication.py.in":
                self._nameCamel + "Application.py.in",
            "tests/test_TemplateApplication.py.in":
                "tests/test_" + self._nameCamel + "Application.py.in",
            "custom_python/test_python_application.cpp.in":
                "custom_python/" + self._nameLower + "_python_application.cpp.in"
        }

        # Rename the files
        for key, value in filesToRename.items():
            os.rename(
                os.path.join(appdir, key),
                os.path.join(appdir, value)
            )

        # Add the sources to the CMakeLists
        self._addClassToSources(
            entityType='Elements',
            fromContainer=self._elements)
        self._addClassToSources(
            entityType='Conditions',
            fromContainer=self._conditions)
        self._addClassToSources(
            entityType='Processes',
            fromContainer=self._processes)

        # Replace the tokens in the app files
        for root, subfolder, files in os.walk(appdir):

            # This is important for some versions of python
            if '.svn' in subfolder:
                subfolder.remove('.svn')

            if '.git' in subfolder:
                subfolder.remove('.git')

            for f in files:
                src = os.path.join(root, f)
                fileName = src.split(".")
                dst = '.'.join(src.split(".")[0:len(fileName) - 1])

                self._applyTemplateRulesToFile(src, dst, self.rules)

                # Here we iterate over existing files, so we have to remove the original
                os.remove(src)

        # Generate the additional classes
        self._generateClassFrom(
            entityType='Elements',
            fromContainer=self._elements)
        self._generateClassFrom(
            entityType='Conditions',
            fromContainer=self._conditions)
        self._generateClassFrom(
            entityType='Processes',
            fromContainer=self._processes)

        # Add Variables
        self._generateVariablesFrom(
            entityType='Variables',
            fromContainer=self._variables)

        # add it to kratos
        self._addApplicationToCMake()
        # self._addApplicationToAppList()

    # Interal goes here
    def _applyTemplateRulesToFile(self, src, dst, rules):
        ''' Creates a 'dst' file by iterating over all pairs of [token,value]
            in 'rules' and replacing all ocurrences of 'token' by its 'value' in
            the source file 'src'.
        '''

        with open(src, 'r') as srcFile, open(dst, 'w+') as dstFile:
            for line in srcFile:
                for rule in rules:
                    line = line.replace(rule['token'], rule['value'])
                dstFile.write(line)

    def _addClassToSources(self, entityType, fromContainer):
        additionalSourceRule = self.GetRule('@{APP_SOURCE_LIST}')

        for entity in fromContainer:
            # Update the list of sources to be compiled
            additionalSourceRule['value'] += '\n' + ctab + '${CMAKE_CURRENT_SOURCE_DIR}/'
            additionalSourceRule['value'] += self._classTemplatePath[entityType][1] + '/' + entity.nameLower + '.cpp'

    def _generateClassFrom(self, entityType, fromContainer):
        # Generate elements
        root = os.path.dirname(os.path.realpath(__file__)) + "/"

        srcpath = root + "../../templates/" + self._classTemplatePath[entityType][0]
        dstpath = self._appDir + self._nameCamel + 'Application/' + self._classTemplatePath[entityType][1]

        if not os.path.exists(dstpath):
            os.makedirs(dstpath)

        for entity in fromContainer:

            files = self._classTemplateFiles[entityType]

            # Insert additional rules needed from the application
            appBaseNameRule = self.GetRule('@{APP_NAME_LOW}')       # To resolve the include path of the app file
            entity.rules.append(appBaseNameRule)

            for f in files:

                src = os.path.join(srcpath, f)
                ext = f.split(".")[-1]
                dst = os.path.join(dstpath, entity.nameLower + "." + ext)
                bkp = None

                if src == dst:
                    bkp = src + '.bkp'
                    shutil.copyfile(src, bkp)
                    os.remove(src)
                    src = bkp

                self._applyTemplateRulesToFile(src, dst, entity.rules)

                if bkp is not None:
                    os.remove(src)

    def _generateVariablesFrom(self, entityType, fromContainer):
        ''' Generate the list of variables for the application

            Input
            -----
            - entityType: type
                desc

            - fromContainer: type
                desc

            Rules
            -----
            @{KRATOS_DEFINE_VARIABLE_LIST}:
                The variable list

            @{KRATOS_CREATE_VARIABLE_LIST}:
                The default init values, in case they exists, for every constructor

            @{KRATOS_REGISTER_VARIABLE_LIST}:
                The init values for the copy constructor.

            @{KRATOS_REGISTER_PYTHON_VARIABLE_LIST}:
                The init values for the copy constructor.
        '''

        # Define the variable rules ( variables, defines, etc...)
        defineRule = {'token': '@{KRATOS_APP_DEFINE_VARIABLE_LIST}', 'value': ''}
        createRule = {'token': '@{KRATOS_APP_CREATE_VARIABLE_LIST}', 'value': ''}
        registerRule = {'token': '@{KRATOS_APP_REGISTER_VARIABLE_LIST}', 'value': ''}
        registerPythonRule = {'token': '@{KRATOS_APP_REGISTER_PYTHON_VARIABLE_LIST}', 'value': ''}

        root = self._appDir + self._nameCamel + 'Application/'

        srcpath = root

        for entity in fromContainer:
            defineRule['value'] += entity.defineString
            createRule['value'] += entity.createString
            registerRule['value'] += entity.registerString
            registerPythonRule['value'] += entity.registerPythonString

        files = self._classTemplateFiles[entityType]
        for f in files:

            src = os.path.join(srcpath, f)
            dst = src
            bkp = None

            if src == dst:
                bkp = src + '.bkp'
                shutil.copyfile(src, bkp)
                os.remove(src)
                src = bkp

            self._applyTemplateRulesToFile(src, dst, [
                defineRule,
                createRule,
                registerRule,
                registerPythonRule
            ])

            if bkp is not None:
                os.remove(src)

    def _addGenericClass(self, entities, container, containerType):
        msg = "{} is not an instance of KratosElementClass"

        for entity in entities:
            if isinstance(entity, containerType):
                container.append(entity)
            else:
                raise TypeError(msg.format(entity))

    def _addApplicationToCMake(self):
        # Locate the CMake file
        srcFile = self._appDir + "CMakeLists.txt"
        dstFile = self._appDir + "CMakeLists.txt.tmp"

        msgCount = 0

        # We can check it this way, as is not that large
        with open(srcFile, 'r') as src:
            if self._nameCamel in src.read():
                msg = Formatc([
                    {'color': bcolors.WARNING, 'msg': '[WARNING]'},
                    {'color': bcolors.CYAN, 'msg': ' {}'.format(self._nameCamel+'Application')},
                    {'color': None, 'msg': ' already has an entry in "'},
                    {'color': bcolors.CYAN, 'msg': '{}'.format('CMakeLists.txt')},
                    {'color': None, 'msg': '". Skipping step.'},
                ], sys.stderr)

                print(msg, file=sys.stderr)

                return

        with open(srcFile, 'r') as src, open(dstFile, 'w+') as dst:
            for l in src:

                # Add the applciation to the list message
                if l.strip() == "\")" :
                    newLine = self._nameUpper + '_APPLICATION;\\\n'

                    dst.write(newLine)

                # TODO: Probably is not a good idea to use the end comment as a
                # token for adding this, but we cannot use the end of the file
                # until that same comment is removed

                # Write the add subdirectory clause
                if '# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)' in l:
                    dst.write('if(${' + self._nameUpper + '_APPLICATION} MATCHES ON)\n')
                    dst.write('  add_subdirectory(' + self._nameCamel+'Application)\n')
                    dst.write('endif(${' + self._nameUpper + '_APPLICATION} MATCHES ON)\n')
                    dst.write('\n')

                # Write the old line
                dst.write(l)

        # Replace the old file with the new one
        os.remove(srcFile)
        os.rename(dstFile, srcFile)

    def _addApplicationToAppList(self):

        # Locate the applications lists
        srcFile = self._appDir + "applications_interface.py"
        dstFile = self._appDir + "applications_interface.py.tmp"

        # For this file is just easier to reconstruct the whole thing
        fileStruct = {
            'header': [],
            'importFalseBlock': [],
            'printMsgBlock': [],
            'appDirBlock': [],
            'importValueBlock': [],
            'prepareBlock': [],
            'initializeBlock': [],
            'footer': [],
        }

        currentBlock = 'header'

        # We can check it this way, as is not that large
        with open(srcFile, 'r') as src:
            if self._nameCamel in src.read():
                msg = Formatc([
                    {'color': bcolors.WARNING, 'msg': '[WARNING]'},
                    {'color': bcolors.CYAN, 'msg': ' {}'.format(self._nameCamel+'Application')},
                    {'color': None, 'msg': ' already has an entry in "'},
                    {'color': bcolors.CYAN, 'msg': '{}'.format('applications_interface.py')},
                    {'color': None, 'msg': '". Skipping step.'},
                ], sys.stderr)

                print(msg, file=sys.stderr)

                return

        # First read the file and parse all the info
        with open(srcFile, 'r') as src:
            for l in src:
                # Select the block
                if 'Import_' in l and 'Application = False' in l:
                    currentBlock = 'importFalseBlock'
                if 'print("Applications Available:")' in l and currentBlock == 'importFalseBlock':
                    currentBlock = 'printMsgBlock'
                if 'application_directory = os.path.dirname' in l and currentBlock == 'printMsgBlock':
                    currentBlock = 'appDirBlock'
                if 'def ImportApplications' in l and currentBlock == 'appDirBlock':
                    currentBlock = 'importValueBlock'
                if 'if(Import_' in l and currentBlock == 'importValueBlock':
                    currentBlock = 'prepareBlock'
                if '# dynamic renumbering of variables' in l and currentBlock == 'prepareBlock':
                    currentBlock = 'initializeBlock'
                if '# def ImportApplications(kernel  ):' in l and currentBlock == 'initializeBlock':
                    currentBlock = 'footer'

                # Append the result if its not null
                if l is not '\n' or currentBlock == 'prepareBlock':
                    fileStruct[currentBlock].append(l)

        # Prepare some blocks
        prepareBlockContent = [
            ptab + 'if(Import_{CAMEL}Application):\n',
            ptab * 2 + 'print("importing Kratos{CAMEL}Application ...")\n',
            ptab * 2 + 'sys.path.append(applications_path + \'/{CAMEL}/python_scripts\')\n',
            ptab * 2 + 'sys.path.append(applications_path + \'/{CAMEL}/Linux\')\n',
            ptab * 2 + 'from Kratos{CAMEL}Application import *\n',
            ptab * 2 + '{LOWER}_application = Kratos{CAMEL}Application()\n',
            ptab * 2 + 'kernel.ImportApplication({LOWER}_application)\n',
            ptab * 2 + 'print("Kratos{CAMEL}Application Succesfully imported")\n',
            '\n'
        ]

        prepareBlockContent = [
            s.format(
                CAMEL=self._nameCamel,
                LOWER=self._nameLower
            ) for s in prepareBlockContent]

        # Add our application to the requeired blocks
        fileStruct['importFalseBlock'].append(
            'Import_' + self._nameCamel + 'Application = False\n'
        )

        fileStruct['printMsgBlock'].append(
            'print("Import_' + self._nameCamel + 'Application: False")\n'
        )

        fileStruct['importValueBlock'].append(
            ptab + 'print("Import_{CAMEL}Application: " + str(Import_{CAMEL}Application))\n'.format(
                CAMEL=self._nameCamel)
        )

        fileStruct['prepareBlock'].append(
            prepareBlockContent
        )

        fileStruct['initializeBlock'].append([
            ptab + 'if(Import_' + self._nameCamel + 'Application):\n',
            ptab * 2 + 'kernel.InitializeApplication(' + self._nameLower + '_application)\n'
        ])

        # Write the whole thing down
        with open(dstFile, 'w+') as dst:
            for l in fileStruct['header']:
                dst.write(l)
            dst.write('\n')
            for l in fileStruct['importFalseBlock']:
                dst.write(l)
            dst.write('\n')
            for l in fileStruct['printMsgBlock']:
                dst.write(l)
            dst.write('\n')
            for l in fileStruct['appDirBlock']:
                dst.write(l)
            dst.write('\n')
            for l in fileStruct['importValueBlock']:
                dst.write(l)
            dst.write('\n')
            for b in fileStruct['prepareBlock']:
                for l in b:
                    dst.write(l)
            dst.write('\n')
            for b in fileStruct['initializeBlock']:
                for l in b:
                    dst.write(l)
            dst.write('\n')
            for l in fileStruct['footer']:
                dst.write(l)

        # Replace the old file with the new one
        os.remove(srcFile)
        os.rename(dstFile, srcFile)

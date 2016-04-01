import os
import sys
import shutil

from kratosClasses import KratosElementClass
from kratosClasses import KratosConditionClass
from kratosClasses import KratosProcessClass

from utils.io import ToUpperFromCamel
from utils.io import ToLowerFromCamel
from utils.io import GetApplicationsDirectory


class ApplicationGenerator(object):

    def __init__(self, name):

        self._appDir = GetApplicationsDirectory()

        self._nameCamel = name
        self._nameUpper = ToUpperFromCamel(self._nameCamel)
        self._nameLower = ToLowerFromCamel(self._nameCamel)

        self._elements = []
        self._conditions = []
        self._processes = []
        self._strategies = []
        self._variables = []

        self.__isDerived = False

        # What templates files need to be used for every class
        self._classTemplateFiles = {
            'Elements': ['element_template.cpp', 'element_template.h'],
            'Conditions': ['condition_template.h'],
            'Processes': ['process_template.h']
        }

        # Where they should be placed
        self._classTemplatePath = {
            'Elements': ['classes', 'custom_elements'],
            'Conditions': ['classes', 'custom_conditions'],
            'Processes': ['classes', 'custom_process']
        }

    def GenerateFile(self, src, dst, subsMap):
        with open(src) as srcFile, open(dst, 'w+') as dstFile:
            for l in srcFile:
                for pair in subsMap:
                    l = l.replace(pair['token'], pair['value'])
                dstFile.write(l)

    def AddCustomElements(self, elements):
        self._addGenericClass(elements, self._elements, KratosElementClass)

    def AddCustomConditions(self, conditions):
        self._addGenericClass(conditions, self._conditions, KratosConditionClass)

    def AddCustomProcesses(self, processes):
        self._addGenericClass(processes, self._processes, KratosProcessClass)

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
        ''' Generates the application by copying it from the template app '''

        root = os.path.dirname(os.path.realpath(__file__)) + "/"

        tpldir = root + "../../templates/" + "template_application"
        appdir = self._appDir + self._nameCamel + 'Application'

        # TODO: Catch the exception
        # Create the directory
        shutil.copytree(tpldir, appdir, symlinks=False, ignore=None)

        # Copy all files from the template application
        filesToRename = {
            "template_application.cpp.in":
                self._nameLower+"_application.cpp.in",
            "template_application.h.in":
                self._nameLower+"_application.h.in",
            "TemplateApplication.py.in":
                self._nameCamel+"Application.py.in",
            "tests/test_TemplateApplication.py.in":
                "tests/test_"+self._nameCamel+"Application.py.in",
            "custom_python/test_python_application.cpp.in":
                "custom_python/"+self._nameLower+"_python_application.cpp.in"
        }

        # Define the substitution rules
        rules = [
            {'token': '@{APP_NAME_CAMEL}', 'value': self._nameCamel},
            {'token': '@{APP_NAME_CAPS}', 'value': self._nameUpper},
            {'token': '@{APP_NAME_LOW}', 'value': self._nameLower},
            {'token': '@{APP_DEPEND_LIST}', 'value': ''}
        ]

        # Rename the files
        for key, value in filesToRename.items():
            os.rename(
                os.path.join(appdir, key),
                os.path.join(appdir, value)
            )

        # Replace the tokens in the app files
        for root, subfolder, files in os.walk(appdir):
            for f in files:
                src = os.path.join(root, f)
                fileName = src.split(".")
                dst = '.'.join(src.split(".")[0:len(fileName)-1])

                self._applyTemplateRulesToFile(src, dst, rules)

                # Here we iterate over existing files, so we have to remove
                # the original
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

        # add it to kratos
        self._addApplicationToCMake()
        self._addApplicationToAppList()

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

    def _generateClassFrom(self, entityType, fromContainer):
        # Generate elements
        root = os.path.dirname(os.path.realpath(__file__)) + "/"

        srcpath = root + "../../templates/" + self._classTemplatePath[entityType][0]
        dstpath = self._appDir + self._nameCamel + 'Application/' + self._classTemplatePath[entityType][1]

        if not os.path.exists(dstpath):
            os.makedirs(dstpath)

        for entity in fromContainer:
            files = self._classTemplateFiles[entityType]
            for f in files:

                src = os.path.join(srcpath, f)
                ext = f.split(".")[1]
                dst = os.path.join(dstpath, entity.nameLower + "." + ext)

                self._applyTemplateRulesToFile(src, dst, entity.rules)

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

        with open(srcFile, 'r') as src, open(dstFile, 'w+') as dst:
            for l in src:

                # Skip the first message
                if 'message( " ")' in l and msgCount == 0:
                    msgCount += 1

                # Add the applciation to the list message
                elif 'message( " ")' in l and msgCount == 1:
                    newLine = ''

                    newLine += 'message("' + self._nameUpper + '_APPLICATION'
                    newLine += ('.'*(25-len(self._nameUpper)))
                    newLine += '${' + self._nameUpper + '}")\n'

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
                if 'if(Import_SolidMechanicsApplication):' in l and currentBlock == 'importValueBlock':
                    currentBlock = 'prepareBlock'
                if '# dynamic renumbering of variables' in l:
                    currentBlock = 'initializeBlock'
                if '# def ImportApplications(kernel  ):' in l:
                    currentBlock = 'footer'

                # Append the result if its not null
                if l is not '\n' or currentBlock == 'prepareBlock':
                    fileStruct[currentBlock].append(l)

        # Prepare some blocks
        prepareBlockContent = [
            '\n'
            '\tif(Import_{CAMEL}):\n'
            '\t\tprint("importing Kratos{CAMEL}Application ...")\n'
            '\t\tsys.path.append(applications_path + \'/{CAMEL}/python_scripts\')\n'
            '\t\tsys.path.append(applications_path + \'/{CAMEL}/Linux\')\n'
            '\t\tfrom Kratos{CAMEL}Application import *\n'
            '\t\t{LOWER}_application = Kratos{CAMEL}Application()\n'
            '\t\tkernel.AddApplication({LOWER}_application)\n'
            '\t\tprint("Kratos{CAMEL}Application Succesfully imported")\n'
            '\n'
        ]

        prepareBlockContent = [
            s.format(
                CAMEL=self._nameCamel,
                LOWER=self._nameLower
            ).replace('\t', '    ') for s in prepareBlockContent
        ]

        # Add our application to the requeired blocks
        fileStruct['importFalseBlock'].append(
            'Import_' + self._nameCamel + ' = False\n'
        )

        fileStruct['printMsgBlock'].append(
            'print("Import_' + self._nameCamel + ': False")\n'
        )

        fileStruct['importValueBlock'].append(
            '\tprint("Import_{CAMEL}: " + str(Import_{CAMEL}))\n'.format(
                CAMEL=self._nameCamel
            ).replace('\t', '    ')
        )

        # This line deletes the last \n
        fileStruct['prepareBlock'] = fileStruct['prepareBlock'][:-1]
        fileStruct['prepareBlock'].append(
            prepareBlockContent
        )

        fileStruct['initializeBlock'].append([
            '\tif(Import_'+self._nameCamel+'Application):\n'.replace('\t', '    '),
            '\t\tkernel.InitializeApplication('+self._nameLower+'_application)\n'.replace('\t', '    ')
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


# TODO: Main sould exists in a separate file
def main():

    # Read the application name and generate Camel, Caps and Low
    appCamel = sys.argv[1]

    # Fetch the applications directory
    debugApp = ApplicationGenerator(appCamel)

    # Add test element
    debugApp.AddCustomElements([
        KratosElementClass("CustomTestElement")
    ])

    debugApp.AddCustomConditions([
        KratosConditionClass("CustomTestCondition")
    ])

    debugApp.AddCustomProcesses([
        KratosProcessClass("CustomTestProcessAlpha"),
        KratosProcessClass("CustomTestProcessDelta")
    ])

    debugApp.Generate()

    # Remind the user to change other files
    print("Your application has been generated in: applications/{}Application".format(appCamel))

if __name__ == "__main__":
    main()

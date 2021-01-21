from __future__ import print_function, absolute_import, division

import difflib
import re
import os
import sys

from utils.io import GetApplicationsDirectory
from utils.io import GetKratosDirectory
from utils.io import ToUpperFromCamel, ToLowerFromCamel
from utils.io import bcolors, Formatc, TestCamel

from utils.constants import ctab
from utils.templateRule import TemplateRule


class ClassCreator(TemplateRule):
    def __init__(
        self, name, base=None, template=None,
        members=None, procedures=None, author='KratosAppGenerator'
    ):
        # Check that name is in camelcase, or at least something parseable
        super(ClassCreator, self).__init__()

        if not TestCamel(name):
            msg = Formatc([
                {'color': bcolors.FAIL, 'msg': '[ERROR]'},
                {'color': None, 'msg': ' Name: "'},
                {'color': bcolors.CYAN, 'msg': '{}'.format(name)},
                {'color': None, 'msg': '" is not camelcase'},
            ], sys.stderr)

            raise ValueError(msg)

        # Check name, if a class already exists with that name in the context
        # specified, ask if you want to continue. In case of continuing a
        # warning message will be shown (in red?) and a backup of the class
        # will be created.

        self.author = author
        self.nameCamel = name
        self.nameUpper = ToUpperFromCamel(name)
        self.nameLower = ToLowerFromCamel(name)
        self.baseHead = [': public ' + base, ''][base is None]
        self.base = base
        self.template = template
        self.procedures = procedures

        # List of common rules for all classes

        self.rules += [
            {'token': '@{KRATOS_APP_AUTHOR}', 'value': ', '+self.author},
            {'token': '@{KRATOS_NAME_CAMEL}', 'value': self.nameCamel},
            {'token': '@{KRATOS_NAME_UPPER}', 'value': self.nameUpper},
            {'token': '@{KRATOS_NAME_LOWER}', 'value': self.nameLower},
            {'token': '@{KRATOS_CLASS_BASE}', 'value': self.base},
            {'token': '@{KRATOS_CLASS_BASE_HEADER}', 'value': self.baseHead},
            {'token': '@{KRATOS_SYSTEM_INCLUDES}', 'value': ''},
            {'token': '@{KRATOS_EXTERNAL_INCLUDES}', 'value': ''},
            {'token': '@{KRATOS_PROJECT_INCLUDES}', 'value': ''},
            {'token': '@{KRATOS_CLASS_TEMPLATE}', 'value': ''},
            {'token': '@{KRATOS_INCLUDE_LIST}', 'value': ''},
            {'token': '@{KRATOS_MEMBERS_LIST}', 'value': ''},
            {'token': '@{KRATOS_STATIC_MEMBERS_LIST}', 'value': ''},
            {'token': '@{KRATOS_INIT_MEMBER_LIST}', 'value': ''},
            {'token': '@{KRATOS_CC_INIT_MEMBER_LIST}', 'value': ''},
            {'token': '@{KRATOS_ASIGN_INIT_MEMBER_LIST}', 'value': ''},
            {'token': '@{KRATOS_PROCEDURES_LIST}', 'value': ''}
        ]

        [a, b, c, d] = self._CheckClassNameAvailable(
            self.nameCamel, "Kratos", True
        )

        for m in c.values():
            msg = Formatc([
                {'color': bcolors.FAIL, 'msg': '[ERROR]'},
                {'color': None, 'msg': ' Class with name: "'},
                {'color': bcolors.CYAN, 'msg': '{}'.format(name)},
                {'color': None, 'msg': '" already defined in:'},
            ], sys.stderr)
            print(msg, file=sys.stderr)
            print("  ["+m[0]+"]: \""+m[1].replace("\n", "")+"\""+"\n",
                  file=sys.stderr)

        for m in d.values():
            msg = Formatc([
                {'color': bcolors.WARNING, 'msg': '[WARNING]'},
                {'color': None, 'msg': ' Class with name: "'},
                {'color': bcolors.CYAN, 'msg': '{}'.format(name)},
                {'color': None, 'msg': '" probably ({}%) defined in:'.format(
                    str(m[2]*100))},
            ], sys.stderr)
            print(msg, file=sys.stderr)
            print("  ["+m[0]+"]: \""+m[1].replace("\n", "")+"\""+"\n",
                  file=sys.stderr)

    def AddClassMemberVariables(self, members):
        ''' Adds a list of class member variables to the class and initializes them
            in case it is possible

            Input
            -----
            members: list of 'KratosClasMemeber'
                the list of class members to be add to the class

            Rules
            -----
            @{KRATOS_MEMBERS_LIST}:
                The variable list

            @{KRATOS_INIT_MEMBER_LIST}:
                The default init values, in case they exists, for every constructor

            @{KRATOS_CC_INIT_MEMBER_LIST}:
                The init values for the copy constructor.
        '''

        membersRule = self.GetRule('@{KRATOS_MEMBERS_LIST}')
        staticMembersRule = self.GetRule('@{KRATOS_STATIC_MEMBERS_LIST}')
        initMembersRule = self.GetRule('@{KRATOS_INIT_MEMBER_LIST}')
        ccInitMembersRule = self.GetRule('@{KRATOS_CC_INIT_MEMBER_LIST}')

        for m in members:

            # Add decorators to the variable name based on its type
            name = 'm' + (m.isPtr * 'p') + (m.isRef * 'r') + m.name

            # Create the variable
            if m.isStatic:
                staticMembersRule['value'] += ctab + m.type + ' ' + name + ';\n'
            else:
                membersRule['value'] += ctab + m.type + ' ' + name + ';\n'

                # We will always preffer initialization list to avoid double-initialization (defaul + body) and
                # to keet the VTable clean

                # Generic constructors
                if m.default is not None:
                    if initMembersRule['value'] is not None:
                        initMembersRule['value'] += '\n' + ctab * 3 + ', '
                    initMembersRule['value'] += name + '(' + m.default + ')'

                # Fill the copy constructor
                if ccInitMembersRule['value'] is not None:
                    ccInitMembersRule['value'] += '\n' + ctab * 3 + ', '
                ccInitMembersRule['value'] += name + '(rOther.' + name + ')'

                # Assignment operator
                # TODO

        return self

    def AddProcudures(self, procedures):
        # TODO
        pass

    def _IsClassDefinition(self, line):
        # Match class definition and possible decorators
        expr = re.compile('(^class) ([\w]+) ([\w]*)')
        r = expr.match(line)
        if r:
            return True, r.group(2), r.group(3)
        else:
            return False, 0, 0

    def _GenerateCandidateFiles(self, files, name):
        nocache = [f for f in files if (len(f.split('.')) < 1 or f.split('.')[-1] in ['h', 'hpp', 'cpp'])]
        return [f for f in nocache if difflib.SequenceMatcher(None, name, f).ratio() > 0.4]

    def _CheckClassNameInDir(
        self, where, name,
        appRep, appSim, repFiles, simFiles
    ):

        for root, subfolder, files in os.walk(where):
            for f in self._GenerateCandidateFiles(files, name):
                src = os.path.join(root, f)
                with open(src, 'r') as _file:
                    for l in _file:
                        [found, m1, m2] = self._IsClassDefinition(l)
                        if found:
                            for m in [k.lower() for k in [m1, m2] if k]:
                                seq = difflib.SequenceMatcher(None, name, m)
                                if seq.ratio() == 1 and f not in repFiles:
                                    appRep = True
                                    repFiles[f] = [f, l, 1]
                                elif seq.ratio() > 0.75 and f not in simFiles:
                                    appSim = True
                                    simFiles[f] = [f, l, round(seq.ratio(), 3)]

        return appRep, appSim, repFiles, simFiles

    def _CheckFileNameInDir(
        self, where, name,
        appRep, appSim, repFiles, simFiles
    ):

        for root, subfolder, files in os.walk(where):
            for f in files:
                seq = difflib.SequenceMatcher(None, name, f)
                if seq.ratio() == 1:
                    pass
                elif seq.ratio() > 0.50:
                    pass

    def _CheckClassNameAvailable(
        self, name,
        context="Application",
        searchSimilar=False
    ):

        appRep = False
        appSim = False

        repFiles = {}
        simFiles = {}

        # Search in the application context
        [appRep, appSim, repFiles, simFiles] = self._CheckClassNameInDir(
            GetApplicationsDirectory()+self.nameCamel+"Application",
            self.nameCamel.lower(),
            appRep, appSim,
            repFiles, simFiles
        )

        if context == "Application":
            if searchSimilar:
                return appRep, appSim, repFiles, simFiles
            else:
                return appRep, repFiles

        # Search also in the kratos context:
        [appRep, appSim, repFiles, simFiles] = self._CheckClassNameInDir(
            GetKratosDirectory(),
            self.nameCamel.lower(),
            appRep, appSim,
            repFiles, simFiles
        )

        if context == "Kratos":
            if searchSimilar:
                return appRep, appSim, repFiles, simFiles
            else:
                return appRep, repFiles

        # Search in other apps... "Everywhere":
        for root, dirs, files in os.walk(GetApplicationsDirectory()):
            for subapp in dirs:
                [appRep, appSim, repFiles, simFiles] = self._CheckClassNameInDir(
                    GetApplicationsDirectory() + subapp,
                    self.nameCamel.lower(),
                    appRep, appSim,
                    repFiles, simFiles
                )

        if searchSimilar:
            return appRep, appSim, repFiles, simFiles
        else:
            return appRep, repFiles

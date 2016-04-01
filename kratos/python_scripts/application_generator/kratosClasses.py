import difflib
import re
import os
import sys

from utils.io import GetApplicationsDirectory
from utils.io import GetKratosDirectory
from utils.io import ToUpperFromCamel, ToLowerFromCamel
from utils.io import bcolors


def TestCamel(s):
    expr = re.compile(r'^(?:[A-Z][a-z]+)+$')
    return expr.match(s)


def formatc(stringList, where):

    string = ""

    for sl in stringList:
        if where.isatty() and sl['color']:
            string += sl['color'] + sl['msg'] + bcolors.ENDC
        else:
            string += sl['msg']

    return string


class KratosClassMember(object):
    def __init__(self, name, memberType=None, initValue=None):
        self.name = name
        self.type = memberType
        self.initValue = initValue

        self.isRef = '&' in memberType


class KratosClassProcedure(object):
    def __init__(self, name, returnType=None, parameters=None, template=None):
        self.name = name
        self.type = returnType
        self.template = template
        self.body = None


class KratosClass(object):
    def __init__(
        self, name, base=None, template=None,
        members=None, procedures=None, author='KratosAppGenerator'
    ):

        # Check that name is in camelcase, or at least something parseable
        if not TestCamel(name):
            msg = formatc([
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
        self.base = [': public '+base, ''][base is None]
        self.template = template
        self.procedures = procedures
        self.rules = [
            {'token': '@{KRATOS_APP_AUTHOR}', 'value': self.author},
            {'token': '@{KRATOS_NAME_CAMEL}', 'value': self.nameCamel},
            {'token': '@{KRATOS_NAME_UPPER}', 'value': self.nameUpper},
            {'token': '@{KRATOS_NAME_LOWER}', 'value': self.nameLower},
            {'token': '@{KRATOS_CLASS_BASE}', 'value': self.base},
            {'token': '@{KRATOS_SYSTEM_INCLUDES}', 'value': ''},
            {'token': '@{KRATOS_EXTERNAL_INCLUDES}', 'value': ''},
            {'token': '@{KRATOS_PROJECT_INCLUDES}', 'value': ''},
            {'token': '@{KRATOS_CLASS_TEMPLATE}', 'value': ''},
            {'token': '@{KRATOS_DEPEND_LIST}', 'value': ''}
        ]

        [a, b, c, d] = self.CheckClassNameAvailable(
            self.nameCamel, "Kratos", True
        )

        for m in c.values():
            msg = formatc([
                {'color': bcolors.FAIL, 'msg': '[ERROR]'},
                {'color': None, 'msg': ' Class with name: "'},
                {'color': bcolors.CYAN, 'msg': '{}'.format(name)},
                {'color': None, 'msg': '" already defined in:'},
            ], sys.stderr)
            print(msg, file=sys.stderr)
            print("  ["+m[0]+"]: \""+m[1].replace("\n", "")+"\""+"\n",
                  file=sys.stderr)

        for m in d.values():
            msg = formatc([
                {'color': bcolors.WARNING, 'msg': '[WARNING]'},
                {'color': None, 'msg': ' Class with name: "'},
                {'color': bcolors.CYAN, 'msg': '{}'.format(name)},
                {'color': None, 'msg': '" probably ({}%) defined in:'.format(
                    str(m[2]*100))},
            ], sys.stderr)
            print(msg, file=sys.stderr)
            print("  ["+m[0]+"]: \""+m[1].replace("\n", "")+"\""+"\n",
                  file=sys.stderr)

    def AddClassMember(self, member):
        self.procedures.append(member)

    def AddClassProcudure(self, procedure):
        self.procedures.append(procedure)

    def _IsClassDefinition(self, line):
        # Match class definition and possible decorators
        expr = re.compile('(^class) ([\w]+) ([\w]*)')
        r = expr.match(line)
        if r:
            return True, r.group(2), r.group(3)
        else:
            return False, 0, 0

    def _GenerateCandidateFiles(self, files, name):
        return [f for f in files if
                difflib.SequenceMatcher(None, name, f).ratio() > 0.4]

    def _CheckClassNameInDir(
      self, where, name,
      appRep, appSim, repFiles, simFiles):

        for root, subfolder, files in os.walk(where):
            for f in self._GenerateCandidateFiles(files, name):
                src = os.path.join(root, f)
                with open(src, 'r', encoding='latin1') as _file:
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
      appRep, appSim, repFiles, simFiles):

        for root, subfolder, files in os.walk(where):
            for f in files:
                seq = difflib.SequenceMatcher(None, name, f)
                if seq.ratio() == 1:
                    pass
                elif seq.ratio() > 0.50:
                    pass

    def CheckClassNameAvailable(
      self, name,
      context="Application",
      searchSimilar=False):

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
                    GetApplicationsDirectory()+subapp,
                    self.nameCamel.lower(),
                    appRep, appSim,
                    repFiles, simFiles
                )

        if searchSimilar:
            return appRep, appSim, repFiles, simFiles
        else:
            return appRep, repFiles


class KratosElementClass(KratosClass):
    def __init__(
        self, name, base='Element', template=None,
        members=None, procedures=None, author='KratosAppGenerator'
    ):
        super(KratosElementClass, self).__init__(
            name, base, template, members,
            procedures, author
        )


class KratosConditionClass(KratosClass):
    def __init__(
        self, name, base='Condition', template=None,
        members=None, procedures=None, author='KratosAppGenerator'
    ):
        super(KratosConditionClass, self).__init__(
            name, base, template, members,
            procedures, author
        )


class KratosProcessClass(KratosClass):
    def __init__(
        self, name, base='Process', template=None,
        members=None, procedures=None, author='KratosAppGenerator'
    ):
        super(KratosProcessClass, self).__init__(
            name, base, template, members,
            procedures, author
        )

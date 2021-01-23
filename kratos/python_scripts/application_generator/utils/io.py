from __future__ import print_function, absolute_import, division

import re
import os

# From blender...
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    CYAN = '\033[36m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def TestCamel(s):
    expr = re.compile(r'^(?:[A-Z][a-z]+)+$')
    return expr.match(s)

def RemoveSubstring(s, d):
    return s.replace(d,'')

def Formatc(stringList, where):

    string = ""

    for sl in stringList:
        if where.isatty() and sl['color']:
            string += sl['color'] + sl['msg'] + bcolors.ENDC
        else:
            string += sl['msg']

    return string

def ToUpperFromCamel(appCamel):
    ''' Converts a Camel-Case string into a upercase snake_case string '''

    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', appCamel)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).upper()


def ToLowerFromCamel(appCamel):
    ''' Converts a Camel-Case string into a lowercase snake_case string'''

    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', appCamel)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


def CheckNameAvail(appsdir, appname):
    ''' Checks that the application does not exists

    Raises
    ------
    Error :
        If the application "appname" already exists

    '''

    appDir = appsdir + appname + "Application"

    if os.path.exists(appDir):
        raise Exception('Application "{}" already exsists'.format(appname))


def GetApplicationsDirectory():
    ''' Return the path to the applications directory '''

    return os.path.dirname(os.path.realpath(__file__)) + "/../../../../applications/"


def GetKratosDirectory():
    ''' Return the path to the applications directory '''

    return os.path.dirname(os.path.realpath(__file__)) + "/../../../../kratos/"

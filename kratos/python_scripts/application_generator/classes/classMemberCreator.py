from __future__ import print_function, absolute_import, division

import sys

from utils.io import bcolors, Formatc


class ClassMemberCreator(object):
    def __init__(self, name, vtype=None, default=None):

        # This wil be extended with valid and invalid types
        if vtype is None:
            msg = Formatc([
                {'color': bcolors.FAIL, 'msg': '[ERROR]'},
                {'color': None, 'msg': ' while creating variable: "'},
                {'color': bcolors.CYAN, 'msg': '{}'.format(name)},
                {'color': None, 'msg': '" :No type defined'},
            ], sys.stderr)

            raise ValueError(msg)

        self.isRef = '&' in vtype
        self.isPtr = '*' in vtype
        self.isStatic = 'static' in vtype
        self.default = default

        if self.isRef and self.default is None:
            msg = Formatc([
                {'color': bcolors.WARNING, 'msg': '[WARNING]'},
                {'color': None, 'msg': ' while creating variable: "'},
                {'color': bcolors.CYAN, 'msg': '{}'.format(name)},
                {'color': None, 'msg': '" : Reference variable without default value.'},
            ], sys.stderr)

            print(msg, file=sys.stderr)

        if self.isStatic and self.default is not None:
            msg = Formatc([
                {'color': bcolors.WARNING, 'msg': '[WARNING]'},
                {'color': None, 'msg': ' while creating variable: "'},
                {'color': bcolors.CYAN, 'msg': '{}'.format(name)},
                {'color': None, 'msg': '" : static variable default values are not currently supported.\n'},
                {'color': None, 'msg': 'Variable will be left unitialized'},
            ], sys.stderr)

            print(msg, file=sys.stderr)

            self.default = None

        self.name = name
        self.type = vtype

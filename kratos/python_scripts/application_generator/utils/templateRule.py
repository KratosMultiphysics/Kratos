from __future__ import print_function, absolute_import, division

import sys

from utils.io import bcolors, Formatc

class TemplateRule(object):

    def __init__(self):
        self.rules = []

    def GetRule(self, ruleName):
        try:
            return next(r for r in self.rules if r['token'] == ruleName)
        except:
            msg = Formatc([
                {'color': bcolors.FAIL, 'msg': '[ERROR]'},
                {'color': None, 'msg': ' Rule: "'},
                {'color': bcolors.CYAN, 'msg': '{}'.format(ruleName)},
                {'color': None, 'msg': '" does not exist'},
            ], sys.stderr)

            raise ValueError(msg)

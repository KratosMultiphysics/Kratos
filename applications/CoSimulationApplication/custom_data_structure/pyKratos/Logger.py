from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

class Logger(object):

    @staticmethod
    def PrintInfo(label, *args):
        print(label, " ".join(map(str,args)))

    @staticmethod
    def PrintWarning(label, *args):
        print(label, " ".join(map(str,args)))
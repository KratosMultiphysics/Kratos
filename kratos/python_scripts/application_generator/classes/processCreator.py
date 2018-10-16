from __future__ import print_function, absolute_import, division

from classes.classCreator import ClassCreator

# from utils.constants import ctab


class ProcessCreator(ClassCreator):
    ''' TO BE IMPLEMENTED '''

    def __init__(
        self, name, base='Process', template=None,
        members=None, procedures=None, author='KratosAppGenerator'
    ):
        super(ProcessCreator, self).__init__(
            name, base, template, members,
            procedures, author
        )

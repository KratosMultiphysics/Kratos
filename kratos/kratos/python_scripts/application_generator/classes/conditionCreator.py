from classes.classCreator import ClassCreator

# from utils.constants import ctab


class ConditionCreator(ClassCreator):
    def __init__(
        self, name, base='Condition', template=None,
        members=None, procedures=None, author='KratosAppGenerator'
    ):
        super(ConditionCreator, self).__init__(
            name, base, template, members,
            procedures, author
        )

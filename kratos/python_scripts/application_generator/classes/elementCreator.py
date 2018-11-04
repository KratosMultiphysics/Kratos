from __future__ import print_function, absolute_import, division

from classes.classCreator import ClassCreator
from utils.constants import ctab


class ElementCreator(ClassCreator):
    ''' Extends KratosClassCreator for a Kratos Element. It expects a element_template.cpp and
        element_template.h as base file.

        Input
        -----
        - name: string
            name of the class

        - base: string
            name of the base element, if any. Defaults to 'Element'.

        - template:: instance of KratosTemplate
            not yet implemented.

        - members: list of KratosClassMember
            class memebers to be added.

        - procedures: list of KratosClassProcedures
            class procedures to be added. NYI

        Output
        ------
        - Instance of a KratosElementCreator

        Rules
        -----
        - @{KRATOS_CLASS_LOCAL_FLAGS}:
            The flag list

        - @{KRATOS_ELEMENT_LIST_DOFS}:
            The dof list

        - @{KRATOS_ELEMENT_ECUATION_ID_DOFS}:
            The dof list for the ecuation id's
    '''

    def __init__(
        self, name, base='Element', template=None,
        members=None, procedures=None, author='KratosAppGenerator'
    ):
        super(ElementCreator, self).__init__(
            name, base, template, members,
            procedures, author
        )

        # Register additional rules for elements
        self.rules += [
            {'token': '@{KRATOS_CLASS_LOCAL_FLAGS}', 'value': ''},
            {'token': '@{KRATOS_ELEMENT_LIST_DOFS}', 'value': ''},
            {'token': '@{KRATOS_ELEMENT_ECUATION_ID_DOFS}', 'value': ''},
            {'token': '@{KRATOS_CLASS_BASE_DIR}', 'value': 'custom_elements'}
        ]

    def AddFlags(self, flagList):
        ''' Adds a list of flags to the element and initializes

            Input
            -----
            - flagList: list of strings
                the list of falg names to be added

            Output
            ------
            - returns itself to allow chaining

            Rules
            -----
            - @{KRATOS_CLASS_LOCAL_FLAGS}:
                The flag list
        '''

        localFlagsValue = ''
        localFlagDefinition = ctab + 'KRATOS_DEFINE_LOCAL_FLAG({});\n'

        for flagName in flagList:
            upperName = flagName.upper()
            if upperName != flagName:
                print('Flag variable name appears to be in incorrect format')
            localFlagsValue += localFlagDefinition.format(upperName)

        # Locate the flag rule and replace the value
        flagRule = self.GetRule('@{KRATOS_CLASS_LOCAL_FLAGS}')
        flagRule['value'] += localFlagsValue

        return self

    def AddDofs(self, dofList):
        ''' Adds a list of dofs to the element and initializes. It is guaranteed that
            the order of the list will be respected.

            Input
            -----
            - dofList: list of strings
                the list of dof names to be added

            Output
            ------
            - returns itself to allow chaining

            Rules
            -----
            - @{KRATOS_ELEMENT_LIST_DOFS}:
                The dof list

            - @{KRATOS_ELEMENT_ECUATION_ID_DOFS}:
                The dof list for the ecuation id's
        '''

        # Dof are used to fill a couple of macros. Dofs can only be replaced, not modified (add/remove)
        # List of dofs
        localDofListValue = ''
        localDofListDefinition = [
            ctab * 1 + 'for (unsigned int i = 0; i < number_of_nodes; i++)\n',
            ctab * 2 + 'rElementalDofList[i] = GetGeometry()[i].pGetDof({});\n',
            '\n'
        ]

        # EcuationID's
        localDofEidValue = ''
        localDofEidDefinition = [
            ctab * 1 + 'for (unsigned int i = 0; i < number_of_nodes; i++)\n',
            ctab * 2 + 'rResult[i] = GetGeometry()[i].GetDof({}).EquationId();\n',
            '\n'
        ]

        # Assemble the code
        for dofName in dofList:
            upperName = dofName.upper()
            if upperName != dofName:
                print('Dof name appears to be in incorrect format')
            for l in localDofListDefinition:
                localDofListValue += l.format(upperName)
            for l in localDofEidDefinition:
                localDofEidValue += l.format(upperName)

        # Locate the rules and replace the values
        dofListRule = self.GetRule('@{KRATOS_ELEMENT_LIST_DOFS}')
        dofEidRule = self.GetRule('@{KRATOS_ELEMENT_ECUATION_ID_DOFS}')

        dofListRule['value'] = localDofListValue
        dofEidRule['value'] = localDofEidValue

        return self

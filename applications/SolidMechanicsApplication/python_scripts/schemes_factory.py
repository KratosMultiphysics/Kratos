## This script collects the available schemes to be used in the SolidMechanicsApplication

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# SolutionScheme class
class SolutionScheme:

    def __init__(self, custom_settings, dofs_list):

        default_settings = KratosMultiphysics.Parameters("""
        {
           "solution_type": "Dynamic",
  	   "analysis_type": "Non-linear",
           "time_integration": "Implicit",
           "integration_method": "Newmark",
           "time_integration_order": 1,
           "buffer_size": 2
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.dofs = []
        for i in range(0, dofs_list.size() ):
            self.dofs.append(dofs_list[i].GetString())

        # add default DISPLACEMENT dof
        if( len(self.dofs) == 0 or (len(self.dofs) == 1 and self.dofs[0] =="ROTATION") ):
            self.dofs.append('DISPLACEMENT')

    def GetVariables(self):

        self._set_variables_and_dofs()

        nodal_variables = self.nodal_variables + self.dof_variables + self.dof_reactions + self.dof_derivatives

        return nodal_variables

    def GetDofsAndReactions(self):

        if not hasattr(self, '_dof_variables'):
            self._set_variables_and_dofs()

        return self.dof_variables, self.dof_reactions

    def GetSolutionScheme(self):

        vector_integration_methods = []
        scalar_integration_methods = []

        for dof in self.dofs:

            ##check if variable type is a vector
            kratos_variable = KratosMultiphysics.KratosGlobals.GetVariable(dof)
            if( isinstance(kratos_variable,KratosMultiphysics.Array1DVariable3) ):

                integration_method_name   = self._get_integration_method_name(dof)
                vector_integration_method = None

                variable_names = self._get_integration_method_variables(dof)
                variables = []
                for variable in variable_names:
                    variables = variables + [KratosMultiphysics.KratosGlobals.GetVariable(variable)]

                integration_method = None
                if( len(variables) == 4 ):
                    vector_integration_method = getattr(KratosSolid, integration_method_name+'VectorIntegration')
                    integration_method = vector_integration_method(variables[0],variables[1],variables[2],variables[3])
                elif( len(variables) == 1 ):
                    if(integration_method_name.find("Step") != -1):
                        vector_integration_method = getattr(KratosSolid, 'StaticStepVectorIntegration')
                    else:
                        vector_integration_method = getattr(KratosSolid, 'StaticVectorIntegration')
                    integration_method = vector_integration_method(variables[0])
                else:
                    raise Exception('len(variables) = ' + str(len(variables)))

                if(integration_method_name.find("Step") != -1):
                    step_variable_name = 'STEP_'+dof
                    integration_method.SetStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(step_variable_name))

                #print("Integration", integration_method)
                vector_integration_methods.append(integration_method)

            elif( isinstance(kratos_variable,KratosMultiphysics.DoubleVariable) ):

                integration_method_name   = self._get_integration_method_name(dof)
                scalar_integration_method = None

                variable_names = self._get_integration_method_variables(dof)
                variables = []
                for variable in variable_names:
                    variables = variables + [KratosMultiphysics.KratosGlobals.GetVariable(variable)]

                integration_method = None
                if( len(variables) == 4 ):
                    scalar_integration_method = getattr(KratosSolid, integration_method_name+'ScalarIntegration')
                    integration_method = scalar_integration_method(variables[0],variables[1],variables[2],variables[3])
                elif( len(variables) == 1 ):
                    if(integration_method_name.find("Step") != -1):
                        scalar_integration_method = getattr(KratosSolid,'StaticStepScalarIntegration')
                    else:
                        scalar_integration_method = getattr(KratosSolid,'StaticScalarIntegration')
                    integration_method = scalar_integration_method(variables[0])
                else:
                    raise Exception('len(variables) = ' + str(len(variables)))

                if(integration_method_name.find("Step") != -1):
                    step_variable_name = 'STEP_'+dof
                    integration_method.SetStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(step_variable_name))

                #print("Integration", integration_method)
                scalar_integration_methods.append(integration_method)

            else:
                raise Exception("DOF is incorrect. Must be a three-component vector or a scalar.")

        solution_scheme = None
        if(self.settings["solution_type"].GetString() == "Dynamic"):
            if(self.settings["time_integration"].GetString() == "Implicit"):
                if(len(vector_integration_methods)):
                    if(len(scalar_integration_methods)):
                        solution_scheme = KratosSolid.DynamicScheme(vector_integration_methods,scalar_integration_methods)
                    else:
                        solution_scheme = KratosSolid.DynamicScheme(vector_integration_methods)
                elif(len(scalar_integration_methods)):
                    solution_scheme = KratosSolid.DynamicScheme(vector_integration_methods,scalar_integration_methods)
                else:
                    print("WARNING: no integration methods")

        elif(self.settings["solution_type"].GetString() == "Static" or self.settings["solution_type"].GetString() == "Quasi-static"):
            if(len(vector_integration_methods)):
                if(len(scalar_integration_methods)):
                    solution_scheme = KratosSolid.StaticScheme(vector_integration_methods,scalar_integration_methods)
                else:
                    solution_scheme = KratosSolid.StaticScheme(vector_integration_methods)
            elif(len(scalar_integration_methods)):
                solution_scheme = KratosSolid.StaticScheme(vector_integration_methods,scalar_integration_methods)
            else:
                print("WARNING: no integration methods")

        return solution_scheme

    #
    def GetComponentIntegrationMethods(self):

        # set solution scheme and integration method dictionary
        integration_methods = {}

        for dof in self.dofs:

            integration_method_name = self._get_integration_method_name(dof)
            variable_names = self._get_integration_method_variables(dof)

            ##check if variable type is a vector
            kratos_variable = KratosMultiphysics.KratosGlobals.GetVariable(dof)
            if( isinstance(kratos_variable,KratosMultiphysics.Array1DVariable3) ):

                component_integration_method = None

                variable_components = ['_X','_Y','_Z']
                for component in variable_components:
                    variables = []
                    for variable in variable_names:
                        variables = variables + [KratosMultiphysics.KratosGlobals.GetVariable(variable+component)]

                    integration_method = None
                    if( len(variables) == 4 ):
                        component_integration_method = getattr(KratosSolid, integration_method_name+'ComponentIntegration')
                        integration_method = component_integration_method(variables[0],variables[1],variables[2],variables[3])
                    elif( len(variables) == 1 ):
                        if(integration_method_name.find("Step") != -1):
                            component_integration_method = getattr(KratosSolid, 'StaticStepComponentIntegration')
                        else:
                            component_integration_method = getattr(KratosSolid, 'StaticComponentIntegration')
                        integration_method = component_integration_method(variables[0])
                    else:
                        raise Exception('len(variables) = ' + str(len(variables)))

                    if(integration_method_name.find("Step") != -1):
                        step_variable_name = 'STEP_'+dof+component
                        integration_method.SetStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(step_variable_name))

                    integration_methods.update({variables[0].Name(): integration_method})
                    if( len(variables) == 4 ):
                        integration_methods.update({variables[1].Name(): integration_method})
                        integration_methods.update({variables[2].Name(): integration_method})


        #print("Component time integration methods ", integration_methods)

        return integration_methods


    #
    def GetScalarIntegrationMethods(self):

        # set solution scheme and integration method dictionary
        integration_methods = {}

        for dof in self.dofs:

            integration_method_name = self._get_integration_method_name(dof)
            variable_names = self._get_integration_method_variables(dof)

            ##check if variable type is a scalar
            kratos_variable = KratosMultiphysics.KratosGlobals.GetVariable(dof)
            if( isinstance(kratos_variable,KratosMultiphysics.DoubleVariable) ):

                scalar_integration_method = None

                variables = []
                for variable in variable_names:
                    variables = variables + [KratosMultiphysics.KratosGlobals.GetVariable(variable)]

                integration_method = None
                if( len(variables) == 4 ):
                    scalar_integration_method = getattr(KratosSolid, integration_method_name+'ScalarIntegration')
                    integration_method = scalar_integration_method(variables[0],variables[1],variables[2],variables[3])
                elif( len(variables) == 1 ):
                    if(integration_method_name.find("Step") != -1):
                        scalar_integration_method = getattr(KratosSolid, 'StaticStepScalarIntegration')
                    else:
                        scalar_integration_method = getattr(KratosSolid, 'StaticScalarIntegration')
                    integration_method = scalar_integration_method(variables[0])
                else:
                    raise Exception('len(variables) = ' + str(len(variables)))

                if(integration_method_name.find("Step") != -1):
                    step_variable_name = 'STEP_'+dof
                    integration_method.SetStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(step_variable_name))

                integration_methods.update({variables[0].Name(): integration_method})

        return integration_methods

    #
    def GetBufferSize(self):
        buffer_size = self.settings["buffer_size"].GetInt()
        time_integration_order = self.settings["time_integration_order"].GetInt()
        if( buffer_size <= time_integration_order ):
            buffer_size = time_integration_order + 1
        return buffer_size

    #
    def GetTimeIntegrationOrder(self):
        time_integration_order = self.settings["time_integration_order"].GetInt()
        return time_integration_order

    #
    def _get_integration_method_name(self, dof):

        integration_method_name = self.settings["integration_method"].GetString()

        if integration_method_name.find("Step") != -1:
            if dof == 'ROTATION':
                integration_method_name  = integration_method_name+'Rotation'
            elif dof != 'DISPLACEMENT':
                start = 0
                end   = integration_method_name.find("Step")
                integration_method_name = integration_method_name[start:end]
                #print("::[---Scheme Factory--]:: Step Method Changed ("+integration_method_name+":"+dof+")")

        return integration_method_name

    #
    def _get_integration_method_variables(self, dof):

        variables = []
        if(self.settings["solution_type"].GetString() == "Dynamic" ):
            if(dof == 'DISPLACEMENT' or dof == 'VELOCITY' or dof == 'ACCELERATION'):
                variables = variables + ['DISPLACEMENT','VELOCITY','ACCELERATION',dof]
            elif(dof == 'ROTATION' or dof == 'ANGULAR_VELOCITY' or dof == 'ANGULAR_ACCELERATION'):
                variables = variables + ['ROTATION','ANGULAR_VELOCITY','ANGULAR_ACCELERATION',dof]
            elif(dof == 'WATER_DISPLACEMENT'):
                variables = variables + ['WATER_DISPLACEMENT','WATER_VELOCITY','WATER_ACCELERATION',dof]
            elif(dof == 'WATER_PRESSURE'):
                variables = variables + ['WATER_PRESSURE','WATER_PRESSURE_VELOCITY','WATER_PRESSURE_ACCELERATION',dof]
            elif(dof == 'PRESSURE'):
                variables = variables + ['PRESSURE','PRESSURE_VELOCITY','PRESSURE_ACCELERATION',dof]
            elif(dof == 'FLUID_PRESSURE'):
                variables = variables + ['FLUID_PRESSURE','FLUID_PRESSURE_VELOCITY','FLUID_PRESSURE_ACCELERATION',dof]
            else:
                variables = variables + [dof]
        else:
            variables = variables + [dof]

        return variables

    #
    def _set_variables_and_dofs(self):

        # Variables and Dofs settings
        self.nodal_variables = []
        self.dof_variables   = []
        self.dof_reactions   = []
        self.dof_derivatives = []

        # Add displacement variables
        if self._check_input_dof("DISPLACEMENT"):
            self.dof_variables = self.dof_variables + ['DISPLACEMENT']
            self.dof_reactions = self.dof_reactions + ['DISPLACEMENT_REACTION']

            # Add dynamic variables
            if(self.settings["solution_type"].GetString() == "Dynamic"):
                self.dof_derivatives = self.dof_derivatives + ['VELOCITY','ACCELERATION']

        if self._check_input_dof("VELOCITY"):
            # Add specific variables for the problem (velocity dofs)
            self.dof_variables = self.dof_variables + ['VELOCITY']
            self.dof_reactions = self.dof_reactions + ['VELOCITY_REACTION']

            # Add dynamic variables
            self.dof_derivatives = self.dof_derivatives + ['DISPLACEMENT','ACCELERATION']


        # Add rotational variables
        if self._check_input_dof("ROTATION"):
            # Add specific variables for the problem (rotation dofs)
            self.dof_variables = self.dof_variables + ['ROTATION']
            self.dof_reactions = self.dof_reactions + ['ROTATION_REACTION']

            if(self.settings["solution_type"].GetString() == "Dynamic"):
                self.dof_derivatives = self.dof_derivatives + ['ANGULAR_VELOCITY','ANGULAR_ACCELERATION']
            # Add large rotation variables
            self.nodal_variables = self.nodal_variables + ['STEP_DISPLACEMENT','STEP_ROTATION']

        # Add pressure variables
        if self._check_input_dof("PRESSURE"):
            # Add specific variables for the problem (pressure dofs)
            self.dof_variables = self.dof_variables + ['PRESSURE']
            self.dof_reactions = self.dof_reactions + ['PRESSURE_REACTION']

            if(self.settings["solution_type"].GetString() == "Dynamic"):
                self.dof_derivatives = self.dof_derivatives + ['PRESSURE_VELOCITY','PRESSURE_ACCELERATION']

        # Add fluid pressure variables
        if self._check_input_dof("FLUID_PRESSURE"):
            # Add specific variables for the problem (pressure dofs)
            self.dof_variables = self.dof_variables + ['FLUID_PRESSURE']
            self.dof_reactions = self.dof_reactions + ['FLUID_PRESSURE_REACTION']

            self.dof_derivatives = self.dof_derivatives + ['FLUID_PRESSURE_VELOCITY','FLUID_PRESSURE_ACCELERATION']

        # Ass thermal variables
        if self._check_input_dof("TEMPERATURE"):
            self.dof_variables = self.dof_variables + ['TEMPERATURE']
            self.dof_reactions = self.dof_reactions + ['TEMPERATURE_REACTION']

        # Add contat variables
        if self._check_input_dof("LAGRANGE_MULTIPLIER"):
            # Add specific variables for the problem (contact dofs)
            self.dof_variables = self.dof_variables + ['LAGRANGE_MULTIPLIER_NORMAL']
            self.dof_reactions = self.dof_reactions + ['LAGRANGE_MULTIPLIER_NORMAL_REACTION']

        # Add water displacement variables
        if self._check_input_dof("WATER_DISPLACEMENT"):
            self.dof_variables = self.dof_variables + ['WATER_DISPLACEMENT','WATER_VELOCITY','WATER_ACCELERATION']
            self.dof_reactions = self.dof_reactions + ['WATER_DISPLACEMENT_REACTION','NOT_DEFINED','NOT_DEFINED']

        # Add water pressure variables
        if self._check_input_dof("WATER_PRESSURE"):
            self.dof_variables = self.dof_variables + ['WATER_PRESSURE', 'WATER_PRESSURE_VELOCITY','WATER_PRESSURE_ACCELERATION']
            self.dof_reactions = self.dof_reactions + ['REACTION_WATER_PRESSURE', 'NOT_DEFINED', 'NOT_DEFINED']

        # Add jacobian variables
        if self._check_input_dof("JACOBIAN"):
            self.dof_variables = self.dof_variables + ['JACOBIAN']
            self.dof_reactions = self.dof_reactions + ['REACTION_JACOBIAN']


    def _check_input_dof(self, variable):

        for i in self.dofs:
            if i == variable:
                return True

        #temporary default
        if(variable == 'DISPLACEMENT' and len(self.dofs) == 0 ):
            return True

        return False

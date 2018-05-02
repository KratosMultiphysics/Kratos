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

    def GetSolutionScheme(self):

        scheme_integration_methods = []

        for dof in self.dofs:

            integration_method_name   = self._get_integration_method_name(dof)
            vector_integration_method = getattr(KratosSolid, integration_method_name+'VectorIntegration')

            variables = self._get_integration_method_variables(dof)
            scheme_integration_method = None
            if( len(variables) == 4 ):
                scheme_integration_method = vector_integration_method(variables[0],variables[1],variables[2],variables[3])
            elif( len(variables) == 1 ):
                scheme_integration_method = vector_integration_method(variables[0])
            else:
                raise Exception('len(variables) = ' + str(len(variables)))
            
            if(integration_method_name.find("Step") != -1):
                step_variable_name = 'STEP_'+dof
                scheme_integration_method.SetStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(step_variable_name))

            print("Integration", scheme_integration_method)
            scheme_integration_methods.append(scheme_integration_method)
            
        solution_scheme = None
        if(self.settings["solution_type"].GetString() == "Dynamic"):

            if(self.settings["time_integration"].GetString() == "Implicit"):
                solution_scheme = KratosSolid.DynamicScheme(scheme_integration_methods)

        elif(self.settings["solution_type"].GetString() == "Static" or self.settings["solution_type"].GetString() == "Quasi-static"):
            solution_scheme = KratosSolid.StaticScheme(scheme_integration_methods)

        return solution_scheme

    #
    def GetIntegrationMethods(self):

        # set solution scheme and integration method dictionary
        integration_methods = {}

        for dof in self.dofs:

            integration_method_name      = self._get_integration_method_name(dof)
            component_integration_method = getattr(KratosSolid, integration_method_name+'ComponentIntegration')
            
            integration_methods.update({dof: component_integration_method()})
           
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

        if(dof == 'ROTATION' and integration_method_name.find("Step") != -1):
            integration_method_name  = integration_method_name+'Rotation'

        return integration_method_name

    #
    def _get_integration_method_variables(self, dof):

        variables = []
        if(self.settings["solution_type"].GetString() == "Dynamic" ):
            if(dof == 'DISPLACEMENT' or dof == 'VELOCITY' or dof == 'ACCELERATION'):
                variables = [KratosMultiphysics.DISPLACEMENT,KratosMultiphysics.VELOCITY,KratosMultiphysics.ACCELERATION]
                if(dof == 'DISPLACEMENT'):
                    variables = variables + [KratosMultiphysics.DISPLACEMENT]
                elif(dof == 'VELOCITY'):
                    variables = variables + [KratosMultiphysics.VELOCITY]
                elif(dof == 'ACCELERATION'):
                    variables = variables + [KratosMultiphysics.ACCELERATION]

            elif(dof == 'ROTATION' or dof == 'ANGULAR_VELOCITY' or dof == 'ANGULAR_ACCELERATION'):
                variables = [KratosMultiphysics.ROTATION,KratosMultiphysics.ANGULAR_VELOCITY,KratosMultiphysics.ANGULAR_ACCELERATION]
                if(dof == 'ROTATION'):
                    variables = variables + [KratosMultiphysics.ROTATION]
                elif(dof == 'ANGULAR_VELOCITY'):
                    variables = variables + [KratosMultiphysics.ANGULAR_VELOCITY]
                elif(dof == 'ANGULAR_ACCELERATION'):
                    variables = variables + [KratosMultiphysics.ANGULAR_ACCELERATION]
            else:
                variables = [KratosMultiphysics.KratosGlobals.GetVariable(dof)]
        else:
            variables = [KratosMultiphysics.KratosGlobals.GetVariable(dof)]

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
            self.dof_reactions = self.dof_reactions + ['REACTION']

            # Add dynamic variables
            self.dof_derivatives = self.dof_derivatives + ['VELOCITY','ACCELERATION']

        if self._check_input_dof("VELOCITY"):
            # Add specific variables for the problem (velocity dofs)
            self.dof_variables = self.dof_variables + ['VELOCITY']
            self.dof_reactions = self.dof_reactions + ['NOT_DEFINED']

            # Add dynamic variables
            self.dof_derivatives = self.dof_derivatives + ['DISPLACEMENT','ACCELERATION']


        # Add rotational variables
        if self._check_input_dof("ROTATION"):
            # Add specific variables for the problem (rotation dofs)
            self.dof_variables = self.dof_variables + ['ROTATION']
            self.dof_reactions = self.dof_reactions + ['TORQUE']

            self.dof_derivatives = self.dof_derivatives + ['ANGULAR_VELOCITY','ANGULAR_ACCELERATION']
            # Add large rotation variables
            self.nodal_variables = self.nodal_variables + ['STEP_DISPLACEMENT','STEP_ROTATION','DELTA_ROTATION']

        # Add pressure variables
        if self._check_input_dof("PRESSURE"):
            # Add specific variables for the problem (pressure dofs)
            self.dof_variables = self.dof_variables + ['PRESSURE']
            self.dof_reactions = self.dof_reactions + ['PRESSURE_REACTION']

        # Add contat variables
        if self._check_input_dof("LAGRANGE_MULTIPLIER"):
            # Add specific variables for the problem (contact dofs)
            self.dof_variables = self.dof_variables + ['LAGRANGE_MULTIPLIER_NORMAL']
            self.dof_reactions = self.dof_reactions + ['LAGRANGE_MULTIPLIER_NORMAL_REACTION']

        # Add water displacement variables
        if self._check_input_dof("WATER_DISPLACEMENT"):
            self.dof_variables = self.dof_variables + ['WATER_DISPLACEMENT','WATER_VELOCITY','WATER_ACCELERATION']
            self.dof_reactions = self.dof_reactions + ['WATER_DISPLACEMENT_REACTION','WATER_VELOCITY_REACTION','WATER_ACCELERATION_REACTION']

        # Add water pressure variables
        if self._check_input_dof("WATER_PRESSURE"):
            self.dof_variables = self.dof_variables + ['WATER_PRESSURE', 'WATER_PRESSURE_VELOCITY','WATER_PRESSURE_ACCELERATION']
            self.dof_reactions = self.dof_reactions + ['REACTION_WATER_PRESSURE', 'WATER_PRESSURE_VELOCITY_REACTION', 'WATER_PRESSURE_ACCELERATION_REACTION']

        # Add jacobian variables
        if self._check_input_dof("JACOBIAN"):
            self.dof_variables = self.dof_variables + ['JACOBIAN']
            self.dof_reactions = self.dof_reactions + ['REACTION_JACOBIAN']


    def _check_input_dof(self, variable):

        for i in self.dofs:
            if i == variable:
                return True

        #temporary default
        if(variable == 'DISPLACEMENT' and len(dofs) == 0 ):
            return True

        return False

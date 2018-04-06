from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_solver as BaseSolver

def CreateSolver(custom_settings):
    return ImplicitMechanicalSolver(custom_settings)

class ImplicitMechanicalSolver(BaseSolver.MechanicalSolver):
    """The solid mechanics implicit dynamic solver.

    This class creates the mechanical solvers for implicit dynamic analysis.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    See solid_mechanics_solver.py for more information.
    """
    def __init__(self, custom_settings):

        # Set defaults and validate custom settings.
        ##TODO : solving_strategy_settings must be time_integration_settings (GiD interface changes needed)
        implicit_solver_settings = KratosMultiphysics.Parameters("""
        {
            "solving_strategy_settings":{
                "bossak_factor" :-0.3,
                "dynamic_factor": 1.0,
                "lumped_mass_matrix" : true,
                "consistent_mass_matrix" : false,
                "rayleigh_damping": false,
                "rayleigh_alpha": 0.0,
                "rayleigh_beta" : 0.0
            }
        }
        """)

        # Validate and transfer settings
        if( custom_settings.Has("solving_strategy_settings") ):
            self._validate_and_transfer_matching_settings(custom_settings["solving_strategy_settings"], implicit_solver_settings["solving_strategy_settings"])
        self.implicit_solver_settings = implicit_solver_settings["solving_strategy_settings"]

        # Construct the base solver.
        super(ImplicitMechanicalSolver, self).__init__(custom_settings)

        print("::[Implicit_Scheme]:: "+self.time_integration_settings["integration_method"].GetString()+" Scheme Ready")


    def GetVariables(self):

        nodal_variables = super(ImplicitMechanicalSolver, self).GetVariables()

        return nodal_variables

    #### Solver internal methods ####

    def _create_solution_scheme(self):

        integration_method = self.time_integration_settings["integration_method"].GetString()

        if( self.implicit_solver_settings["rayleigh_damping"].GetBool() == True ):
            self.process_info[KratosSolid.RAYLEIGH_ALPHA] = self.implicit_solver_settings["rayleigh_alpha"].GetDouble()
            self.process_info[KratosSolid.RAYLEIGH_BETA]  = self.implicit_solver_settings["rayleigh_beta"].GetDouble()
        else:
            self.process_info[KratosSolid.RAYLEIGH_ALPHA] = 0.0
            self.process_info[KratosSolid.RAYLEIGH_BETA]  = 0.0

        # compute mass lumped matrix
        if( self.implicit_solver_settings["lumped_mass_matrix"].GetBool() == True ):
            self.process_info[KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX] = True
        else:
            # compute consistent dynamic tangent/mass matrix
            if( self.implicit_solver_settings["consistent_mass_matrix"].GetBool() == True ):
                self.process_info[KratosMultiphysics.COMPUTE_DYNAMIC_TANGENT] = True

        # set bossak factor
        bossak_factor = self.implicit_solver_settings["bossak_factor"].GetDouble()
        self.process_info[KratosMultiphysics.BOSSAK_ALPHA] = bossak_factor;

        # set solution scheme and integration method dictionary
        self.integration_methods = {}
        if(integration_method == "Newmark"):
            self.integration_methods.update({'DISPLACEMENT': KratosSolid.NewmarkMethod(),
                                             'ROTATION': KratosSolid.NewmarkMethod()}) #shells
            mechanical_scheme = KratosSolid.DisplacementNewmarkScheme()
        elif(integration_method == "Bossak"):
            self.integration_methods.update({'DISPLACEMENT': KratosSolid.BossakMethod(),
                                             'ROTATION': KratosSolid.BossakMethod()}) #shells
            mechanical_scheme = KratosSolid.DisplacementBossakScheme()
        elif(integration_method == "Simo"):
            self.integration_methods.update({'DISPLACEMENT': KratosSolid.SimoMethod(),
                                             'ROTATION': KratosSolid.SimoMethod()}) #shells
            mechanical_scheme = KratosSolid.DisplacementSimoScheme()
        elif(integration_method == "RotationNewmark"):
            self.integration_methods.update({'DISPLACEMENT': KratosSolid.NewmarkStepMethod(),
                                             'ROTATION': KratosSolid.NewmarkStepRotationMethod()}) #beams
            mechanical_scheme = KratosSolid.DisplacementRotationNewmarkScheme()
        elif(integration_method == "RotationBossak"):
            self.integration_methods.update({'DISPLACEMENT': KratosSolid.BossakStepMethod(),
                                             'ROTATION': KratosSolid.BossakStepRotationMethod()}) #beams
            mechanical_scheme = KratosSolid.DisplacementRotationBossakScheme()
        elif(integration_method == "RotationSimo"):
            self.integration_methods.update({'DISPLACEMENT': KratosSolid.SimoStepMethod(),
                                             'ROTATION': KratosSolid.SimoStepRotationMethod()}) #beams
            mechanical_scheme = KratosSolid.DisplacementRotationSimoScheme()
        elif(integration_method == "RotationEMC"):
            self.integration_methods.update({'DISPLACEMENT': KratosSolid.EmcStepMethod(),
                                             'ROTATION': KratosSolid.EmcStepRotationMethod()}) #shells
            mechanical_scheme = KratosSolid.DisplacementRotationEmcScheme()
        else:
            raise Exception("Unsupported integration_method: " + integration_method)

        # set integration parameters
        self._set_time_integration_methods()

        return mechanical_scheme

    def _set_time_integration_methods(self):

        # assign an default integration method (static) for all not set dofs
        for i in range(0, self.settings["dofs"].size() ):
            if( not (self.settings["dofs"][i].GetString() in self.integration_methods.keys()) ):
                self.integration_methods.update({self.settings["dofs"][i].GetString() : KratosSolid.StaticMethod()})

        # first: calculate parameters (only once permitted) (set is included)
        #self.integration_methods['DISPLACEMENT'].CalculateParameters(self.process_info) #calculate
        # second: for the same method the parameters (already calculated)
        #self.integration_methods['ROTATION'].SetParameters(self.process_info) #set parameters
        # third... somehow it must be a integration method for each dof supplied in "dofs"=[]

        # calculate method parameters and set to process info internally
        main_dof = next(iter(self.integration_methods))
        self.integration_methods[main_dof].CalculateParameters(self.process_info)

        # add to integration methods container and set to process_info for processes acces
        integration_methods_container = KratosSolid.TimeIntegrationMethodsContainer()
        for dof, method in self.integration_methods.items():
            method.SetParameters(self.process_info) #set same parameters to all methods from process_info values
            integration_methods_container.Set(dof,method)
        integration_methods_container.AddToProcessInfo(KratosSolid.TIME_INTEGRATION_METHODS, integration_methods_container, self.process_info)

        #print(integration_methods_container)

    def _create_mechanical_solver(self):
        if(self.solving_strategy_settings["line_search"].GetBool() == True):
            mechanical_solver = self._create_line_search_strategy()
        else:
            mechanical_solver = self._create_newton_raphson_strategy()
        return mechanical_solver


    def _create_line_search_strategy(self):
        mechanical_scheme = self._get_solution_scheme()
        #linear_solver = self._get_linear_solver()
        convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()

        options = KratosMultiphysics.Flags()
        options.Set(KratosSolid.SolverLocalFlags.COMPUTE_REACTIONS, self.solving_strategy_settings["compute_reactions"].GetBool())
        options.Set(KratosSolid.SolverLocalFlags.REFORM_DOFS, self.solving_strategy_settings["reform_dofs_at_each_step"].GetBool())
        options.Set(KratosSolid.SolverLocalFlags.MOVE_MESH, self.solving_strategy_settings["move_mesh_flag"].GetBool())

        return KratosSolid.LineSearchStrategy(self.model_part, mechanical_scheme, builder_and_solver, convergence_criterion,
                                              options, self.solving_strategy_settings["max_iteration"].GetInt())

        #return KratosSolid.ResidualBasedNewtonRaphsonLineSearchStrategy(self.model_part,
        #                                                                mechanical_scheme,
        #                                                                linear_solver,
        #                                                                convergence_criterion,
        #                                                                builder_and_solver,
        #                                                                self.solving_strategy_settings["max_iteration"].GetInt(),
        #                                                                self.solving_strategy_settings["compute_reactions"].GetBool(),
        #                                                                self.solving_strategy_settings["reform_dofs_at_each_step"].GetBool(),
        #                                                                self.solving_strategy_settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_strategy(self):
        mechanical_scheme = self._get_solution_scheme()
        #linear_solver = self._get_linear_solver()
        convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()

        options = KratosMultiphysics.Flags()
        options.Set(KratosSolid.SolverLocalFlags.COMPUTE_REACTIONS, self.solving_strategy_settings["compute_reactions"].GetBool())
        options.Set(KratosSolid.SolverLocalFlags.REFORM_DOFS, self.solving_strategy_settings["reform_dofs_at_each_step"].GetBool())
        options.Set(KratosSolid.SolverLocalFlags.IMPLEX, self.solving_strategy_settings["implex"].GetBool())

        return KratosSolid.NewtonRaphsonStrategy(self.model_part, mechanical_scheme, builder_and_solver, convergence_criterion,
                                                 options, self.solving_strategy_settings["max_iteration"].GetInt())


        #return KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.model_part,
        #                                                             mechanical_scheme,
        #                                                             linear_solver,
        #                                                             convergence_criterion,
        #                                                             builder_and_solver,
        #                                                             self.solving_strategy_settings["max_iteration"].GetInt(),
        #                                                             self.solving_strategy_settings["compute_reactions"].GetBool(),
        #                                                             self.solving_strategy_settings["reform_dofs_at_each_step"].GetBool(),
        #                                                             self.solving_strategy_settings["move_mesh_flag"].GetBool())

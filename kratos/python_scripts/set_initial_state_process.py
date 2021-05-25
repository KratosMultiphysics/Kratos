import KratosMultiphysics


def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return SetInitialStateProcess(Model, settings["Parameters"])


def param_to_str(param):
    """Convert input to string type.

    Keyword arguments:
    param -- Accepted types are Parameters double, Parameters string, int, float, str.

    Returns:
    string expression of the input
    """
    if type(param) is KratosMultiphysics.Parameters:
        if param.IsNumber():
            sparam = str(param.GetDouble())
        elif param.IsString():
            sparam = param.GetString()
    elif type(param) in [int, float, str]:
        sparam = str(param)
    else:
        msg = "SetInitialStateProcess: Value must be scalar or string"
        raise Exception(msg)
    return sparam


def str_to_function(expr):
    """Convert text to a mathematical function.

    Input string must be a valid mathematical expression (e.g."1.2", "3.4 + t",
    "4.5 * t"). It returns generated functions (e.g. f(1.2), f(3.4+t), f(4.5*t)),
    where 't' will be later evaluated as the current simulation time.

    Keyword arguments:
    expr -- strings with valid expression

    Returns:
    Generated functions rendered from the input
    """
    function = KratosMultiphysics.GenericFunctionUtility(expr)
    if function.DependsOnSpace():
        msg = "SetInitialStateProcess: 'x', 'y' and 'z' variables not supported yet"
        raise Exception(msg)
    return function


class SetInitialStateProcess(KratosMultiphysics.Process):
    """This process sets a given value for a certain flag in all the nodes of a submodelpart.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings):
        """The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        # The value can be a double or a string (function)
        default_settings = KratosMultiphysics.Parameters(
            """
        {
            "help"            : "This sets the initial conditions in terms of imposed strain, stress or deformation gradient",
            "mesh_id"         : 0,
            "model_part_name" : "please_specify_model_part_name",
            "dimension"       : 3,
            "imposed_strain_magnitude": "1.0",
            "imposed_strain"  : ["0","0","0",0,0,0],
            "imposed_stress_magnitude": "1.0",
            "imposed_stress"  : [0,0,0,0,0,0],
            "imposed_deformation_gradient": [[1,0,0],[0,1,0],[0,0,1]],
            "interval"        : [0.0, 1e30]
        }
        """
        )

        # assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.dimension = settings["dimension"].GetInt()

        # init strain
        m = param_to_str(settings["imposed_strain_magnitude"])
        v = [param_to_str(x) for x in settings["imposed_strain"]]
        self.strain_functions = [str_to_function("{}*{}".format(m, x)) for x in v]
        nr_comps = len(self.strain_functions)
        self.imposed_strain = KratosMultiphysics.Vector(nr_comps)

        # init stress
        m = param_to_str(settings["imposed_stress_magnitude"])
        v = [param_to_str(x) for x in settings["imposed_stress"]]
        self.stress_functions = [str_to_function("{}*{}".format(m, x)) for x in v]
        nr_comps = len(self.stress_functions)
        self.imposed_stress = KratosMultiphysics.Vector(nr_comps)

        # init deformation gradient
        aux_matrix = settings["imposed_deformation_gradient"]
        self.deformation_functions = []
        for row in aux_matrix:
            v = [param_to_str(x) for x in row]
            aux_vector = [str_to_function(x) for x in v]
            self.deformation_functions.append(aux_vector)
        nrows = aux_matrix.size()
        ncols = aux_matrix[0].size()
        self.imposed_deformation_gradient = KratosMultiphysics.Matrix(nrows, ncols)

    def ExecuteInitializeSolutionStep(self):
        """This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval.IsInInterval(current_time):

            # strain
            nr_comps = len(self.strain_functions)
            for i in range(nr_comps):
                self.imposed_strain[i] = self.strain_functions[i].CallFunction(
                    0, 0, 0, current_time, 0, 0, 0
                )

            # stress
            nr_comps = len(self.stress_functions)
            for i in range(nr_comps):
                self.imposed_stress[i] = self.stress_functions[i].CallFunction(
                    0, 0, 0, current_time, 0, 0, 0
                )

            # deformation gradient
            nr_rows = len(self.deformation_functions)
            nr_cols = len(self.deformation_functions[0])
            for r in range(nr_rows):
                for c in range(nr_cols):
                    self.imposed_deformation_gradient[
                        r, c
                    ] = self.deformation_functions[r][c].CallFunction(
                        0, 0, 0, current_time, 0, 0, 0
                    )

            self.SetInitialState()

    def SetInitialState(self):
        """This method creates the c++ process and sets each entity

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if self.dimension == 3:
            KratosMultiphysics.SetInitialStateProcess3D(
                self.model_part,
                self.imposed_strain,
                self.imposed_stress,
                self.imposed_deformation_gradient,
            ).ExecuteInitializeSolutionStep()
        else:  # 2D case
            KratosMultiphysics.SetInitialStateProcess2D(
                self.model_part,
                self.imposed_strain,
                self.imposed_stress,
                self.imposed_deformation_gradient,
            ).ExecuteInitializeSolutionStep()

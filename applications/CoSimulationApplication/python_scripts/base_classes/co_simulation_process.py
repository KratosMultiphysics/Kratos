import KratosMultiphysics

def Factory(settings, Model):
    return CoSimulationProcess(Model, settings["Parameters"])

class CoSimulationProcess(object):
    """
    This class defines the base CoSimulationProcess class. The rest of CoSimulationProcess must be defined from this class
    """
    def __init__(self, Model, settings):
        """ Constructor of CoSimulationProcess object
        
        It is intended to be called from the constructor
        of deriving classes:
        super().__init__(Model, settings)

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- The Model to be used
        settings -- The cosimulation process settings
        """
        self.Model = Model
        settings.AddMissingParameters(self.GetDefaultParameters())
        self.settings = settings

    def ExecuteInitialize(self):
        """This function it is supposed to not be called from outside. It is defined in case by error is called as an standard process
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        KratosMultiphysics.Logger.PrintWarning("CoSimulationProcess", "You are calling a CoSimulationProcess, these methods are intended to do nothing")

    def ExecuteBeforeCoSimulationInitialize(self):
        """This function is executed before ExecuteInitialize of the cosimulated solvers
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteAfterCoSimulationInitialize(self):
        """This function is executed after ExecuteInitialize of the cosimulated solvers
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteBeforeSolutionLoop(self):
        """This function it is supposed to not be called from outside. It is defined in case by error is called as an standard process
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        KratosMultiphysics.Logger.PrintWarning("CoSimulationProcess", "You are calling a CoSimulationProcess, these methods are intended to do nothing")

    def ExecuteBeforeCoSimulationBeforeSolutionLoop(self):
        """This function is executed before ExecuteBeforeSolutionLoop of the cosimulated solvers
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteAfterCoSimulationBeforeSolutionLoop(self):
        """This function is executed after ExecuteBeforeSolutionLoop of the cosimulated solvers
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteInitializeSolutionStep(self):
        """This function it is supposed to not be called from outside. It is defined in case by error is called as an standard process
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        KratosMultiphysics.Logger.PrintWarning("CoSimulationProcess", "You are calling a CoSimulationProcess, these methods are intended to do nothing")

    def ExecuteBeforeCoSimulationInitializeSolutionStep(self):
        """This function is executed before ExecuteInitializeSolutionStep of the cosimulated solvers
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteAfterCoSimulationInitializeSolutionStep(self):
        """This function is executed after ExecuteInitializeSolutionStep of the cosimulated solvers
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalizeSolutionStep(self):
        """This function it is supposed to not be called from outside. It is defined in case by error is called as an standard process
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        KratosMultiphysics.Logger.PrintWarning("CoSimulationProcess", "You are calling a CoSimulationProcess, these methods are intended to do nothing")

    def ExecuteBeforeCoSimulationFinalizeSolutionStep(self):
        """This function is executed before ExecuteFinalizeSolutionStep of the cosimulated solvers
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteAfterCoSimulationFinalizeSolutionStep(self):
        """This function is executed after ExecuteFinalizeSolutionStep of the cosimulated solvers
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteBeforeOutputStep(self):
        """This function it is supposed to not be called from outside. It is defined in case by error is called as an standard process
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        KratosMultiphysics.Logger.PrintWarning("CoSimulationProcess", "You are calling a CoSimulationProcess, these methods are intended to do nothing")

    def ExecuteAfterOutputStep(self):
        """This function it is supposed to not be called from outside. It is defined in case by error is called as an standard process
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        KratosMultiphysics.Logger.PrintWarning("CoSimulationProcess", "You are calling a CoSimulationProcess, these methods are intended to do nothing")

    def ExecuteBeforeCoSimulationOutputSolutionStep(self):
        """This function is executed before OutputSolutionStep of the cosimulated solvers
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteAfterCoSimulationOutputSolutionStep(self):
        """This function is executed after OutputSolutionStep of the cosimulated solvers
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalize(self):
        """This function it is supposed to not be called from outside. It is defined in case by error is called as an standard process
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        KratosMultiphysics.Logger.PrintWarning("CoSimulationProcess", "You are calling a CoSimulationProcess, these methods are intended to do nothing")

    def ExecuteBeforeCoSimulationExecuteFinalize(self):
        """This function is executed before ExecuteFinalize of the cosimulated solvers
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteAfterCoSimulationExecuteFinalize(self):
        """This function is executed after ExecuteFinalize of the cosimulated solvers
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def Check(self):
        """This function checks the CoSimulationProcess
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def Clear(self):
        """This function clears the CoSimulationProcess
        
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def GetDefaultParameters(self):
        """This function provides the default settings of the CoSimulationProcess

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        return KratosMultiphysics.Parameters(""" {} """)

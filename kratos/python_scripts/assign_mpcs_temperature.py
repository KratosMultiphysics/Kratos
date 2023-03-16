import KratosMultiphysics
import numpy as np
import KratosMultiphysics.MeshingApplication as KratosMA

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignMPCsTemperature(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignMPCsTemperature(KratosMultiphysics.Process):
    """This class is a dummy-process that shows how the functions that can be implemented
    in order to customize the behavior

    Public member variables:
    model -- the container of the different model parts.
    settings -- Kratos parameters containing process settings.
    """

    def __init__(self, model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        KratosMultiphysics.Process.__init__(self) # calling the baseclass constructor

        default_settings = KratosMultiphysics.Parameters("""{
            "computing_model_part_name" : "ComputingModelPartName",
            "slave_model_part_name" : "SlaveModelPartName",
            "master_model_part_name" : "MasterModelPartName",
            "average_temperature": false
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        settings.ValidateAndAssignDefaults(default_settings)
        
        self.model = model
        self.settings = settings

    def ExecuteInitialize(self):
        self.computing_model_part_name = self.settings["computing_model_part_name"].GetString()
        self.slave_model_part_name = self.settings["slave_model_part_name"].GetString()
        self.master_model_part_name = self.settings["master_model_part_name"].GetString()
        self.computing_model_part = self.model.GetModelPart(self.computing_model_part_name)
        self.slave_model_part = self.model.GetModelPart(self.slave_model_part_name)
        self.master_model_part = self.model.GetModelPart(self.master_model_part_name)
        self.assign_mpcs_utility = KratosMultiphysics.AssignMPCsToNeighboursUtility(self.master_model_part.Nodes)
        average_temperature_flag = self.settings["average_temperature"].GetBool()
        if average_temperature_flag:
            self.NewMultiPointConstraint()
        else:
            self.AssignMPCsToNodes()
    
    def AssignMPCsToNodes(self):
        tol = 0.01
        self.assign_mpcs_utility.AssignMPCsToNodes(self.slave_model_part.Nodes,tol,self.computing_model_part, KratosMultiphysics.TEMPERATURE)
    
    def NewMultiPointConstraint(self):
        number_of_master_nodes = self.slave_model_part.NumberOfNodes()
        weights_vector = KratosMultiphysics.Matrix(1,number_of_master_nodes)
        weight = 1/number_of_master_nodes
        slave_dof_container = KratosMultiphysics.Vector(1)
        master_dofs_container = KratosMultiphysics.Vector(number_of_master_nodes)
        master_dofs_container = []
        constant_vector = KratosMultiphysics.Vector(number_of_master_nodes)
        i = 0
        for node in self.slave_model_part.Nodes:
            constant_vector[i] = 0
            weights_vector[0,i] = weight
            master_dofs_container.append(node.GetDof(KratosMultiphysics.TEMPERATURE))
            i += 1
        counter = 1
        for node in self.master_model_part.Nodes:
            slave_dof_container = []
            slave_dof_container.append(node.GetDof(KratosMultiphysics.TEMPERATURE))
            self.computing_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", counter, master_dofs_container, slave_dof_container, weights_vector, constant_vector)
            counter += 1
        

    def Check(self):
        """ This method verifies that the input is correct

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed just before the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteBeforeOutputStep(self):
        """ This method is executed before writing the output (if output
        is being written in this step)

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteAfterOutputStep(self):
        """ This method is executed after writing the output (if output
        is being written in this step)

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalize(self):
        """ This method is executed after the computations, at the end of the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass
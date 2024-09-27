import KratosMultiphysics as KM
import numpy as np
from pathlib import Path
import KratosMultiphysics.time_based_ascii_file_writer_utility as AsciiWriter

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return LocalCoordinateSystemDefinition(Model, settings["Parameters"])

class LocalCoordinateSystemDefinition(KM.Process):

    """
    This process sets a local coordinate system for each element in the model_part. 
    The Z axis of the local system is aligned with the through-thickness direction 
    (calculated using the normal vector of the surface elements), and the X axis is set to (1, 0, 0). 
    The Y axis is computed as the cross product of Z and X.


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
        KM.Process.__init__(self)

        # The value can be a double or a string (function)
        default_settings = KM.Parameters("""{
            "model_part_name"              : "please_specify_model_part_name"
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        # Initialize model_part  
        self.model_part = Model[settings["model_part_name"].GetString()]

        # Normal calculation utility
        self.normal_calculation_utils = KM.NormalCalculationUtils()
        


    def ExecuteInitialize(self):
        """This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.normal_calculation_utils.CalculateUnitNormalsNonHistorical(self.model_part, True, KM.NORMAL)
        
        # We go thorough all elements to set the local axis
        for element in self.model_part.Elements:
            for node in element.GetGeometry():
                if node.Has(KM.NORMAL):
                    local_z = node.GetValue(KM.NORMAL)
                    local_z = np.array(local_z)
                    local_z /= np.linalg.norm(local_z)

                    global_x = np.array([1.0, 0.0, 0.0])
                    
                    dot_product = np.dot(global_x, local_z)
                    if abs(dot_product) > 0.999:  # Nearly collinear case (dot product close to ±1)
                        # If collinear, switch to an alternative axis, e.g., [0, 1, 0]
                        alternative_x = np.array([0.0, 1.0, 0.0])
                        projection_x = alternative_x - np.dot(alternative_x, local_z) * local_z
                    else:
                        # If not collinear, project global X onto the plane orthogonal to Z
                        projection_x = global_x - np.dot(global_x, local_z) * local_z
                    
                    local_x = projection_x / np.linalg.norm(projection_x)

                    local_y = np.cross(local_z, local_x)
                    local_y /= np.linalg.norm(local_y)
                    
            element.SetValue(KM.LOCAL_AXIS_1, local_x)  # X direction
            element.SetValue(KM.LOCAL_AXIS_2, local_y)  # Y direction
            element.SetValue(KM.LOCAL_AXIS_3, local_z)  # Z direction (through-thickness)
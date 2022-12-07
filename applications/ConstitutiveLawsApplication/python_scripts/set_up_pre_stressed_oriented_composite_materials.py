import KratosMultiphysics as KM
import KratosMultiphysics.ConstitutiveLawsApplication as CLApp
from pathlib import Path
import math

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string"
        )
    return SetUpPreStressedOrientedCompositeMaterials(Model, settings["Parameters"])

class SetUpPreStressedOrientedCompositeMaterials(KM.Process):

    """This process sets a proper orientation of the local axes of the elements intersected by line elements (steel tendons). Besides it also computes and sets a volumetric participation of steel within the concrete FE as well as an indicated pre-stressing strain. It also creates a submodelpart for each steel tendon intersected FE.

    Format of the file:
    # Intersection points (D is the diameter of the tendon and Ep the imposed pre-stressing strain):

    Begin Tendon - Hexahedra intersection: Tendon_Inf     D = 1.0e-3    Ep = 0.00001
	    4807	0.01000	     0.05000	   0.01500	    0.02000	     0.05000   0.01500
        ...
    End Tendon - Hexahedra intersection

    # FE inside each tendon:

    Begin Tendon - Hexahedra: Tendon_Inf
        4807    4812    4817    4822    4827    4832    4837    4842    4847    4852   ...
    End Tendon - Hexahedra


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
        default_settings = KM.Parameters(
            """
        {
            "help" : "This sets the initial conditions in terms of imposed pre-stressing strain, local axes and % participation of fiber",
            "model_part_name" : "please_specify_model_part_name",
            "intersection_file_name" : "please_include_directory_and_full_name_with_extension",
            "minimum_fiber_participation_threshold" : 1.0e-7,
            "local_axis_colineal_tolerance" : 1.0e-7,
            "echo_level" : 0
        }
        """
        )
        settings.ValidateAndAssignDefaults(default_settings)

        self.intersection_file_name = settings["intersection_file_name"].GetString()
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.fiber_participation_threshold = settings["minimum_fiber_participation_threshold"].GetDouble()
        self.local_axis_colineal_tol = settings["local_axis_colineal_tolerance"].GetDouble()
        self.echo = settings["echo_level"].GetInt()

    def ExecuteInitializeSolutionStep(self):
        """This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if self.model_part.ProcessInfo[KM.STEP] == 1:
            if self.echo > 0:
                KM.Logger.PrintInfo("SetUpPreStressedOrientedCompositeMaterials ", "Reading intersections file " + self.intersection_file_name + "...")
            intersections_file = open(self.intersection_file_name, "r")
            lines = intersections_file.readlines()

            phi = 0.0; Ep = 0.0 # Diameter and pre-stressing strain
            intersection_block = False
            for line in lines:
                # We loop over the whole file...
                split_line = line.split()
                if len(split_line) > 0: # We skip empty lines
                    if line.find("intersection:") != -1: # a new tendon intersection block starts
                        intersection_block = True
                        tendon_name = split_line[5]
                        # We create a submodel for each tendon
                        self.model_part.CreateSubModelPart(tendon_name)
                        phi = float(split_line[8])   # Diameter of the tendon
                        Ep  = float(split_line[11])  # Imposed pre-stressing strain
                        if self.echo > 0:
                            KM.Logger.PrintInfo("SetUpPreStressedOrientedCompositeMaterials", "Reading block of ", tendon_name + " with and imposed strain of " + str(Ep) + " and a diameter of " + str(phi))
                    elif line.find("End") != -1 and line.find("intersection") != -1:
                        intersection_block = False
                    
                    if intersection_block and line.find("Begin") == -1:
                        id_elem = int(split_line[0])
                        elem = self.model_part.GetElement(id_elem)
                        self.model_part.GetSubModelPart(tendon_name).AddElement(elem, 0) # We add the element to the tendon submodelpart

                        # Here we apply the imposed strain
                        elem.Initialize(self.model_part.ProcessInfo) # necessary to initialize the element first...
                        array_bool = elem.CalculateOnIntegrationPoints(CLApp.IS_PRESTRESSED, self.model_part.ProcessInfo)
                        bool_vect = []
                        for index in range(len(array_bool)):
                            bool_vect.append(True)
                        elem.SetValue(CLApp.SERIAL_PARALLEL_IMPOSED_STRAIN, Ep)
                        elem.SetValuesOnIntegrationPoints(CLApp.IS_PRESTRESSED, bool_vect, self.model_part.ProcessInfo)

                        # Here we set a proper volumetric participation of the fiber
                        elem_volume = elem.GetGeometry().DomainSize()
                        intersection_vector = self.CalculateIntersectionVectors(split_line)
                        tendon_volume = self.Norm2(intersection_vector) * math.pi * phi**2 / 4
                        kf = tendon_volume / (elem_volume)
                        kf = self.CheckFiberParticipation(kf)
                        kf_vect = []
                        for index in range(len(array_bool)):
                            kf_vect.append(kf)
                        elem.SetValuesOnIntegrationPoints(CLApp.FIBER_VOLUMETRIC_PARTICIPATION, kf_vect, self.model_part.ProcessInfo)

                        # Here we set the local axes... (Gram–Schmidt)
                        axis_1 = self.CalculateIntersectionVectors(split_line)
                        if self.Norm2(axis_1) > 0.0:
                            axis_1 = axis_1 / self.Norm2(axis_1)
                            elem.SetValue(KM.LOCAL_AXIS_1, axis_1)
                            axis_2 = KM.Vector(3)

                            axis_2[0] = 0.0
                            axis_2[1] = 1.0
                            axis_2[2] = 0.0

                            if self.Norm2(axis_2-axis_1) <= self.local_axis_colineal_tol:
                                axis_2[0] = 0.0
                                axis_2[1] = 0.0
                                axis_2[2] = 1.0
                            # here we apply (Gram–Schmidt)
                            axis_2 = axis_2 - self.InnerProd(axis_1, axis_2) / self.InnerProd(axis_1, axis_1) * axis_1
                            axis_2 = axis_2 / self.Norm2(axis_2)
                            elem.SetValue(KM.LOCAL_AXIS_2, axis_2)
                        else: # we avoid that starting/ending element
                            self.model_part.GetSubModelPart(tendon_name).RemoveElement(elem, 0)

            intersections_file.close()

    def CalculateIntersectionVectors(self, Line):
        coords_1 = KM.Vector(3)
        coords_1[0] = float(Line[1])
        coords_1[1] = float(Line[2])
        coords_1[2] = float(Line[3])

        coords_2 = KM.Vector(3)
        coords_2[0] = float(Line[4])
        coords_2[1] = float(Line[5])
        coords_2[2] = float(Line[6])

        intersection_vector = KM.Vector(3)
        intersection_vector[0] = coords_2[0] - coords_1[0]
        intersection_vector[1] = coords_2[1] - coords_1[1]
        intersection_vector[2] = coords_2[2] - coords_1[2]

        return intersection_vector
    
    def Norm2(self, Vector):
        return math.sqrt(Vector[0]**2 + Vector[1]**2+ Vector[2]**2)
    
    def InnerProd(self, Vector1, Vector2):
        return Vector1[0]*Vector2[0] + Vector1[1]*Vector2[1] + Vector1[2]*Vector2[2]
    
    def CheckFiberParticipation(self, kf):
        if kf <= 0.0:
            kf = self.fiber_participation_threshold
        elif kf > 1.0:
            kf = 0.99999999
        return kf
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IGA

import math
import numpy

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetDirectorsProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class SetDirectorsProcess(KratosMultiphysics.Process):
    """
    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"      : "please_specify_model_part_name"
        }
        """)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.director_utilities = IGA.DirectorUtilities(self.model_part, settings)

    def ExecuteInitialize(self):
        director = KratosMultiphysics.Vector(3) 
        director[0] = 0
        director[1] = 0
        director[2] = 1

        dirTang = KratosMultiphysics.Matrix(3,2) 
        dirTang[0,0] = 1
        dirTang[1,0] = 0
        dirTang[2,0] = 0

        dirTang[0,1] = 0
        dirTang[1,1] = 1
        dirTang[2,1] = 0
        for node in self.model_part.Nodes:
            node.SetValue(IGA.DIRECTOR, director)
            node.SetValue(IGA.DIRECTORTANGENTSPACE, dirTang)

        self.director_utilities.ComputeDirectors()
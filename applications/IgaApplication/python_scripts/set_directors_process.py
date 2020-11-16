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


    def TangentSpaceFromStereographicProjection(director):
        BLA = KratosMultiphysics.Matrix(3,2) 
        y = KratosMultiphysics.Vector(2) 
        st =  numpy.sign(director[2])
        s  = 1/(1+st*(director[2]))
        y[0] = director[0] * s;
        y[1] = director[1] * s;
        ys1 = y[0]*y[0];
        ys2 = y[1]*y[1];
        s2 = 2*(1+ys1+ys2);

        BLA[0,0] = s2-4*ys1
        BLA[1,0] = -4*y[0]*y[1]
        BLA[2,0] = -st*4*y[0]

        BLA[0,1] = -4*y[0]*y[1]
        BLA[1,1] = s2-4*ys2
        BLA[2,1] = -st*4*y[1]

        normcol0 = 1/( math.sqrt(BLA[0,0] *BLA[0,0] + BLA[1,0] *  BLA[1,0] + BLA[2,0] * BLA[2,0]))
        normcol1 = 1/( math.sqrt(BLA[0,1] *BLA[0,1] + BLA[1,1] *  BLA[1,1] + BLA[2,1] * BLA[2,1]))

        for i in range(3):
            BLA[i,0] =  BLA[i,0] /normcol0
            BLA[i,1] =  BLA[i,1] /normcol1

        return BLA

    def FinalizeNonLinearIteration(self):
        for node in self.model_part.Nodes:
            director = node.GetValue(IGA.DIRECTOR)
            inc2d3 = node.GetSolutionStepValue(IGA.DIRECTORINC)
            inc2d = KratosMultiphysics.Vector(2) 
            inc2d[0] =inc2d3[0]
            inc2d[1] =inc2d3[1]
            print(inc2d)
            print(director)
            BLA = node.GetValue(IGA.DIRECTORTANGENTSPACE)
            print(BLA)

            inc3d = BLA * inc2d
            print(inc3d)
            director = director + inc3d
            director = director / math.sqrt(director[0] * director[0] + director[1] * director[1] + director[2] * director[2])
            print("director",director)
      
            node.SetValue(IGA.DIRECTOR, director)
            node.SetValue(IGA.DIRECTORTANGENTSPACE, SetDirectorsProcess.TangentSpaceFromStereographicProjection(director))
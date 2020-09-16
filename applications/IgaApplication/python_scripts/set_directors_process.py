from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics

import math

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

     dirTang[0,1] = 1
     dirTang[1,1] = 0
     dirTang[2,1] = 0
     for node in self.model_part.Nodes:
         node.SetValue(DIRECTOR, director)
         node.SetValue(DIRECTORTANGENTSPACE, dirTang)

#       for node in self.model_part.Nodes:
#          
#          node.SetValue(DIRECTOR, director);
#
#
#        for element in self.model_part.Elements:
#            A3 = element.GetGeometry().Normal(point_number); //this makes only sense if the geometry is undeformed
#            A3 = A3 / sqrt((A3 * A3));
#            //Vector m_Nvec = trans(row(m_N, point_number));
#            rRightHandSideMatrix = rRightHandSideMatrix + outer_prod(trans(row(m_N, point_number)),  trans(A3)     ) ;

#            rLeftHandSideMatrix = rLeftHandSideMatrix + outer_prod(row(m_N,point_number), row(m_N, point_number));


    def ExecuteInitializeSolutionStep(self):
        for node in self.model_part.Nodes:
            director = node.GetValue(DIRECTOR)
            inc2d = node.GetValue(DIRECTORINC)

            BLA = node.GetValue(DIRECTORTANGENTSPACE)

            inc3d = BLA * inc2d

            director = director + inc3d
            director = director / sqrt(director * director)

            node.SetValue(DIRECTOR, director)
            # node.SetValue(DIRECTORTANGENTSPACE, director)

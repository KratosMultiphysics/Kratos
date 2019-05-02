from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPFEMSolid
import KratosMultiphysics.ContactMechanicsApplication as KratosContact


import os.path


def Factory( custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CavityProcess(Model, custom_settings["Parameters"])


class CavityProcess(KratosMultiphysics.Process):

    #
    def __init__(self, model_part, custom_settings):

        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
        }
        """)

        self.model_part = model_part['Main_Domain']

        # initialize figure path
        problem_path = os.getcwd()
        self.csv_path = os.path.join(problem_path, "cavity.csv")
        csv_file = open(self.csv_path, "a")
        line = " 0.0 , 0.0 , 0.0 , 0.0  \n"
        csv_file.write(line)
        csv_file.close()


    def ExecuteFinalizeSolutionStep(self):
        
        time = self._GetStepTime()
        
        SomeProcessInfo = KratosMultiphysics.ProcessInfo()
        import math
        Force = 0.0
        for node in self.model_part.GetNodes(0): 
            x = node.X0;
            y = node.Y0
            r = x*x + y*y;
            r = math.sqrt(r)
            xDir = x;
            yDir = y;
            if ( abs(r) > 1e-12):
                xDir = x/r;
                yDir = y/r;
            if ( r < 1.00001):
                res = node.GetSolutionStepValue(KratosSolid.DISPLACEMENT_REACTION)
                Force = Force + res[0]*xDir + res[1]*yDir;

        line = str(time) + " , " + str(Force) + " \n"

        csv_file = open(self.csv_path, "a")
        csv_file.write(line)
        csv_file.close()

    def _AddToLine(self, line, vector):

        for i in vector:
            line = line + str(i) + " , "

        return line

    def  _GetStepTime(self):

        return self.model_part.ProcessInfo[KratosMultiphysics.TIME]


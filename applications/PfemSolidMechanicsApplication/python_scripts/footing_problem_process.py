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
    return FootingProblemProcess(Model, custom_settings["Parameters"])


class FootingProblemProcess(KratosMultiphysics.Process):

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
        self.csv_path = os.path.join(problem_path, "footing.csv")
        csv_file = open(self.csv_path, "a")
        line = " 0.0 , 0.0 , 0.0 , 0.0  \n"
        csv_file.write(line)
        csv_file.close()


    def ExecuteFinalizeSolutionStep(self):
        
        time = self._GetStepTime()
        
        SomeProcessInfo = KratosMultiphysics.ProcessInfo()

        Force = 0.0*KratosMultiphysics.Array3()
        for node in self.model_part.GetNodes(0): 
            if ( abs(node.Y0) < 1e-15):
                if (node.X0 < 1.0001):
                    res = node.GetSolutionStepValue(KratosSolid.DISPLACEMENT_REACTION)
                    Force = Force + res

        line = str(time) + " , "

        line = self._AddToLine( line, Force)
        line = line + " \n"

        csv_file = open(self.csv_path, "a")
        csv_file.write(line)
        csv_file.close()

    def _AddToLine(self, line, vector):

        for i in vector:
            line = line + str(i) + " , "

        return line

    def  _GetStepTime(self):

        return self.model_part.ProcessInfo[KratosMultiphysics.TIME]


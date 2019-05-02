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
    return GaussPointUtility(Model, custom_settings["Parameters"])


class GaussPointUtility(KratosMultiphysics.Process):

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
        self.csv_path = os.path.join(problem_path, "GaussPoint.csv")


    def ExecuteFinalizeSolutionStep(self):
        
        time = self._GetStepTime()
        
        SomeProcessInfo = KratosMultiphysics.ProcessInfo()

        elems = self.model_part.GetElements(0)
        for elem in elems: 
            cauchy_vector = elem.GetValuesOnIntegrationPoints( KratosMultiphysics.CAUCHY_STRESS_TENSOR, SomeProcessInfo)
            epsi_vector = elem.GetValuesOnIntegrationPoints( KratosMultiphysics.GREEN_LAGRANGE_STRAIN_TENSOR, SomeProcessInfo)
            break

        line = str(time) + " , "

        line = self._AddToLine( line, cauchy_vector[0])
        line = self._AddToLine( line, epsi_vector[0])
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


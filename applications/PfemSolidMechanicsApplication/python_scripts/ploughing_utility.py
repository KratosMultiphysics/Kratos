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
    return PloughingUtility(Model, custom_settings["Parameters"])


class PloughingUtility(KratosMultiphysics.Process):

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
        self.csv_path = os.path.join(problem_path, "Forces_Ploughing.csv")


    def ExecuteBeforeOutputStep(self):

        delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]


        for node in self.model_part.GetNodes(0):
            delta_disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            delta_disp = delta_disp -  node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 1)
            for ii in range(0,2):
                delta_disp[ii] = delta_disp[ii] / delta_time
            if (node.SolutionStepsDataHas(KratosMultiphysics.VELOCITY)):
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, delta_disp)
            
    def ExecuteAfterOutputStep(self):

        for node in self.model_part.GetNodes(0):
            if ( node.SolutionStepsDataHas( KratosMultiphysics.CONTACT_FORCE)):
                CF = node.GetSolutionStepValue(KratosMultiphysics.CONTACT_FORCE)
                CF = 0.0*CF
                node.SetSolutionStepValue( KratosMultiphysics.CONTACT_FORCE, CF)
            if ( node.SolutionStepsDataHas( KratosMultiphysics.CONTACT_NORMAL)):
                CF = node.GetSolutionStepValue(KratosMultiphysics.CONTACT_NORMAL)
                CF = 0.0*CF
                node.SetSolutionStepValue( KratosMultiphysics.CONTACT_NORMAL, CF)

            if ( node.SolutionStepsDataHas( KratosContact.CONTACT_STRESS)):
                CF = node.GetSolutionStepValue( KratosContact.CONTACT_STRESS)
                CF = 0.0*CF
                node.SetSolutionStepValue( KratosContact.CONTACT_STRESS, CF)


            if ( node.SolutionStepsDataHas( KratosContact.EFFECTIVE_CONTACT_FORCE)):
                CF = node.GetSolutionStepValue( KratosContact.EFFECTIVE_CONTACT_FORCE)
                CF = 0.0*CF
                node.SetSolutionStepValue( KratosContact.EFFECTIVE_CONTACT_FORCE, CF)

            if ( node.SolutionStepsDataHas( KratosContact.EFFECTIVE_CONTACT_STRESS)):
                CF = node.GetSolutionStepValue( KratosContact.EFFECTIVE_CONTACT_STRESS)
                CF = 0.0*CF
                node.SetSolutionStepValue( KratosContact.EFFECTIVE_CONTACT_STRESS, CF)


    def ExecuteFinalizeSolutionStep(self):

        time = self._GetStepTime()

        ContactForce = KratosMultiphysics.Array3();
        for ii in range(0,3):
            ContactForce[ii] = 0

        for node in self.model_part.GetNodes(0):
            Force = node.GetSolutionStepValue( KratosMultiphysics.CONTACT_FORCE)
            ContactForce = ContactForce + Force


        line = str(time) + " , " + str(ContactForce[0]) + " , " + str(ContactForce[1]) + " \n "

        csv_file = open(self.csv_path, "a")
        csv_file.write(line)
        csv_file.close()



    def  _GetStepTime(self):

        return self.model_part.ProcessInfo[KratosMultiphysics.TIME]


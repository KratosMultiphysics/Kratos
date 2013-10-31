#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *
CheckForPreviousImport()

class ConstitutiveLawUtility:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part  = model_part
        self.domain_size = domain_size 
        

    #######################################################################
    def Initialize(self):
        self.SetConstitutiveLaw();

    #######################################################################   
    def SetConstitutiveLaw(self):
        if(self.domain_size == 2):
            for prop in self.model_part.Properties:
                ConstitutiveLawName=prop.GetValue(CONSTITUTIVE_LAW_NAME)
                if(ConstitutiveLawName == "LinearElasticPlaneStrain"):
                    prop.SetValue(CONSTITUTIVE_LAW, LinearElasticPlaneStrain2DLaw() )
                    print "Linear Elastic Plane Strain 2D model selected"
                elif(ConstitutiveLawName == "LinearElasticPlaneStress"):
                    prop.SetValue(CONSTITUTIVE_LAW, LinearElasticPlaneStress2DLaw() )
                    print "Linear Elastic Plane Stress 2D model selected"
                elif(ConstitutiveLawName == "HyperElasticPlaneStrain"):
                    prop.SetValue(CONSTITUTIVE_LAW, HyperElasticPlaneStrain2DLaw() )
                    print "Hyperelastic Plane Strain 2D model selected"
                elif(ConstitutiveLawName == "HyperElasticPlaneStrainUP"):
                    prop.SetValue(CONSTITUTIVE_LAW, HyperElasticUPPlaneStrain2DLaw() )
                    print "Hyperelastic UP Plane Strain 2D model selected"
                elif(ConstitutiveLawName == "LinearElasticAxisym"):
                    prop.SetValue(CONSTITUTIVE_LAW, LinearElasticAxisym2DLaw() )
                    print "Linear Elastic Axisym 2D model selected"
                elif(ConstitutiveLawName == "HyperElasticAxisym"):
                    prop.SetValue(CONSTITUTIVE_LAW, HyperElasticAxisym2DLaw() )
                    print "Hyper Elastic Axisym 2D model selected"
                elif(ConstitutiveLawName == "HyperElasticAxisymUP"):
                    prop.SetValue(CONSTITUTIVE_LAW, HyperElasticUPAxisym2DLaw() )
                    print "Hyper Elastic UP Axisym 2D model selected"
                elif(ConstitutiveLawName == "HyperElasticPlasticJ2PlaneStrain"):
                    prop.SetValue(CONSTITUTIVE_LAW, HyperElasticPlasticJ2PlaneStrain2DLaw() )
                    print "Hyper Elastic Plastic J2 Plane Strain 2D model selected"
                elif(ConstitutiveLawName == "HenckyElasticMatsuokaPlasticPlaneStrain"):
                    prop.SetValue(CONSTITUTIVE_LAW, HenckyMatsuokaPlasticPlaneStrain2DLaw() )
                    print "Hencky Elastic Matsuoka Plastic Plane Strain 2D model selected"
                elif(ConstitutiveLawName == "HenckyElasticMatsuokaPlasticAxisym"):
                    prop.SetValue(CONSTITUTIVE_LAW, HenckyMatsuokaPlasticAxisym2DLaw() )
                    print "Hencky Elastic Matsuoka Plastic Axisym 2D model selected"
                elif(ConstitutiveLawName == "NonLinearHenckyCamClayPlasticPlaneStrain2DLaw"):
                    prop.SetValue(CONSTITUTIVE_LAW, NonLinearHenckyCamClayPlasticPlaneStrain2DLaw() )
                    print "Hencky Elastic Matsuoka Plastic Axisym 2D model selected"
                else:
                    print "ERROR: CONSTITUTIVE_LAW 2D not defined properly"
                    
            else:
                for prop in self.model_part.Properties:
                    ConstitutiveLawName=prop.GetValue(CONSTITUTIVE_LAW_NAME);
                    if(ConstitutiveLawName == "LinearElastic"):
                        prop.SetValue(CONSTITUTIVE_LAW, LinearElastic3DLaw() )
                        print "Linear Elastic 3D model selected"
                    elif(ConstitutiveLawName == "HyperElastic"):
                        prop.SetValue(CONSTITUTIVE_LAW, HyperElastic3DLaw() )
                        print "Hyperelastic 3D model selected"
                    elif(ConstitutiveLawName == "HyperElasticUP"):
                        prop.SetValue(CONSTITUTIVE_LAW, HyperElasticUP3DLaw() )
                        print "Hyperelastic UP 3D model selected"
                    elif(ConstitutiveLawName == "HyperElasticPlasticJ2"):
                        prop.SetValue(CONSTITUTIVE_LAW, HyperElasticPlasticJ23DLaw() )
                        print "Hyper Elastic Plastic J2 3D model selected"
                    else:
                        print "ERROR: CONSTITUTIVE_LAW 3D not defined properly"




    #######################################################################   



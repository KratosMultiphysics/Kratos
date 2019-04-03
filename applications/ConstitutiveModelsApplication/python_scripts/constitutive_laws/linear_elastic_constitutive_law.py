from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterial

class LinearElasticLaw:

    def __init__(self):
        #print(" python: elastic law ")
        self.constitutive_law = KratosMaterial.Linear3DLaw()
    #
    def CallLaw(self):
        print(" python: call the python law ")

    #
    def GetLawFeatures(self, LawFeatures):
        #print(" python: law features ")
        Options = LawFeatures.GetOptions()
        Options.Set(KratosMultiphysics.ConstitutiveLaw.THREE_DIMENSIONAL_LAW, True)
        Options.Set(KratosMultiphysics.ConstitutiveLaw.INFINITESIMAL_STRAINS, True)
        Options.Set(KratosMultiphysics.ConstitutiveLaw.ISOTROPIC, True)

        #LawFeatures.SetOptions(Options)

        LawFeatures.SetStrainMeasure(KratosMultiphysics.StrainMeasure.StrainMeasure_Infinitesimal)
        LawFeatures.SetStrainMeasure(KratosMultiphysics.StrainMeasure.StrainMeasure_Deformation_Gradient)

        LawFeatures.SetStrainSize(self.GetStrainSize())
        LawFeatures.SetSpaceDimension(self.WorkingSpaceDimension())

    #
    def WorkingSpaceDimension(self):
        return 3
    #
    def GetStrainSize(self):
        return 6

    #
    def CalculateMaterialResponsePK2(self, LawParameters):
        #print(" python: PK2 ")
        self.constitutive_law.CalculateMaterialResponsePK2(LawParameters)

    #
    def CalculateMaterialResponseKirchhoff(self, LawParameters):
        #print(" python: Kirchhoff ")

        LawOptions = LawParameters.GetOptions()
        Properties = LawParameters.GetMaterialProperties()

        # be careful only assignation of the scalar components of these tensors is allowed
        # any local copy will be lost out of the python scope
        self.StrainVector = LawParameters.GetStrainVector()
        self.StressVector = LawParameters.GetStressVector()
        self.ConstitutiveMatrix = LawParameters.GetConstitutiveMatrix()

        young_modulus = Properties.GetValue(KratosMultiphysics.YOUNG_MODULUS)
        poisson_coefficient = Properties.GetValue(KratosMultiphysics.POISSON_RATIO)

        if(LawOptions.Is(KratosMultiphysics.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN) == True):
            deformation_gradient = LawParameters.GetDeformationGradientF()
            self.CalculateStrainMatrix(deformation_gradient)

        if(LawOptions.Is(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRESS) == True):
            #if(LawOptions.Is(KratosMultiphysics.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR) == True):
            self.CalculateLinearElasticMatrix(young_modulus,poisson_coefficient)
            self.CalculateStress()

        '''
        if(LawOptions.Is(KratosMultiphysics.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR) == False):
            print(" strain",self.StrainVector)
            print(" constitutive",self.ConstitutiveMatrix)
            print(" stress",self.StressVector)
        '''

        '''
        self.constitutive_law.CalculateMaterialResponseKirchhoff(LawParameters)
        StrainVector = LawParameters.GetStrainVector()
        StressVector = LawParameters.GetStressVector()
        ConstitutiveMatrix = LawParameters.GetConstitutiveMatrix()
        print(" STRAIN",StrainVector)
        print(" CONSTITUTIVE", ConstitutiveMatrix)
        print(" STRESS",StressVector)
        '''

    #
    def CalculateStrainMatrix(self, DeformationGradient):

        DeformationGradientT = KratosMultiphysics.Matrix(3,3)
        for i in range(0,3):
            for j in range(0,3):
                DeformationGradientT[i,j] = DeformationGradient[j,i]

        RightCauchyGreen = DeformationGradientT * DeformationGradient

        self.StrainVector = Vector(6)
        for i in range(0,3):
            self.StrainVector[i] = 0.5 * (RightCauchyGreen[i,i] - 1.0)

        self.StrainVector[3] = RightCauchyGreen[0,1]
        self.StrainVector[4] = RightCauchyGreen[1,2]
        self.StrainVector[5] = RightCauchyGreen[0,2]

        #return self.StrainVector

    #
    def CalculateLinearElasticMatrix(self, YoungModulus, PoissonCoefficient):

        for i in range(0,6):
            for j in range(0,6):
                self.ConstitutiveMatrix[i,j] = 0

        self.ConstitutiveMatrix[0,0] = (YoungModulus*(1.0-PoissonCoefficient)/((1.0+PoissonCoefficient)*(1.0-2.0*PoissonCoefficient)));
        self.ConstitutiveMatrix[1,1] = self.ConstitutiveMatrix[0,0]
        self.ConstitutiveMatrix[2,2] = self.ConstitutiveMatrix[0,0]

        self.ConstitutiveMatrix[3,3] = self.ConstitutiveMatrix[0,0] * (1.0-2.0*PoissonCoefficient)/(2.0*(1.0-PoissonCoefficient));
        self.ConstitutiveMatrix[4,4] = self.ConstitutiveMatrix[3,3]
        self.ConstitutiveMatrix[5,5] = self.ConstitutiveMatrix[3,3]

        self.ConstitutiveMatrix[0,1] = self.ConstitutiveMatrix[0,0] * PoissonCoefficient/(1.0-PoissonCoefficient);
        self.ConstitutiveMatrix[1,0] = self.ConstitutiveMatrix[0,1]

        self.ConstitutiveMatrix[0,2] = self.ConstitutiveMatrix[0,1]
        self.ConstitutiveMatrix[2,0] = self.ConstitutiveMatrix[0,1]

        self.ConstitutiveMatrix[1,2] = self.ConstitutiveMatrix[0,1]
        self.ConstitutiveMatrix[2,1] = self.ConstitutiveMatrix[0,1]

        #return self.ConstitutiveMatrix

    #
    def CalculateStress(self):
        #self.StressVector.resize(6)
        #for i in range(0,6):
        #    self.StressVector[i] = 0

        stress_vector = self.ConstitutiveMatrix * self.StrainVector

        for i in range(0,6):
            self.StressVector[i] = stress_vector[i]


        #for i in range(0,6):
        #    for j in range(0,6):
        #        self.StressVector[i] += self.ConstitutiveMatrix[i,j] * self.StrainVector[j]



        #return self.StressVector
    #
    def Check(self, Properties, Geometry, Info):

        check = 1
        if( Properties.Has(KratosMultiphysics.YOUNG_MODULUS) ):
            E = Properties.GetValue(KratosMultiphysics.YOUNG_MODULUS)
            if( E <= 0.0 ):
                check = 0
                print(" YOUNG_MODULUS has an invalid value:", E)
        else:
            check = 0
            print(" YOUNG_MODULUS has Key zero ")

        if( Properties.Has(KratosMultiphysics.YOUNG_MODULUS) ):
            nu = Properties.GetValue(KratosMultiphysics.POISSON_RATIO)
            if( (nu > 0.499 and nu < 0.501) or (nu < - 0.999 and nu > -1.01) ):
                check = 0
                print(" POISSON_RATIO has an invalid value:", E)
        else:
            check = 0
            print(" POISSON_RATIO has Key zero ")

        if( Properties.Has(KratosMultiphysics.DENSITY) ):
            rho = Properties.GetValue(KratosMultiphysics.DENSITY)
            if( rho <= 0.0 ):
                check = 0
                print(" DENSITY has an invalid value:", rho)
        else:
            check = 0
            print(" DENSITY has Key zero ")

        return check

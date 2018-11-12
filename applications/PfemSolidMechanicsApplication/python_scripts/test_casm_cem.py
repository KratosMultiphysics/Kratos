from __future__ import print_function, absolute_import, division

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication     as KratosSolid
#import KratosMultiphysics.ExternalSolversApplication    as KratosSolvers
import KratosMultiphysics.PfemApplication           as KratosPfem
#import KratosMultiphysics.ContactMechanicsApplication   as KratosContact
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPfemSolid
#import KratosMultiphysics.PfemFluidDynamicsApplication  as KratosPfemFluid

import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestCASMCementedModel(KratosUnittest.TestCase):

    def __init__(self, *args, **kwargs):
        try:
            import importlib
            import math as math
        except ImportError as e:
            self.skipTest("Missing python libraries ( importlib and math)")

        super(KratosUnittest.TestCase, self).__init__(*args, **kwargs)

        
    #oedometer
    def test_IsotropicCompression(self):
        import math
        import numpy as np
        import matplotlib.pyplot as plt

        ID = 'IsoC4Mh12G0s1000_'
        NumberIncrements = 1000
        F11max = 0.8
        IncrementalF = self._set_identity_matrix()

        IncrementalF[0,0] = F11max**(1/NumberIncrements)
        IncrementalF[1,1] = F11max**(1/NumberIncrements)
        IncrementalF[2,2] = F11max**(1/NumberIncrements)
        
        #initialize plot
        fig = plt.figure('iso comp', figsize=(10,5))
        xx = np.linspace(0, 200, 10)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.plot(xx, xx*0.9, color='red')
        ax1.set_xlabel('p')
        ax1.set_ylabel('q')
        ax2.set_xlabel('eps_v')
        ax2.set_ylabel('log(p)')
        
        #reference test without bounding
        self._create_material_model_and_law()
        self.material_law.SetPlasticVariables(-80, 0)
        pp, qq, epsVol, FF01, FF11 = self._compute_strain_driven_problem(IncrementalF, NumberIncrements, ID+str(0))
        ax1.plot(pp, qq)
        ax2.plot(np.log10(pp), epsVol)
        

        if(1==1):
            #variation of bonding parameters
            for pt in np.arange(1.0,5.0,1.0):
                #create and initialize element, gauss point & constitutive law
                self._create_material_model_and_law()
                self.properties.SetValue( KratosMultiphysics.INITIAL_BONDING, pt)
                #self.properties.SetValue( KratosMultiphysics.DEGRADATION_RATE_COMPRESSION, pt)
                #self.properties.SetValue( KratosMultiphysics.DEGRADATION_RATE_SHEAR, pt)
                self.parameters.SetMaterialProperties( self.properties )
                self.material_law.SetPlasticVariables(-80, pt)
                #stress integration
                pp, qq, epsVol, FF01, FF11 = self._compute_strain_driven_problem(IncrementalF, NumberIncrements, ID+str(pt))
                ax1.plot(pp, qq)
                ax2.plot(np.log10(pp), epsVol)

        #
        #plt.show()


    #oedometer
    def _test_OedometricLoading(self):
        import math
        import numpy as np
        import matplotlib.pyplot as plt

        ID = 'OedA4Mh12G0longs1000_'
        NumberIncrements = 1000
        F11max = 0.7
        IncrementalF = self._set_identity_matrix()

        #IncrementalF[0,0] = F11max**(1/NumberIncrements)
        IncrementalF[1,1] = F11max**(1/NumberIncrements)
        #IncrementalF[2,2] = F11max**(1/NumberIncrements)
        
        #initialize plot
        fig = plt.figure('oedometer', figsize=(10,5))
        xx = np.linspace(0, 200, 10)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.plot(xx, xx*0.9, color='red')
        ax1.set_xlabel('p')
        ax1.set_ylabel('q')
        ax2.set_xlabel('eps_v')
        ax2.set_ylabel('log(p)')
        
        #reference test without bounding
        self._create_material_model_and_law()
        self.material_law.SetPlasticVariables(-80, 0)
        pp, qq, epsVol, FF01, FF11 = self._compute_strain_driven_problem(IncrementalF, NumberIncrements, ID+str(0))
        ax1.plot(pp, qq)
        ax2.plot(np.log10(pp), epsVol)
        

        if(1==1):
            #variation of bonding parameters
            for pt in np.arange(1.0,5.0,1.0):
                #create and initialize element, gauss point & constitutive law
                self._create_material_model_and_law()
                self.properties.SetValue( KratosMultiphysics.INITIAL_BONDING, pt)
                #self.properties.SetValue( KratosMultiphysics.DEGRADATION_RATE_COMPRESSION, pt)
                self.parameters.SetMaterialProperties( self.properties )
                self.material_law.SetPlasticVariables(-80, pt)
                #stress integration
                pp, qq, epsVol, FF01, FF11 = self._compute_strain_driven_problem(IncrementalF, NumberIncrements, ID+str(pt))
                ax1.plot(pp, qq)
                ax2.plot(np.log10(pp), epsVol)

        #
        plt.show()


    #simple shear
    def _test_SimpleShear(self):
        import math
        import numpy as np
        import matplotlib.pyplot as plt
        
        ID = 'SSA4Yh12G0s1000_'
        NumberIncrements = 1000
        F01max = 1.0
        IncrementalF = self._set_identity_matrix()
        IncrementalF[0,1] = F01max/NumberIncrements
        #IncrementalF[1,0] = F01max/NumberIncrements
        
        #initialize plot
        fig = plt.figure('simple shear', figsize=(10,5))
        xx = np.linspace(0, 200, 10)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        ax1.plot(xx, xx*0.9, color='red')
        ax1.set_xlabel('p')
        ax1.set_ylabel('q')
        ax2.set_xlabel('F12')
        ax2.set_ylabel('q')
        
        #reference test without bounding
        self._create_material_model_and_law()
        self.material_law.SetPlasticVariables(-80, 0)
        pp, qq, epsVol, FF01, FF11 = self._compute_strain_driven_problem(IncrementalF, NumberIncrements, ID+str(0))
        #plt.subplot(1,2,1)
        ax1.plot(pp, qq)
        #plt.subplot(1,2,2)
        ax2.plot(FF01, qq)
        #plt.show()
        
        #write to output file
        #self._OpenOutputFile(self)
        #self._WriteThisToFile(self, , stress, strain)
        
        if(1==1):
            #variation of bounding parameters
            for pt in np.arange(1.0,5.0,1.0):
                #create and initialize element, gauss point & constitutive law
                self._create_material_model_and_law()
                self.properties.SetValue( KratosMultiphysics.INITIAL_BONDING, pt)
                #self.properties.SetValue( KratosMultiphysics.DEGRADATION_RATE_COMPRESSION, pt)
                self.parameters.SetMaterialProperties( self.properties )
                self.material_law.SetPlasticVariables(-80, pt)

                #stress integration
                pp, qq, epsVol, FF01, FF11 = self._compute_strain_driven_problem(IncrementalF, NumberIncrements, ID+str(pt))
                ax1.plot(pp, qq)
                ax2.plot(FF01, qq)
                #plt.show()
                #write to output file
                

        #
        plt.show()


    #
    def _compute_strain_driven_problem(self, IncrF, nIncr, ID):
    
        import numpy as np
        
        print(" ")
        print("*** COMPUTING STRAIN DRIVEN PROBLEM *** F_incr: ", IncrF)

        #set initial values
        self.parameters.SetDeformationGradientF( self.F )
        self.parameters.SetDeterminantF( self.detF )
        #self.material_law.SetPlasticVariables(-60, 3.0)
        
            #print( "  F: ", self.parameters.GetDeformationGradientF(), " det(F): ", self.parameters.GetDeterminantF() )
            #print( "  Initial state parameters set ")
            #print( " " )
        
            #print(" ")
            #print("+++++++++++++++++++++++++++++")
            #print("+++ Compute Initial state +++")
            #print("+++++++++++++++++++++++++++++")   
                 
        #initial step
        self.material_law.CalculateMaterialResponseKirchhoff( self.parameters )
        self.material_law.FinalizeMaterialResponseKirchhoff( self.parameters )
        self.material_law.FinalizeSolutionStep( self.properties, self.geometry, self.N, self.model_part.ProcessInfo )
        p0 = self.material_law.GetPreconPressure()
        bb = self.material_law.GetBonding()
        
        #initialize output
        pp = 1.1*np.arange(nIncr+1)
        qq = 1.1*np.arange(nIncr+1)
        epsVol = 1.1*np.arange(nIncr+1)
        #epsDev = 1.1*np.arange(nIncr+1)
        FF01 = 1.1*np.arange(nIncr+1)
        FF11 = 1.1*np.arange(nIncr+1)
        
        
        stress = self.parameters.GetStressVector()
        self.strain = self.parameters.GetStrainVector()
        strain = self._ComputeHenckyStrainFromF(self.F)
        
        self.stress = self.parameters.GetStressVector()
        pp[0], qq[0] = self._calculate_invariants()
        epsVol[0] = strain[0] + strain[1] + strain[2]
        epsVol[0] /= -3.0
        FF01[0] = self.F[0,1]
        FF11[0] = self.F[1,1]
        
        #write to output file
        self._OpenOutputFile(ID)
        self._WriteThisToFile(0, pp[0], qq[0], p0, bb, epsVol[0], FF01[0], FF11[0], stress, strain)
        
        #loop over nIncr
        for step in range(1, nIncr+1):
        
            #print(" ")
            #print("+++++++++++++++++++++++++++++++++++++++++++++++++")
            #print("+++ applying deformation increment ", step, " +++")
            #print("+++++++++++++++++++++++++++++++++++++++++++++++++")
        
            #add increment
            self.F = IncrF * self.F
            self.detF = self._compute_determinant( self.F )
            
            #compute increment
            self.parameters.SetDeformationGradientF( self.F )
            self.parameters.SetDeterminantF( self.detF )
            self.material_law.CalculateMaterialResponseKirchhoff( self.parameters )
            self.material_law.FinalizeMaterialResponseKirchhoff( self.parameters )
            self.material_law.FinalizeSolutionStep( self.properties, self.geometry, self.N, self.model_part.ProcessInfo )
            p0 = self.material_law.GetPreconPressure()
            bb = self.material_law.GetBonding()
            
            #calculate invariants from self.stress
            self.stress = self.parameters.GetStressVector()
            strain = self._ComputeHenckyStrainFromF( self.F )
            pp[step], qq[step] = self._calculate_invariants()
            epsVol[step] = strain[0] + strain[1] + strain[2]
            epsVol[step] /= -3.0
            FF01[step] = self.F[0,1]
            FF11[step] = self.F[1,1]
            
            #
            stress = self.parameters.GetStressVector()
            
            #write to output file
            #self._OpenOutputFile()
            self._WriteThisToFile(step, pp[step], qq[step], p0, bb, epsVol[step], FF01[step], FF11[step], stress, strain)
          
        return pp, qq, epsVol, FF01, FF11    
            
        
    #
    def _test_triaxial_compression(self):
    
        #import stuff
        import numpy as np
        import matplotlib.pyplot as plt
        
        #create constitutive law
        self._create_material_model_and_law()
        #set initial F, detF and hardening parameters
        self.parameters.SetDeformationGradientF( self.F )
        self.parameters.SetDeterminantF( self.detF )
        self.material_law.SetPlasticVariables(-80, 0.0)
        #calculate initial step
        self.material_law.CalculateMaterialResponseCauchy( self.parameters )
        self.material_law.FinalizeMaterialResponseCauchy( self.parameters )
        self.material_law.FinalizeSolutionStep( self.properties, self.geometry, self.N, self.model_part.ProcessInfo )
        
        #
        stress = self.parameters.GetStressVector()
        self.strain = self.parameters.GetStrainVector()
        strain = self._ComputeHenckyStrainFromF(self.F)
        
        #create output files
        
        
        #BC for the test
        sigma_h = -80;
        nLoadingSteps = 200
        FinalAxialDeformation = 0.35
        
        #strain vector (Hencky????)
        XX = 0.0*self.parameters.GetStrainVector()
        print('XX = ',XX)

        
        #initialize output
        VolStrain = 1.1*np.arange(nLoadingSteps+1)
        DevStrain = 1.1*np.arange(nLoadingSteps+1)
        VolStress = 1.1*np.arange(nLoadingSteps+1)
        DevStress = 1.1*np.arange(nLoadingSteps+1)
        #set initial volumetric & deviatoric strain
        VolStrain[0] = strain[0] + strain[1] + strain[2]
        VolStrain[0] *= -1.0
        DevStrain[0] = strain[1] - strain[0]
        DevStrain[0] *= -2.0/3.0
        #set initial volumetric & deviatoric stress
        VolStress[0] = stress[0] + stress[1] + stress[2]
        VolStress[0] /= -3.0
        DevStress[0] = stress[1] - stress[0]
        DevStress[0] *= -1.0
        
        
        #stepwise increase of axial strain
        for tt in range(1, nLoadingSteps+1):
        
            #total axial strain
            verticalStrain = -tt/nLoadingSteps * FinalAxialDeformation
            
            for iteracio in range(0,100):
                #get current stress and strain
                stress = self.parameters.GetStressVector()
                strain = self.parameters.GetStrainVector()
                
                #compute the residual and its norm
                residual = 0.0*stress + stress
                
                residual[0] = stress[0] - sigma_h
                verticalF = self.F[1,1]
                residual[1] = verticalF - verticalStrain
                residual[2] = stress[2] - sigma_h
                
                #check stopping criteria
                norm = 0
                for i in range (0,4):
                    norm += residual[i]*residual[i]
                if (norm < 1e-12):
                    break
                
                #construct system matrix
                C = self.parameters.GetConstitutiveMatrix()
                
                print(stress)
                print(C)
                print(strain)
                
                SystemMatrix = np.zeros((4,4))
                for i in range(0,4):
                    for j in range(0,4):
                        SystemMatrix[i,j] = C[i,j]
                #introduce axial strain condition
                for i in range(0,4):
                    SystemMatrix[1,i] = 0.0
                SystemMatrix[1,1] = 1.0
                #solve dx = C_inv * residual
                dx = np.linalg.solve(SystemMatrix, residual)
                #update Hencky strain vector
                for i in range(0,4):
                    XX[i] = XX[i] - dx[i]
                #set deformation gradient
                self.F = self._set_identity_matrix()
                for j in range(0,3):
                    self.F[j,j] = self.F[j,j] + 1.0*XX[j]
                self.detF = self._compute_determinant(self.F)
                
                print('XX = ',XX)
                print('dX = ',dx)
                
                #compute constitutive response based on updated F
                print(self.F)
                print(self.detF)
                self.parameters.SetDeformationGradientF( self.F )
                self.parameters.SetDeterminantF( self.detF )
                self.material_law.CalculateMaterialResponseCauchy( self.parameters )
                
                
            #finalize solution step for appropriate constitutive matrix
            self.material_law.FinalizeMaterialResponseCauchy( self.parameters )
            self.material_law.FinalizeSolutionStep( self.properties, self.geometry, self.N, self.model_part.ProcessInfo )
            stress = self.parameters.GetStressVector()
            self.strain = self.parameters.GetStrainVector()
            strain = self._ComputeHenckyStrainFromF( self.F )
            
            #write to file
            
            #compute output
            VolStrain[tt] = strain[0] + strain[1] + strain[2]
            VolStrain[tt] *= -1.0
            DevStrain[tt] = strain[1] - strain[0]
            DevStrain[tt] *= -2.0/3.0
            VolStress[tt] = stress[0] + stress[1] + stress[2]
            VolStress[tt] /= -3.0
            DevStress[tt] = stress[1] - stress[0]
            DevStress[tt] *= -1.0
            
            
        #plots q-eps_q, q-p
        plt.subplot(2,2,2)
        plt.plot( DevStrain, DevStress )
        plt.subplot(2,2,1)
        plt.plot( VolStress, DevStress )
        
        plt.subplot(2,2,3)
        plt.plot( VolStress, VolStrain )
        plt.subplot(2,2,4)
        plt.plot( DevStrain, VolStrain )
        plt.show(block=False)
        
        

    #
    def _create_material_model_and_law(self):

        settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "MaterialDomain",
            "properties_id"   : 1,
            "material_name"   : "soil",
            "constitutive_law": {
                "law_name"   : "KratosMultiphysics.PfemSolidMechanicsApplication.BorjaHenckyCasmCemPlasticAxisym2DLaw",
                "model_name" : "KratosMultiphysics.PfemSolidMechanicsApplication.GensNovaModel"
            },
            "variables": {
                "KratosMultiphysics.DENSITY": 1.7,
                "KratosMultiphysics.YOUNG_MODULUS": 0.0,
                "KratosMultiphysics.NORMAL_COMPRESSION_SLOPE": 0.085,
                "KratosMultiphysics.SWELLING_SLOPE": 0.0078,
                "KratosMultiphysics.INITIAL_SHEAR_MODULUS": 500.0,
                "KratosMultiphysics.ALPHA_SHEAR": 0.0,
                "KratosMultiphysics.PRE_CONSOLIDATION_STRESS": 0.8e+02,
                "KratosMultiphysics.OVER_CONSOLIDATION_RATIO": 4.0,
                "KratosMultiphysics.CRITICAL_STATE_LINE": 0.9,
                "KratosMultiphysics.INTERNAL_FRICTION_ANGLE": 0.0,
                "KratosMultiphysics.SPACING_RATIO": 1.5,
                "KratosMultiphysics.SHAPE_PARAMETER": 4.0,
                "KratosMultiphysics.INITIAL_BONDING": 0.0,
                "KratosMultiphysics.DEGRADATION_THRESHOLD": 0.0,
                "KratosMultiphysics.DEGRADATION_RATE_COMPRESSION": 10.0,
                "KratosMultiphysics.DEGRADATION_RATE_SHEAR": 10.0,
                "KratosMultiphysics.ALPHA_TENSILE": 0.0,
                "KratosMultiphysics.PLASTIC_DEVIATORIC_STRAIN_HARDENING": 0.0,
                "KratosMultiphysics.DENSITY_WATER": 1.0,
                "KratosMultiphysics.K0": 0.7,
                "KratosMultiphysics.THICKNESS": 1.0
            },
            "element_type": "Triangle2D3",
            "nodes" : [ [0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0] ],
            "strain": {
                "deformation_gradient" : [ [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0] ],
                "jacobian": 1.0
            },
            "echo_level" : 0

        }
        """)

        #create model_part
        self.model_part = KratosMultiphysics.ModelPart(settings["model_part_name"].GetString())
        self.echo_level = settings["echo_level"].GetInt()

        #add nodes to model_part
        self.number_of_nodes = settings["nodes"].size()
        self.nodes = []
        for ii in range (0, self.number_of_nodes):
            node = self.model_part.CreateNewNode(ii, settings["nodes"][ii][0].GetDouble(), settings["nodes"][ii][1].GetDouble(), settings["nodes"][ii][2].GetDouble())
            self.nodes.append(node)

        #create element
        self.geometry = KratosMultiphysics.Geometry()
        self.dimension = 3
        if(settings["element_type"].GetString() == "Tetrahedra3D4"):
            if( self.number_of_nodes != 4 ):
                print(" number of nodes:",self.number_of_nodes," do not matches geometry :", settings["element_type"].GetString() )
            else:
                self.geometry = KratosMultiphysics.Tetrahedra3D4(self.nodes[0],self.nodes[1],self.nodes[2],self.nodes[3])

        if(settings["element_type"].GetString() == "Triangle2D3"):
            if( self.number_of_nodes != 3 ):
                print(" number of nodes:",self.number_of_nodes," do not matches geometry :", settings["element_type"].GetString() )
            else:
                self.geometry  = KratosMultiphysics.Triangle2D3(self.nodes[0],self.nodes[1],self.nodes[2])
                self.dimension = 2

        #create and set material properties
        self.properties = self.model_part.Properties[settings["properties_id"].GetInt()]
        self.variables = settings["variables"]
        for key, value in self.variables.items():
            variable = self._GetItemFromModule(key)
            if(value.IsDouble()):
                self.properties.SetValue(variable, value.GetDouble())
            elif(value.IsArray()):
                vector_value = KratosMultiphysics.Vector(value.size())
                for i in range(0, value.size() ):
                    vector_value[i] = value[i].GetDouble()
                self.properties.SetValue(variable, vector_value)

        #set strain
        self.F = KratosMultiphysics.Matrix(3,3)
        self.strain_measure = settings["strain"]["deformation_gradient"]

        for i in range(0,3):
            for j in range(0,3):
                self.F[i,j] = self.strain_measure[i][j].GetDouble()

        self.detF = settings["strain"]["jacobian"].GetDouble()

        #initialize element parameters
        self.N = KratosMultiphysics.Vector(self.number_of_nodes)
        self.DN_DX = KratosMultiphysics.Matrix(self.number_of_nodes, self.dimension)
        
        #create constitutive law
        self.material_law = self._GetItemFromModule( settings["constitutive_law"]["law_name"].GetString() )()
        self.material_law.InitializeMaterial(self.properties, self.geometry, self.N)
        
        #check dimension
        if(self.material_law.WorkingSpaceDimension() != self.dimension):
            raise Exception( "mismatch between the ConstitutiveLaw dimension and the dimension of the space")

        #set calculation flags
        self.options = KratosMultiphysics.Flags()
        self.options.Set(KratosMultiphysics.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, True)
        self.options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRESS, True)
        self.options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)

        #set calculation variables
        self.stress_vector          = KratosMultiphysics.Vector( self.material_law.GetStrainSize() )
        self.strain_vector          = KratosMultiphysics.Vector( self.material_law.GetStrainSize() )
        self.constitutive_matrix    = KratosMultiphysics.Matrix( self.material_law.GetStrainSize(), self.material_law.GetStrainSize() )
        
        self.parameters = KratosMultiphysics.ConstitutiveLawParameters()
        
        self.parameters.SetOptions( self.options )
        self.parameters.SetDeformationGradientF( self.F )
        self.parameters.SetDeterminantF( self.detF )
        self.parameters.SetStrainVector( self.strain_vector )
        self.parameters.SetStressVector( self.stress_vector )
        self.parameters.SetConstitutiveMatrix( self.constitutive_matrix )
        self.parameters.SetShapeFunctionsValues( self.N )
        self.parameters.SetShapeFunctionsDerivatives( self.DN_DX )
        self.parameters.SetProcessInfo( self.model_part.ProcessInfo )
        self.parameters.SetMaterialProperties( self.properties )
        self.parameters.SetElementGeometry( self.geometry )


    #
    def _GetItemFromModule(self,my_string):

        import importlib 

        splitted = my_string.split(".")
        if(len(splitted) == 0):
            raise Exception("something wrong. Trying to split the string "+my_string)
        if(len(splitted) == 1):
            return eval(my_string)
        else:
            module_name = ""
            for i in range(len(splitted)-1):
                module_name += splitted[i] 
                if i != len(splitted)-2:
                    module_name += "."

            module = importlib.import_module(module_name)
            return getattr(module,splitted[-1])
    
    #        
    def _OpenOutputFile(self, ID):
        import os.path

        #return;
        problem_path = os.getcwd()
        self.csv_path = os.path.join(problem_path, "GaussPoint" + str(ID) + ".csv")

        csv_file = open(self.csv_path, "w")
        csv_file.close()
        
    #
    def _WriteThisToFile(self, t, pp, qq, p0, bb, epsVol, FF01, FF11, stress, strain):

        line = str(t) + " " + str(pp) + " " + str(qq) + " " + str(p0) + " " + str(bb) + " " + str(epsVol) + " " + str(FF01) + " " + str(FF11) + " "
        for st in stress:
            line = line + str(st) + " "

        for st in strain:
            line = line + str(st) + " "

        line = line[:-2]
        line = line + "\n"
        csv_file = open(self.csv_path, "a")
        csv_file.write(line)
        csv_file.close()

    #
    def _ComputeHenckyStrainFromF(self, F):

        import numpy as np
        #print('############# F before: ', F)
        HenckyF = F + F
        #print('############# F after: ', F)
        #print('############# Hencky after: ', HenckyF)
        #print('############# self.strain before: ', self.strain)
        strain = self.strain + self.strain * 0.0
        #print('############# self.strain after: ', self.strain)
        #print('############# strain: ', strain)

        for i in range(0,3):
            HenckyF[i,i] = F[i,i]*F[i,i]
            HenckyF[i,i] = np.log(HenckyF[i,i])/2.0
            strain[i] = HenckyF[i,i]

        return strain

    #
    def _calculate_invariants(self):

        import math

        if(self.material_law.GetStrainSize() == 6):
            #Compute invariants
            Pressure = -(self.stress[0] + self.stress[1] + self.stress[2])/3.0 #Geotech sign convention
            dev = self.stress + 0.0*self.stress
            for i in range(0,3):
                dev[i] = dev[i] + Pressure

            J2 = 0
            for i in range(0,3):
                J2 = J2 + dev[i] *dev[i]
            for i in range(3,6):
                J2 = J2 + 2.0*dev[i] *dev[i]
            J2 = math.sqrt(J2*0.5)
            J2 = math.sqrt(3.0) * J2
            DeviatoricQ = J2
        else:
            #Compute invariants
            Pressure = -(self.stress[0] + self.stress[1] + self.stress[2])/3.0 #Geotech sign convention
            dev = self.stress + 0.0*self.stress
            for i in range(0,3):
                dev[i] = dev[i] + Pressure

            J2 = 0
            for i in range(0,3):
                J2 = J2 + dev[i] *dev[i]
            for i in range(3,4):
                J2 = J2 + 2.0*dev[i] *dev[i]
            J2 = math.sqrt(J2*0.5)
            J2 = math.sqrt(3.0) * J2
            DeviatoricQ = J2
        return Pressure, DeviatoricQ


    def _set_identity_matrix(self):
        identity = KratosMultiphysics.Matrix(3,3)
        for i in range(0,3):
            for j in range(0,3):
                if ( i == j ):
                    identity[i,j] = 1.0
                else:
                    identity[i,j] = 0.0

        return identity
        
    def _compute_determinant(self, A):

        det = 0

        det = det + A[0,0]*A[1,1]*A[2,2]
        det = det + A[1,0]*A[2,1]*A[0,2]
        det = det + A[2,0]*A[0,1]*A[1,2]

        det = det - A[0,2]*A[1,1]*A[2,0]
        det = det - A[1,2]*A[2,1]*A[0,0]
        det = det - A[2,2]*A[0,1]*A[1,0]

        return det

        


if __name__ == '__main__':
    KratosUnittest.main()


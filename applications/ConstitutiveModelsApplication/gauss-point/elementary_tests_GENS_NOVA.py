from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterialModels
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestModifiedCamClayModel(KratosUnittest.TestCase):

    def __init__(self, *args, **kwargs):
        try:
            import importlib
            import math as math
        except ImportError as e:
            self.skipTest("Missing python libraries ( importlib and math)")

        self.size_parametric_analysis = 0;

        super(KratosUnittest.TestCase, self).__init__(*args, **kwargs)



    def test_OedometricLoading(self):
        import math
        import numpy as np
        
        NumberIncrements = 500
        IncrementalF = self._set_identity_matrix()
        IncrementalF[0,0] = 0.999

        self._create_material_model_and_law()
        for case in range(0, self.size_parametric_analysis):
            self._create_material_model_and_law()
            self.properties.SetValue(self.parametric_analysis_variable,  self.parametric_analysis_values[case])
            self.parameters.SetMaterialProperties( self.properties )
    
            self._compute_strain_driven_problem(IncrementalF, NumberIncrements)
            Pressure, DeviatoricQ = self._calculate_invariants()

        import matplotlib.pyplot as plt
        plt.show()


    # Triaxial test at gauss point
    def test_TriaxialLoading(self):

        import math
        import numpy as np

        self._create_material_model_and_law()
        for case in range(0, self.size_parametric_analysis):
            self._create_material_model_and_law()
            
            self.properties.SetValue(self.parametric_analysis_variable,  self.parametric_analysis_values[case])
            self.parameters.SetMaterialProperties( self.properties )

            self._test_TriaxialCompression()
    
        import matplotlib.pyplot as plt
        plt.show(block=True)

    # driver of the triaxial compression
    def _test_TriaxialCompression(self):
    
        #self._create_material_model_and_law()
        
        self.F = self._set_identity_matrix()
        self.detF = 1.0
        self.parameters.SetDeformationGradientF( self.F )
        self.parameters.SetDeterminantF( self.detF )
        self.material_law.CalculateMaterialResponseCauchy( self.parameters )
        self.material_law.FinalizeMaterialResponseCauchy( self.parameters )
        self.material_law.FinalizeSolutionStep( self.properties, self.geometry, self.N, self.model_part.ProcessInfo )


        import numpy as np
        import matplotlib.pyplot as plt

        #isotropic compression
        #
        #
        XX = 0.0*self.parameters.GetStrainVector()
        
        sigma_v = -100.0
        k0 = 1.0
        sigma_h = sigma_v * k0

        for iteracio in range(0, 200):
            stress = self.parameters.GetStressVector()

            residual = 0.0*stress + stress;
            
            residual[0] = stress[0] - sigma_h
            residual[2] = stress[2] - sigma_h
            residual[1] = stress[1] - sigma_v
            
            #Compute the norm of the residual
            norm = 0
            for i in range(0,6):
                norm += residual[i]*residual[i]

            if (norm < 1e-9):
                break
            if ( iteracio > 195):
                disp('The initial state is not converging')
                pNode


            C = self.parameters.GetConstitutiveMatrix()
            SystemMatrix = np.zeros((6,6));
            for i in range(0,6):
                for j in range(0,6):
                    SystemMatrix[i,j] = C[i,j];
            dx = np.linalg.solve(SystemMatrix, residual)
            for i in range(0,6):
                XX[i] = XX[i] - dx[i];
            self.F = self._set_identity_matrix()
            for j in range(0,3):
                self.F[j,j] = self.F[j,j] + 1.0*XX[j]
            self.detF = self._compute_determinant(self.F)

            #Compute
            self.parameters.SetDeformationGradientF( self.F )
            self.parameters.SetDeterminantF( self.detF )
            self.material_law.CalculateMaterialResponseCauchy( self.parameters)


        self.material_law.FinalizeMaterialResponseCauchy( self.parameters )
        self.material_law.FinalizeSolutionStep( self.properties, self.geometry, self.N, self.model_part.ProcessInfo )
        stress = self.parameters.GetStressVector();
        self.strain = self.parameters.GetStrainVector();
        strain = self.ComputeStrainFromF(self.F)

        
        self._OpenOutputFile()
        self._WriteThisToFile(0, stress, strain)

        isotropicLoadingStrain = self.F[1,1]

        # Second part
        nLoadingSteps = 400
        FinalAxialDeformation = 0.35



        VolStrain = 1.1*np.arange( float(nLoadingSteps+1) )
        DevStrain = 1.1*np.arange( float(nLoadingSteps+1) )
        VolStress = 1.1*np.arange( float(nLoadingSteps+1) )
        DevStress = 1.1*np.arange( float(nLoadingSteps+1) )

        VolStrain[0] = strain[0] + strain[1] + strain[2]
        VolStrain[0] *= -1.0
        DevStrain[0] = strain[1] - strain[0]
        DevStrain[0] *= -1.0

        VolStress[0] = stress[0] + stress[1] + stress[2]
        VolStress[0] /= -3.0
        DevStress[0] = stress[1]-stress[0]
        DevStress[0] *= -1.0




        for t in range(1, nLoadingSteps+1):

            verticalStrain = isotropicLoadingStrain - float(t)/float(nLoadingSteps) * FinalAxialDeformation
            for iteracio in range(0, 100):

                stress = self.parameters.GetStressVector()
                strain = self.parameters.GetStrainVector();

                residual = 0.0*stress + stress;
            
                residual[0] = stress[0] - sigma_h
                residual[2] = stress[2] - sigma_h
                verticalF = self.F[1,1]
                residual[1] = verticalF - verticalStrain

                norm = 0
                for i in range(0,6):
                    norm += residual[i]*residual[i]
                if (norm < 1e-12):
                    break
                if (iteracio > 90):
                    print('In the triaxial loading. Not converging')
                    pNode

                C = self.parameters.GetConstitutiveMatrix()
                SystemMatrix = np.zeros((6,6));
                for i in range(0,6):
                    for j in range(0,6):
                        SystemMatrix[i,j] = C[i,j];
                # partial inversion of the system
                for i in range(0,6): 
                    SystemMatrix[1,i] = 0.0
                SystemMatrix[1,1] = 1.0;
                dx = np.linalg.solve(SystemMatrix, residual)
                for i in range(0,6):
                    XX[i] = XX[i] - dx[i];
                self.F = self._set_identity_matrix()
                for j in range(0,3):
                    self.F[j,j] = self.F[j,j] + 1.0*XX[j]
                self.detF = self._compute_determinant(self.F)

                #Compute
                
                self.parameters.SetDeformationGradientF( self.F )
                self.parameters.SetDeterminantF( self.detF )
                self.material_law.CalculateMaterialResponseCauchy( self.parameters)


            self.material_law.FinalizeMaterialResponseCauchy( self.parameters )
            self.material_law.FinalizeSolutionStep( self.properties, self.geometry, self.N, self.model_part.ProcessInfo )
            stress = self.parameters.GetStressVector();
            self.strain = self.parameters.GetStrainVector();
            strain = self.ComputeStrainFromF(self.F)
        
            self._WriteThisToFile(t, stress, strain)

            VolStrain[t] = strain[0] + strain[1] + strain[2]
            VolStrain[t] *= -1.0
            DevStrain[t] = strain[1] - strain[0]
            DevStrain[t] *= -1.0

            VolStress[t] = stress[0] + stress[1] + stress[2]
            VolStress[t] /= -3.0
            DevStress[t] = stress[1]-stress[0]
            DevStress[t] *= -1.0

        plt.subplot(2,2,1)
        plt.plot(VolStress, DevStress)
        plt.xlabel('pressure')
        plt.ylabel('DevStress')

        plt.subplot(2,2,2)
        plt.plot(DevStrain, DevStress )
        plt.xlabel('DevStrain')
        plt.ylabel('DevStress')

        
        plt.subplot(2,2,3)
        plt.plot(VolStress, VolStrain)
        plt.xlabel('pressure')
        plt.ylabel('VolStrain')

        plt.subplot(2,2,4)
        plt.plot(DevStrain, VolStrain)
        plt.xlabel('DevStrain')
        plt.ylabel('VolStrain')
        plt.show(block=False)


    def _compute_strain_driven_problem(self, IncrF, nIncr):
    
        self.parameters.SetDeformationGradientF(self.F)
        self.parameters.SetDeterminantF(self.detF)
        self.material_law.CalculateMaterialResponseKirchhoff(self.parameters)
        self.material_law.FinalizeMaterialResponseKirchhoff(self.parameters)
        self.material_law.FinalizeSolutionStep(self.properties, self.geometry, self.N, self.model_part.ProcessInfo)

        import numpy as np 
        pp = 1.1*np.arange( nIncr)
        qq = 1.1*np.arange( nIncr)

        self._OpenOutputFile()
        stress = self.parameters.GetStressVector()
        self.strain = self.parameters.GetStrainVector()
        strain = self.ComputeStrainFromF(self.F)
        self._WriteThisToFile(0, stress, strain)

        for step in range(0, nIncr):

            #Actualize strain
            self.F = IncrF * self.F
            self.detF = self._compute_determinant(self.F)

            #Compute
            self.parameters.SetDeformationGradientF( self.F )
            self.parameters.SetDeterminantF( self.detF )
            self.material_law.CalculateMaterialResponseKirchhoff( self.parameters )
            self.material_law.FinalizeMaterialResponseKirchhoff( self.parameters )
            self.material_law.FinalizeSolutionStep( self.properties, self.geometry, self.N, self.model_part.ProcessInfo )

            self.stress = self.parameters.GetStressVector()
            pp[step], qq[step] = self._calculate_invariants()
    
            stress = self.parameters.GetStressVector()
            strain = self.ComputeStrainFromF(self.F)
            self._WriteThisToFile(step, stress, strain)

        if (nIncr > 2):
           
           import matplotlib.pyplot as plt
           plt.figure(1)
           plt.plot(pp, qq)
           plt.show(block=False)

    def _compute_determinant(self, A):

        det = 0

        det = det + A[0,0]*A[1,1]*A[2,2]
        det = det + A[1,0]*A[2,1]*A[0,2]
        det = det + A[2,0]*A[0,1]*A[1,2]

        det = det - A[0,2]*A[1,1]*A[2,0]
        det = det - A[1,2]*A[2,1]*A[0,0]
        det = det - A[2,2]*A[0,1]*A[1,0]

        return det


    def _set_identity_matrix(self):
        identity = KratosMultiphysics.Matrix(3,3)
        for i in range(0,3):
            for j in range(0,3):
                if ( i == j ):
                    identity[i,j] = 1.0
                else:
                    identity[i,j] = 0.0

        return identity

    def _multiply_matrix_by_number(self, A, alpha):
        for i in range(0,3):
            for j in range(0,3):
                A[i,j] = alpha*A[i,j]
        return A

    def _create_material_model_and_law(self):

        parameter_file = open("ProjectParameters.json", 'r')
        settings = KratosMultiphysics.Parameters(parameter_file.read() )
        self.model_part = KratosMultiphysics.ModelPart(settings["model_part_name"].GetString())
        self.echo_level = settings["echo_level"].GetInt()

        #read nodes
        self.number_of_nodes = settings["nodes"].size()
        self.nodes = [] #self.model_part.GetNodes()
        for i in range(0, self.number_of_nodes):
            node = self.model_part.CreateNewNode(i,settings["nodes"][i][0].GetDouble(),settings["nodes"][i][1].GetDouble(),settings["nodes"][i][2].GetDouble())
            self.nodes.append(node)

        self.geometry  = KratosMultiphysics.Geometry()
        self.dimension = 3
        if( settings["element_type"].GetString() == "Tetrahedra3D4"):
            if( self.number_of_nodes != 4 ):
                print(" number of nodes:",self.number_of_nodes," do not matches geometry :", settings["element_type"].GetString() )
            else:
                self.geometry = KratosMultiphysics.Tetrahedra3D4(self.nodes[0],self.nodes[1],self.nodes[2],self.nodes[3])

                
        if( settings["element_type"].GetString() == "Triangle2D3"):
            if( self.number_of_nodes != 4 ):
                print(" number of nodes:",self.number_of_nodes," do not matches geometry :", settings["element_type"].GetString() )
            else:
                self.geometry  = KratosMultiphysics.Triangle2D3(self.nodes[0],self.nodes[1],self.nodes[2])
                self.dimension = 2
                
        #material properties
        self.properties = self.model_part.Properties[settings["properties_id"].GetInt()]

        #read variables
        self.variables = settings["variables"]
        for key, value in self.variables.items():
            variable = self._GetItemFromModule(key)
            if( value.IsDouble() ):
                self.properties.SetValue(variable, value.GetDouble())
            elif( value.IsArray() ):
                vector_value = KratosMultiphysics.Vector(value.size())
                for i in range(0, value.size() ):
                    vector_value[i] = value[i].GetDouble()
                self.properties.SetValue(variable, vector_value)

        
        #create constitutive law
        self.material_model = self._GetItemFromModule( settings["constitutive_law"]["model_name"].GetString())()
        self.material_law   = self._GetItemFromModule( settings["constitutive_law"]["law_name"].GetString())(self.material_model)

        #check dimension
        if(self.material_law.WorkingSpaceDimension() != self.dimension):
            raise Exception( "mismatch between the ConstitutiveLaw dimension and the dimension of the space")

        #set strain
        self.F = KratosMultiphysics.Matrix(3,3)
        self.strain_measure = settings["strain"]["deformation_gradient"]
       
        for i in range(0,3):
            for j in range(0,3):
                self.F[i,j] = self.strain_measure[i][j].GetDouble()
        
        self.detF = settings["strain"]["jacobian"].GetDouble()        
        
        #element parameters
        self.N     = KratosMultiphysics.Vector(self.number_of_nodes)
        self.DN_DX = KratosMultiphysics.Matrix(self.number_of_nodes, self.dimension)
        
        #set calculation flags
        self.options = KratosMultiphysics.Flags()
        self.options.Set(KratosMultiphysics.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, True)
        self.options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRESS, True)
        self.options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)

        #set calculation variables
        self.stress_vector       = KratosMultiphysics.Vector(self.material_law.GetStrainSize())
        self.strain_vector       = KratosMultiphysics.Vector(self.material_law.GetStrainSize())
        self.constitutive_matrix = KratosMultiphysics.Matrix(self.material_law.GetStrainSize(),self.material_law.GetStrainSize())

                
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

        if ( self.size_parametric_analysis == 0):
            self.parametric_analysis_variable = KratosMultiphysics.KratosGlobals.GetVariable( settings["parametric_analysis_variable"].GetString() )  
            self.size_parametric_analysis = settings["parametric_analysis_values"].size()
            self.parametric_analysis_values = [] #self.model_part.GetNodes()
            for i in range(0, self.size_parametric_analysis):
                value = settings["parametric_analysis_values"][i].GetDouble()
                self.parametric_analysis_values.append(value)


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

    def _calculate_invariants(self):

        import math

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
        return Pressure, DeviatoricQ

    def _OpenOutputFile(self):
        import os.path

        problem_path = os.getcwd()
        self.csv_path = os.path.join(problem_path, "GaussPoint.csv")

        csv_file = open(self.csv_path, "w")
        csv_file.close()

    def ComputeStrainFromF(self, F):

        import numpy as np
        HenckyF = F + F;
        strain = self.strain + self.strain * 0.0;

        for i in range(0,3):
            HenckyF[i,i] = F[i,i]*F[i,i]
            HenckyF[i,i] = np.log(HenckyF[i,i])/2.0
            strain[i] = HenckyF[i,i]



        return strain



    def _WriteThisToFile(self, t, stress, strain):

        import numpy as np
        strainV = strain[0]+strain[1]+strain[2];
        J = np.exp(strainV);
        stress = J * stress;

        line = str(t) + " , "
        for st in stress:
            line = line + str(st) + " , "

        for st in strain:
            line = line + str(st) + " , "

        ps = 0;
        ps = self.material_law.GetValue( KratosMultiphysics.ConstitutiveModelsApplication.PS, ps)
        pt = 0;
        pt = self.material_law.GetValue( KratosMultiphysics.ConstitutiveModelsApplication.PT, pt)
        pm = 0;
        pm = self.material_law.GetValue( KratosMultiphysics.ConstitutiveModelsApplication.PM, pm)

        line = line + str(ps) + " , " + str(pt) + " , " + str(pm)

        line = line + " \n"
        csv_file = open(self.csv_path, "a")
        csv_file.write(line)
        csv_file.close()
 

if __name__ == '__main__':
    KratosUnittest.main()

from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterialModels
import KratosMultiphysics.UmatApplication as KratosUmatModels
import KratosMultiphysics.KratosUnittest as KratosUnittest
import time as timer

import numpy 



class TestParametricAnalysisOnConstitutiveModel(KratosUnittest.TestCase):

    def __init__(self, *args, **kwargs):
        try:
            import importlib
            import math as math
        except ImportError as e:
            self.skipTest("Missing python libraries ( importlib and math)")

        self.size_parametric_analysis = 0;

        super(KratosUnittest.TestCase, self).__init__(*args, **kwargs)

        self.WriteToAFile = True
        self.IMPLEX = False
        self._create_material_model_and_law()



    # (Drained) oedometric compression
    def test_OedometricLoading(self):
        
        nLoadingSteps = 10
        FinalDeformation = 0.1

        for case in range(0, self.size_parametric_analysis ):
            self.current_case = case 

            self._create_material_model_and_law()
            nLoadingSteps = self._SetParametricAnalysisVariable(case, nLoadingSteps)

            self._setInitialStressState()
            self._OpenOutputFile('oedometer.csv')

            self.strain_vector = self.parameters.GetStrainVector()
            delta_strain = 0.0*self.strain_vector
            delta_strain[1] = -FinalDeformation / float(nLoadingSteps)
            self._compute_strain_driven_problem(delta_strain, nLoadingSteps)

    # Isotropic compression
    def _test_IsotropicLoading(self):
        
        nLoadingSteps = 1000
        FinalDeformation = 0.1

        SavedInitialStress = self.initial_stress_state
        for i in range(1,3):
            self.initial_stress_state[i] = self.initial_stress_state[0]

        for case in range(0, self.size_parametric_analysis ):
            self.current_case = case 

            self._create_material_model_and_law()
            nLoadingSteps = self._SetParametricAnalysisVariable(case, nLoadingSteps)

            for i in range(1,3):
                self.initial_stress_state[i] = self.initial_stress_state[0]

            self._setInitialStressState()
            self._OpenOutputFile('isotropic.csv')

            self.strain_vector = self.parameters.GetStrainVector()
            delta_strain = 0.0*self.strain_vector
            for ii in range(0,3):
                delta_strain[ii] = -FinalDeformation / float(nLoadingSteps)
            self._compute_strain_driven_problem(delta_strain, nLoadingSteps)

        self.initial_stress_state = SavedInitialStress

    # undrained triaxial compression
    def test_TriaxialLoading_Undrained(self):
        
        nLoadingSteps = 1000
        FinalDeformation = 0.2

        for case in range(0, self.size_parametric_analysis ):
            self.current_case = case 

            self._create_material_model_and_law()
            nLoadingSteps = self._SetParametricAnalysisVariable(case, nLoadingSteps)

            IncrementalF = self._set_identity_matrix()
            IncrementalF[0,1] = (FinalDeformation/float(nLoadingSteps) )

            self._setInitialStressState()
            self._OpenOutputFile('undrained_triaxial.csv')

            self.strain_vector = self.parameters.GetStrainVector()
            delta_strain = 0.0*self.strain_vector
            delta_strain[4] = -FinalDeformation / float(nLoadingSteps)
            self._compute_strain_driven_problem(delta_strain, nLoadingSteps)


    # Triaxial test at gauss point
    def test_TriaxialLoading_Drained(self):


        FinalDeformation = 0.14
        nLoadingSteps = 1000

        for case in range(0, self.size_parametric_analysis ):
            self.current_case = case

            self._create_material_model_and_law()
            nLoadingSteps = self._SetParametricAnalysisVariable(case, nLoadingSteps)

            self._setInitialStressState()
            self._OpenOutputFile('drained_triaxial.csv')

            self._test_TriaxialCompression(FinalDeformation, nLoadingSteps)


    # Triaxial test at gauss point
    def test_TriaxialLoadingExtension_Drained(self):


        FinalDeformation = -0.14
        nLoadingSteps = 1000

        for case in range(0, self.size_parametric_analysis ):
            self.current_case = case

            self._create_material_model_and_law()
            nLoadingSteps = self._SetParametricAnalysisVariable(case, nLoadingSteps)

            self._setInitialStressState()
            self._OpenOutputFile('drained_triaxial_ext.csv')

            self._test_TriaxialCompression(FinalDeformation, nLoadingSteps)


    # driver of the triaxial compression
    def _test_TriaxialCompression(self, FinalAxialDeformation, nLoadingSteps):
   
        Tolerance = 1e-12

        self.model_part.ProcessInfo.SetValue( KratosMultiphysics.IMPLEX, False)

        stress = self.parameters.GetStressVector()
        self.strain_vector = self.parameters.GetStrainVector()
        strain = self.strain_vector
        self._WriteThisToFile(0, stress, strain)

        isotropicLoadingStrain = self.strain_vector[1]

        self.model_part.ProcessInfo.SetValue( KratosMultiphysics.IMPLEX, self.IMPLEX)

        Penalty = 10000.0;
        sigma_h = self.initial_stress_state[0]
        residual = 0.0*stress + stress;
        SystemMatrix = numpy.zeros((6,6));
        XX = numpy.zeros((6,1))
        C = self.parameters.GetConstitutiveMatrix()
        self.parameters.GetStrainVector( self.strain_vector )
        self.material_law.CalculateMaterialResponseCauchy( self.parameters)


        for t in range(1, nLoadingSteps+1):

            verticalStrainImposed = isotropicLoadingStrain - float(t)/float(nLoadingSteps) * FinalAxialDeformation

            for iteracio in range(0, 100):

                residual[0] = stress[0] - sigma_h*(1.0-2.0*float(t)/float(nLoadingSteps))
                residual[2] = stress[2] - sigma_h*(1.0-2.0*float(t)/float(nLoadingSteps))
                verticalSS = self.strain_vector[1]
                residual[1] = Penalty*(verticalSS - verticalStrainImposed)

                norm = 0.0
                for i in range(0,6):
                    norm += residual[i]*residual[i]

                print('Iteration '+ str(iteracio)+ ' residual '+ str(norm) )
                if (norm < Tolerance):
                    break
                #if (iteracio > 100):
                #    print('In the triaxial loading. Not converging')
                #    pNode

                for i in range(0,6):
                    for j in range(0,6):
                        SystemMatrix[i,j] = C[i,j];

                # partial inversion of the system
                for i in range(0,6): 
                    SystemMatrix[1,i] = 0.0
                SystemMatrix[1,1] = Penalty;

                dx = numpy.linalg.solve(SystemMatrix, residual)


                for iii in range(0,6):
                    self.strain_vector[iii] = self.strain_vector[iii] - dx[iii]

                #Compute
                self.parameters.GetStrainVector( self.strain_vector )
                self.material_law.CalculateMaterialResponseCauchy( self.parameters)


            self.model_part.ProcessInfo.SetValue( KratosMultiphysics.IMPLEX, False)
            self.material_law.FinalizeMaterialResponseCauchy( self.parameters )
            self.model_part.ProcessInfo.SetValue( KratosMultiphysics.IMPLEX, self.IMPLEX)


            if ( self.WriteToAFile == True):
                stress = self.parameters.GetStressVector();
                self.strain_vector = self.parameters.GetStrainVector();
                strain = strain
                self._WriteThisToFile(t, stress, strain)
       
        if ( self.WriteToAFile == False):
            stress = self.parameters.GetStressVector();
            self.strain_vector = self.parameters.GetStrainVector();
            strain = self.strain_vector
            self._WriteThisToFile(t, stress, strain)




    def _setInitialStressState(self):

        self.strain_vector = self.parameters.GetStrainVector()
        self.strain_vector = 0.0*self.strain_vector
        self.parameters.SetStrainVector(self.strain_vector)
        self.material_law.CalculateMaterialResponseCauchy( self.parameters )
        self.material_law.FinalizeMaterialResponseCauchy( self.parameters )


        XX = 0.0*self.parameters.GetStrainVector()
        
        sigma_v = self.initial_stress_state[1]
        sigma_h = self.initial_stress_state[0]
        sigma_h2 = self.initial_stress_state[2]


        someProcessInfo = KratosMultiphysics.ProcessInfo()
        self.material_law.SetValue( KratosMultiphysics.ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS, self.initial_stress_state, someProcessInfo)

        self.material_law.CalculateMaterialResponseCauchy( self.parameters)


        self.material_law.FinalizeMaterialResponseKirchhoff( self.parameters )

        

    def _compute_strain_driven_problem(self, delta_strain, nIncr):
    
        self.parameters.SetDeformationGradientF(self.F)
        self.parameters.SetDeterminantF(self.detF)
        self.material_law.CalculateMaterialResponseKirchhoff(self.parameters)
        self.material_law.FinalizeMaterialResponseKirchhoff(self.parameters)


        stress = self.parameters.GetStressVector()
        self.strain_vector = self.parameters.GetStrainVector()
        strain = self.strain_vector
        self._WriteThisToFile(0, stress, strain)

        for step in range(1, nIncr+1):

            #Actualize strain
            self.strain_vector = self.strain_vector + delta_strain

            #Compute
            self.parameters.SetStrainVector( self.strain_vector )

            self.material_law.CalculateMaterialResponseCauchy( self.parameters )
            self.material_law.FinalizeMaterialResponseCauchy( self.parameters )

            if ( self.WriteToAFile == False):
                continue;

            self.stress = self.parameters.GetStressVector()
            stress = self.parameters.GetStressVector()
            strain = self.strain_vector
            self._WriteThisToFile(step, stress, strain)
    

        if ( self.WriteToAFile == False):
            stress = self.parameters.GetStressVector()
            strain = self.ComputeStrainFromF(self.F)
            self._WriteThisToFile(step, stress, strain)


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

    def _create_material_model_and_law(self):

        parameter_file = open("ProjectParameters.json", 'r')
        settings = KratosMultiphysics.Parameters(parameter_file.read() )
        parameter_file.close()
        model = KratosMultiphysics.Model()
        self.model_part = model.CreateModelPart(settings["model_part_name"].GetString())
        #self.model_part = KratosMultiphysics.ModelPart(settings["model_part_name"].GetString())
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
            self.parametric_analysis_values = [] 
            for i in range(0, self.size_parametric_analysis):
                value = settings["parametric_analysis_values"][i].GetDouble()
                self.parametric_analysis_values.append(value)

        self.initial_stress_state = [];
        for i in range(0, 3):
            value = settings["initial_stress_state"][i].GetDouble()
            self.initial_stress_state.append(value)
        self.IMPLEX = settings["use_implex"].GetBool()
        self.NameOfSeries = settings["name_computational_series"].GetString()


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


    def _SetParametricAnalysisVariable(self, case, nLoadingSteps = 1000):

        if ( self.parametric_analysis_variable == KratosMultiphysics.KratosGlobals.GetVariable('STEP') ):
            nLoadingSteps = int(self.parametric_analysis_values[case])
        else:
            self.properties.SetValue(self.parametric_analysis_variable,  self.parametric_analysis_values[case])
            self.parameters.SetMaterialProperties( self.properties )

        return nLoadingSteps

    def ComputeStrainFromF(self, F):

        #HenckyF = F + F;
        strain = self.strain_vector_1234212 + self.strain * 0.0;

        #HenckyF = np.matrix( np.zeros(3,3))
        HenckyF = numpy.zeros( (3,3) )

        for i in range(0,3):
            for j in range(0,3):
                for k in range(0,3):
                    HenckyF[i,k] = HenckyF[i,k] + F[i,j]*F[k,j];


        values,vectors = numpy.linalg.eig( HenckyF)

        valuesMatrix = numpy.zeros( (3,3) )
        for i in range(0,3):
            valuesMatrix[i,i] = numpy.log( values[i] )/2.0


        Hencky = numpy.zeros( (3,3) )

        Hencky = numpy.matmul( vectors, valuesMatrix)
        Hencky = numpy.matmul( Hencky, vectors.transpose() )

        strain = numpy.zeros( 6)

        for i in range(0,3):
            strain[i] = Hencky[i,i]

        strain[3] = 2.0*Hencky[0,1]
        strain[4] = 2.0*Hencky[0,2]
        strain[5] = 2.0*Hencky[1,2]

        return strain



    def _OpenOutputFile(self, FileName):


        import os.path

        FileName = str(self.current_case)+ "-"+ self.NameOfSeries+ "-" + FileName
        problem_path = os.getcwd()
        self.csv_path = os.path.join(problem_path, FileName)

        csv_file = open(self.csv_path, "w")
        csv_file.close()

        self.t0 = timer.clock()


    def _WriteThisToFile(self, t, stress, strain):

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

        time = timer.clock() - self.t0
    
        line = line + " , " +  str(time)

        line = line + " \n"
        csv_file = open(self.csv_path, "a")
        csv_file.write(line)
        csv_file.close()
 

if __name__ == '__main__':
    KratosUnittest.main()

import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterialModels
import importlib


def CreateProcess(custom_settings):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return TestConstitutiveModelProcess(custom_settings)

class TestConstitutiveModelProcess(KratosMultiphysics.Process):
    def __init__(self, custom_settings ):

        KratosMultiphysics.Process.__init__(self)
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
	    "model_part_name" : "MaterialDomain",
	    "properties_id"   : 1,
            "material_name"   : "steel",
	    "constitutive_law": {
                "law_name"   : "LargeStrain3DLaw",
		"model_name" : "SaintVenantKirchhoff"
            },
	    "variables": {},
	    "tables": {},
	    "element_type": "Tetrahedra3D4N",
            "nodes" : [],
            "strain": {
		"deformation_gradient" : [ [1,0,0], [0,1,0], [0,0,1] ],
		"jacobian": 1.0
            },
            "echo_level" : 0
        }
        """)

        #overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #build model part and element
        self.model_part = KratosMultiphysics.ModelPart(self.settings["model_part_name"].GetString())
        self.echo_level = self.settings["echo_level"].GetInt()
        
        #read nodes
        self.number_of_nodes = self.settings["nodes"].size()
        self.nodes = [] #self.model_part.GetNodes()
        for i in range(0, self.number_of_nodes):
            node = self.model_part.CreateNewNode(i,self.settings["nodes"][i][0].GetDouble(),self.settings["nodes"][i][1].GetDouble(),self.settings["nodes"][i][2].GetDouble())
            self.nodes.append(node)

        self.geometry  = KratosMultiphysics.Geometry()
        self.dimension = 3
        if( self.settings["element_type"].GetString() == "Tetrahedra3D4"):
            if( self.number_of_nodes != 4 ):
                print(" number of nodes:",self.number_of_nodes," do not matches geometry :", self.settings["element_type"].GetString() )
            else:
                self.geometry = KratosMultiphysics.Tetrahedra3D4(self.nodes[0],self.nodes[1],self.nodes[2],self.nodes[3])
                print(" geometry ",self.geometry)
                
        if( self.settings["element_type"].GetString() == "Triangle2D3"):
            if( self.number_of_nodes != 4 ):
                print(" number of nodes:",self.number_of_nodes," do not matches geometry :", self.settings["element_type"].GetString() )
            else:
                self.geometry  = KratosMultiphysics.Triangle2D3(self.nodes[0],self.nodes[1],self.nodes[2])
                self.dimension = 2
                
        #material properties
        self.properties = self.model_part.Properties[self.settings["properties_id"].GetInt()]

        self.properties.SetValue(KratosMultiphysics.YOUNG_MODULUS,100)
        #read variables
        self.variables = self.settings["variables"]
        for key, value in self.variables.items():
            variable = self._GetItemFromModule(key)
            if( value.IsDouble() ):
                self.properties.SetValue(variable, value.GetDouble())
            elif( value.IsArray() ):
                vector_value = KratosMultiphysics.Vector(value.size())
                for i in range(0, value.size() ):
                    vector_value[i] = value[i].GetDouble()
                self.properties.SetValue(variable, vector_value)


        #read table
        self.tables  = self.settings["tables"]
        for key, table in self.tables.items():
            table_name = key
            input_variable  = self._GetItemFromModule(table["input_variable"].GetString())
            output_variable = self._GetItemFromModule(table["output_variable"].GetString())

            new_table = KratosMultiphysics.PiecewiseLinearTable()
            for i in range(0, table["data"].size() ):
                new_table.AddRow(table["data"][i][0].GetDouble(), table["data"][i][1].GetDouble())
                
            self.properties.SetTable(input_variable,output_variable,new_table)

        
        #create constitutive law
        self.material_model = self._GetItemFromModule( self.settings["constitutive_law"]["model_name"].GetString())()
        self.material_law   = self._GetItemFromModule( self.settings["constitutive_law"]["law_name"].GetString())(self.material_model)

        #check dimension
        if(self.material_law.WorkingSpaceDimension() != self.dimension):
            raise Exception( "mismatch between the ConstitutiveLaw dimension and the dimension of the space")

        #set strain
        self.F = KratosMultiphysics.Matrix(3,3)
        self.strain_measure = self.settings["strain"]["deformation_gradient"]
        
        self.F[0,0] = self.strain_measure[0][0].GetDouble()
        self.F[0,1] = self.strain_measure[0][1].GetDouble()
        self.F[0,2] = self.strain_measure[0][2].GetDouble()
        self.F[1,0] = self.strain_measure[1][0].GetDouble()
        self.F[1,1] = self.strain_measure[1][1].GetDouble()
        self.F[1,2] = self.strain_measure[1][2].GetDouble()
        self.F[2,0] = self.strain_measure[2][0].GetDouble()
        self.F[2,1] = self.strain_measure[2][1].GetDouble()
        self.F[2,2] = self.strain_measure[2][2].GetDouble()
        
        self.detF = self.settings["strain"]["jacobian"].GetDouble()        
        
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
        self.strain_vector = 0.0 * self.strain_vector
        self.parameters.SetStrainVector( self.strain_vector )
        self.parameters.SetStressVector( self.stress_vector )
        self.parameters.SetConstitutiveMatrix( self.constitutive_matrix )
        self.parameters.SetShapeFunctionsValues( self.N )
        self.parameters.SetShapeFunctionsDerivatives( self.DN_DX )
        self.parameters.SetProcessInfo( self.model_part.ProcessInfo )
        self.parameters.SetMaterialProperties( self.properties )
        self.parameters.SetElementGeometry( self.geometry )

        
    #
    def ExecuteInitialize(self):
        pass

    #
    def Execute(self):

        #check parameters
        self.parameters.CheckAllParameters()
        self.parameters.CheckMechanicalVariables()
        self.parameters.CheckShapeFunctions()
        
        self.CalculateMaterialResponse()

    #
    def ExecuteFinalize(self):

        self.echo_level = 1
        self.CalculateMaterialResponse()

        

    #
    def CalculateMaterialResponse(self):
    
        
        self.material_law.FinalizeMaterialResponseCauchy( self.parameters )
        self.material_law.FinalizeSolutionStep( self.properties, self.geometry, self.N, self.model_part.ProcessInfo )


        import numpy as np
        import matplotlib.pyplot as plt

        #isotropic compression

        XX = self.parameters.GetStrainVector()
        
        sigma_v = -1500.0
        k0 = 0.5
        sigma_h = sigma_v * k0

        for iteracio in range(0, 100):
            stress = self.parameters.GetStressVector()

            residual = 0.0*stress + stress;
            
            residual[0] = stress[0] - sigma_h
            residual[2] = stress[2] - sigma_h
            residual[1] = stress[1] - sigma_v
            
            #Compute the norm of the residual
            norm = 0
            for i in range(0,6):
                norm += residual[i]*residual[i]
            #print('iteracio')
            #print(iteracio)
            #print('norm')
            #print(norm)
            if (norm < 1e-9):
                break
            if ( iteracio > 95):
                pNode


            C = self.parameters.GetConstitutiveMatrix()
            SystemMatrix = np.zeros((6,6));
            for i in range(0,6):
                for j in range(0,6):
                    SystemMatrix[i,j] = C[i,j];
            dx = np.linalg.solve(SystemMatrix, residual)
            for i in range(0,6):
                XX[i] = XX[i] - dx[i];
                self.strain_vector[i] = XX[i]
                self.parameters.SetStrainVector(self.strain_vector)
            self.material_law.CalculateMaterialResponseCauchy( self.parameters)


        self.material_law.FinalizeMaterialResponseCauchy( self.parameters )
        self.material_law.FinalizeSolutionStep( self.properties, self.geometry, self.N, self.model_part.ProcessInfo )
        stress = self.parameters.GetStressVector();
        strain = self.parameters.GetStrainVector();

        isotropicLoadingStrain = strain[1]
        
        self._OpenOutputFile()
        self._WriteThisToFile(0, stress, strain)

        # Second part
        nLoadingSteps = 1000
        FinalAxialDeformation = 0.5

        VolStrain = np.arange( float(nLoadingSteps+1) )
        DevStrain = np.arange( float(nLoadingSteps+1) )
        VolStress = np.arange( float(nLoadingSteps+1) )
        DevStress = np.arange( float(nLoadingSteps+1) )

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
                residual[1] = strain[1] - verticalStrain
            
                #Compute the norm of the residual
                norm = 0
                for i in range(0,6):
                    norm += residual[i]*residual[i]
                #print('iteracio')
                #print(iteracio)
                #print('norm')
                #print(norm)
                if (norm < 1e-12):
                    break
                if (iteracio > 90):
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
                    self.strain_vector[i] = XX[i]
                    self.parameters.SetStrainVector(self.strain_vector)
                self.material_law.CalculateMaterialResponseCauchy( self.parameters)


            self.material_law.FinalizeMaterialResponseCauchy( self.parameters )
            self.material_law.FinalizeSolutionStep( self.properties, self.geometry, self.N, self.model_part.ProcessInfo )
            stress = self.parameters.GetStressVector();
            strain = self.parameters.GetStrainVector();
        
            self._WriteThisToFile(t, stress, strain)

            #print('Stress')
            #print(stress)
            #print('Strain')
            #print(strain)

            VolStrain[t] = strain[0] + strain[1] + strain[2]
            VolStrain[t] *= -1.0
            DevStrain[t] = strain[1] - strain[0]
            DevStrain[t] *= -1.0

            VolStress[t] = stress[0] + stress[1] + stress[2]
            VolStress[t] /= -3.0
            DevStress[t] = stress[1]-stress[0]
            DevStress[t] *= -1.0

        plt.subplot(2,2,2)
        plt.plot(DevStrain, DevStress, 'ro-')
        plt.subplot(2,2,1)
        plt.plot(VolStress, DevStress, 'ro-')
        
        plt.subplot(2,2,3)
        plt.plot(VolStress, VolStrain, 'ro-')
        plt.subplot(2,2,4)
        plt.plot(DevStrain, VolStrain, 'ro-')
        plt.show()
        plt.show(block=False)

        
    def _GetItemFromModule(self,my_string):
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

    def _OpenOutputFile(self):
        import os.path

        problem_path = os.getcwd()
        self.csv_path = os.path.join(problem_path, "GaussPoint.csv")

        csv_file = open(self.csv_path, "w")
        csv_file.close()

    def _WriteThisToFile(self, t, stress, strain):

        line = str(t) + " , "
        for st in stress:
            line = line + str(st) + " , "

        for st in strain:
            line = line + str(st) + " , "

        line = line[:-2]
        line = line + "\n"
        csv_file = open(self.csv_path, "a")
        csv_file.write(line)
        csv_file.close()
 

import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterials
import importlib

def Factory(custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignMaterialsProcess(Model, custom_settings["Parameters"])

class AssignMaterialsProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
	    "model_part_name" : "MaterialDomain",
	    "properties_id"   : 1,
            "material_name"   : "steel",
	    "constitutive_law": {
                "name"   : "KratosMultiphysics.ConstitutiveModelsApplication.LargeStrain3DLaw.LinearElasticModel"
            },
	    "variables": {},
	    "tables": {},
            "echo_level" : 0
        }
        """)

        #overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #build model part and element
        self.model_part = Model[self.settings["model_part_name"].GetString()]
        self.echo_level = self.settings["echo_level"].GetInt()
        self.material_name =  self.settings["material_name"].GetString()

        #material properties
        self.main_model_part = self.model_part.GetRootModelPart()
        self.properties      = self.main_model_part.Properties[self.settings["properties_id"].GetInt()]
        #read variables
        self.variables = self.settings["variables"]
        for key, value in self.variables.items():
            try:
                my_key = "KratosMultiphysics."+key
                variable = self._GetItemFromModule(my_key)
            except:
                my_key = "KratosMultiphysics.SolidMechanicsApplication."+key
                variable = self._GetItemFromModule(my_key)

            if( value.IsDouble() ):
                self.properties.SetValue(variable, value.GetDouble())
            elif( value.IsInt() ):
                self.properties.SetValue(variable, value.GetInt())
            elif( value.IsArray() ):
                vector_value = KratosMultiphysics.Vector(value.size())
                for i in range(0, value.size() ):
                    vector_value[i] = value[i].GetDouble()
                self.properties.SetValue(variable, vector_value)


        #read table
        self.tables = self.settings["tables"]
        properties_layout = KratosMaterials.PropertiesLayout()
        number_of_tables = 0
        for key, table in self.tables.items():
            number_of_tables += 1
            table_name = key
            print(" table",key)
            if table.Has("table_file_name"):
                import os
                problem_path = os.getcwd()
                table_path = os.path.join(problem_path, table["table_file_name"].GetString() )
                import csv
                with open(table_path, 'r') as table_file:
                    reader = csv.DictReader(table_file)
                    new_table = KratosMultiphysics.PiecewiseLinearTable()
                    input_variable_name  = reader.fieldnames[0]
                    output_variable_name = reader.fieldnames[1]
                    input_variable  = self._GetItemFromModule("KratosMultiphysics."+str(input_variable_name))
                    output_variable = self._GetItemFromModule("KratosMultiphysics."+str(output_variable_name))
                    for line in reader:
                        new_table.AddRow(float(line[input_variable_name]), float(line[output_variable_name]))
                    self.properties.SetTable(input_variable,output_variable,new_table)
                    properties_layout.RegisterTable(input_variable,output_variable)
            else:
                input_variable  = self._GetItemFromModule(table["input_variable"].GetString())
                output_variable = self._GetItemFromModule(table["output_variable"].GetString())
                new_table = KratosMultiphysics.PiecewiseLinearTable()
                for i in range(0, table["data"].size() ):
                    new_table.AddRow(table["data"][i][0].GetDouble(), table["data"][i][1].GetDouble())
                self.properties.SetTable(input_variable,output_variable,new_table)
                properties_layout.RegisterTable(input_variable,output_variable)

        #set properties layout and table keys and arguments
        if number_of_tables > 0:
            properties_layout.AddToProperties(KratosMaterials.PROPERTIES_LAYOUT, properties_layout, self.properties)

        #create constitutive law
        self.material_law = self._GetLawFromModule(self.settings["constitutive_law"]["name"].GetString())

        self.properties.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, self.material_law.Clone())

        #self.model_part.SetProperties(self.main_model_part.GetProperties())
        self.model_part.AddProperties(self.properties)

        splitted_law_name = (self.settings["constitutive_law"]["name"].GetString()).split(".")

        print("::[------Material-----]:: "+self.material_name+" [Model: "+splitted_law_name[len(splitted_law_name)-1]+"]")

    #
    def ExecuteInitialize(self):
        pass

    #
    def Execute(self):

        self._AssignMaterialProperties()

    #
    def ExecuteFinalize(self):
        pass

    #
    def _AssignMaterialProperties(self):

        # Check dimension
        self.dimension = self.model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]

        if(self.material_law.WorkingSpaceDimension() != self.dimension):
            #feature flags
            self.features =KratosMultiphysics.ConstitutiveLawFeatures()
            self.material_law.GetLawFeatures(self.features)
            self.feature_options = self.features.GetOptions()
            if( self.feature_options.IsNot(KratosMultiphysics.ConstitutiveLaw.PLANE_STRESS_LAW) ):
                raise Exception("mismatch between the ConstitutiveLaw dimension and the dimension of the space")

        # Assign properties to the model_part elements
        for Element in self.model_part.Elements:
            Element.Properties = self.properties

        # Assign properties to the model_part conditions
        for Condition in self.model_part.Conditions:
            Condition.Properties = self.properties


    def _GetLawFromModule(self,my_string):
        splitted = my_string.split(".")
        if(len(splitted) == 0):
            raise Exception("something wrong. Trying to split the string "+my_string)
        if(len(splitted) == 1):
            material_law = "KratosMaterials."+my_string+"()"
            return eval(material_law)
        elif(len(splitted) == 2):
            material_law = "KratosMaterials."+splitted[0]+"(KratosMaterials."+splitted[1]+"())"
            return eval(material_law)
        elif(len(splitted) == 3):
            module_name = ""
            for i in range(len(splitted)-1):
                module_name += splitted[i]
                if i != len(splitted)-2:
                    module_name += "."
            module = importlib.import_module(module_name)
            material_law = module_name+"."+splitted[-1]+"()"
            return eval(material_law)
        elif(len(splitted) == 4):
            module_name = ""
            for i in range(len(splitted)-2):
                module_name += splitted[i]
                if i != len(splitted)-3:
                    module_name += "."
            module = importlib.import_module(module_name)
            material_law = module_name+"."+splitted[-2]+"("+module_name+"."+splitted[-1]+"())"
            return eval(material_law)
        else:
            raise Exception("something wrong .. Trying to split the string "+my_string)



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
            try:
                return getattr(module,splitted[-1])
            except:
                raise

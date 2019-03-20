import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterialModels
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
        self.properties = self.model_part.Properties[self.settings["properties_id"].GetInt()]

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
        self.material_law = self._GetLawFromModule(self.settings["constitutive_law"]["name"].GetString())

        self.properties.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, self.material_law.Clone())


    #
    def ExecuteInitialize(self):
        pass

    #
    def Execute(self):

        self._AssignMaterialProperties()

        print(" Material ", self.material_name, " assigned " )

    #
    def ExecuteFinalize(self):
        pass


    #
    def _AssignMaterialProperties(self):

        # Check dimension
        self.dimension = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        if(self.material_law.WorkingSpaceDimension() != self.dimension):
            raise Exception( "mismatch between the ConstitutiveLaw dimension and the dimension of the space")


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
            material_law = "KratosMaterialModels."+my_string+"()"
            return eval(material_law)
        elif(len(splitted) == 2):
            material_law = "KratosMaterialModels."+splitted[0]+"(KratosMaterialModels."+splitted[1]+"())"
            return eval(material_law)
        elif(len(splitted) == 3):
            module_name = ""
            for i in range(len(splitted)-1):
                module_name += splitted[i]
                if i != len(splitted)-2:
                    module_name += "."
            module = importlib.import_module(module_name)
            return getattr(module,splitted[-1])
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
            return getattr(module,splitted[-1])



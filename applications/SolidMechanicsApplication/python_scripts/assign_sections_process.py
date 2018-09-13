import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import importlib
import math

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignSectionsProcess(Model, custom_settings["Parameters"])

class AssignSectionsProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
	    "model_part_name" : "MaterialDomain",
	    "properties_id"   : 1,
            "material_name"   : "steel",
	    "section_type"    : "Rectangular",
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
            try:
                my_key = "KratosMultiphysics."+key
                variable = self._GetItemFromModule(my_key)
            except:
                my_key = "KratosMultiphysics.SolidMechanicsApplication."+key
                variable = self._GetItemFromModule(my_key)

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


        #set section characteristics
        self.properties = self._AddSectionProperties( self.properties, self.variables, self.settings["section_type"].GetString() )


    #
    def ExecuteInitialize(self):
        pass

    #
    def Execute(self):

        self._AssignMaterialProperties()
        print("::[--Section_Assigned-]::", self.settings["section_type"].GetString())

    #
    def ExecuteFinalize(self):
        pass


    #
    def _AssignMaterialProperties(self):

        # Assign properties to the model_part elements
        for Element in self.model_part.Elements:
            Element.Properties = self.properties

        # Assign properties to the model_part conditions
        for Condition in self.model_part.Conditions:
            Condition.Properties = self.properties


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


    def _AddSectionProperties(self,properties,variables,section_type):

        if( (section_type == "IPN") or (section_type == "IPE") or (section_type == "HEB") or (section_type == "HEA") or (section_type == "HEM") or (section_type == "UPN") ):

            size = str(int(self._GetScalarVariableValue("SECTION_SIZE")))

            import os
            csv_file = os.path.dirname(__file__) + '/beam_profiles.csv'

            if( ( (section_type == "HEB") or (section_type == "HEA") or (section_type == "HEM") ) and ( int(self._GetScalarVariableValue("SECTION_SIZE")) < 100 ) ):
                raise Exception(" HEB, HEA, HEM profiles range is [100, 600] ")

            section_properties = self._SearchCSVProperties(csv_file, section_type, size)

            area = float(section_properties["A(m2)"])
            inertia_matrix = KratosMultiphysics.Matrix(2, 2)
            inertia_matrix[0, 0] = float(section_properties["Iz(m4)"])  # z is the horizontal axis of the section
            inertia_matrix[0, 1] = float(section_properties["Iz(m4)"]) + float(section_properties["Iy(m4)"])
            inertia_matrix[1, 0] = float(section_properties["Iz(m4)"]) + float(section_properties["Iy(m4)"])
            inertia_matrix[1, 1] = float(section_properties["Iy(m4)"])  # y is the vertical axis of the section

            properties.SetValue(KratosMultiphysics.LOCAL_INERTIA_TENSOR, inertia_matrix)
            properties.SetValue(KratosSolid.CROSS_SECTION_AREA, area)
            mean_radius = float(size) * 0.5
            properties.SetValue(KratosSolid.CROSS_SECTION_RADIUS, mean_radius)
            sides = 4
            properties.SetValue(KratosSolid.CROSS_SECTION_SIDES, sides)
            return properties

        elif( section_type == "Circular" ):

            diameter = self._GetScalarVariableValue("DIAMETER")
            radius = diameter*0.5

            circular_area = 3.14 * (radius ** 2)
            circular_inertia = (3.14 * (radius ** 4)) * 0.25
            circular_inertia_polar = circular_inertia + circular_inertia
            circular_module = (3.14 * (radius ** 3)) * 0.25
            circular_turning_radius = (((3.14 * (radius ** 4)) * 0.25) / (3.14 *(radius**2)))**(0.5)
            inertia_matrix = KratosMultiphysics.Matrix(2, 2)
            inertia_matrix[0, 0] = circular_inertia
            inertia_matrix[0, 1] = circular_inertia_polar
            inertia_matrix[1, 0] = circular_inertia_polar
            inertia_matrix[1, 1] = circular_inertia

            properties.SetValue(KratosMultiphysics.LOCAL_INERTIA_TENSOR, inertia_matrix)
            properties.SetValue(KratosSolid.CROSS_SECTION_AREA, circular_area)
            properties.SetValue(KratosSolid.CROSS_SECTION_RADIUS, radius)
            sides = 25
            properties.SetValue(KratosSolid.CROSS_SECTION_SIDES, sides)
            return properties

        elif( section_type == "Rectangular" ):

            height_square = self._GetScalarVariableValue("SECTION_HEIGHT")
            base_square   = self._GetScalarVariableValue("SECTION_WIDTH")

            square_area = base_square * height_square
            square_inertia_z = (base_square * height_square ** 3) / 12.0
            square_inertia_y = (height_square * base_square ** 3) / 12.0
            square_inertia_polar = square_inertia_z + square_inertia_y
            square_module_z = (base_square * height_square ** 2) / 6.0
            square_module_y = (height_square * base_square ** 2) / 6.0
            square_turning_radius_z = ((((base_square * height_square ** 3) / 12.0) / (base_square * height_square)) ** (0.5))
            square_turning_radius_y = ((((height_square * base_square ** 3) / 12.0) / (base_square * height_square)) ** (0.5))
            inertia_matrix = KratosMultiphysics.Matrix(2, 2)
            inertia_matrix[0, 0] = square_inertia_z
            inertia_matrix[0, 1] = square_inertia_polar
            inertia_matrix[1, 0] = square_inertia_polar
            inertia_matrix[1, 1] = square_inertia_y

            properties.SetValue(KratosMultiphysics.LOCAL_INERTIA_TENSOR, inertia_matrix)
            properties.SetValue(KratosSolid.CROSS_SECTION_AREA, square_area)
            mean_radius = math.sqrt(square_area)
            properties.SetValue(KratosSolid.CROSS_SECTION_RADIUS, mean_radius)
            sides = 4
            properties.SetValue(KratosSolid.CROSS_SECTION_SIDES, sides)
            return properties

        elif( section_type == "Tubular" ):

            diameter_D = self._GetScalarVariableValue("DIAMETER")
            thickness  = self._GetScalarVariableValue("THICKNESS")

            radius = diameter_D * 0.5

            diameter_d = diameter_D - 2*thickness

            distance_a = (diameter_D ** 2) - (diameter_d ** 2)
            distance_b = (diameter_D ** 2) + (diameter_d ** 2)

            circular_area = 3.14 * (distance_a ** 2) * 0.25

            # for thin tubes
            #circular_inertia = (3.14 * (diameter_D ** 3) * thickness) * 0.25

            # for thick tubes
            circular_inertia = (3.14 * (radius ** 4 - (radius-thickness) ** 4) ) * 0.25

            circular_inertia_polar = circular_inertia + circular_inertia

            circular_module = (3.14 * (diameter_D ** 2) * thickness) * 0.25
            circular_turning_radius = math.sqrt(distance_b) * 0.25

            inertia_matrix = KratosMultiphysics.Matrix(2, 2)
            inertia_matrix[0, 0] = circular_inertia
            inertia_matrix[0, 1] = circular_inertia_polar
            inertia_matrix[1, 0] = circular_inertia_polar
            inertia_matrix[1, 1] = circular_inertia

            properties.SetValue(KratosMultiphysics.LOCAL_INERTIA_TENSOR, inertia_matrix)
            properties.SetValue(KratosSolid.CROSS_SECTION_AREA, circular_area)
            properties.SetValue(KratosSolid.CROSS_SECTION_RADIUS, radius)
            sides = 25
            properties.SetValue(KratosSolid.CROSS_SECTION_SIDES, sides)
            return properties

        elif( section_type == "UserDefined" ):

            area = self._GetScalarVariableValue("CROSS_SECTION_AREA")
            inertia_z = self._GetScalarVariableValue("INERTIA_X")
            inertia_y = self._GetScalarVariableValue("INERTIA_Y")

            radius = math.sqrt(area)
            inertia_polar = inertia_z + inertia_y

            inertia_matrix = KratosMultiphysics.Matrix(2, 2)
            inertia_matrix[0, 0] = inertia_z
            inertia_matrix[0, 1] = inertia_polar
            inertia_matrix[1, 0] = inertia_polar
            inertia_matrix[1, 1] = inertia_y

            properties.SetValue(KratosMultiphysics.LOCAL_INERTIA_TENSOR, inertia_matrix)
            properties.SetValue(KratosSolid.CROSS_SECTION_AREA, area)
            properties.SetValue(KratosSolid.CROSS_SECTION_RADIUS, radius)
            sides = 25
            properties.SetValue(KratosSolid.CROSS_SECTION_SIDES, sides)

            return properties

        elif( section_type == "UserParameters" ):

            ConstitutiveMatrix = KratosMultiphysics.Matrix(6,6)

            for i in range(0,6):
                for j in range(0,6):
                    ConstitutiveMatrix[i,j] = 0

            ConstitutiveMatrix[0,0] = self._GetScalarVariableValue("SHEARxREDUCED_AREA")  # GAy
            ConstitutiveMatrix[1,1] = self._GetScalarVariableValue("SHEARxREDUCED_AREA")  # GAz
            ConstitutiveMatrix[2,2] = self._GetScalarVariableValue("YOUNGxAREA")  # EA

            ConstitutiveMatrix[3,3] = self._GetScalarVariableValue("YOUNGxINERTIA_X")  # EIy
            ConstitutiveMatrix[4,4] = self._GetScalarVariableValue("YOUNGxINERTIA_Y")  # EIz
            ConstitutiveMatrix[5,5] = self._GetScalarVariableValue("SHEARxPOLAR_INERTIA")  # GJ

            properties.SetValue(KratosMultiphysics.LOCAL_CONSTITUTIVE_MATRIX, ConstitutiveMatrix)
            return properties

    def _GetScalarVariableValue(self,variable_name):

        variable = None
        for key, value in self.variables.items():
            if( key == variable_name ):
                variable = value.GetDouble()
                break

        if( variable == None ):
            raise Exception( "variable", variable_name, "not found in section properties")

        return variable


    def _SearchCSVProperties(self,file_name,shape,size):
        separator = ";"

        print(" shape", shape, "size", size)

        with open(file_name, "rU") as f:
            header = f.readline()
            header = header.rstrip("\n")
            header = header.split(separator)
            line = f.readline()

            while line != "":
                line = line.rstrip("\n")
                line = line.split(separator)
                if (line[0] == shape and line[1] == size):
                    result = dict(list(zip(header, line)))
                    break
                line = f.readline()
                print(" line",line)

        if line == "":
            raise Exception(" Section",shape," size", size," not found in section definition files" )

        return result

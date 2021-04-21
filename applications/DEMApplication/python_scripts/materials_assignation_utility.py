import KratosMultiphysics as Kratos
from KratosMultiphysics.DEMApplication import *

class MaterialsAssignationUtility():

    def __init__(self, model, spheres_model_part, DEM_material_parameters):
        self.model = model
        self.spheres_model_part = spheres_model_part
        self.DEM_material_parameters = DEM_material_parameters

    def AssignMaterialParametersToProperties(self):

        materials_parameters = self.DEM_material_parameters
        list_of_materials = materials_parameters["materials"]
        list_of_material_relations = materials_parameters["material_relations"]

        for material in list_of_materials:
            material_id = material["material_id"].GetInt()
            if  not self.spheres_model_part.HasProperties(material_id):
                self.spheres_model_part.AddProperties(Kratos.Properties(material_id))

            properties_of_model_part_with_this_id = self.spheres_model_part.GetProperties()[material_id]
            properties = material["properties"]
            properties_of_model_part_with_this_id[Kratos.PARTICLE_MATERIAL] = material_id
            self.SetInputMaterialPropertiesIntoModelPartProperties(properties, properties_of_model_part_with_this_id)

            for material_relation in list_of_material_relations:
                subprops = None
                material_ids_list = material_relation["material_ids_list"].GetVector()
                if material_id == material_ids_list[0]:
                    index_of_the_other_material = int(material_ids_list[1])
                    subprops = Kratos.Properties(index_of_the_other_material)
                elif material_id == material_ids_list[1]:
                    index_of_the_other_material = int(material_ids_list[0])
                    subprops = Kratos.Properties(index_of_the_other_material)

                if subprops:
                    contact_properties = material_relation["properties"]
                    self.SetInputMaterialRelationPropertiesIntoModelPartSubProperties(contact_properties, subprops)
                    properties_of_model_part_with_this_id.AddSubProperties(subprops)


    def SetInputMaterialPropertiesIntoModelPartProperties(self, input_properties, mp_properties):
        if input_properties.Has("PARTICLE_DENSITY"):
            mp_properties[PARTICLE_DENSITY] = input_properties["PARTICLE_DENSITY"].GetDouble()
        mp_properties[Kratos.YOUNG_MODULUS] = input_properties["YOUNG_MODULUS"].GetDouble()
        mp_properties[Kratos.POISSON_RATIO] = input_properties["POISSON_RATIO"].GetDouble()
        if input_properties.Has("COMPUTE_WEAR"):
            mp_properties[COMPUTE_WEAR] = input_properties["COMPUTE_WEAR"].GetBool()
        else:
            mp_properties[COMPUTE_WEAR] = False


    def SetInputMaterialRelationPropertiesIntoModelPartSubProperties(self, contact_properties, subprops):
        subprops[COEFFICIENT_OF_RESTITUTION] = contact_properties["COEFFICIENT_OF_RESTITUTION"].GetDouble()
        subprops[STATIC_FRICTION] = contact_properties["STATIC_FRICTION"].GetDouble()
        subprops[DYNAMIC_FRICTION] = contact_properties["DYNAMIC_FRICTION"].GetDouble()
        subprops[FRICTION_DECAY] = contact_properties["FRICTION_DECAY"].GetDouble()
        subprops[ROLLING_FRICTION] = contact_properties["ROLLING_FRICTION"].GetDouble()
        subprops[ROLLING_FRICTION_WITH_WALLS] = contact_properties["ROLLING_FRICTION_WITH_WALLS"].GetDouble()
        if contact_properties.Has("SEVERITY_OF_WEAR"):
            subprops[SEVERITY_OF_WEAR] = contact_properties["SEVERITY_OF_WEAR"].GetDouble()
        if contact_properties.Has("IMPACT_WEAR_SEVERITY"):
            subprops[IMPACT_WEAR_SEVERITY] = contact_properties["IMPACT_WEAR_SEVERITY"].GetDouble()
        if contact_properties.Has("BRINELL_HARDNESS"):
            subprops[BRINELL_HARDNESS] = contact_properties["BRINELL_HARDNESS"].GetDouble()
        if contact_properties.Has("CONICAL_DAMAGE_CONTACT_RADIUS"):
            subprops[CONICAL_DAMAGE_CONTACT_RADIUS] = contact_properties["CONICAL_DAMAGE_CONTACT_RADIUS"].GetDouble()
        if contact_properties.Has("CONICAL_DAMAGE_MAX_STRESS"):
            subprops[CONICAL_DAMAGE_MAX_STRESS] = contact_properties["CONICAL_DAMAGE_MAX_STRESS"].GetDouble()
        if contact_properties.Has("CONICAL_DAMAGE_ALPHA"):
            subprops[CONICAL_DAMAGE_ALPHA] = contact_properties["CONICAL_DAMAGE_ALPHA"].GetDouble()
        if contact_properties.Has("CONICAL_DAMAGE_GAMMA"):
            subprops[CONICAL_DAMAGE_GAMMA] = contact_properties["CONICAL_DAMAGE_GAMMA"].GetDouble()
        if contact_properties.Has("LEVEL_OF_FOULING"):
            subprops[LEVEL_OF_FOULING] = contact_properties["LEVEL_OF_FOULING"].GetDouble()
        if contact_properties.Has("PARTICLE_COHESION"):
            subprops[PARTICLE_COHESION] = contact_properties["PARTICLE_COHESION"].GetDouble()
        if contact_properties.Has("INITIAL_COHESION"):
            subprops[INITIAL_COHESION] = contact_properties["INITIAL_COHESION"].GetDouble()
        if contact_properties.Has("AMOUNT_OF_COHESION_FROM_STRESS"):
            subprops[AMOUNT_OF_COHESION_FROM_STRESS] = contact_properties["AMOUNT_OF_COHESION_FROM_STRESS"].GetDouble()

        subprops[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME] = contact_properties["DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME"].GetString()
        discontinuum_constitutive_law_instance = globals().get(subprops[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME])()
        discontinuum_constitutive_law_instance.SetConstitutiveLawInProperties(subprops, True)

        if contact_properties.Has("continuum_contact_law_parameters"):
            subprops[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME] = contact_properties["continuum_contact_law_parameters"]["DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME"].GetString()
            continuum_constitutive_law_instance = globals().get(subprops[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME])()
            continuum_constitutive_law_instance.SetConstitutiveLawInPropertiesWithParameters(subprops, contact_properties["continuum_contact_law_parameters"], True)


    def AssignPropertiesToEntities(self):
        materials_parameters = self.DEM_material_parameters
        list_of_materials = materials_parameters["materials"]
        material_assignation_table = materials_parameters["material_assignation_table"]

        for pair in material_assignation_table:
            submodelpart_name_in_assignation_table = pair[0].GetString()
            material_name_in_assignation_table = pair[1].GetString()
            submodelpart = None
            for material in list_of_materials:
                material_name_in_materials_list = material["material_name"].GetString()
                if material_name_in_assignation_table == material_name_in_materials_list:
                    submodelpart = self.model.GetModelPart(submodelpart_name_in_assignation_table)
                    material_id = material["material_id"].GetInt()
                    props = self.spheres_model_part.GetProperties()[material_id]
            for element in submodelpart.Elements:
                element.Properties = props
            for condition in submodelpart.Conditions:
                condition.Properties = props
from jsonschema import validate, ValidationError
import json

class Schema_Validation:
    def __init__(self, config_data):
        self.config_data = self.param_to_dict(config_data)
        
        self.create_schemas()
        self.create_type_schemas()
        self.validate_structural_element_input()

    def param_to_dict(self, params):
        return json.loads(params.WriteJsonString())

    def create_schemas(self):
        """Add templates for the configuration file here.
        """

        self.panel_schema = {
            "type": "object",
            "properties": {
                "type": {"const": "Panel"},
                "submodelpart": {"type": "string"},
                #"panel_origin_node": {"type": "integer"},
                #"corner_node_x" : {"type": "integer"},
                #"corner_node_y" : {"type": "integer"},
                "boundary_conditions" : {"type": "array", "items" : {"type" : "number"}},
                "analysis_methods"  :   {"type": "array", "items": {"type": "string"}}
            },
            "required": ["type", "submodelpart", "analysis_methods"]
        }

        self.column_schema = {
            "type": "object",
            "properties": {
                "type": {"const": "Column"},
                "submodelpart": {"type": "string"},
                "euler_case": {"type": "array", "items" : {"type" : "number"}},
                "analysis_methods": {"type": "array", "items": {"type": "string"}}
            },
            "required": ["type", "submodelpart", "analysis_methods"]
        }

    def create_type_schemas(self):
        self.type_schemas = {
            "Panel": self.panel_schema,
            "Column": self.column_schema
        }

    def validate_structural_element_input(self):
        structural_elements = self.config_data.get("Structural_Elements", [])
        for i, struct_elem in enumerate(structural_elements):
            element_type = struct_elem.get("type")
            schema = self.type_schemas.get(element_type)
            if not schema:
                raise ValueError(f"Unknown element type '{element_type}' in element #{i+1}")
            try:
                validate(instance=struct_elem, schema=schema)
            except ValidationError as e:
                raise ValidationError(f"Validation error in element #{i+1} ({element_type}): {e.message}")
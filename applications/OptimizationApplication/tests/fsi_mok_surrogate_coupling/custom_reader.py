import KratosMultiphysics

def ReadMaterials(model,materialSettings):
    settings = materialSettings["settings"]
    read_mat_util = KratosMultiphysics.ReadMaterialsUtility(model)
    read_mat_util.ReadMaterials(settings)
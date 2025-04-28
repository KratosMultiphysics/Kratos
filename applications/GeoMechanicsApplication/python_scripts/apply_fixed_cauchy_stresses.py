import KratosMultiphysics as Core

def Factory(settings, model):
    if not isinstance(settings, Core.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")

    model_part = model[settings["Parameters"]["model_part_name"].GetString()]
    return ApplyFixedCauchyStressesProcess(model_part, settings["Parameters"])


def read_plx_stress_table_txt(file_name):
    """
    Read the .txt table file and return the data as a dictionary
    with stress point numbers as keys.
    """
    with open(file_name, 'r', encoding="utf-8-sig") as file:
        lines = file.readlines()

    # Extract the header
    header = [variable.strip() for variable in lines[0].strip().split("\t")]

    # Extract the data
    data = {}

    plx_2_kratos_map = {1: 3, 2:1, 3:2}
    data["header"] = header
    for line in lines[1:]:
        values = line.strip().split("\t")
        if values[0].startswith("Clus."):
            element_number = int(values[0].split()[-1])
        local_no = plx_2_kratos_map[int(values[2])]

        if element_number not in data:
            data[element_number] = {}
        if local_no not in data[element_number]:
            data[element_number][local_no] = []

        data[element_number][local_no] = [1000*float(value) for value in values[5:9]]

    return data

class ApplyFixedCauchyStressesProcess(Core.Process):
    def __init__(self, model_part, settings):
        Core.Process.__init__(self)
        self.model_part = model_part
        self.settings = settings

    def ExecuteInitialize(self):
        data = read_plx_stress_table_txt(self.settings["file_name"].GetString())
        for element in self.model_part.Elements:
            plx_stress_vectors = data[element.Id]
            print(f"Element {element.Id} stress vectors: {plx_stress_vectors}")
            kratos_stress_vectors = [Core.Vector(plx_stress_vectors[i]) for i in range(1,4)]
            print(f"Element {element.Id} kratos stress vectors: {kratos_stress_vectors}")
            element.Initialize(self.model_part.ProcessInfo)
            element.SetValuesOnIntegrationPoints(Core.CAUCHY_STRESS_VECTOR, kratos_stress_vectors , 3, self.model_part.ProcessInfo)






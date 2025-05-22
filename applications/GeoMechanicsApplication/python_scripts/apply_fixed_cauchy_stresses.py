import KratosMultiphysics as Core

def Factory(settings, model): # noqa: N802
    if not isinstance(settings, Core.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")

    model_part = model[settings["Parameters"]["model_part_name"].GetString()]
    return ApplyFixedCauchyStressesProcess(model_part, settings["Parameters"])

class ApplyFixedCauchyStressesProcess(Core.Process):
    def __init__(self, model_part, settings):
        Core.Process.__init__(self)
        self.model_part = model_part
        self.settings = settings


    def ExecuteInitialize(self):
        fixed_stress_data = {
            1: {
                1: [-100.0, -100.0, -100.0, 0.0, 0.0, 0.0],
                2: [-100.0, -100.0, -100.0, 0.0, 0.0, 0.0],
                3: [-100.0, -100.0, -100.0, 0.0, 0.0, 0.0],
            },
            2: {
                1: [-100.0, -100.0, -100.0, 0.0, 0.0, 0.0],
                2: [-100.0, -100.0, -100.0, 0.0, 0.0, 0.0],
                3: [-100.0, -100.0, -100.0, 0.0, 0.0, 0.0]
            }
        }
        for element in self.model_part.Elements:
            kratos_stress_vectors = fixed_stress_data[element.Id]
            stress_vectors = [Core.Vector(kratos_stress_vectors[i]) for i in range(1,4)]
            element.Initialize(self.model_part.ProcessInfo)
            element.SetValuesOnIntegrationPoints(Core.CAUCHY_STRESS_VECTOR, stress_vectors, 3, self.model_part.ProcessInfo)

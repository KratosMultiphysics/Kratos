import KratosMultiphysics
from KratosMultiphysics import SLIP, NORMAL, TIME, DELTA_TIME, MOMENTUM

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication import decelerating_wall_process


class TestDeceleratingWallProcess(KratosUnittest.TestCase):

    @classmethod
    def _Rotate(cls, tangential, normal, cos, sin):
        vec = KratosMultiphysics.Vector(3)
        vec[0] = tangential*cos - normal*sin
        vec[1] = tangential*sin + normal*cos
        vec[2] = 0.0
        return vec

    @classmethod
    def _GenerateDomain(cls):
        model = KratosMultiphysics.Model()
        root_model_part = model.CreateModelPart("root")
        child_model_part = root_model_part.CreateSubModelPart("child")

        root_model_part.AddNodalSolutionStepVariable(MOMENTUM)
        root_model_part.AddNodalSolutionStepVariable(NORMAL)

        initial_mom = cls._Rotate(3.0, 5.0, 4/5, 3/5)
        child_model_part.CreateNewNode(1, 0.0, 0.0, 0.0).SetSolutionStepValue(MOMENTUM, initial_mom)
        child_model_part.CreateNewNode(2, 4.0, 3.0, 0.0).SetSolutionStepValue(MOMENTUM, initial_mom)
        child_model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], child_model_part.GetProperties()[0])

        root_model_part.ProcessInfo.SetValue(TIME, 0.0)
        root_model_part.ProcessInfo.SetValue(DELTA_TIME, 0.05)

        return model

    @classmethod
    def _SimulateStep(cls, modelpart, process):
        new_t = modelpart.ProcessInfo[TIME] + modelpart.ProcessInfo[DELTA_TIME]
        modelpart.CloneTimeStep(new_t)
        process.ExecuteInitializeSolutionStep()
        process.ExecuteFinalizeSolutionStep()
        process.ExecuteBeforeOutputStep()
        process.ExecuteAfterOutputStep()

    def testDeceleratingWallProcess(self):
        model = self._GenerateDomain()
        root_model_part = model["root"]
        child_model_part = model["root.child"]

        settings = KratosMultiphysics.Parameters("""
        {
            "Parameters": {
                "model_part_name": "root.child",
                "period": 0.95
            }
        }""")

        process = decelerating_wall_process.Factory(settings, model)
        process.Check()
        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()

        step = 0
        while True:
            self._SimulateStep(root_model_part, process)
            step += 1

            if(root_model_part.ProcessInfo[TIME] > 0.95):
                break

            self.assertFalse(child_model_part.GetNode(1).Is(SLIP), "Slip enabled too soon (time-step #{})".format(step))
            self.assertFalse(child_model_part.GetNode(2).Is(SLIP), "Slip enabled too soon (time-step #{})".format(step))

        self.assertTrue(child_model_part.GetNode(1).Is(SLIP))
        self.assertTrue(child_model_part.GetNode(2).Is(SLIP))

        expected = self._Rotate(3.0, 5.0e-3, 4/5, 3/5)

        self.assertVectorAlmostEqual(child_model_part.GetNode(1).GetSolutionStepValue(MOMENTUM), expected, 2)
        self.assertVectorAlmostEqual(child_model_part.GetNode(2).GetSolutionStepValue(MOMENTUM), expected, 2)


if __name__ == "__main__":
    KratosUnittest.main()

import KratosMultiphysics
from KratosMultiphysics import SLIP, NORMAL, TIME, DELTA_TIME, VELOCITY

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication import compressible_slip_wall_process


class TestCompressibleSlipWallProcess(KratosUnittest.TestCase):

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

        root_model_part.AddNodalSolutionStepVariable(VELOCITY)
        root_model_part.AddNodalSolutionStepVariable(NORMAL)

        initial_mom = cls._Rotate(3.0, 5.0, 4/5, 3/5)
        child_model_part.CreateNewNode(1, 0.0, 0.0, 0.0).SetSolutionStepValue(VELOCITY, initial_mom)
        child_model_part.CreateNewNode(2, 4.0, 3.0, 0.0).SetSolutionStepValue(VELOCITY, initial_mom)
        child_model_part.CreateNewCondition("LineCondition2D2N", 1, [1, 2], child_model_part.GetProperties()[0])

        root_model_part.ProcessInfo.SetValue(TIME, 0.0)
        root_model_part.ProcessInfo.SetValue(DELTA_TIME, 0.1)

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
                "interval" : [1, 2.5],
                "rampup_time": 0.95,
                "variable_name" : "VELOCITY"
            }
        }""")

        process = compressible_slip_wall_process.Factory(settings, model)
        process.Check()
        process.ExecuteInitialize()
        process.ExecuteBeforeSolutionLoop()

        stages_reached = {
            process.AWAITING: False,
            process.RAMP_UP: False,
            process.STEADY: False,
            process.FINISHED: False
        }

        for step in range(27):
            self._SimulateStep(root_model_part, process)
            time = root_model_part.ProcessInfo[TIME]

            node_1 = child_model_part.GetNode(1)
            node_2 = child_model_part.GetNode(2)

            stages_reached[process._stage] = True

            if(time <= 1.0):
                self.assertTrue(process._stage == process.AWAITING)
                self.assertTrue(node_1.IsNot(SLIP))
                self.assertTrue(node_2.IsNot(SLIP))

            elif(time <= 1.95):
                self.assertTrue(process._stage == process.RAMP_UP)
                self.assertTrue(node_1.IsNot(SLIP))
                self.assertTrue(node_2.IsNot(SLIP))

                decay_time = time - 1.0
                normal = 5.0 * process.decay_constant**(-decay_time/0.95)
                expected = self._Rotate(3.0, normal, 4/5, 3/5)
                self.assertVectorAlmostEqual(node_1.GetSolutionStepValue(VELOCITY), expected, 2)
                self.assertVectorAlmostEqual(node_2.GetSolutionStepValue(VELOCITY), expected, 2)

            elif(time <= 2.5):
                self.assertTrue(process._stage == process.STEADY)
                self.assertTrue(node_1.Is(SLIP))
                self.assertTrue(node_2.Is(SLIP))

            else:
                self.assertTrue(process._stage == process.FINISHED)
                self.assertTrue(node_1.IsNot(SLIP))
                self.assertTrue(node_2.IsNot(SLIP))

        for (stage, reached) in stages_reached.items():
            self.assertTrue(reached, msg="Failed to reach stage {}".format(stage))


if __name__ == "__main__":
    KratosUnittest.main()

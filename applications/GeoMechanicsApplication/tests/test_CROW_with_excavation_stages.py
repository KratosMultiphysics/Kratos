import os
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers
import test_helper



class KratosGeoMechanicsCrowValidationWithStages:
    def __init__(self):
        # fmt: off
        self.stages_info = {
            "initial_stage":      {"end_time": -1.0, "base_name": "1_Initial_stage"},
            "null_step":          {"end_time":  0.0, "base_name": "2_Null_step"},
            "wall_installation":  {"end_time":  1.0, "base_name": "3_Sheetpile_installation"},
            "first_excavation":   {"end_time":  2.0, "base_name": "4_First_excavation"},
            "strut_installation": {"end_time":  3.0, "base_name": "5_Anchor_installation"},
            "second_excavation":  {"end_time":  4.0, "base_name": "6_Second_excavation"},
            "third_excavation":   {"end_time":  5.0, "base_name": "7_Third_excavation"},
        }
        # fmt: on

    def run_simulation(self):
        project_path = (
                Path(__file__).resolve().parent
                / "crow_validation"
                / "with_excavation_stages"
        )

        with context_managers.set_cwd_to(project_path):
            model = Kratos.Model()
            for stage in self.stages_info.values():
                with open(f"{stage['base_name']}.json", "r") as f:
                    stage_parameters = Kratos.Parameters(f.read())
                analysis.GeoMechanicsAnalysis(model, stage_parameters).Run()


if __name__ == "__main__":
    model_runner = KratosGeoMechanicsCrowValidationWithStages()
    model_runner.run_simulation()
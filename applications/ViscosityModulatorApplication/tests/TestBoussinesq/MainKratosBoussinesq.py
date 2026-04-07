from BoussinesqGidCase.AnalysisStageWithFlush import CreateAnalysisStageWithFlushInstance
import KratosMultiphysics
import importlib

import scripts.utils as utils
import json
import os


def change_peclet(project_parameters, config):
    peclet = config["peclet"]
    buoyancy_materials_file_path = project_parameters["solver_settings"]["scalar_solver_settings"]["material_import_settings"]["materials_filename"].GetString()
    with open(buoyancy_materials_file_path, 'r') as parameter_file:
        data = json.load(parameter_file)
        data["properties"][0]["Material"]["Variables"]["CONDUCTIVITY"] = 1.0 / peclet
    with open(buoyancy_materials_file_path, 'w') as parameter_file:
        json.dump(data, parameter_file, indent=4)


def change_reynolds(project_parameters, config):
    reynolds = config["reynolds"]
    fluid_materials_file_path = project_parameters["solver_settings"]["fluid_solver_settings"]["material_import_settings"]["materials_filename"].GetString()
    with open(fluid_materials_file_path, 'r') as parameter_file:
        data = json.load(parameter_file)
        data["properties"][0]["Material"]["Variables"]["DYNAMIC_VISCOSITY"] = 1.0 / reynolds
    with open(buoyancy_materials_file_path, 'w') as parameter_file:
        json.dump(data, parameter_file, indent=4)


def run(config):
    file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(file_path)
    params_file = os.path.join(current_dir, "ProjectParametersBoussinesqCoupled.json")

    with open(params_file, 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

        utils.select_mesh(config, os.path.join(current_dir, "BoussinesqGidCase.mdpa"))
        change_peclet(parameters, config)
        utils.change_output_name_in_params(parameters, config)

    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

    global_model = KratosMultiphysics.Model()
    simulation = CreateAnalysisStageWithFlushInstance(analysis_stage_class, global_model, parameters, config)
    simulation.Run()


if __name__ == "__main__":
    run()

import argparse
import atexit
import sys


def _clean_up_api_object(api_object):
    del api_object


def run_flow_analysis(
    working_dir, project_parameters_filename, critical_head_boundary_model_part_name
) -> int:
    import KratosMultiphysics.GeoMechanicsApplication as GeoMechanicsApplication

    flow_api = GeoMechanicsApplication.CustomWorkflowFactory.CreateKratosGeoFlow()
    atexit.register(_clean_up_api_object, flow_api)

    critical_head_info = GeoMechanicsApplication.KratosExecuteCriticalHeadInfo(
        0.0, 0.0, 0.0
    )

    no_logging = lambda msg: None
    no_progress_reporting = lambda fraction_done: None
    no_progress_message = lambda msg: None
    do_not_cancel = lambda: False
    callbacks = GeoMechanicsApplication.KratosExecuteCallBackFunctions(
        no_logging, no_progress_reporting, no_progress_message, do_not_cancel
    )

    return flow_api.ExecuteFlowAnalysis(
        working_dir,
        project_parameters_filename,
        critical_head_info,
        critical_head_boundary_model_part_name,
        callbacks,
    )


def _main() -> int:
    parser = argparse.ArgumentParser(
        description="Run GeoMechanics backward erosion piping analysis"
    )
    parser.add_argument(
        "project_parameters_filename",
        metavar="PROJECT_PARAMETERS_FILENAME",
        help="'ProjectParameters.json' filename",
    )
    parser.add_argument(
        "-w", "--working-dir", help="working directory of the analysis run", default="."
    )
    parser.add_argument(
        "--critical-head-boundary-model-part-name",
        help="Model part name of the critical head boundary; may be empty",
        default="",
    )

    args = parser.parse_args()

    return run_flow_analysis(
        args.working_dir,
        args.project_parameters_filename,
        args.critical_head_boundary_model_part_name,
    )


if __name__ == "__main__":
    sys.exit(_main())

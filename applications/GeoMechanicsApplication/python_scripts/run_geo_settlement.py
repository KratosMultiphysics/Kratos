import argparse
import atexit
import sys

import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers


def _clean_up_api_object(api_object):
    del api_object


def run_stages(working_dir, project_parameters_filenames) -> int:
    import KratosMultiphysics.GeoMechanicsApplication as GeoMechanicsApplication

    no_logging = lambda msg: None
    no_progress_reporting = lambda fraction_done: None
    no_progress_message = lambda msg: None
    do_not_cancel = lambda: False

    settlement_api = (
        GeoMechanicsApplication.CustomWorkflowFactory.CreateKratosGeoSettlement()
    )
    atexit.register(_clean_up_api_object, settlement_api)

    with context_managers.set_cwd_to(working_dir):
        for filename in project_parameters_filenames:
            status = settlement_api.RunStage(
                working_dir,
                filename,
                no_logging,
                no_progress_reporting,
                no_progress_message,
                do_not_cancel,
            )
            if status != 0:
                print(f"Running analysis stage failed: status = {status}")
                return status

    return 0


def _main() -> int:
    parser = argparse.ArgumentParser(description="Run GeoMechanics settlement analysis")
    parser.add_argument(
        "-w", "--working-dir", help="working directory of the analysis run", default="."
    )
    parser.add_argument(
        "project_parameters_filenames",
        metavar="PROJECT_PARAMETERS_FILENAME",
        help="'ProjectParameters.json' filenames, one for each stage",
        nargs="+",
    )

    args = parser.parse_args()

    return run_stages(args.working_dir, args.project_parameters_filenames)


if __name__ == "__main__":
    sys.exit(_main())

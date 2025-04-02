import argparse
import sys


def clean_up_api_object(api_object):
    del api_object


def _main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-w", "--working-dir", help="working directory of the app", default="."
    )
    parser.add_argument(
        "project_parameters_file",
        help="'ProjectParameters.json' files of the stages",
        nargs="+",
    )

    args = parser.parse_args()

    import KratosMultiphysics.GeoMechanicsApplication as GeoMechanicsApplication

    no_logging = lambda msg: None
    no_progress_reporting = lambda fraction_done: None
    no_progress_message = lambda msg: None
    do_not_cancel = lambda: False

    import atexit

    settlement_api = (
        GeoMechanicsApplication.CustomWorkflowFactory.CreateKratosGeoSettlement()
    )
    atexit.register(clean_up_api_object, settlement_api)

    for project_parameters_filename in args.project_parameters_file:
        status = settlement_api.RunStage(
            args.working_dir,
            project_parameters_filename,
            no_logging,
            no_progress_reporting,
            no_progress_message,
            do_not_cancel,
        )
        if status != 0:
            print(f"Running analysis stage failed: status = {status}")
            return status

    return 0


if __name__ == "__main__":
    sys.exit(_main())

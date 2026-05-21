"""A script for creating XDMF files for results stored in HDF5.

license: HDF5Application/license.txt
"""


from argparse import ArgumentParser


from KratosMultiphysics.HDF5Application.xdmf_utils import WriteMultifileTemporalAnalysisToXdmf
from KratosMultiphysics.HDF5Application.xdmf_utils import WriteSinglefileTemporalAnalysisToXdmf


def main():
    """Parse the command line arguments and write the corresponding XDMF file."""
    parser = ArgumentParser(
        description="Write an XDMF file for post-processing results in HDF5.")
    parser.add_argument(dest="file_name", metavar="<filename>",
                        help="path to an HDF5 file for which XDMF metadata should be written")
    parser.add_argument("-t", "--type", dest="type", metavar="<type>",
                        choices=['single', 'multiple'], default="multiple", help="type of HDF5 file")
    parser.add_argument("-a", "--analysis", dest="analysis", metavar="<analysis>",
                        choices=['static', 'temporal'], default="temporal", help="type of analysis")
    parser.add_argument("-m", "--mesh-path", dest="mesh_path", metavar="<mesh-path>",
                        default="/ModelData", help="internal HDF5 file path to the mesh")
    parser.add_argument("-r", "--results-path", dest="results_path", metavar="<results-path>",
                        default="/ResultsData", help="internal HDF5 file path to the results")
    parser.add_argument("--require-results",
                        dest = "require_results",
                        action = "store_const",
                        default = False,
                        const = True,
                        help = "Ignore outputs that have mesh data but lack results.")
    print('\nCreate XDMF:')
    args = parser.parse_args()
    if args.type == "multiple" and args.analysis == "temporal":
        WriteMultifileTemporalAnalysisToXdmf(args.file_name,
                                             args.mesh_path,
                                             args.results_path)
    elif args.type == "single" and args.analysis == "temporal":
        WriteSinglefileTemporalAnalysisToXdmf(args.file_name,
                                              args.mesh_path,
                                              args.results_path,
                                              require_results = args.require_results)
    else:
        raise RuntimeError("Unsupported command line options.")


if __name__ == "__main__":
    main()

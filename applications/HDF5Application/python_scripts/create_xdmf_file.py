"""A script for creating XDMF files for results stored in HDF5.

license: HDF5Application/license.txt
"""


from argparse import ArgumentParser


"""A script for creating XDMF files for results stored in HDF5.

license: HDF5Application/license.txt
"""


from argparse import ArgumentParser

from KratosMultiphysics.HDF5Application.xdmf_utils import WriteDataSetsToXdmf
from KratosMultiphysics.HDF5Application.xdmf_utils import IdentifyPattern
from KratosMultiphysics.HDF5Application.core.xdmf_dataset_generator import SingleMeshMultiFileSameDatasetsGenerator
from KratosMultiphysics.HDF5Application.core.xdmf_dataset_generator import SingleFileDatasetsGenerator
from KratosMultiphysics.HDF5Application.core.xdmf_dataset_generator import MultiFileDatasetsGenerator
from KratosMultiphysics.HDF5Application.core.xdmf_dataset_generator import HasTags
from KratosMultiphysics.HDF5Application.core.xdmf_dataset_generator import GetDataSetPatterns


def main():
    """Parse the command line arguments and write the corresponding XDMF file."""
    parser = ArgumentParser(description="Write an XDMF file for post-processing results in HDF5.")

    parser.add_argument("dataset_pattern",
                        metavar="<dataset_pattern>",
                        help="The pattern of the dataset with prefix. Eg: test_<step>:/ResultsData, test.h5:/Results_<step>")
    parser.add_argument("-0",
                        "--output_name",
                        dest="output_name",
                        type=str,
                        default="output.xdmf",
                        help="Output filen name. (default = output.xdmf if a pattern is given.)")

    parser.add_argument("-p",
                        "--time_pos",
                        dest="temporal_tag_position",
                        type=int,
                        default=0,
                        help="The temporal value tag position. (default = 0)")
    parser.add_argument("-n",
                        "--no-identify",
                        dest="identify_pattern",
                        action="store_false",
                        default=True,
                        help="Use the provided <dataset_pattern> to identify a pattern string (default = on)")
    parser.add_argument("-m",
                        "--mesh-type",
                        dest="mesh_type",
                        choices=["single", "multiple"],
                        default="single",
                        help="Specify whether there is a \"single\" mesh or \"multiple\" meshes. (default = \"single\")")

    print('\nCreate XDMF:')
    args = parser.parse_args()
    dataset_pattern:str = args.dataset_pattern
    temporal_tag_position:int = args.temporal_tag_position
    mesh_type:str = args.mesh_type
    identify_pattern:bool = args.identify_pattern
    output_name:str = args.output_name

    # by default the dataset_pattern if not a pattern then, the pattern will be identified.
    # if it is already a pattern, then nothing is done.
    if dataset_pattern.find("<time>") != -1 or dataset_pattern.find("<step>") != -1:
        identify_pattern = False

    if identify_pattern:
        output_name = dataset_pattern.split(":")[0][:-2] + "xdmf"
        dataset_pattern, tag_type_dict = IdentifyPattern(dataset_pattern)
    else:
        tag_type_dict = {"<step>" : int, "<time>": float}

    # now check where the tags are
    h5_file_name, _ = GetDataSetPatterns(dataset_pattern)

    if HasTags(h5_file_name, tag_type_dict):
        if mesh_type == "single":
            dataset_generator = SingleMeshMultiFileSameDatasetsGenerator(dataset_pattern, temporal_tag_position, tag_type_dict)
        else:
            dataset_generator = MultiFileDatasetsGenerator(dataset_pattern, temporal_tag_position, tag_type_dict)
    else:
        dataset_generator = SingleFileDatasetsGenerator(dataset_pattern, temporal_tag_position, tag_type_dict)

    WriteDataSetsToXdmf(dataset_generator, output_name)

if __name__ == "__main__":
    main()

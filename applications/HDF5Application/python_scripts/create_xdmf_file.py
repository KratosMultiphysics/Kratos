"""A script for creating XDMF files for results stored in HDF5.

license: HDF5Application/license.txt
"""

"""A script for creating XDMF files for results stored in HDF5.

license: HDF5Application/license.txt
"""


from argparse import ArgumentParser

import KratosMultiphysics as Kratos
from KratosMultiphysics.HDF5Application.xdmf_utils import WriteDataSetsToXdmf
from KratosMultiphysics.HDF5Application.xdmf_utils import WriteMeshToXdmf
from KratosMultiphysics.HDF5Application.xdmf_utils import IdentifyPattern
from KratosMultiphysics.HDF5Application.core.dataset_generator import SingleMeshMultiFileSameDatasetsGenerator
from KratosMultiphysics.HDF5Application.core.dataset_generator import GenericDatasetsGenerator
from KratosMultiphysics.HDF5Application.core.dataset_generator import HasTags
from KratosMultiphysics.HDF5Application.core.dataset_generator import GetDataSetPatterns

def CreateXDMFFile(
        dataset_pattern: str,
        temporal_tag_position: int = 0,
        is_single_file: bool = False,
        is_mesh_only: bool = False,
        identify_pattern: bool = True,
        dynamic_mesh: bool = False,
        output_xdmf_file_name: str = "output.xdmf"):

    # by default the dataset_pattern if not a pattern then, the pattern will be identified.
    # if it is already a pattern, then nothing is done.
    if dataset_pattern.find("<time>") != -1 or dataset_pattern.find("<step>") != -1:
        identify_pattern = False

    h5_file_name, dataset_prefix = GetDataSetPatterns(dataset_pattern)

    if is_mesh_only:
        WriteMeshToXdmf(h5_file_name, output_xdmf_file_name, dataset_prefix)
    else:
        tag_type_dict = {"<step>": int, "<time>": float}
        if identify_pattern:
            output_xdmf_file_name = h5_file_name[:-2] + "xdmf"
            if is_single_file:
                dataset_prefix, tag_type_dict = IdentifyPattern(dataset_prefix)
                if not HasTags(dataset_prefix, tag_type_dict):
                    raise RuntimeError(f"Single file xdmf creation requires tags to be present in the dataset prefix.")
            else:
                h5_file_name, tag_type_dict = IdentifyPattern(h5_file_name)

        dataset_pattern = f"{h5_file_name}:{dataset_prefix}"

        if not dynamic_mesh and HasTags(h5_file_name, tag_type_dict) and not HasTags(dataset_prefix, tag_type_dict):
            Kratos.Logger.PrintInfo("XDMF", "Using SingleMeshMultiFileSameDatasetsGenerator.")
            dataset_generator = SingleMeshMultiFileSameDatasetsGenerator(dataset_pattern, temporal_tag_position, tag_type_dict)
        else:
            Kratos.Logger.PrintInfo("XDMF", "Using GenericDatasetsGenerator.")
            dataset_generator = GenericDatasetsGenerator(dataset_pattern, temporal_tag_position, tag_type_dict)

        WriteDataSetsToXdmf(dataset_generator, output_xdmf_file_name)

def main():
    """Parse the command line arguments and write the corresponding XDMF file."""
    parser = ArgumentParser(description="Write an XDMF file for post-processing results in HDF5.")

    parser.add_argument("dataset_pattern",
                        metavar="<dataset_pattern>",
                        help="The pattern of the dataset with prefix. Eg: test_<step>:/ResultsData, test.h5:/Results_<step>")
    parser.add_argument("-o",
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
    parser.add_argument("-i",
                        "--no-identify",
                        dest="identify_pattern",
                        action="store_false",
                        help="Use the provided <dataset_pattern> to identify a pattern string (default = on)")
    parser.add_argument("-m",
                        "--mesh-only",
                        dest="mesh_only",
                        action="store_true",
                        help="Find and write datasets (default = off)")
    parser.add_argument("-s",
                        "--single-file",
                        dest="single_file",
                        action="store_true",
                        help="Specify whether all data is present in one file. (default = off")
    parser.add_argument("-d",
                        "--dynamic-mesh",
                        dest="dynamic_mesh",
                        action="store_true",
                        help="Specify whether multiple meshes are present. (default = off")

    print('\nCreate XDMF:')
    args = parser.parse_args()
    dataset_pattern:str = args.dataset_pattern
    output_name:str = args.output_name
    temporal_tag_position:int = args.temporal_tag_position
    identify_pattern:bool = args.identify_pattern
    mesh_only:bool = args.mesh_only
    single_file: bool = args.single_file
    dynamic_mesh:bool = args.dynamic_mesh

    CreateXDMFFile(dataset_pattern, temporal_tag_position, single_file, mesh_only, identify_pattern, dynamic_mesh, output_name)

if __name__ == "__main__":
    main()
